from gpkit import Variable, Model, parse_variables, units
from gpkit import SignomialsEnabled, SignomialEquality, Vectorize
from gpkit.constraints.tight import Tight
from gpkit.constraints.bounded import Bounded

from relaxations import relaxed_constants, post_process

import numpy as np

Tight.reltol = 1e-2


class SRM(Model):

    """ basic SRM model

    Variables
    ---------
    radius                            [m]          radius of section
    l                                 [m]          total length
    l_sec                             [m]          section length
    dt                                [s]          time step
    P_end                             [Pa]         endcap pressure
    mdot_out                          [kg/s]       mass flow rate out
    T_t_out                           [K]          stagnation temperature out
    T_out                             [K]          static temperature out
    u_out                             [m/s]        flow velocity out
    P_out                             [Pa]         static pressure at outlet
    rho_out                           [kg/m^3]     density out
    l_b_max              3            [-]          maximum burn length factor
    k_A_max              1.5          [-]          maximum area ratio
    R                    287          [J/kg/K]     gas constant of air
    T_amb                273          [K]          ambient temperature
    P_amb                1e6          [Pa]         ambient pressure
    r_c                  5.606        [mm/s]       burn rate coefficient
    r_k                  0.05         [1/(m/s)]    erosive burn rate coefficient
    rho_p                1700         [kg/m^3]     propellant density
    k_comb_p             1.23e6       [J/kg]       propellant specific heat of combustion
    c_p                  1000         [J/kg/K]     specific heat of combustion products

    Variables of length n
    ---------------------
    k_A                               [-]          area ratio
    A_in                              [m^2]        area in
    A_out                             [m^2]        area out
    A_avg                             [m^2]        average area
    A_b                               [m^2]        burn area
    A_p_in                            [m^2]        initial propellant area
    A_p_out                           [m^2]        final propellant area
    l_b                               [m]          avg burn length
    mdot                              [kg/s]       mass flow rate
    rho                               [kg/m^3]     density
    P                                 [Pa]         static pressure
    P_chamb                           [Pa]         section static pressure
    V_chamb                           [m^3]        section volume
    T_t                               [K]          stagnation temperature
    T                                 [K]          static temperature
    u                                 [m/s]        flow velocity
    r                                 [mm/s]       burn rate
    q                                 [kg/s]       rate of generation of products

    Upper Unbounded
    ---------------
    radius, T_out

    Lower Unbounded
    ---------------
    dt, A_p_out, mdot_out, P_out

    LaTex Strings
    -------------
    radius                   \mathrm{r}
    l                        \mathrm{l}
    l_sec                    l_{sec}
    dt                       \delta t
    l_b_max                  l_{b,\mathrm{max}}
    R                        R
    T_amb                    T_{\mathrm{amb}}
    rho_p                    \rho_p
    k_comb_p                 k_{p,\mathrm{comb}}
    k_A_max                  k_{A,\mathrm{max}}
    A_in                     A_{\mathrm{in}}
    A_out                    A_{\mathrm{out}}
    A_avg                    A_{avg}
    A_slack                  A_{slack}
    A_p_in                   A_{p,\mathrm{in}}
    A_p_out                  A_{p,\mathrm{out}}
    mdot_out                 \dot{m}_{\mathrm{out}}
    mdot                     \dot{m}
    rho_out                  \rho_{\mathrm{out}}
    P_out                    P_{\mathrm{out}}
    P_chamb                  P_{\mathrm{chamb}}
    V_chamb                  V_{\mathrm{chamb}}
    T_t_out                  T_{t_{\mathrm{out}}}
    T_out                    T_{\mathrm{out}}
    u_out                    u_{\mathrm{out}}

    """

    def setup(self, n):
        # where n is the number of sections
        exec parse_variables(SRM.__doc__)
        constraints = [
                    # Bounding area ratio
                    k_A[:] >= 1/k_A_max,
                    k_A[:] <= k_A_max,
                    # Mass flow
                    q[0] == mdot[0], #massCons
                    # Momentum cons
                    Tight([P_end*A_in[0] >= P[0]*A_out[0] + q[0]*u[0]/V_chamb[0]*l_sec*A_out[0]], name='momCons', printwarning=True),
                    # Setting section lengths,
                    n*l_sec == l,
                    # Exit quantities
                    u_out == u[n-1],
                    T_t_out == T_t[n-1],
                    T_out == T[n-1],
                    mdot_out == mdot[n-1],
                    rho_out == rho[n-1],
                    P_out == P[n-1],
         ]
        with SignomialsEnabled():
            constraints += [
                # Inlet temperatures
                Tight([T_t[0]*mdot[0] <= q[0]*T_amb + q[0]*k_comb_p/c_p], name='energyCons' , printwarning=True),
                T_t >= T_amb,
                # Chamber pressure
                P_chamb[0]**2 == P_end*P[0],
                # # Burn rate
                Tight([r[0] >= r_c * (P_chamb[0]/1e6*units('1/Pa'))** 0.35 * (1 + 0.5*r_k*u[0])], name='burnRate', printwarning=True),

            ]

        for i in range(1,n):
            constraints += [
                # Coupling
                A_in[i]    == A_out[i-1],
                # Pressure increase
                P[i] <= P[i-1],
                # Chamber pressure
                P_chamb[i]**2 == P[i-1]*P[i],
            ]

            with SignomialsEnabled():
                constraints += [
                # Flow acceleration (conservation of momentum)
                # Note: assumes constant rate of burn through the chamber
                Tight([P[i-1]*A_in[i] + rho[i-1]*u[i-1]**2*A_in[i] >= P[i]*A_out[i] + q[i]*u[i]/V_chamb[i]*l_sec*A_out[i] +
                    rho[i-1]*u[i-1]*u[i]*A_out[i]], name='momCons', printwarning=True),
                # Temperatures
                Tight([T_t[i]*mdot[i] <= mdot[i-1]*T_t[i-1] + q[i]*T_amb + q[i]*k_comb_p/c_p], name='energyCons', printwarning=True), #
                # Mass flows
                Tight([mdot[i-1] + q[i] >= mdot[i]], name='massCons', printwarning=True),
                # # Burn rate (Saint-Robert's Law, coefficients taken for Space Shuttle SRM)
                Tight([r[i] >= r_c * (P_chamb[i]/1e6*units('1/Pa'))**0.35 * (1 + 0.5*r_k*(u[i]+u[i-1]))], name='burnRate', printwarning=True)
                ]

        # for i in range(n):
        constraints += [
                # Volume of chamber
                V_chamb == A_avg*l_sec,
                # Area ratio
                A_in / A_out == k_A,
                # Mass flow rate
                rho * u * A_out == mdot,
                # Product generation rate
                q == rho_p * A_b * r,
                # Burn area
                A_b == l_b * l,
                # Conservation of fuel
                Tight([A_p_out + r*l_b*dt <= A_p_in], name='fuel', printwarning=True),
                # Ideal gas law
                P == rho*R*T,
                # Making sure burn surface length is feasible
                l_b >= 2*np.pi**0.5*(A_avg)**0.5,
                l_b <= l_b_max*2*np.pi**0.5*(A_avg)**0.5,
                # Constraining flow speed
                u**2 <= 1.4*R*T,
                ]
        with SignomialsEnabled():
            constraints += [
                    # Taking averages (with a slack variable)
                    A_in*A_out == A_avg**2,
                    # Stagnation quantities
                    Tight([T_t <= (T + u**2/(2*c_p))], name='Tt', printwarning=True),
                    # Constraining areas
                    SignomialEquality(A_avg + A_p_in, np.pi*radius**2),
                ]
        return constraints

if __name__ == "__main__":
    n = 5
    m = SRM(n)
    radius = 20*units('cm')
    m.substitutions.update({
        m.k_A_max         :5,
        m.radius          :radius,
        # m.A_p_in         :0.1*np.ones(n)*np.pi*radius**2,
        m.A_p_out         :1e-20*np.ones(n)*np.pi*radius**2,
        m.mdot_out        :5*units('kg/s'),
        m.T_t_out         :1500*units('K'),
        m.P_out           :2e7*units('Pa'),
        m.u_out           :50*units('m/s'),
        m.dt              :0.25*units('s'),
    })
    m.cost = np.sum(m.A_p_in)*m.l
    m_relax = relaxed_constants(m)
    sol = m_relax.localsolve(reltol = 1e-5, verbosity=4)
