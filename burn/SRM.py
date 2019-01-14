from gpkit import Variable, Model, parse_variables, units
from gpkit import SignomialsEnabled, SignomialEquality, Vectorize
from gpkit.constraints.tight import Tight
from gpkit.constraints.bounded import Bounded

from relaxed_constants import relaxed_constants, post_process

import numpy as np

Tight.reltol = 1e-3


class SRM(Model):

    """ basic SRM model

    Variables
    ---------
    radius                            [m]          radius of section
    l                                 [m]          total length
    l_sec                             [m]          section length
    dt                                [s]          time step
    l_b_max              3            [-]          maximum burn length factor
    k_A_max              5            [-]          maximum area ratio
    R                    287          [J/kg/K]     gas constant of air
    T_amb                273          [K]          ambient temperature
    r_c                  5.606        [mm/s]       burn rate coefficient
    r_k                  0.05        [1/(m/s)]    erosive burn rate coefficient
    rho_p                1700         [kg/m^3]     propellant density
    k_comb_p             1.23e6       [J/kg]       propellant specific heat of combustion
    c_p                  1000         [J/kg/K]     specific heat of combustion products

    Variables of length n
    ---------------------
    k_A                               [-]          area ratio
    A_in                              [m^2]        area in
    A_out                             [m^2]        area out
    A_avg                             [m^2]        average area
    A_slack                           [-]          area slack variable
    A_b                               [m^2]        burn area
    A_p_in                            [m^2]        initial propellant area
    A_p_out                           [m^2]        final propellant area
    l_b                               [m]          avg burn length
    mdot_out                          [kg/s]       mass flow rate out
    rho_out                           [kg/m^3]     density out
    P_out                             [Pa]         static pressure at outlet
    dP                                [Pa]         pressure increase
    P_chamb                           [Pa]         section static pressure
    V_chamb                           [m^3]        section volume
    T_t_out                           [K]          stagnation temperature out
    T_out                             [K]          static temperature out
    u_out                             [m/s]        velocity out
    r                                 [mm/s]       burn rate
    q                                 [kg/s]       rate of generation of products

    Upper Unbounded
    ---------------
    A_slack

    Lower Unbounded
    ---------------
    dt, A_p_out

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
    rho_out                  \rho_{\mathrm{out}}
    P_out                    P_{\mathrm{out}}
    dP                       \delta P
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
                    # Flow acceleration (conservation of momentum)
                    # Note: assumes constant rate of burn through the chamber
                    dP[0] == q[0]*u_out[0]/V_chamb[0]*l_sec,
                    # Mass flow
                    q[0] == mdot_out[0],
                    # Setting section lengths,
                    n*l_sec    == l,
                    # Getting geometric average of propellant remaining
                    # A_p_avg**n >= np.prod(A_p_out),
         ]
        with SignomialsEnabled():
            constraints += [
                # Inlet temperatures
                T_t_out[0]*mdot_out[0] <= q[0]*T_amb + q[0]*k_comb_p/c_p,
                T_t_out >= T_amb,
                T_t_out[0] >= T_out + u_out[0]**2/(2*c_p),
                # Tight([A_slack[0]*T_t_out[0]*mdot_out[0] >= q[0]*T_amb + q[0]*k_comb_p/c_p]),
                # Tight([A_slack[0]*P_chamb[0] >= P_out[0] + 0.25*rho_out[0]*u_out[0]**2]),
                SignomialEquality(P_chamb[0], P_out[0] + 0.25*rho_out[0]*u_out[0]**2),
                # # Burn rate
                # Tight([r[0] >= r_c * (P_chamb[0]/1e6*units('1/Pa'))** 0.35 * (1 + 0.5*r_k*u_out[0])]),

            ]

        for i in range(1,n):
            constraints += [
                # Coupling
                mdot_out[i] >= mdot_out[i-1],
                A_in[i]    == A_out[i-1],
                # Pressure increase
                P_out[i] + dP[i] <= P_out[i-1],
                # Chamber pressure
                P_chamb[i]**2 == P_out[i-1]*P_out[i],
                # Constraining pressures
                P_chamb[i] <= P_chamb[i-1],
            ]
            with SignomialsEnabled():
                constraints += [
                # Flow acceleration (conservation of momentum)
                # Note: assumes constant rate of burn through the chamber
                dP[i] + rho_out[i-1]*u_out[i-1]**2 >= q[i]*u_out[i]/V_chamb[i]*l_sec +
                    rho_out[i-1]*u_out[i-1]*u_out[i],
                # Temperatures
                Tight([T_t_out[i]*mdot_out[i] <= mdot_out[i-1]*T_t_out[i-1] + q[i]*T_amb + q[i]*k_comb_p/c_p]), #
                # Mass flows
                Tight([mdot_out[i-1] + q[i] >= mdot_out[i]]),
                # # Burn rate (Saint-Robert's Law, coefficients taken for Space Shuttle SRM)
                # Tight([r[i] >= r_c * (P_chamb[i]/1e6*units('1/Pa'))** 0.35 * (1 + 0.5*r_k*(u_out[i]+u_out[i-1]))]),
                ]

        # for i in range(n):
        constraints += [
                # Volume of chamber
                V_chamb == A_avg*l_sec,
                # Area ratio
                A_in / A_out == k_A,
                # Mass flow rate
                rho_out * u_out * A_out == mdot_out,
                # Product generation rate
                q == rho_p * A_b * r,
                A_b == l_b * l,
                A_p_out + r*l_b*dt <= A_p_in,
                # Ideal gas law
                P_out == rho_out*R*T_out,
                # Making sure burn surface length is feasible
                l_b >= 2*np.pi**0.5*(A_avg)**0.5,
                l_b <= l_b_max*2*np.pi**0.5*(A_avg)**0.5,
                # Constraining flow speed
                u_out**2 <= 1.4*R*T_out,
                ]
        with SignomialsEnabled():
            constraints += [
                    # Taking averages (with a slack variable)
                    A_in*A_out == A_avg**2,
                    A_slack >= 1,
                    # Burn rate (Saint-Robert's Law, coefficients taken for Space Shuttle SRM)
                    Tight([r == r_c * (P_chamb/1e6*units('1/Pa'))** 0.35]), # * (1 + 0.5*r_k*(u_in+u_out))]),
                    # Stagnation quantities
                    Tight([T_t_out <= T_out + u_out**2/(2*c_p)]),
                    # Constraining areas
                    Tight([A_avg + A_p_in <= A_slack*np.pi*radius**2]),
                    Tight([(A_avg + A_p_in)*A_slack >= np.pi*radius**2]),
                ]
        return constraints

if __name__ == "__main__":
    n = 5
    m = SRM(n)
    radius = 10*units('cm')
    m.substitutions.update({
        m.k_A_max         :5,
        m.radius          :radius,
        m.l               :200*units('cm'),
        m.A_p_in          :0.1*np.ones(n)*np.pi*radius**2,
        m.A_p_out         :0.05*np.ones(n)*np.pi*radius**2,
        # m.k_A         :np.ones(n),
        # m.A_in            :0.9*np.pi*radius**2,
        # m.A_out           :0.9*np.pi*radius**2,
        m.dt              :0.25*units('s'),
    })
    m.cost = 1/(m.mdot_out[n-1]*m.T_t_out[n-1])*np.prod(m.A_slack**3)
    m = Model(m.cost, Bounded(m), m.substitutions)
    sol = m.localsolve(reltol = 1e-3)
