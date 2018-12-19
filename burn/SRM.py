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
    l_b_max                           [-]          maximum burn length factor
    R                    287          [J/kg/K]     gas constant of air
    T_amb                273          [K]          ambient temperature
    r_c                  5.606        [mm/s]       burn rate coefficient
    r_k                  0.05         [1/(m/s)]    erosive burn rate coefficient
    rho_p                1700         [kg/m^3]     propellant density
    k_comb_p             1.23e6       [J/kg]       propellant specific heat of combustion
    c_p                  1000         [J/kg/K]     specific heat of combustion products
    A_p_avg                           [m^2]        average sectional propellant

    Variables of length n
    ---------------------
    A_ratio                           [-]          area ratio
    A_in                              [m^2]        area in
    A_out                             [m^2]        area out
    A_avg                             [m^2]        average area
    A_b                               [m^2]        burn area
    A_p_in                            [m^2]        initial propellant area
    A_p_out                           [m^2]        final propellant area
    l_b                               [m]          avg burn length
    mdot_in                           [kg/s]       mass flow rate in
    mdot_out                          [kg/s]       mass flow rate out
    rho_in                            [kg/m^3]     density in
    rho_out                           [kg/m^3]     density out
    P_t_in                            [Pa]         stagnation pressure in
    P_t_out                           [Pa]         stagnation pressure out
    P_in                              [Pa]         static pressure at inlet
    P_out                             [Pa]         static pressure at outlet
    dP                                [Pa]         pressure increase
    P_chamb                           [Pa]         section static pressure
    V_chamb                           [m^3]        section volume
    T_t_in                            [K]          stagnation temperature in
    T_t_out                           [K]          stagnation temperature out
    T_in                              [K]          static temperature in
    T_out                             [K]          static temperature out
    u_in                              [m/s]        velocity in
    u_out                             [m/s]        velocity out
    r                                 [mm/s]       burn rate
    q                                 [kg/s]       rate of generation of products
    """

    def setup(self, n):
        # where n is the number of sections
        exec parse_variables(SRM.__doc__)
        constraints = [
                    # Initializing inlet constraints
                    mdot_in[0] == 1e-20*units('kg/s'),
                    u_in[0]    == 1e-20*units('m/s'),
                    rho_in[0]  == 1e-20*units('kg/(m^3)'),
                    P_in[0]    == P_t_in,
                    T_t_in[0]  == T_amb,
                    T_in[0]    == T_amb,
                    # Flow acceleration (conservation of momentum)
                    # Note: assumes constant rate of burn through the chamber
                    dP[0] == q[0]*u_out[0]/V_chamb[0]*(2./3.*l_sec),
                    # Mass flow
                    q[0] == mdot_out[0],
                    # Setting section lengths,
                    n*l_sec    == l,
                    # Getting geometric average of propellant remaining
                    A_p_avg**n >= np.prod(A_p_out),
         ]
        with SignomialsEnabled():
            constraints += [
                # Inlet temperatures
                Tight([T_t_out[0]*mdot_out[0] <= q[0]*T_amb + q[0]*k_comb_p/c_p]),
            ]

        for i in range(1,n):
            constraints += [
                # Mass flow rate
                rho_in * u_in * A_in == mdot_in,
                # Coupling
                mdot_in[i] == mdot_out[i-1],
                P_in[i]    == P_out[i-1],
                P_t_in[i]  == P_t_out[i-1],
                T_in[i]    == T_out[i-1],
                T_t_in[i]  == T_t_out[i-1],
                rho_in[i]  == rho_out[i-1],
                A_in[i]    == A_out[i-1],
                u_in[i]    == u_out[i-1],
            ]
            with SignomialsEnabled():
                constraints += [
                # Flow acceleration (conservation of momentum)
                # Note: assumes constant rate of burn through the chamber
                dP[i] + rho_in[i]*u_in[i]**2 >= q[i]*u_out[i]/V_chamb[i]*(2./3.*l_sec) +
                    rho_in[i]*u_in[i]*u_out[i],
                # Temperatures
                Tight([T_t_out[i]*mdot_out[i] <= mdot_in[i]*T_t_in[i] + q[i]*T_amb + q[i]*k_comb_p/c_p]),
                # Mass flows
                Tight([mdot_in[i] + q[i] >= mdot_out[i]]),
                ]

        # for i in range(n):
        constraints += [
                # Taking averages,
                A_in*A_out == A_avg**2,
                # Volume of chamber
                V_chamb == A_avg*l_sec,
                # Area ratio
                A_in / A_out == A_ratio,
                # Mass flow rate
                rho_out * u_out * A_out == mdot_out,
                # Product generation rate
                q == rho_p * A_b * r,
                A_b == l_b * l,
                A_p_out + r*l_b*dt <= A_p_in,
                # Ideal gas law
                P_out == rho_out*R*T_out,
                # Pressure increase
                P_out + dP <= P_in,
                # Making sure burn surface length is feasible
                l_b >= 2*np.pi**0.5*(A_avg)**0.5,
                l_b <= l_b_max*2*np.pi**0.5*(A_avg)**0.5,
                # Chamber pressure
                P_chamb**2 == P_in*P_out,
                # Constraining flow speed
                u_out**2 <= 1.4*R*T_out,
                ]
        with SignomialsEnabled():
            constraints += [
                    # Burn rate (Saint-Robert's Law, coefficients taken for Space Shuttle SRM)
                    Tight([r == r_c * (P_chamb/1e6*units('1/Pa'))]), #** 0.35 * (1 + 0.5*r_k*(u_in+u_out))]),
                    # Mass flows
                    mdot_out >= mdot_in,
                    # Stagnation quantities
                    Tight([P_t_out <= P_out + 0.5*rho_out*u_out**2]),
                    Tight([T_t_out <= T_out + u_out**2/(2*c_p)]),
                    # Constraining areas
                    SignomialEquality(A_avg + A_p_in, np.pi*radius**2),
                ]
        return constraints

if __name__ == "__main__":
    n = 3
    m = SRM(n)
    radius = 10*units('cm')
    m.substitutions.update({
        m.radius          :radius,
        m.l_b_max         :3,
        # m.P_t_out[-1]     :2000*units('kPa'),
        m.l               :200*units('cm'),
        m.A_p_in          :0.1*np.ones(n)*np.pi*radius**2,
        m.A_p_out         :0.05*np.ones(n)*np.pi*radius**2,
        # m.A_ratio         :np.ones(n),
        # m.A_in            :0.9*np.pi*radius**2,
        # m.A_out           :0.9*np.pi*radius**2,
        m.dt              :0.01*units('s'),
    })
    m.cost = 1/(m.mdot_out[-1]*m.T_t_out[-1]*m.P_t_out[-1])
    m = Model(m.cost, Bounded(m), m.substitutions)
    m_relax = relaxed_constants(m)
    sol = m_relax.localsolve(reltol = 1e-3)
