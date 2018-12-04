from gpkit import Variable, Model, parse_variables, units
from gpkit import SignomialsEnabled, SignomialEquality, Vectorize
from gpkit.constraints.tight import Tight

from relaxed_constants import relaxed_constants, post_process

import numpy as np


class Section(Model):
    """ basic SRM section model

    Variables
    -------------
    radius                            [m]          radius of section
    Aratio                            [-]          area ratio
    A_in                              [m^2]        area in
    A_out                             [m^2]        area out
    A_avg                             [m^2]        average area
    A_b                               [m^2]        burn area
    A_p_in                            [m^2]        initial propellant area
    A_p_out                           [m^2]        final propellant area
    l_b                               [m]          avg burn length
    l_b_max                           [-]          maximum burn length factor
    l                                 [m]          section length
    mdot_in                           [kg/s]       mass flow rate in
    mdot_out                          [kg/s]       mass flow rate out
    rho_in                            [kg/m^3]     density in
    rho_out                           [kg/m^3]     density out
    P_t_in                            [Pa]         stagnation pressure in
    P_t_out                           [Pa]         stagnation pressure out
    P_in                              [Pa]         static pressure at inlet
    P_out                             [Pa]         static pressure at outlet
    dP                                [Pa]         pressure increase
    P_chamb                           [Pa]         chamber static pressure
    V_chamb                           [m^3]        chamber volume
    T_t_in                            [K]          stagnation temperature in
    T_t_out                           [K]          stagnation temperature out
    T_in                              [K]          static temperature in
    T_out                             [K]          static temperature out
    u_in                              [m/s]        velocity in
    u_out                             [m/s]        velocity out
    u_avg                             [m/s]        average velocity
    r                                 [mm/s]       burn rate
    q                                 [kg/s]       rate of generation of products
    dt                                [s]          time step
    R                    287          [J/kg/K]     gas constant of air
    T_amb                273          [K]          ambient temperature
    r_c                  5.606        [mm/s]       burn rate coefficient
    r_k                  0.05         [1/(m/s)]    erosive burn rate coefficient
    rho_p                1700         [kg/m^3]     propellant density
    k_comb_p             1.23e6       [J/kg]       propellant specific heat of combustion
    c_p                  1000         [J/kg/K]     specific heat of combustion products
    """

    def setup(self):
        exec parse_variables(Section.__doc__)
        constraints = [
            # Taking averages,
            u_in*u_out == u_avg**2,
            A_in*A_out == A_avg**2,
            # Volume of chamber
            V_chamb == A_avg*l,
            # Constraining areas
            Tight([A_avg + A_p_in <= np.pi*radius**2]),
            # Area ratio
            A_in / A_out == Aratio,
            # Mass flow rate
            rho_in * u_in * A_in == mdot_in,
            rho_out * u_out * A_out == mdot_out,
            # Chamber pressure
            P_chamb**2 == P_in*P_out,
            # Product generation rate
            q == rho_p * A_b * r,
            A_b == l_b * l,
            A_p_out + r*l_b*dt <= A_p_in,
            # Stagnation quantities
            P_t_in >= P_in + 0.5*rho_in*u_in**2,
            T_t_in >= T_in + u_in**2/(2*c_p),
            # Ideal gas law
            P_out == rho_out*R*T_out,
            # Pressure increase
            P_out + dP <= P_in,
            # Making sure burn surface length is feasible
            l_b >= 2*np.pi**0.5*(A_avg)**0.5,
            l_b <= l_b_max*2*np.pi**0.5*(A_avg)**0.5,

        ]
        with SignomialsEnabled():
            constraints += [
                # Flow acceleration (conservation of momentum)
                # Note: assumes constant rate of burn through the chamber
                dP >= q*u_out/V_chamb*(2./3.*l) + rho_in*u_in*(u_out - u_in),
                u_out >= u_in,
                # Burn rate (Saint-Robert's Law, coefficients taken for Space Shuttle SRM)
                SignomialEquality(r, r_c * (P_chamb/1e6*units('1/Pa')) ** 0.35 * (1 + r_k*u_avg)),
                # Mass flows
                SignomialEquality(mdot_in + q, mdot_out),
                mdot_out >= mdot_in,
                # Temperatures
                Tight([T_t_out*mdot_out <= mdot_in*T_t_in + q*T_amb + q*k_comb_p/c_p]),
                # Stagnation quantities
                Tight([P_t_out <= P_out + 0.5*rho_out*u_out**2]),
                Tight([T_t_out <= T_out + u_out**2/(2*c_p)]),

            ]
        return constraints

if __name__ == "__main__":
    m = Section()
    radius = 10*units('cm')
    m.substitutions.update({
        m.radius          :radius,
        m.l_b_max         :3,
        m.P_t_in          :1000*units('kPa'),
        m.l               :10*units('cm'),
        m.T_t_in          :700*units('K'),
        m.T_in            :690*units('K'),
        m.mdot_in         :0.5*units('kg/s'),
        m.A_p_in          :0.1*np.pi*radius**2,
        m.A_p_out         :0.08*units('m^2'),
        m.A_in            :0.9*np.pi*radius**2,
        m.A_out           :0.9*np.pi*radius**2,
        m.dt              :0.01*units('s'),
    })
    m.cost = 1/(m.mdot_out*m.dt*m.T_t_out*m.P_t_out)
    m_relax = relaxed_constants(m)
    sol = m_relax.localsolve()
