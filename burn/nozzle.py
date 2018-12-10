from gpkit import Variable, Model, parse_variables, units
from gpkit import SignomialsEnabled, SignomialEquality, Vectorize

from relaxed_constants import relaxed_constants, post_process

class Nozzle(Model):
    """ Nozzle model

    Variables
    -------------
    A_ratio                            [-]        nozzle area ratio
    A_star                             [m^2]      area of throat
    A_e                               [m^2]      nozzle exit area

    Upper Unbounded
    ---------------
    A_star, A_e, A_ratio

    Lower Unbounded
    ---------------
    A_star, A_e, A_ratio

    LaTex Strings
    -------------
    A_star        A_{star}
    A_ratio       A_{ratio}

    """
    def setup(self):
        exec parse_variables(Nozzle.__doc__)
        constraints = [
            # Area ratio
            A_e/A_star == A_ratio,
            ]
        return constraints

class NozzlePerformance(Model):
    """ Nozzle isentropic expansion model

    Variables
    ---------
    c_star                             [m/s]      characteristic velocity
    rho_star                           [kg/m^3]   density at throat
    u_star                             [m/s]      velocity at throat
    a_star                             [m/s]      speed of sound at throat
    P_star                             [Pa]       pressure at throat
    T_star                             [K]        temperature at throat
    P_t                               [Pa]       stagnation pressure
    T_t                               [K]        stagnation temperature
    mdot                              [kg/s]     mass flow rate
    T                                 [N]        thrust
    c_T                                [-]        thrust coefficient
    rho_e                             [kg/m^3]   density at exit
    u_e                               [m/s]      horizontal jet velocity
    a_e                               [m/s]      speed of sound at exit
    P_e                               [Pa]       exit static pressure
    T_e                               [K]        exit static temperature
    M_e                               [-]        jet Mach number
    chokeConstant          1.281      [-]        choked flow constant
    P_atm                  101000     [Pa]       atmospheric pressure
    c_p                    1000       [J/kg/K]   specific heat of air (constant pressure)
    c_v                    718        [J/kg/K]   specific heat of air (constant volume)
    R                      287        [J/kg/K]   gas constant of air
    gm1                    0.4        [-]        gamma-1
    gp1                    2.4        [-]        gamma+1
    g                      1.4        [-]        gamma

    Upper Unbounded
    ---------------
    rho_e, mdot, rho_star, P_e, M_e, P_star, c_star, P_t, T_t

    Lower Unbounded
    ---------------
    a_star, u_star, T_star, T_e, a_e, c_T

    LaTex Strings
    -------------
    c_star              c^*
    rho_star            \\rho^*
    u_star              u^*
    a_star              a^*
    P_star              P^*
    T_star              T^*
    mdot                \\dot{m}
    rho_e               \\rho_e
    chokeConstant       \\bar{c}
    P_atm               P_{atm}
    gm1                 (\\gamma-1)
    gp1                 (\\gamma+1)
    g                   \\gamma

    """
    def setup(self,nozzle):
        exec parse_variables(NozzlePerformance.__doc__)
        constraints = [
            # Mass flow rate
            mdot == rho_star*u_star*nozzle.A_star,
            mdot == rho_e*u_e*nozzle.A_e,
            # Characteristic velocity
            c_star == P_t*nozzle.A_star/mdot,
            c_star**2 == 1/1.4*(2.4/2)**(2.4/.4)*R*T_t,
            # Thrust coefficient
            T == mdot*c_star*c_T,
            # Universal gas law
            P_e == rho_e*R*T_e,
            P_star == rho_star*R*T_star,
            # Stagnation pressure and static pressure at throat
            (P_t/P_star) == (2.4/2)**(1.4/0.4),
            # Choked throat relations
            mdot*c_p**0.5*T_t**0.5/nozzle.A_star/P_t == chokeConstant,
            u_star/a_star == 1,
            # Stagnation temperature to velocity at exit
            2*c_p*T_e + u_e**2 <= 2*c_p*T_t,
            # Exit Mach number
            M_e == u_e/a_e,
            M_e >= 1,
            # Speed of sound
            a_e**2 == g*R*T_e,
            a_star**2 == g*R*T_star,
            # Exit pressure
            P_e == rho_e*R*T_e,
            (T_e/T_t)**(1.4/.4) == P_e/P_t,
            ]
        with SignomialsEnabled():
            constraints += [
                # Thrust
                T <= mdot*u_e + (P_e - P_atm)*nozzle.A_e,
                # Speed of sound at throat
                SignomialEquality(T_star + u_star**2/(2*c_p),T_t),
            ]
        return constraints

class Rocket(Model):
    """
    Basic integrated Rocket model

    Variables
    ---------

    """

    def setup(self, nt):
        exec parse_variables(Rocket.__doc__)
        self.nt = nt
        self.nozzle = Nozzle()
        with Vectorize(nt):
            self.nozzlePerformance = NozzlePerformance(self.nozzle)
        constraints = []

        return constraints, self.nozzle, self.nozzlePerformance

def test_Nozzle():
    print "Testing Nozzle..."
    m = Rocket(1)
    m.substitutions.update({
                            'T_t':         1800*units('K'),
                            'P_t':         5000*units('kPa'),
                            'mdot':        1*units('kg/s'),
                            })
    m.cost = 1/sum(m.nozzlePerformance.T)
    m_relax = relaxed_constants(m, None, [m.nozzlePerformance.chokeConstant])
    sol = m_relax.localsolve(verbosity=2)
    post_process(sol)
    print sol.table()


if __name__ == "__main__":
    test_Nozzle()
