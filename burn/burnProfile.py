from gpkit import Variable, Model, parse_variables, units
from gpkit import SignomialsEnabled, SignomialEquality

class Nozzle(Model):
    """
    Defining isentropic rocket nozzle model

    Variables
    -------------
    Aratio                            [-]        nozzle area ratio
    A_e                               [m^2]      nozzle exit area
    rhostar                           [kg/m^3]   density at throat
    ustar                             [m/s]      velocity at throat
    Astar                             [m^2]      area of throat
    pstar                             [Pa]       pressure at throat
    Tstar                             [K]        temperature at throat
    p_t                               [Pa]       stagnation pressure
    T_t                               [K]        stagnation temperature
    mdot                              [kg/s]     mass flow rate through nozzle
    T                                 [N]        thrust
    chokeConstant          1.281      [-]        choked flow constant
    rho_e                             [kg/m^3]   density at exit
    p_e                               [Pa]       exit static pressure
    V_e                               [m/s]      horizontal jet velocity
    M_e                               [-]        jet Mach number
    T_e                               [K]        exit static temperature
    a_e                               [m/s]      speed of sound at exit
    p_atm                  101000     [Pa]       atmospheric pressure
    c_p                    1000       [J/kg/K]   specific heat of air (constant pressure)
    c_v                    718        [J/kg/K]   specific heat of air (constant volume)
    R                      287        [J/kg/K]   gas constant of air
    gm1                    0.4        [-]        gamma-1
    gp1                    2.4        [-]        gamma+1
    g                      1.4        [-]        gamma
    dummy                  1.         [-]        unit
    """
    def setup(self):
        exec parse_variables(Nozzle.__doc__)
        constraints = [
            # Mass flow rate
            mdot == rhostar*ustar*Astar,
            mdot == rho_e*V_e*A_e,
            # Universal gas law
            pstar == rhostar*R*Tstar,
            p_e == rho_e*R*T_e,
            # Area ratio
            A_e/Astar == Aratio,
            # Stagnation pressure and static pressure at throat
            (p_t/pstar) == (2.4/2)**(1.4/0.4),
            # Choked throat relations
            mdot*c_p**0.5*T_t**0.5/Astar/p_t == chokeConstant,
            # Stagnation temperature to velocity at exit
            2*c_p*T_e + V_e**2 <= 2*c_p*T_t,
            # Exit Mach number
            M_e == V_e/a_e,
            # Exit speed of sound
            a_e**2 == g*R*T_e,
            ]
        with SignomialsEnabled():
            constraints += [
                # Thrust
                T <= mdot*V_e + (p_e - p_atm)*A_e,
            ]
        return constraints

class Rocket(Model):
    """
    Basic integrated Rocket model

    Scalar variables
    -------------
    rshell                            [cm]       radius of shell
    lshell                            [cm]       length of shell



    Variables of length nt
    -------------

    Variables of length nt-1
    -------------

    """

    def setup(self, nt):
        exec parse_variables(Rocket.__doc__)
        self.nt = nt
        constraints = []

        return constraints

def test_Nozzle():
    print "Testing Nozzle..."
    m = Nozzle()
    m.substitutions.update({'Aratio':      5,
                            'T_t':         1000*units('K'),
                            'p_t':         500*units('kPa')})
    m.cost = m.T
    sol = m.localsolve(verbosity = 4)

if __name__ == "__main__":
    test_Nozzle()