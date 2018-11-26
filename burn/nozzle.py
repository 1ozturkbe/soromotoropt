from gpkit import Variable, Model, parse_variables, units
from gpkit import SignomialsEnabled, SignomialEquality, Vectorize

from relaxed_constants import relaxed_constants, post_process

class Nozzle(Model):
    """ Nozzle model

    Variables
    -------------
    Aratio                            [-]        nozzle area ratio
    Astar                             [m^2]      area of throat
    A_e                               [m^2]      nozzle exit area
    """
    def setup(self):
        exec parse_variables(Nozzle.__doc__)
        constraints = [
            # Area ratio
            A_e/Astar == Aratio,
            ]
        return constraints

class NozzlePerformance(Model):
    """ Nozzle isentropic expansion model

    Variables
    ---------
    cstar                             [m/s]      characteristic velocity
    rhostar                           [kg/m^3]   density at throat
    ustar                             [m/s]      velocity at throat
    astar                             [m/s]      speed of sound at throat
    pstar                             [Pa]       pressure at throat
    Tstar                             [K]        temperature at throat
    p_t                               [Pa]       stagnation pressure
    T_t                               [K]        stagnation temperature
    mdot                              [kg/s]     mass flow rate
    T                                 [N]        thrust
    cT                                [-]        thrust coefficient
    rho_e                             [kg/m^3]   density at exit
    u_e                               [m/s]      horizontal jet velocity
    a_e                               [m/s]      speed of sound at exit
    p_e                               [Pa]       exit static pressure
    T_e                               [K]        exit static temperature
    M_e                               [-]        jet Mach number
    chokeConstant          1.281      [-]        choked flow constant
    p_atm                  101000     [Pa]       atmospheric pressure
    c_p                    1000       [J/kg/K]   specific heat of air (constant pressure)
    c_v                    718        [J/kg/K]   specific heat of air (constant volume)
    R                      287        [J/kg/K]   gas constant of air
    gm1                    0.4        [-]        gamma-1
    gp1                    2.4        [-]        gamma+1
    g                      1.4        [-]        gamma
    """
    def setup(self,nozzle):
        exec parse_variables(NozzlePerformance.__doc__)
        constraints = [
            # Mass flow rate
            mdot == rhostar*ustar*nozzle.Astar,
            mdot == rho_e*u_e*nozzle.A_e,
            # Characteristic velocity
            cstar == p_t*nozzle.Astar/mdot,
            cstar**2 == 1/1.4*(2.4/2)**(2.4/.4)*R*T_t,
            # Thrust coefficient
            T == mdot*cstar*cT,
            # Universal gas law
            p_e == rho_e*R*T_e,
            # Stagnation pressure and static pressure at throat
            (p_t/pstar) == (2.4/2)**(1.4/0.4),
            # Choked throat relations
            mdot*c_p**0.5*T_t**0.5/nozzle.Astar/p_t == chokeConstant,
            # Stagnation temperature to velocity at exit
            2*c_p*T_e + u_e**2 <= 2*c_p*T_t,
            # Exit Mach number
            M_e == u_e/a_e,
            M_e >= 1,
            # Speed of sound
            a_e**2 == g*R*T_e,
            # Exit pressure
            (T_e/T_t)**(1.4/.4) == p_e/p_t,
            ]
        with SignomialsEnabled():
            constraints += [
                # Thrust
                T <= mdot*u_e + (p_e - p_atm)*nozzle.A_e,
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
    m.substitutions.update({'Aratio':      10,
                            'T_t':         1800*units('K'),
                            'p_t':         5000*units('kPa'),
                            'mdot':        1*units('kg/s'),
                            'rhostar':     1.5*units('kg/m^3')})
    m.cost = 1/sum(m.nozzlePerformance.T)
    m_relax = relaxed_constants(m, None, [m.nozzlePerformance.chokeConstant])
    sol = m_relax.localsolve(verbosity=2)
    post_process(sol)

if __name__ == "__main__":
    pass
