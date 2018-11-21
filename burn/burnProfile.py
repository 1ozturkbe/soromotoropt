from gpkit import Variable, Model, parse_variables



class Nozzle(Model):
    """
    Defining the most basic rocket nozzle model

    Scalar variables
    -------------
    rhostar                           [kg/m^3]   density at throat
    ustar                             [m/s]      velocity at throat
    astar                             [m^2]      area of throat
    cstar                             [m/s]      characteristic velocity

    p_t                               [Pa]       stagnation pressure
    mdot                              [kg/s]     mass flow rate through nozzle
    T                                 [N]        thrust

    gm1                    0.3        [-]        gamma-1
    gp1                    2.3        [-]        gamma+1
    g                      1.3        [-]        gamma

    """
    def setup(self,nt):
        exec parse_variables(Nozzle.__doc__)
        self.nt = nt
        constraints = [
            cstar == p_t*astar/mdot,
            mdot == rhostar*ustar*astar,
            (p_t/pstar) == (gp1/2)**(g/gm1),
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

    def setup(self):
        exec parse_variables(Rocket.__doc__)

        constraints = []

        return constraints

def test_Nozzle():
    print "Testing Nozzle..."

if __name__ == "__main__":
    test_Nozzle()
