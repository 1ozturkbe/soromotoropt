from gpkit import Variable, Model, parse_variables, units
from gpkit import SignomialsEnabled, SignomialEquality, Vectorize

from nozzle import Nozzle, NozzlePerformance
from basicSRM import Section

from relaxed_constants import relaxed_constants, post_process

import numpy as np

class Rocket(Model):
    """
    Basic 1-section, many timestep Rocket model

    Variables
    ---------
    r                    [m]     shell radius
    l                    [m]     total length

    """

    def setup(self, nt):
        exec parse_variables(Rocket.__doc__)
        self.nt = nt
        self.nozzle = Nozzle()
        with Vectorize(nt):
            self.nozzlePerformance = NozzlePerformance(self.nozzle)
            self.section = Section()

        constraints = [
            # All fuel burnt
            self.section.A_p_out[nt-1] == 1e-20*units('m^2'),
        ]

        for i in range(nt-1):
            constraints += [
                # Decreasing fuel
                self.section.A_p_out[i] == self.section.A_p_in[i+1],
                # Equal time segments
                self.section.dt[i] == self.section.dt[i+1],
            ]

        for i in range(nt):
            constraints += [
                # Rocket length is section length,
                self.section.l[i] == l,
                # Setting inlet conditions
                self.section.P_in[i] == self.section.P_t_in[i],
                self.section.mdot_in[i] == 1e-20*units('kg/s'),
                self.section.T_t_in[i] == self.section.T_amb[i],
                # Matching nozzle and section conditions
                self.nozzlePerformance.mdot[i] == self.section.mdot_out[i],
                # Same area ratio (TODO: relax)
                self.section.A_in[i] == self.section.A_out[i],

        ]

        with SignomialsEnabled():
            constraints += [

            ]

        return constraints, self.nozzle, self.nozzlePerformance

if __name__ == "__main__":
    nt = 10
    m = Rocket(nt)
    radius = 0.1*units('m')
    length = 5*units('m')
    m.substitutions.update({
        m.l                       :length,
        m.section.A_p_in[0]       :0.5*np.pi*radius**2,
    })
    m.cost = 1/np.prod(m.nozzlePerformance.T)
    m_relax = relaxed_constants(m)
    sol = m_relax.localsolve()
