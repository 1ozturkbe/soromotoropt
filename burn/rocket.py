from gpkit import Variable, Model, parse_variables, units
from gpkit import SignomialsEnabled, SignomialEquality, Vectorize
from gpkit.constraints.bounded import Bounded
from gpkit.constraints.tight import Tight

from nozzle import Nozzle, NozzlePerformance
from basicSRM import Section

from relaxed_constants import relaxed_constants, post_process

import numpy as np

Tight.reltol = 1e-3

class Rocket(Model):
    """
    Basic 1-section, many timestep Rocket model

    Variables
    ---------
    t_T                  [s]     total burn time
    c_TV                 [-]     thrust variation coefficient
    r                    [m]     shell radius
    l                    [m]     total length
    P_max                [Pa]    maximum chamber pressure
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
            # Limiting nozzle size
            self.nozzle.A_e <= np.pi*r**2,
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
                # Thrust variation coefficient
                c_TV >= np.prod(self.nozzlePerformance.T)**(1./nt)/self.nozzlePerformance.T[i],
                # Rocket length is section length,
                self.section.radius[i] == r,
                self.section.l[i] == l,
                # Setting inlet conditions
                self.section.P_in[i] == self.section.P_t_in[i],
                self.section.mdot_in[i] == 1e-20*units('kg/s'),
                self.section.u_in[i] == 1e-20*units('m/s'),
                self.section.T_t_in[i] == self.section.T_amb[i],
                # self.section.T_in[i] == self.section.T_amb[i],
                # Matching nozzle and section conditions
                self.nozzlePerformance.mdot[i] == self.section.mdot_out[i],
                self.nozzlePerformance.P_t[i] == self.section.P_t_out[i],
                self.nozzlePerformance.T_t[i] == self.section.T_t_out[i],
                # self.nozzlePerformance.u_star[i] >= self.section.u_out[i],
                # self.nozzlePerformance.P_star[i] <= self.section.P_out[i],
                # Same area ratio (TODO: relax)
                self.section.A_in[i] == self.section.A_out[i],
                # Dumb bounds (cutting planes)
                self.section.A_p_in[i] >= self.section.A_p_out[i],
                # Maximum chamber pressure
                self.section.P_chamb <= P_max,

        ]

        with SignomialsEnabled():
            constraints += [

            ]

        return constraints, self.nozzle, self.nozzlePerformance, self.section

if __name__ == "__main__":
    nt = 10
    m = Rocket(nt)
    radius = 1*units('m')
    length = 0.2*units('m')
    m.substitutions.update({
        m.c_TV                                       :2,
        # m.nozzle.A_ratio                              :5,
        m.l                                          :length,
        m.r                                          :radius,
        m.P_max                                      :10.**8*units('Pa'),
        m.t_T                                        :100*units('s'),
        m.section.A_p_in[0]                          :0.5*np.pi*radius**2,
        m.section.T_amb                              :300*np.ones(nt)*units('K'),
    })
    m.cost = 1/np.prod(m.nozzlePerformance.T*m.section.dt)
    m = Model(m.cost, Bounded(m), m.substitutions)
    m_relax = relaxed_constants(m)
    sol = m_relax.localsolve()
    post_process(sol)
