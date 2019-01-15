from gpkit import Variable, Model, parse_variables, units
from gpkit import SignomialsEnabled, SignomialEquality, Vectorize
from gpkit.constraints.bounded import Bounded
from gpkit.constraints.tight import Tight

from nozzle import Nozzle, NozzlePerformance
from SRM import SRM

from relaxed_constants import relaxed_constants, post_process

import numpy as np

Tight.reltol = 1e-3

class Rocket(Model):
    """
    Basic rocket model

    Variables
    ---------
    t_T                  [s]     total burn time
    r                    [m]     shell radius
    l                    [m]     total length
    P_max                [Pa]    maximum chamber pressure

    Variables of length nsections
    -----------------------------
    A_fuel               [m^2]   initial fuel areas

    Variables of length nt
    ----------------------
    T                    [N]     thrust
    c_T                  [-]     thrust coefficient

    Variables of length (nsections, nt)
    -----------------------------------
    A_slack              [-]     slack variables

    Lower Unbounded
    ---------------
    t_T, T, c_T

    Upper Unbounded
    ---------------
    P_max, A_slack

    """

    def setup(self, nt, nsections):
        exec parse_variables(Rocket.__doc__)
        self.nt = nt
        self.nsections = nsections
        self.nozzle = Nozzle()
        with Vectorize(nt):
            self.nozzlePerformance = NozzlePerformance(self.nozzle)
            self.section = SRM(nsections)

        constraints = [
            # Limiting nozzle size to cross-sectional area
            self.nozzle.A_e <= np.pi*r**2,
            # Equal time segments
            self.section.dt == t_T/nt,
            # All fuel is consumed
            self.section.A_p_out[:,-1] >= 1e-20*np.ones(nsections)*np.pi*r**2,
            A_fuel == self.section.A_p_in[:,0],
            T == self.nozzlePerformance.T,
            c_T == self.nozzlePerformance.c_T,
            A_slack == self.section.A_slack,
        ]

        for i in range(nt-1):
            constraints += [
                # Decreasing fuel
                self.section.A_p_out[:, i] == self.section.A_p_in[:, i+1],
                # Decreasing sectional area
                self.section.A_in[:,i] <= self.section.A_in[:,i+1],
                self.section.A_out[:,i] <= self.section.A_out[:,i+1],
            ]

        for i in range(nt):
            constraints += [
                # Rocket length is section length,
                self.section.radius[i] == r,
                self.section.l[i] == l,
                # Matching nozzle and section conditions
                self.nozzlePerformance.mdot[i] == self.section.mdot_out[nsections-1, i],
                self.nozzlePerformance.T_t[i] == self.section.T_t_out[nsections-1, i],
                # Maximum chamber pressure
                self.section.P_chamb[:, i] <= P_max,
                self.nozzlePerformance.P_star[i] <= P_max,
            ]
            with SignomialsEnabled():
                constraints += [
                # Matching nozzle stagnation pressure
                self.nozzlePerformance.P_t[i] <= self.section.P_out[nsections-1, i] +
                                0.5*self.section.rho_out[nsections-1, i]*self.section.u_out[nsections-1, i]**2,
                ]

        return constraints, self.nozzle, self.nozzlePerformance, self.section

if __name__ == "__main__":
    nt = 2
    nsections = 5
    m = Rocket(nt, nsections)
    radius = 0.1*units('m')
    length = 2*units('m')
    m.substitutions.update({
        m.nozzle.k_A                                 :5,
        m.t_T                                        :1.5*units('s'),
        # m.l                                          :length,
        # m.r                                          :radius,
        m.P_max                                      :2*10.**8*units('Pa'),
        m.section.l_b_max                            :3*np.ones(nt),
        # m.nozzlePerformance.T                        :np.linspace(2e5,2e5,nt)*units('N'),
        m.A_fuel                                     :0.1*np.ones(nsections)*np.pi*radius**2,
    })
    for i in range(nt):
        m.substitutions.update()

    # m.cost = np.prod(m.section.A_slack**3)*np.sum(m.A_fuel)*m.l/nsections #np.prod(m.nozzlePerformance.T**-1)
    m.cost = np.prod(m.section.A_slack**3)*np.prod(m.nozzlePerformance.T**-1)
    # m = Model(m.cost, Bounded(m), m.substitutions)
    # m_relax = relaxed_constants(m)
    sol = m.localsolve(verbosity=4)
    post_process(sol)
