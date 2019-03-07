import numpy as np
from gpkit.small_scripts import mag
from gpkit import units

from pint import UnitRegistry


def allocate_fuel(sol, model, relaxDict, nt, nx):
    rho_p = sol('rho_p')
    # Assuming accelerant is the same density for now
    rho_a = rho_p
    # Beta matrix
    beta = np.array((nt,nx,3))
    # Porosity of fuel
    porosity = relaxDict['massCons']
    # Accelerant to propellant ratio
    alpha_a2p = relaxDict['burnRate']
    # Filler to (accelerant+propellant) ratio
    alpha_f2p = relaxDict['energyCons']*(np.ones((nt,nx))-porosity)
    # Trying to figure out the effect of energy conservation
    # on momentum conservation
    E1o2 = (2*(np.ones((nt,nx))+relaxDict['energyCons'])/(np.ones((nt,nx))-porosity))**0.5
    dv1o2 = E1o2/(np.ones((nt,nx))-porosity)

    T_amb = sol(m.section.T_amb)
    T = sol(m.section.T)
    u = sol(m.section.u)
    c_p = sol(m.section.c_p)
    T_t = [[] for i in range(nx)]
    for i in range(nt):
        for j in range(nx):
            T_t[j].append((T[j,i] + u[j,i]**2/(2*c_p[i])))
    mdot = sol(m.section.mdot)
    q = sol(m.section.q)
    k_comb_p = sol(m.section.k_comb_p)

    # Computing relative burn rate...
    # Tight([T_t[i]*mdot[i] <= mdot[i-1]*T_t[i-1] + q[i]*T_amb + q[i]*k_comb_p/c_p], name='energyCons', printwarning=True), #
    q_comb = np.zeros((nt,nx))
    for i in range(nt):
        for j in range(nx):
            if j >= 1:
                q_comb[i,j] = np.round(mag(((T_t[j][i])*mdot[j,i] - q[j,i]*T_amb[i] - T_t[j][i]*mdot[j,i-1])*c_p[i]/k_comb_p[i]),decimals=2)
            else:
                q_comb[i,j] = np.round(mag((T_t[j][i]*mdot[j,i] - q[j,i]*T_amb[i])*c_p[i]/k_comb_p[i]),decimals=2)

    # k_comb_rel = np.zeros((nt,nx))
    # for i in range(nt):
    #     k_comb_rel[i,:] = k_comb_act[i,:]/k_comb_p[i]
