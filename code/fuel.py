import numpy as np
from gpkit.small_scripts import mag

def allocate_fuel(sol, m, relaxDict, nt, nx):
    """
    :param sol: rocket solution
    :param m: rocket model
    :param relaxDict: dictionary of relaxations of model
    :param nt: time steps
    :param nx: spatial discretization
    :return: mass proportions of propellant, accelerant and filler material
    """
    # Mapping solution values
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

    # Porosity of fuel
    porosity = relaxDict['massCons']

    # Computing relative burn rate...
    q_comb = [[] for i in range(nx)]
    for i in range(nt):
        for j in range(nx):
            if j >= 1:
                q_comb[j].append(np.round(mag(((T_t[j][i])*mdot[j,i] - q[j,i]*T_amb[i] - T_t[j][i]*mdot[j,i-1])*c_p[i]/k_comb_p[i]),decimals=2))
            else:
               q_comb[j].append(np.round(mag((T_t[j][i]*mdot[j,i] - q[j,i]*T_amb[i])*c_p[i]/k_comb_p[i]),decimals=2))

    # Since this isn't working properly, we hack and offset q_comb
    q_min = np.min(q_comb)
    qratmin = np.min(mag(q_comb/q))
    q_comb = q_comb-(-0.1+qratmin)*mag(q)
    qrat = mag(q_comb/q) # Sum of propellant+accelerant ratios
    beta_f = np.ones((nx,nt))-qrat-porosity # Filler ratio

    # Use burn rate relaxation to obtain beta_a and beta_p
    beta_p = 1/(np.ones((nx,nt))+relaxDict['burnRate'])*qrat
    beta_a = np.ones((nx,nt))-beta_p-beta_f
    return beta_p, beta_a, beta_f, porosity
