import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d

from gpkit.small_scripts import mag


def draw_2D(sol, m, vectorvar, title):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    nsections = len(m.section.A_p_in)
    nt = len(m.section.A_p_in[0])
    x = np.linspace(1, nt, nt)
    y = np.linspace(1, nsections, nsections)
    X, Y = np.meshgrid(x,y)
    Z = sol(vectorvar)
    ax.contour3D(X, Y, Z, 50, cmap='binary')
    plt.xlabel('Time step')
    plt.ylabel('Axial coordinate')
    plt.title(title)
    plt.show()

def draw_2D_bar(sol, m, vectorvar, title):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    nsections = len(m.section.A_p_in)
    nt = len(m.section.A_p_in[0])
    x = np.linspace(1, nsections, nsections)
    y = np.linspace(1, nt, nt)
    X, Y = np.meshgrid(x,y)
    Z = sol(vectorvar)
    for i in range(nt):
        color = [0.5, 0., 1.0*i/nt]
        ax.bar(x, mag(Z[:,i]), zs = i, zdir='y', alpha=0.5, color = color)
    plt.xlabel('Axial coordinate')
    plt.ylabel('Time step')
    plt.title(title)
    plt.show()


def draw1D(sol, m, var, title):
    fig = plt.figure()
    ax = plt.axes()
    nsections = len(m.section.A_p_in)
    nt = len(m.section.A_p_in[0])
    if len(var) == nt:
        n = nt
        plt.xlabel('Time step')
    elif len(var) == nsections:
        n = nsections
        plt.xlabel('Axial coordinate')
    else:
        print 'Warning: axis not recognized'
    x = np.linspace(1, n, n)
    plt.plot(x, sol(var))
    plt.title(title)
    plt.grid()
    plt.show()

def draw_fuel(sol,m):
    fig = plt.figure()
    ax = plt.axes()
    nsections = len(m.section.A_p_in)
    nt = len(m.section.A_p_in[0])
    z = sol(m.section.A_p_in)
    fig, ax = plt.subplots(1,nsections)
    for i in range(nsections):
        ax[i].pie(z[i]/sum(z[i]), labels = np.linspace(1,nt,nt), radius = mag(z[i,0])/max(z.flat))
        ax[i].axis('equal')
    plt.show()

if __name__ == "__main__":
    draw_2D(sol, m, m.section.A_in, 'Section inlet area evolution')
    draw_2D(sol, m, m.section.A_in, 'Section outlet area evolution')
    draw_2D(sol, m, m.section.A_p_in, 'Section propellant area evolution')
    draw_2D(sol, m, m.section.r, 'Section burn rate evolution')
    draw_2D(sol, m, m.section.mdot_out, 'Section mass flow evolution')
    draw_2D(sol, m, m.section.T_t_out, 'Section stagnation temperature evolution')
    draw_2D_bar(sol, m, m.section.A_p_in, 'Section propellant evolution')

    draw1D(sol, m, m.T_target, 'Thrust evolution')
