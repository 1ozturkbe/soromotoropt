import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d

from gpkit.small_scripts import mag


def draw_2D(sol, m, vectorvar):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    nsections = len(m.section.A_p_in)
    nt = len(m.section.A_p_in[0])
    x = np.linspace(1, nsections, nsections)
    y = np.linspace(1, nt, nt)
    X, Y = np.meshgrid(x,y)
    Z = sol(vectorvar)
    ax.contour3D(X, Y, Z, 50, cmap='binary')
    plt.show()

def draw1D(sol, m, var):
    fig = plt.figures()
    ax = plt.axes(projection='2d')
    x = np.linspace(1,len(var), len(var))
    plt.plot(x, mag(sol(var)))
    plt.show()

if __name__ == "__main__":
