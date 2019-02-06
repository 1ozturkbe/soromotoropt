# This file includes methods to calculate the dimensions of regular pointed stars that can
# represent the circumference and area distributions given by the optimization model

from gpkit import Model, Variable, Vectorize, SignomialsEnabled, SignomialEquality
from gpkit.constraints.tight import Tight
import numpy as np

class Star(Model):
    def setup(self, nPoints):
        r_i = Variable('r_i', 'm', 'internal radius')
        r_o = Variable('r_o', 'm', 'outer radius')
        r_ii= Variable('r_{ii}', 'm', 'radius of biggest internal circle')
        A = Variable('A', 'm^2', 'total internal area')
        C = Variable('C', 'm', 'total circumference')
        S_e = Variable('S', 'm', 'semi-perimeter of outer triangles')
        l_e = Variable('l_e', 'm', 'internal polyhedron edge length')
        h_e = Variable('h_e', 'm', 'external polyhedron height')
        e = Variable('e', 'm', 'external triangle external edge length')
        A_e = Variable('A_e', 'm^2', 'area of external triangle')
        cos_ai = Variable('cos(a_i)', np.cos(np.pi/nPoints),  '-', 'cosine of internal triangle half angle')
        sin_ai = Variable('sin(a_i)', np.sin(np.pi/nPoints),  '-', 'sine of internal triangle half angle')

        with SignomialsEnabled():
            constraints = [
                r_ii/r_i == cos_ai,
                l_e/(2*r_i) == sin_ai,
                r_o >= r_ii + h_e,
                Tight([A >= 0.5*nPoints*(r_ii*l_e) + nPoints*A_e], printwarning=True),
                Tight([A_e >= 0.5*(h_e*l_e)], printwarning=True),
                # Heron's formula for outer triangle outer edge
                S_e >= 0.5*(l_e + 2*e),
                Tight([A_e**2 >= S_e*(0.5*l_e)*(0.5*l_e)*(S_e-l_e)], printwarning=True),
                C == 2.*nPoints*e,
            ]
        return constraints


def return_radii(nPoints, sol, mrocket):
    nt = len(sol('A_avg')[0])
    nx = len(sol('A_avg')[:,0])
    with Vectorize(nt):
        with Vectorize(nx):
            m = Star(nPoints)
    linkingConstraints = []
    for i in m.variables_byname('r_o'):
        linkingConstraints+= [i <= sol(mrocket.r)]
    for j in range(nt-1):
            linkingConstraints += [
                m['r_i'][:,j] <= m['r_i'][:,j+1],
                m['r_o'][:,j] <= m['r_o'][:,j+1],
                m['r_{ii}'][:,j] <= m['r_{ii}'][:,j+1],
            ]
    cost = np.prod(m['r_o'])/np.prod(m['r_{ii}']*m['h_e'])
    m = Model(cost, [m, linkingConstraints], m.substitutions)
    m.substitutions.update({'A': sol('A_avg'),
                            'C': sol('l_b')})
    r_sol = m.localsolve(verbosity=4)
    return r_sol

def plot_star(nPoints, r_sol, title='Star'):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    nsections = len(r_sol('A'))
    nt = len(r_sol('A')[0])
    r_o = mag(r_sol('r_o'))
    r_i = mag(r_sol('r_i'))
    x = np.linspace(1, nsections, nsections)
    theta = np.linspace(0, 2*np.pi, 2*nPoints + 1)
    X, Theta = np.meshgrid(x,theta)
    for i in range(nt):
        r_o_s = r_o[:,i]
        r_i_s = r_i[:,i]
        O, Theta = np.meshgrid(r_o_s, theta)
        I, Theta = np.meshgrid(r_i_s, theta)
        Z = O + np.ceil(np.mod(Theta, 2*np.pi/nPoints))*I
        ax.contour3D(X, Z*np.sin(Theta), Z*np.cos(Theta), 50, cmap='binary')
    plt.xlabel('Time step')
    plt.ylabel('Axial coordinate')
    plt.title(title)
    plt.show()

    # This import registers the 3D projection, but is otherwise unused.
# from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
#
# import matplotlib.pyplot as plt
# import numpy as np
#
#
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
#
# # Create the mesh in polar coordinates and compute corresponding Z.
# r = np.linspace(0, 1.25, 50)
# p = np.linspace(0, 2*np.pi, 50)
# R, P = np.meshgrid(r, p)
# Z = ((R**2 - 1)**2)
#
# # Express the mesh in the cartesian system.
# X, Y = R*np.cos(P), R*np.sin(P)
#
# # Plot the surface.
# ax.plot_surface(X, Y, Z, cmap=plt.cm.YlGnBu_r)
#
# # Tweak the limits and add latex math labels.
# ax.set_zlim(0, 1)
# ax.set_xlabel(r'$\phi_\mathrm{real}$')
# ax.set_ylabel(r'$\phi_\mathrm{im}$')
# ax.set_zlabel(r'$V(\phi)$')
#
# plt.show()

if __name__ == "__main__":
    nPoints = 5
    r_sol = return_radii(5, sol, m)
