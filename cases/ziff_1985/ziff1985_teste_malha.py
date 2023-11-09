from numpy import arange, linspace, array, piecewise, zeros, ones, diff, insert
from itertools import cycle
import sys
import os.path as path
import pickle

dir = path.dirname(__file__)
if __name__ == "__main__":
    sys.path.append(path.abspath(path.join(dir, "..\\..")))
from pbe.solvers.moc_ramk import MOCSolution


"""
Case setup based on:

    Ziff, R.M and McGrady, E.D
    "New solutions to the fragmentation equation", J. Phys A, vol. 18, 1985

    We are looking at case 4. from their paper with kernel F(x, y) = x + y.
    This corresponds to choosing:
        beta = 1.0/y,
        Gamma = y^2,
    for our kernels.

    Teste com mesma malha utilizada pelo bettersize
"""

grids = [74]  # Número de classes utilizadas na discretização
time = arange(0.0, 100.0, 0.01)  # Tempo e passo
v0 = 7.87011e-05/1e9
dv0 = 2.65025e-05
r = 1.445124894

#v0 = 1.0
#r = 1
caso = 3
vmax = 1.0

# Caso 4 funcionou bem
#v0 = 0.0001
#dv0 = v0

pbe_solutions = dict()

for g in grids:
    # This is modelling Dirac's delta
    # initial number density function
    N0 = zeros(g)
    N0[-1] = 1

    if caso == 1:
        xi = v0 * r ** arange(g)
        dxi = dv0 * r ** arange(g)
    elif caso == 2:
        dxi = dv0 * 1.1 ** arange(g)
        dxi = diff(xi)
        dxi1 = insert(dxi, 0, dxi[0])
        dxi2 = insert(dxi, -1, dxi[-1])
        dxi = (dxi1 + dxi2)/2
        dxi[0] = dxi[0]/2

    elif caso == 3:
        xi = v0 + v0 * arange(g)
        xi = linspace(v0, vmax, g, endpoint=True)

    elif caso == 4:
        xi = v0 * r ** arange(g)

    pbe_solutions[g] = MOCSolution(
        g,  # number of classes
        time,
        xi=xi,
        N0=N0,
        beta=lambda x, y: 1.0 / y,  # DDSD
        gamma=lambda x: x**2,  # breakage rate
    )

pbe_solutions["vmax"] = xi[-1]  # Max volume
pbe_solutions["v0"] = xi[-1]  # Initial volume
pbe_solutions["time"] = time
pbe_solutions["grid"] = grids

with open(path.join(dir, "ziff1985_teste_malha.pickle"), "wb") as f:
    pickle.dump(pbe_solutions, f)
