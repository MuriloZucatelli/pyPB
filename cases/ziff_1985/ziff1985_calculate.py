from numpy import arange, linspace, array, piecewise
from itertools import cycle
import sys
import os.path as path
import pickle

dir = path.dirname(__file__)
if __name__ == "__main__":
    sys.path.append(path.abspath(path.join(dir, "..\\..")))
from pbe.solvers.moc import MOCSolution


"""
Case setup based on:

    Ziff, R.M and McGrady, E.D
    "New solutions to the fragmentation equation", J. Phys A, vol. 18, 1985

    We are looking at case 4. from their paper with kernel F(x, y) = x + y.
    This corresponds to choosing:
        beta = 1.0/y,
        Gamma = y^2,
    for our kernels.
"""

grids = [10, 40, 160]  # Número de classes utilizadas na discretização
grids = [10]
time = arange(0.0, 1000.0, 0.01)  # Tempo e passo
# N0 = 1
vmax = 1.0  # Volume máximo
v0 = 1.0
pbe_solutions = dict()

for g in grids:
    threshold = vmax / g / 2

    # This is modelling Dirac's delta
    # initial number density function
    def n0_init(x):
        return piecewise(
            x, [x < v0 - threshold, x >= v0 - threshold], [0, g / vmax]
        )  # 0 ou g/vmax   função delta de dirac

    def n0_init2(x):
        return piecewise(
            x, [abs(v0 - x) < threshold, v0 - x >= threshold], [g / vmax, 0]
        )  # 0 ou g/vmax   função delta de dirac generica

    pbe_solutions[g] = MOCSolution(
        g,  # number of classes
        time,
        vmax / g,  # dxi
        n0=n0_init2,
        beta=lambda x, y: 1.0 / y,  # DDSD   era 2.0 / y
        gamma=lambda x: x**2,  # breakage rate
    )

pbe_solutions["vmax"] = vmax
pbe_solutions["v0"] = v0
pbe_solutions["time"] = time
pbe_solutions["grid"] = grids

with open(path.join(dir, "ziff1985.pickle"), "wb") as f:
    pickle.dump(pbe_solutions, f)
