from numpy import arange, linspace, exp, piecewise, where, array, delete, zeros
from itertools import cycle
import sys
import os.path as path

dir = path.dirname(__file__)
if __name__ == "__main__":
    sys.path.append(path.abspath(path.join(dir, "..\\..")))
from pbe.solvers.moc_ramk import MOCSolution
from pbe.setup.helpers import plt_config2
from tests.test_moc import scott_total_number_solution3
from tests.test_moc import scott_pbe_solution3
import matplotlib.pyplot as plt

plt_config2(relative_fig_width=0.7)

g = 3
time = arange(0.0, 100, 0.001)
C = 0.1  # constante de coalescencia

N0 = zeros(g)
N0[-1] = 2  # initial Number of droplets
sol = dict()

xi = array([1, 3, 5])



sol[g] = MOCSolution(g, time, xi=xi, N0=N0, Q=lambda x, y: C)
