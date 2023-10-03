from numpy import arange, linspace, exp, piecewise
from itertools import cycle
import sys
import os.path as path

dir = path.dirname(__file__)
if __name__ == "__main__":
    sys.path.append(path.abspath(path.join(dir, "..\\")))
from pbe.solvers.moc import MOCSolution
from tests.test_moc import scott_total_number_solution3
from tests.test_moc import scott_pbe_solution3
import matplotlib.pyplot as plt

"""
Case setup based on:

    Scott, W.T.
    "Analytical studies in cloud droplet coalescence", I. J. Atmos. Sci.
    vol. 25, 1968

    We are looking at constant coalescence kernels (case III in the paper).

    NOTE: We use a simplified version of the "Gaussian-like" distribution of
    initial IC.
"""

grids = [20, 40, 80, 160]  # number os classes
time = arange(0.0, 10, 0.001)
vmax = 2e1  # max volume
C = 0.1  # constante de coalescencia
N0 = 2  # initial Number of droplets
v0 = 0.5  # initial volume

pbe_solutions = dict()


# Distribuição inicial Gaussiana
# Number density function
def n0_init(v):
    # return (N0 / v0) * v/v
    return (N0 / v0) * (v / v0) * exp(-v / v0)


import matplotlib.pyplot as plt
import numpy as np

plt.plot(
    np.arange(v0, vmax, vmax / 160),
    n0_init(np.arange(v0, vmax, vmax / 160)),
    label="Ninit",
)
for g in grids:
    pbe_solutions[g] = MOCSolution(g, time, vmax / g, n0=n0_init, Q=lambda x, y: C)
plt.legend()
plt.show()
totals = dict((n, pbe_solutions[n].total_numbers) for n in pbe_solutions)
volume = dict((n, pbe_solutions[n].total_volume) for n in pbe_solutions)
print(volume)
v = linspace(0, vmax, 100)
Na = scott_total_number_solution3(time, C=C, N0=N0)

fig = plt.figure()
ax = fig.gca()
linestyles = cycle(["-", "--", ":"])
for n in sorted(totals):
    ax.plot(
        time,
        totals[n] / totals[n][0],
        linestyle=next(linestyles),
        label="MOC with N={0}".format(n),
    )
ax.plot(time, Na / N0, "-k", linewidth=2, label="Analytical")
ax.legend(loc="lower right", shadow=True)
ax.set_xlabel("t")
ax.set_ylabel("N/N0")
plt.show()


fig = plt.figure()
ax = fig.gca()
markers = cycle(["o", "s", "v", "*", ".", ","])
v = linspace(0, vmax, 100)

for n in sorted(pbe_solutions):
    ax.semilogy(
        pbe_solutions[n].xi,
        pbe_solutions[n].number_density[-1],
        "+",
        marker=next(markers),
        label="MOC with N={0}".format(n),
    )
ax.semilogy(
    v,
    scott_pbe_solution3(v, time[-1], C=C, xi0=2.0 * v0, N0=N0),
    "-k",
    linewidth=2,
    label="Analytical $t=\infty$",
)
ax.legend(loc="lower left", shadow=True)
ax.set_xlabel("Particle volume")
ax.set_ylabel("Number density function")
plt.show()
