from numpy import arange, linspace, piecewise
from itertools import cycle
import sys
import os.path as path

dir = path.dirname(__file__)
if __name__ == "__main__":
    sys.path.append(path.abspath(path.join(dir, "..\\..")))
from pbe.solvers.moc import MOCSolution
from tests.test_moc import blatz_and_tobolsky_pbe_solution
import matplotlib.pyplot as plt
from pbe.setup.helpers import plt_config2

"""
Case setup based on:

    Blatz, B.J and Tobolsky, V.
    "Note on the Kinetic of systems manifesting simulaneous
    polymerization-depolymerization phenomena"

    B&T equations capture the evolution of polymer sized from a distribution of
    monomers. The kernels give constant coefficients for breakup and
    coalescence. Special care needs to be taken as the case is discrete.

    NOTE: B&T equations are normalised with respect to initial number of
    monomers N0. Because of the difference with our formulation the coalescence
    rates have to be divided by N0 to keep the problems equivalent. This is
    because density function appears twice.
"""
plt_config2(relative_fig_width=0.7)

time = arange(0.0, 20, 0.005)
# Initial number of monomers
N0 = 10000
grids = [30]
kc = 1.0
kb = 0.25
vmax = 100

sol = dict()


def N0Init(x):
    return piecewise(x, [x < 2, x > 2], [N0, 0])


for g in grids:
    sol[g] = MOCSolution(
        g,
        time,
        1.0,
        N0=N0Init,
        # Dividing coalescence coefficient by the number of monomers to make
        # formulations equivalent
        Q=lambda x, y: kc / N0,
        # Guard has to be imlemented in order to avoid division by zero
        beta=lambda x, y: 1.0 / max([y - 1.0, 1e-6]),
        gamma=lambda x: kb * (x - 1.0),
    )


fig = plt.figure()
ax = fig.gca()
markers = cycle(["o", "s", "v", "*", ".", ","])
v = linspace(0.9, vmax, 10000)

ax.set_ylim([1e-10, 10e3])
for n in sorted(sol):
    ax.loglog(
        sol[n].xi,
        sol[n].N[-1],
        marker=next(markers),
        ms=6,
        label="MOC with N={0}".format(n),
    )
ax.loglog(
    v,
    N0 * blatz_and_tobolsky_pbe_solution(v, time[-1], kc, kb),
    "-k",
    linewidth=2,
    label="Analytical $t=\infty$",
)
ax.legend(loc="upper right", shadow=True)
ax.set_xlabel("Volume")
ax.set_ylabel("N/N0")
plt.show()
