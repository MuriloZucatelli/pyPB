from numpy import arange, linspace, piecewise, geomspace, zeros
from itertools import cycle
import sys
import os.path as path
import seaborn as sns
dir = path.dirname(__file__)
if __name__ == "__main__":
    sys.path.append(path.abspath(path.join(dir, "..\\..")))
from pbe.solvers.moc_ramk import MOCSolution as MOCramk
from pbe.solvers.moc import MOCSolution as MOChidy
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
cor = sns.color_palette("husl", n_colors=5)
plt_config2(relative_fig_width=0.7)
time = arange(0.0, 1000, 0.1)
# Initial number of monomers
N0 = 10000
grids = [50]
kc = 1.0
kb = 0.25
vmax = 100

sol_ramk = dict()
sol_hidy = dict()


def N0Init(x):
    return piecewise(x, [x < 2, x > 2], [N0, 0])


def beta(v, vl, vast, g):
    """Função quantidade de partículas na quebra

    Args:
        v (_type_): v
        vl (_type_): v'
        vast (_type_): v* volume unitário
        g (_type_): grid
    """
    for i in g:
        if v - vast * i == 0:
            return 1.0 / (vl - 1.0)
        else:
            return 0


for g in grids:
    xi = linspace(1, vmax, g, endpoint=True)
    N = zeros(g)
    N[0] = N0
    sol_ramk[g] = MOCramk(
        g,
        time,
        dxi=1.0,
        N0=N,
        # Dividing coalescence coefficient by the number of monomers to make
        # formulations equivalent
        Q=lambda x, y: kc / N0,
        # Guard has to be imlemented in order to avoid division by zero
        beta=lambda x, y: 1.0 / max([y - 1.0, 1e-10]),  #
        gamma=lambda x: kb * (x - 1.0),
    )

for g in grids:
    xi = linspace(1, vmax, g, endpoint=True)
    N = zeros(g)
    N[0] = N0
    sol_hidy[g] = MOChidy(
        g,
        time,
        dxi=1.0,
        N0=N,
        # Dividing coalescence coefficient by the number of monomers to make
        # formulations equivalent
        Q=lambda x, y: kc / N0,
        # Guard has to be imlemented in order to avoid division by zero
        beta=lambda x, y: 1.0 / max([y - 1.0, 1e-10]),  #
        gamma=lambda x: kb * (x - 1.0),
    )

print(sol_ramk[g].initial_total_volume, sol_ramk[g].total_volume)

fig = plt.figure()
ax = fig.gca()
markers = cycle(["v", "s", "o", "*", ".", ","])
cor1 = cycle(cor)
cor2 = cycle(cor)
v = linspace(0.9, vmax, 10000)

ax.set_ylim([1e-10, 10e3])
for n in sorted(sol_ramk):
    ax.loglog(
        sol_ramk[n].xi,
        sol_ramk[n].N[-1],
        marker=next(markers),
        color=next(cor1),
        ms=9,
        mfc="none",
        mec=next(cor2),
        mew=1.5,
        label="MOC Ramk N={0}".format(n),
    )

for n in sorted(sol_hidy):
    ax.loglog(
        sol_hidy[n].xi,
        sol_hidy[n].N[-1],
        marker=next(markers),
        color=next(cor1),
        ms=8,
        mfc="none",
        mec=next(cor2),
        mew=1.5,
        label="MOC Hidy N={0}".format(n),
    )
ax.loglog(
    v,
    N0 * blatz_and_tobolsky_pbe_solution(v, time[-1], kc, kb),
    "-k",
    linewidth=2,
    label="Analytical $t=\infty$",
)
ax.legend(loc="upper right", shadow=True)
ax.set_xlabel("Volume do polímero")
ax.set_ylabel("N/N0")
plt.show()

fig.savefig(
        path.join(dir, "blatz_x.pdf"),
        backend="pgf",
        bbox_inches="tight",
        pad_inches=0.05,
    )