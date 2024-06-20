from numpy import arange, linspace, exp, geomspace, piecewise, where, zeros, diff
from itertools import cycle
import sys
import os.path as path

dir = path.dirname(__file__)
if __name__ == "__main__":
    sys.path.append(path.abspath(path.join(dir, "..\\..")))
from pbe.solvers.moc_ramk import MOCSolution as MOCramk
from pbe.solvers.moc import MOCSolution as MOChidy
from pbe.setup.helpers import plt_config2
from tests.test_moc import scott_total_number_solution3
from tests.test_moc import scott_pbe_solution3
import matplotlib.pyplot as plt

plt_config2(relative_fig_width=0.7)
"""
Case setup based on:

    Scott, W.T.
    "Analytical studies in cloud droplet coalescence", I. J. Atmos. Sci.
    vol. 25, 1968

    We are looking at constant coalescence kernels (case III in the paper).

    NOTE: We use a simplified version of the "Gaussian-like" distribution of
    initial IC.
"""


v0 = 7.87011e-05 / 1e9
r = 1.445124894
r = 1.2

grids = [10, 40]  # number os classes
time = arange(0.0, 100, 0.01)
vmax = 50  # max volume
C = 0.1  # constante de coalescencia
N0 = 2  # initial Number of droplets
v0 = 0.5  # initial volume
malha = 3
sol_ramk = dict()
sol_hidy = dict()

# Distribuição inicial Gaussiana
# Number density function


def n0_init(v):
    # return (N0 / v0) * v/v
    return (N0 / v0) * (v / v0) * exp(-v / v0)


# TODO: dxi esta definido errado no inicio e no fim, o valor de metade do intervalo 
#       não está dando certo aqui mas deu certo no caso de pura quebra, averiguar...
for g in grids:
    if malha == 1:
        xi = v0 * r ** arange(g)
    elif malha == 2:
        xi = linspace(v0, vmax, g, endpoint=True)
        xi = (vmax / g) + (vmax / g) * arange(g)
    elif malha == 3:
        xi = geomspace(1e-4 * v0, vmax, g, endpoint=True)
    sol_ramk[g] = MOCramk(g, time, xi=xi, n0=n0_init, Q=lambda x, y: C)

    xi = linspace(v0, vmax, g, endpoint=True)
    sol_hidy[g] = MOChidy(g, time, dxi=vmax / g, n0=n0_init, Q=lambda x, y: C)
    # pbe_solutions[g] = MOCSolution(g, time, xi=xi, N0=N, Q=lambda x, y: C)

totals_ramk = dict((n, sol_ramk[n].total_numbers) for n in sol_ramk)
totals_hidy = dict((n, sol_hidy[n].total_numbers) for n in sol_ramk)
volume = dict((n, sol_ramk[n].phase_fraction) for n in sol_ramk)

v = linspace(0, vmax, 200)
Na = scott_total_number_solution3(time, C=C, N0=N0)


def plot_n0_init():
    fig = plt.figure()
    ax = fig.gca()
    ax.plot(
        arange(v0, vmax, vmax / 160),
        n0_init(arange(v0, vmax, vmax / 160)),
        label="Ninit",
    )
    ax.grid()
    ax.legend()
    fig.savefig(
        path.join(dir, "plot_n0_init_teste.pdf"),
        backend="pgf",
        bbox_inches="tight",
        pad_inches=0.05,
    )


def total_numbers():
    fig = plt.figure()
    ax = fig.gca()
    linestyles = cycle(["-", "--", ":"])
    for n in sorted(totals_ramk):
        ax.plot(
            time,
            totals_ramk[n] / totals_ramk[n][0],
            linestyle=next(linestyles),
            markevery=0.05,
            label="MOC Ramk. M={0}".format(n),
        )
        ax.plot(
            time,
            totals_hidy[n] / totals_hidy[n][0],
            linestyle=next(linestyles),
            markevery=0.05,
            label="MOC Hidy M={0}".format(n),
        )
    ax.plot(time, Na / N0, "--k", linewidth=1.5, label="Analítico")
    ax.legend(loc="best", shadow=True)
    ax.set_xlabel("Tempo [s]")
    ax.set_ylabel("N/N0")
    ax.text(
        0.79,
        0.63,
        "C = {:.1f}".format(C),
        fontsize=11,
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
    )
    ax.grid()
    plt.show()
    fig.savefig(
        path.join(dir, "scott1968_N.pdf"),
        backend="pgf",
        bbox_inches="tight",
        pad_inches=0.05,
    )
    return fig


"""
    Densidade numérica
"""


def densi_n():
    fig = plt.figure()
    ax = fig.gca()
    markers = cycle(["o", "s", "v", "*", ".", ","])

    for n in sorted(sol_ramk):
        ax.loglog(
            sol_ramk[n].xi,
            sol_ramk[n].number_density[-1],
            ls="",
            marker=next(markers),
            label="MOC Ramk. M={0}".format(n),
        )

    for n in sorted(sol_hidy):
        ax.loglog(
            sol_hidy[n].xi,
            sol_hidy[n].number_density[-1],
            ls="",
            marker=next(markers),
            label="MOC Hidy M={0}".format(n),
        )

    ax.loglog(
        v,
        scott_pbe_solution3(v, time[-1], C=C, xi0=2.0 * v0, N0=N0),
        ":k",
        linewidth=1.5,
        label="Analítico {:.0f} s".format(time[-1]),
    )
    ax.text(
        0.29,
        0.79,
        r"N/N0 $\approx$ " + "{:.3f}".format(totals_ramk[n][-1] / totals_ramk[n][0]),
        fontsize=11,
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
    )
    ax.text(
        0.46,
        0.42,
        "C = {:.1f}".format(C),
        fontsize=11,
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
    )
    ax.legend(loc="best", shadow=True)
    ax.set_xlabel("Volume [mm³]")
    ax.set_ylabel("Densidade numérica [-]")
    ax.grid()
    fig.savefig(
        path.join(dir, "scott1968_density.pdf"),
        backend="pgf",
        bbox_inches="tight",
        pad_inches=0.05,
    )

    return fig


def densi_n_t():
    t2 = where(time == 10)[0].squeeze()
    fig = plt.figure()
    ax = fig.gca()
    markers = cycle(["o", "s", "v", "*", ".", ","])

    n = grids[-1]
    ax.loglog(
        sol_ramk[n].xi,
        sol_ramk[n].number_density[-1],
        ls="",
        marker=next(markers),
        label="MOC Ramk. t={1:.0f} s".format(n, time[-1]),
    )

    ax.loglog(
        sol_ramk[n].xi,
        sol_ramk[n].number_density[-1],
        ls="",
        marker=next(markers),
        label="MOC Ramk. t={1:.0f} s".format(n, time[t2]),
    )

    ax.loglog(
        sol_hidy[n].xi,
        sol_hidy[n].number_density[t2],
        ls="",
        marker=next(markers),
        label="MOC Hidy t={1:.0f} s".format(n, time[-1]),
    )

    ax.loglog(
        sol_hidy[n].xi,
        sol_hidy[n].number_density[t2],
        ls="",
        marker=next(markers),
        label="MOC Hidy t={1:.0f} s".format(n, time[t2]),
    )
    ax.loglog(
        v,
        scott_pbe_solution3(v, time[-1], C=C, xi0=2.0 * v0, N0=N0),
        "--k",
        linewidth=1.5,
        label="Analítico {:.0f} s".format(time[-1]),
    )
    ax.loglog(
        v,
        scott_pbe_solution3(v, time[t2], C=C, xi0=2.0 * v0, N0=N0),
        ":k",
        linewidth=1.5,
        label="Analítico {:.0f} s".format(time[t2]),
    )
    ax.text(
        0.6,
        0.9,
        r"N/N0 $\approx$ " + "{:.3f}".format(totals_ramk[n][-1] / totals_ramk[n][0]),
        fontsize=11,
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
    )
    ax.text(
        0.55,
        0.72,
        r"N/N0 $\approx$ " + "{:.1f}".format(totals_ramk[n][t2] / totals_ramk[n][0]),
        fontsize=11,
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
    )
    ax.text(
        0.485,
        0.45,
        "C = {:.1f}".format(C),
        fontsize=11,
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
    )
    ax.text(
        0.48,
        0.35,
        "M = {:.0f}".format(grids[-1]),
        fontsize=11,
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
    )
    ax.legend(loc="best", shadow=True)
    ax.set_xlabel("Volume [mm³]")
    ax.set_ylabel("Densidade numérica [-]")
    ax.grid()
    fig.savefig(
        path.join(dir, "scott1968_n_t.pdf"),
        backend="pgf",
        bbox_inches="tight",
        pad_inches=0.05,
    )

    return fig


# plot_n0_init()
fig = total_numbers()
fig = densi_n()
#fig = densi_n_t()
plt.show()
