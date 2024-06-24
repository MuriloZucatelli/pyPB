from numpy import (
    arange,
    linspace,
    array,
    zeros,
    isclose,
    diff,
    insert,
    where,
    piecewise,
    geomspace,
)
from itertools import cycle
from sys import path as sph
import os.path as path
import matplotlib.pyplot as plt
from matplotlib import ticker

dir = path.dirname(__file__)
if __name__ == "__main__":
    sph.insert(1, path.abspath(path.join(dir, "..\\..")))
from pbe.solvers.moc_ramk import MOCSolution as MOCramk
from pbe.solvers.moc import MOCSolution as MOChidy
from pbe.setup.helpers import plt_config2
from testes.test_moc import ziff_total_number_solution, ziff_pbe_solution

"""
Comparative between Hidy and Ramkrishna

    Ziff, R.M and McGrady, E.D
    "New solutions to the fragmentation equation", J. Phys A, vol. 18, 1985

    We are looking at case 4. from their paper with kernel F(x, y) = x + y.
    This corresponds to choosing:
        beta = 1.0/y,
        Gamma = y^2,
    for our kernels.

"""
save = False
grids = [10, 50]  # Número de classes utilizadas na discretização
time = arange(0.0, 100.0, 0.01)  # Tempo e passo
v0 = 7.87011e-05 / 1e9
dv0 = 2.65025e-05
r = 1.445124894

# v0 = 0.05
# v0 = 10.0
# r = 1.2
malha = 3
# vmax = 0.1
vmax = 1.0
v0 = 1.0
# Caso 4 funcionou bem
# v0 = 0.0001
# dv0 = v0

sol_ramk = dict()
sol_hidy = dict()

for g in grids:
    threshold = vmax / g / 2

    def n0_init2(x):
        return piecewise(
            x, [abs(v0 - x) < threshold, v0 - x >= threshold], [g / vmax, 0]
        )  # 0 ou g/vmax   função delta de dirac generica

    malha = 3
    if malha == 2:
        xi = linspace(1e-5, vmax, g, endpoint=True)

    elif malha == 3:
        xi = v0 + v0 * arange(g)
        xi = geomspace(1e-3, vmax, g, endpoint=True)

    elif malha == 5:
        xi = (vmax / g) + (vmax / g) * arange(g)

    sol_ramk[g] = MOCramk(
        g,  # number of classes
        time,
        xi=xi,
        n0=n0_init2,
        # N0=N0,
        beta=lambda x, y, i: 1.0 / y,  # DDSD
        gamma=lambda x: x**2,  # breakage rate
    )

    xi = linspace(1e-5, vmax, g, endpoint=True)
    sol_hidy[g] = MOChidy(
        g,  # number of classes
        time,
        xi=xi,
        n0=n0_init2,
        # N0=N0,
        beta=lambda x, y: 1.0 / y,  # DDSD
        gamma=lambda x: x**2,  # breakage rate
    )


# set_plt_params()
"""
    plotss
"""


vmax = xi[-1]  # Max volume
v0 = xi[-1]  # Initial volume
totals_ramk = dict((n, sol_ramk[n].total_numbers) for n in grids)
totals_hidy = dict((n, sol_hidy[n].total_numbers) for n in grids)
volume = dict((n, sol_ramk[n].phase_fraction) for n in grids)


def total_numbers():
    plt_config2(relative_fig_width=0.7)
    v = linspace(0, v0, 100)
    Na = array([ziff_total_number_solution(v, t, v0) for t in time])

    fig = plt.figure()
    ax = fig.gca()
    linestyles = cycle(["-", "--", ":"])
    markers = cycle(["o", "s", "v", "*", ".", ","])
    # Iterando sobre cada resultado para o numero de classes n
    for n in sorted(totals_ramk):
        ax.loglog(
            time,
            totals_ramk[n]
            / totals_ramk[n][
                0
            ],  # Total numbers of droplets / Initial total n drop. t=0
            linestyle=next(linestyles),
            marker=next(markers),
            markevery=0.05,
            label="MOC Ramkrishna M={0}".format(n),
        )

    for n in sorted(totals_hidy):
        ax.loglog(
            time,
            totals_hidy[n]
            / totals_hidy[n][
                0
            ],  # Total numbers of droplets / Initial total n drop. t=0
            linestyle=next(linestyles),
            marker=next(markers),
            markevery=0.05,
            label="MOC Hidy  M={0}".format(n),
        )

    ax.loglog(time, Na, "--k", linewidth=1.5, label="Analítico")
    # ax.set_ylim(top=100)
    ax.legend(loc="best", shadow=True)
    ax.set_xlabel("Tempo [s]")
    ax.set_ylabel("N/N0")
    ax.grid()
    y_formatter = ticker.ScalarFormatter(useOffset=False)
    y_formatter.set_scientific(False)
    # ax.yaxis.set_major_formatter(y_formatter)
    if save:
        fig.savefig(
            path.join(dir, "ziff_N.pdf"),
            backend="pgf",
            bbox_inches="tight",
            pad_inches=0.05,
        )
    return fig


"""
    Função densidade numerica
"""


def densidade_n():
    plt_config2(relative_fig_width=0.7)
    fig = plt.figure()
    ax = fig.gca()
    markers = cycle(["o", "s", "v", "*", ".", ","])
    v = linspace(0, v0, 10000)

    for n in sorted(grids):
        ax.semilogx(
            sol_ramk[n].xi,  # volume pivot
            sol_ramk[n].number_density[-1],  # last time
            marker=next(markers),
            label="Ramk. M={0}".format(n),
        )

    for n in sorted(grids):
        ax.semilogx(
            sol_hidy[n].xi,  # volume pivot
            sol_hidy[n].number_density[-1],  # last time
            marker=next(markers),
            label="Hidy M={0}".format(n),
        )

    ax.semilogx(
        v,
        ziff_pbe_solution(v, time[-1], v0),  # analytical in last time
        "--k",
        linewidth=1.0,
        label="Analítico {:.0f} s".format(time[-1]),
    )

    ax.text(
        0.18,
        0.85,
        r"N/N0 $\approx$ " + "{:.0f}".format(totals_ramk[n][-1] / totals_ramk[n][0]),
        fontsize=11,
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
    )
    ax.legend(loc="best", shadow=True)
    ax.set_xlabel("Volume [mm³]")
    # ax.set_ylabel("Number density function")
    ax.set_ylabel("Densidade numérica [-]")
    ax.grid()
    if save:
        fig.savefig(
            path.join(dir, "ziff_n_x.pdf"),
            backend="pgf",
            bbox_inches="tight",
            pad_inches=0.05,
        )
    plt.show()


"""
    Densidade numerica para tempos diferentes
"""


def densidade_n_t():
    plt_config2(relative_fig_width=0.7)
    time2 = 0.3 * time[-1]
    tol = 0.2 * diff(time).mean()
    loct2 = where(isclose(time, time2, atol=tol))[-1].squeeze()
    fig = plt.figure()
    ax = fig.gca()
    markers = cycle(["o", "s", "v", "*", ".", ","])
    v = linspace(0, v0, 10000)

    n = grids[-1]
    ax.semilogx(
        sol_ramk[n].xi,  # volume pivot
        sol_ramk[n].number_density[-1],  # last time
        marker=next(markers),
        label="Ramk. t={1:.0f} s".format(n, time[-1]),
    )

    n = grids[-1]
    ax.semilogx(
        sol_ramk[n].xi,  # volume pivot
        sol_ramk[n].number_density[loct2],  # last time
        marker=next(markers),
        label="Ramk. t={1:.0f} s".format(n, time[loct2]),
    )

    ax.semilogx(
        sol_hidy[n].xi,  # volume pivot
        sol_hidy[n].number_density[-1],  # last time
        marker=next(markers),
        label="Hidy t={1:.0f} s".format(n, time[-1]),
    )

    n = grids[-1]
    ax.semilogx(
        sol_hidy[n].xi,  # volume pivot
        sol_hidy[n].number_density[loct2],
        marker=next(markers),
        label="Hidy t={1:.0f} s".format(n, time[loct2]),
    )

    ax.semilogx(
        v,
        ziff_pbe_solution(v, time[-1], v0),  # analytical in last time
        "--k",
        linewidth=1.0,
        label="Analítico {:.0f} s".format(time[-1]),
    )

    ax.semilogx(
        v,
        ziff_pbe_solution(v, time[loct2], v0),  # analytical in last time
        ":k",
        linewidth=1.0,
        label="Analítico {:.0f} s".format(time[loct2]),
    )

    ax.text(
        0.20,
        0.87,
        r"N/N0 Ramk. $\approx$ "
        + "{:.0f}".format(totals_ramk[n][-1] / totals_ramk[n][0]),
        fontsize=11,
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
    )

    ax.text(
        0.16,
        0.35,
        r"N/N0 Ramk. $\approx$ "
        + "{:.0f}".format(totals_ramk[n][loct2] / totals_ramk[n][0]),
        fontsize=11,
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
    )

    ax.text(
        0.85,
        0.45,
        r"M = " + "{:.0f}".format(grids[-1]),
        fontsize=11,
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
    )

    ax.legend(loc="best", shadow=False, framealpha=0.65)
    ax.set_xlabel("Volume [mm³]")
    # ax.set_ylabel("Number density function")
    ax.set_ylabel("Densidade numérica [-]")
    ax.grid()
    if save:
        fig.savefig(
            path.join(dir, "ziff_n_x_t.pdf"),
            backend="pgf",
            bbox_inches="tight",
            pad_inches=0.05,
        )
    plt.show()


# call plots
fig = total_numbers()
fig = densidade_n()
# fig = densidade_n_t()
plt.show()
