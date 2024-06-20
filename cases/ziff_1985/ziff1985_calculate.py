from numpy import arange, linspace, array, piecewise, where, isclose, zeros, diff
from itertools import cycle
import sys
import os.path as path
from matplotlib import ticker
import matplotlib.pyplot as plt

dir = path.dirname(__file__)
if __name__ == "__main__":
    sys.path.append(path.abspath(path.join(dir, "..\\..")))
from pbe.solvers.moc import MOCSolution
from pbe.setup.helpers import plt_config2
from tests.test_moc import ziff_total_number_solution, ziff_pbe_solution


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
grids = [150]
time = arange(0.0, 100.0, 0.01)  # Tempo e passo
# N0 = 1
vmax = 1.0  # Volume máximo
v0 = 1.0
pbe_solutions = dict()

for g in grids:
    threshold = vmax / g / 2
    N0 = zeros(g)
    N0[-1] = 1
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
        #n0=n0_init2,
        N0=N0,
        beta=lambda x, y: 1.0 / y,  # DDSD   era 2.0 / y
        gamma=lambda x: x**2,  # breakage rate
    )


# set_plt_params()
"""
    plotss
"""

totals = dict((n, pbe_solutions[n].total_numbers) for n in grids)
volume = dict((n, pbe_solutions[n].total_volume) for n in grids)


def total_numbers():
    plt_config2(relative_fig_width=0.7)
    v = linspace(0, v0, 100)
    Na = array([ziff_total_number_solution(v, t, v0) for t in time])

    fig = plt.figure()
    ax = fig.gca()
    linestyles = cycle(["-", "--", ":"])
    markers = cycle(["o", "s", "v", "*", ".", ","])
    # Iterando sobre cada resultado para o numero de classes n
    for n in sorted(totals):
        ax.loglog(
            time,
            totals[n]
            / totals[n][0],  # Total numbers of droplets / Initial total n drop. t=0
            linestyle=next(linestyles),
            marker=next(markers),
            markevery=0.05,
            label="MOC com M={0}".format(n),
        )
    ax.loglog(time, Na, "--k", linewidth=1.5, label="Analítico")
    ax.set_ylim(top=100)
    ax.legend(loc="best", shadow=True)
    ax.set_xlabel("Tempo [s]")
    ax.set_ylabel("N/N0")
    ax.grid()
    y_formatter = ticker.ScalarFormatter(useOffset=False)
    y_formatter.set_scientific(False)
    # ax.yaxis.set_major_formatter(y_formatter)
    fig.savefig(
        path.join(dir, "total_number_of_droplets_ziff1985.pdf"),
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
            pbe_solutions[n].xi,  # volume pivot
            pbe_solutions[n].number_density[-1],  # last time
            marker=next(markers),
            label="MOC com M={0}".format(n),
        )

    ax.semilogx(
        v,
        ziff_pbe_solution(v, time[-1], v0),  # analytical in last time
        "--k",
        linewidth=1.0,
        label="Analítico t={:.0f} s".format(time[-1]),
    )

    ax.text(
        0.18,
        0.81,
        r"N/N0 $\approx$ " + "{:.0f}".format(totals[n][-1] / totals[n][0]),
        fontsize=11,
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
    )
    ax.legend(loc="lower left", shadow=True)
    ax.set_xlabel("Volume [mm³]")
    # ax.set_ylabel("Number density function")
    ax.set_ylabel("Densidade numérica [-]")
    ax.grid()
    fig.savefig(
        path.join(dir, "number_density_func_ziff1985.pdf"),
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
    tol = 0.5 * diff(time).mean()
    loct2 = where(isclose(time, time2, atol=tol))[-1].squeeze()
    fig = plt.figure()
    ax = fig.gca()
    markers = cycle(["o", "s", "v", "*", ".", ","])
    v = linspace(0, v0, 10000)

    n = grids[-1]
    ax.semilogx(
        pbe_solutions[n].xi,  # volume pivot
        pbe_solutions[n].number_density[-1],  # last time
        marker=next(markers),
        label="MOC com M={0}, t={1:.0f} s".format(n, time[-1]),
    )

    n = grids[-1]
    ax.semilogx(
        pbe_solutions[n].xi,  # volume pivot
        pbe_solutions[n].number_density[loct2],
        marker=next(markers),
        label="MOC com M={0}, t={1:.0f} s".format(n, time[loct2]),
    )

    ax.semilogx(
        v,
        ziff_pbe_solution(v, time[-1], v0),  # analytical in last time
        "--k",
        linewidth=1.0,
        label="Analítico t={:.0f} s".format(time[-1]),
    )

    ax.semilogx(
        v,
        ziff_pbe_solution(v, time[loct2], v0),  # analytical in last time
        ":k",
        linewidth=1.0,
        label="Analítico t={:.0f} s".format(time[loct2]),
    )

    ax.text(
        0.18,
        0.81,
        r"N/N0 $\approx$ " + "{:.0f}".format(totals[n][-1] / totals[n][0]),
        fontsize=11,
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
    )

    ax.text(
        0.16,
        0.21,
        r"N/N0 $\approx$ " + "{:.0f}".format(totals[n][loct2] / totals[n][0]),
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
    fig.savefig(
        path.join(dir, "number_density_t10&t1000_ziff1985.pdf"),
        backend="pgf",
        bbox_inches="tight",
        pad_inches=0.05,
    )
    plt.show()


# call plots
fig = total_numbers()
fig = densidade_n()
fig = densidade_n_t()
plt.show()
