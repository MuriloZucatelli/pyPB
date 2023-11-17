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
    geomspace
)
from itertools import cycle
import sys
import os.path as path
import matplotlib.pyplot as plt
from matplotlib import ticker

dir = path.dirname(__file__)
if __name__ == "__main__":
    sys.path.append(path.abspath(path.join(dir, "..\\..")))
from pbe.solvers.moc_ramk import MOCSolution
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

    Teste com mesma malha utilizada pelo bettersize
"""

grids = [100]  # Número de classes utilizadas na discretização
time = arange(0.0, 100.0, 0.01)  # Tempo e passo
v0 = 7.87011e-05 / 1e9
dv0 = 2.65025e-05
r = 1.445124894

#v0 = 0.05
#v0 = 10.0
#r = 1.2
malha = 3
#vmax = 0.1
vmax = 1.0

# Caso 4 funcionou bem
# v0 = 0.0001
# dv0 = v0

sol = dict()

for g in grids:
    # This is modelling Dirac's delta
    # initial number density function

    #threshold = vmax / g / 2
    N0 = zeros(g)
    N0[-1] = 1

    def n0_init2(x):
        return piecewise(
            x, [abs(v0 - x) < threshold, v0 - x >= threshold], [g / vmax, 0]
        )  # 0 ou g/vmax   função delta de dirac generica

    if malha == 1:
        xi = v0 * r ** arange(g)
        dxi = dv0 * r ** arange(g)

    elif malha == 2:
        xi = linspace(1e-5, vmax, g, endpoint=True)

    elif malha == 3:
        xi = v0 + v0 * arange(g)
        xi = geomspace(1e-3, vmax, g, endpoint=True)

    elif malha == 4:
        xi = v0 * r ** arange(g)

    elif malha == 5:
        xi = (vmax / g) + (vmax / g) * arange(g)

    sol[g] = MOCSolution(
        g,  # number of classes
        time,
        xi=xi,
        #n0=n0_init2,
        N0=N0,
        beta=lambda x, y: 1.0 / y,  # DDSD
        gamma=lambda x: x**2,  # breakage rate
    )


# set_plt_params()
"""
    plotss
"""


vmax = xi[-1]  # Max volume
v0 = xi[-1]  # Initial volume
totals = dict((n, sol[n].total_numbers) for n in grids)
volume = dict((n, sol[n].total_volume) for n in grids)


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
    # ax.set_ylim(top=100)
    ax.legend(loc="best", shadow=True)
    ax.set_xlabel("Tempo [s]")
    ax.set_ylabel("N/N0")
    ax.grid()
    y_formatter = ticker.ScalarFormatter(useOffset=False)
    y_formatter.set_scientific(False)
    # ax.yaxis.set_major_formatter(y_formatter)
    fig.savefig(
        path.join(dir, "total_number_of_droplets_teste_malha.pdf"),
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
            sol[n].xi,  # volume pivot
            sol[n].number_density[-1],  # last time
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
        0.85,
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
        path.join(dir, "number_density_teste_malha.pdf"),
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
        sol[n].xi,  # volume pivot
        sol[n].number_density[-1],  # last time
        marker=next(markers),
        label="MOC com M={0}, t={1:.0f} s".format(n, time[-1]),
    )

    n = grids[-1]
    ax.semilogx(
        sol[n].xi,  # volume pivot
        sol[n].number_density[loct2],
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
        0.85,
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
        path.join(dir, "number_density_teste_malha.pdf"),
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
