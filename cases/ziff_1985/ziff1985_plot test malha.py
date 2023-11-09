from numpy import genfromtxt, linspace, array, where
import pickle
import os
import sys
import os.path as path

dir = path.dirname(__file__)
if __name__ == "__main__":
    sys.path.append(path.abspath(path.join(dir, "..\\..")))
from pbe.setup.helpers import set_plt_params
from pbe.setup.helpers import plt_config2
from matplotlib.pyplot import figure
from itertools import cycle
from tests.test_moc import ziff_total_number_solution, ziff_pbe_solution
import matplotlib.pyplot as plt
from matplotlib import ticker

# set_plt_params()
"""
    plotss
"""


dir = os.path.dirname(__file__)
with open(os.path.join(dir, "ziff1985_teste_malha.pickle"), "rb") as f:
    pbe_solutions = pickle.load(f)

vmax = pbe_solutions["vmax"]
v0 = pbe_solutions["v0"]
time = pbe_solutions["time"]

totals = dict((n, pbe_solutions[n].total_numbers) for n in pbe_solutions["grid"])
volume = dict((n, pbe_solutions[n].total_volume) for n in pbe_solutions["grid"])


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
    #ax.set_ylim(top=100)
    ax.legend(loc="best", shadow=True)
    ax.set_xlabel("Tempo [s]")
    ax.set_ylabel("N/N0")
    ax.grid()
    y_formatter = ticker.ScalarFormatter(useOffset=False)
    y_formatter.set_scientific(False)
    #ax.yaxis.set_major_formatter(y_formatter)
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

    for n in sorted(pbe_solutions["grid"]):
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

def densidade_n_t(t=100):
    plt_config2(relative_fig_width=0.7)
    t2 = where(time == t)[0].squeeze()
    fig = plt.figure()
    ax = fig.gca()
    markers = cycle(["o", "s", "v", "*", ".", ","])
    v = linspace(0, v0, 10000)

    n = pbe_solutions["grid"][-1]
    ax.semilogx(
        pbe_solutions[n].xi,  # volume pivot
        pbe_solutions[n].number_density[-1],  # last time
        marker=next(markers),
        label="MOC com M={0}, t={1:.0f} s".format(n, time[-1]),
    )

    n = pbe_solutions["grid"][-1]
    ax.semilogx(
        pbe_solutions[n].xi,  # volume pivot
        pbe_solutions[n].number_density[t2],
        marker=next(markers),
        label="MOC com M={0}, t={1:.0f} s".format(n, time[t2]),
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
        ziff_pbe_solution(v, time[t2], v0),  # analytical in last time
        ":k",
        linewidth=1.0,
        label="Analítico t={:.0f} s".format(time[t2]),
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
        r"N/N0 $\approx$ " + "{:.0f}".format(totals[n][t2] / totals[n][0]),
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
fig = densidade_n_t(t=100)
plt.show()
