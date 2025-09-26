from os.path import join, abspath, dirname
from os import getcwd
from sys import path as sph
import pickle
from numpy import abs, pi, set_printoptions, arange
from pbe.setup.helpers import plt_config2
from matplotlib import pyplot as plt
from matplotlib import ticker, rcParams
from seaborn import scatterplot, color_palette, move_legend
from pbe.app.dtg_class import Import_flow_DSD2

plt_config2(relative_fig_width=0.7)
set_printoptions(precision=4)

cor = color_palette("Paired", n_colors=5)
rcParams["legend.fontsize"] = 5
y_name_error = r"Erro relativo ANM $\psi$ \% "


def d_max_kolmogorov_length(
    data,
    name,
    dir,
    save=False,
    x="lambdak [um]",
    y="dv99d MV01 [um]",
    style=None,
    hue=None,
    size=None,
):
    plt_config2(relative_fig_width=0.5)
    Altura = rcParams["figure.figsize"][0]
    fig = plt.figure(figsize=[Altura, Altura])
    ax = fig.gca()
    ax.set_xlabel(r"Comprimento de Kolmogorov $\lambda_k$ [$\mu$m]")
    ax.set_ylabel(r"Tamanho da gota $D_{v99}$ [$\mu$m]")
    ax.text(
        0.5,
        0.8,
        r"Subescala inercial" + "\n" + r"d > $\lambda_k$",
        fontsize=10,
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
    )
    ax.text(
        0.8,
        0.55,
        r"Subescala viscosa" + "\n" + r"d < $\lambda_k$",
        fontsize=10,
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
    )

    ax.plot([0, 1], [0, 1], transform=ax.transAxes, color="black", linestyle="dashed")

    scatterplot(
        data=data,
        x=x,
        y=y,
        style=style,
        hue=hue,
        size=size,
        palette=color_palette("bright"),
        ax=ax,
    )
    # ax.set_xlim(10)
    # ax.set_ylim(10)
    low_x, high_x = ax.get_xlim()
    low_y, high_y = ax.get_ylim()
    low = min(low_x, low_y)
    high = max(high_x, high_y)
    ax.set(xlim=(low, high), ylim=(low, high))

    # ax.set_ylim(0, ax.get_xlim()[1])
    ax.set_xmargin(1 / 100)
    # ax.semilogx()
    # ax.legend()
    # move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    # ax.grid()
    y_formatter = ticker.ScalarFormatter(useOffset=False)
    y_formatter.set_scientific(False)
    if save:
        fig.savefig(
            join(dir, f"{name}.pdf"),
            backend="pgf",
            bbox_inches="tight",
            pad_inches=0.05,
        )
    return fig


def d_We_Re(
    data,
    name,
    dir,
    save=False,
    x="",
    y="",
    D=20.95 / 1000,
    style=None,
    hue=None,
    size=None,
    dp_name="MV01",
):
    plt_config2(relative_fig_width=0.6)

    fig = plt.figure()
    ax = fig.gca()
    ax.set_xlabel(r"$We^{4/5} / Re$ [-]")
    ax.set_ylabel(r"$ (\overline{d} / D) We^{3/5}$ [-]")
    ax.text(
        0.5,
        0.8,
        r"Subescala inercial" + "\n" + r"d > $\lambda_k$",
        fontsize=10,
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
    )
    ax.text(
        0.8,
        0.55,
        r"Subescala viscosa" + "\n" + r"d < $\lambda_k$",
        fontsize=10,
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
    )

    data["dDWe"] = data["We_pipe"] ** (4 / 5) / data["Re_pipe"]
    data["We_Re"] = data[f"dv50d {dp_name} [um]"] / (D * 1e6) * data["We_pipe"] ** (3 / 5)

    Re = arange(50, 1000, 1)
    We = (0.0674 * Re) ** (5 / 4)

    ax.plot(Re, We, color="black", linestyle="dashed")

    scatterplot(
        data=data,
        x="Re_pipe",
        y="We_pipe",
        style=style,
        hue=hue,
        size=size,
        palette=color_palette("bright"),
        ax=ax,
    )
    # ax.set_xlim(10)
    # ax.set_ylim(10)
    # low_x, high_x = ax.get_xlim()
    # low_y, high_y = ax.get_ylim()
    # low = min(low_x, low_y)
    # high = max(high_x, high_y)
    # ax.set(xlim=(low, high), ylim=(low, high))

    # ax.set_ylim(0, ax.get_xlim()[1])
    ax.set_xmargin(1 / 100)
    ax.semilogx()
    ax.semilogy()
    # ax.legend()
    # move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    # ax.grid()
    y_formatter = ticker.ScalarFormatter(useOffset=False)
    y_formatter.set_scientific(False)
    if save:
        fig.savefig(
            join(dir, f"{name}.pdf"),
            backend="pgf",
            bbox_inches="tight",
            pad_inches=0.05,
        )
    return fig
