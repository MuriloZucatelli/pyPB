from os.path import join, abspath, dirname
from os import getcwd
from sys import path as sph
import pickle
from numpy import abs, pi, set_printoptions
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


def D43Reduction_x_Epsilon(data, name, dir, save=False):
    fig = plt.figure()
    ax = fig.gca()
    ax.set_xlabel(r"Dissipação de energia cinética média $\overline{\epsilon}$ [m²/s³]")
    ax.set_ylabel(r"Razão $D_{[4,3]}d / D_{[4,3]}a$ ANM")
    ax.text(
        0.3,
        0.9,
        r"$L_{diss}$ = 5D",
        fontsize=10,
        horizontalalignment="left",
        verticalalignment="center",
        transform=ax.transAxes,
    )

    scatterplot(
        data=data,
        x="epsilon",
        y="Red D43",
        style="H2O [%]",
        hue="ANM [%]",
        palette=color_palette("bright"),
        ax=ax,
    )
    ax.hlines(
        y=1.0, xmin=0, xmax=0.9 * ax.get_xlim()[1], linestyles="dashed", linewidth=0.8
    )
    ax.set_xmargin(1 / 100)
    # ax.semilogx()
    ax.legend()
    move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
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


def D95_x_Epsilon(data, name, dir, save=False):
    fig = plt.figure()
    ax = fig.gca()
    ax.set_xlabel(r"Dissipação de energia cinética média $\overline{\epsilon}$ [m²/s³]")
    ax.set_ylabel(r"Diâmetro $D_{90}$ depois ANM $[\mu m]$")
    ax.text(
        0.3,
        0.9,
        r"$L_{diss}$ = 5D",
        fontsize=10,
        horizontalalignment="left",
        verticalalignment="center",
        transform=ax.transAxes,
    )

    scatterplot(
        data=data,
        x="epsilon",
        y="D95",
        style="H2O [%]",
        hue="ANM [%]",
        palette=color_palette("bright"),
        ax=ax,
    )
    # ax.hlines(
    #    y=1.0, xmin=0, xmax=0.9 * ax.get_xlim()[1], linestyles="dashed", linewidth=0.8
    # )
    ax.set_xmargin(1 / 100)
    # ax.semilogx()
    ax.legend()
    move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
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


def D90_d_a_Epsilon(data, name, dir, save=False):
    fig = plt.figure()
    ax = fig.gca()
    ax.set_xlabel(r"Dissipação de energia cinética média $\overline{\epsilon}$ [m²/s³]")
    ax.set_ylabel(r"Diâmetro $D_{90}d$ - $D_{90}a$ ANM $[\mu m]$")
    ax.text(
        0.3,
        0.9,
        r"$L_{diss}$ = 5D",
        fontsize=10,
        horizontalalignment="left",
        verticalalignment="center",
        transform=ax.transAxes,
    )

    scatterplot(
        data=data,
        x="epsilon",
        y="D90d_a",
        style="H2O [%]",
        hue="ANM [%]",
        palette=color_palette("bright"),
        ax=ax,
    )
    # ax.hlines(
    #    y=1.0, xmin=0, xmax=0.9 * ax.get_xlim()[1], linestyles="dashed", linewidth=0.8
    # )
    ax.set_xmargin(1 / 100)
    # ax.semilogx()
    ax.legend()
    move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
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


def error_x_(
    data,
    name,
    dir,
    save=False,
    x="N_emul",
    y="Erro",
    style="Parâmetro",
    hue="ANM [%]",
    size=None,
    xlabel="Teste",
    ylabel=y_name_error,
    annot=r"$L_{diss}$ = 5D",
    fw=0.7,
    xtext=1.05,
):
    """Grafico do Erro (Função objetivo vs x)
    O Erro é calculado considerando
    """
    plt_config2(relative_fig_width=fw)
    fig = plt.figure()
    ax = fig.gca()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.text(
        xtext,
        0.2,  # 0.05 for error_x_test_D43red
        annot,
        fontsize=10,
        horizontalalignment="left",
        verticalalignment="center",
        transform=ax.transAxes,
    )

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
    # ax.semilogx()
    if style is not None or hue is not None or size is not None:
        ax.legend()
        move_legend(ax, "upper left", bbox_to_anchor=(1, 1), fontsize=8)
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


def error_x_D43red(
    data, name, dir, x="Red D43", y="Erro", style="Parâmetro", hue="ANM [%]", save=False
):
    fig = plt.figure()
    ax = fig.gca()
    ax.set_xlabel(r"Razão $D_{[4,3]}d / D_{[4,3]}a$ ANM")
    ax.set_ylabel(y_name_error)
    ax.text(
        1.05,
        0.2,
        r"$L_{diss}$ = 5D",
        fontsize=10,
        horizontalalignment="left",
        verticalalignment="center",
        transform=ax.transAxes,
    )

    scatterplot(
        data=data,
        x=x,
        y=y,
        style=style,
        hue=hue,
        palette=color_palette("bright"),
        ax=ax,
    )
    # ax.semilogx()
    if style is not None or hue is not None:
        ax.legend()
        move_legend(ax, "upper left", bbox_to_anchor=(1, 1), fontsize=8)
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


def DTG_best_worst(data, runs, name, dir, save=False, nrows=2, ws=0.2):

    plt_config2(relative_fig_width=1, page_width=180, golden_mean=0.5)

    hr = [10, 10]
    if nrows == 1:
        hr = [10]
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=2,
        sharey=True,
        gridspec_kw={
            "width_ratios": [10, 10],
            "height_ratios": hr,
            "hspace": 0.2,
            "wspace": ws,
        },
        subplot_kw=(dict(box_aspect=1) if False else None),
    )
    for ax in fig.get_axes():
        pass

    i = 0
    enum = ["a", "b", "c", "d"]
    for ax in fig.get_axes():

        df = data.iloc[i]
        run = df["Parâmetro"]
        N_emul = df["N_emul"]
        m = df["Marco"]

        exp = runs[run]["experiments"]  # Pega o experimento
        IDs = exp.compares[data["Compare"].iloc[0]]

        sol = runs[run][df["idict"]]
        N_esc = sol["N_escoam"]
        ax.plot(
            sol["pbe_sol"].moc.xi_d * 1e6,
            exp.get_DTG(ID=IDs[0], marco=m, teste=N_esc).freq_v,
            label="Antes da válvula",
            dashes=[2, 2],
        )
        ax.plot(
            sol["pbe_sol"].moc.xi_d * 1e6,
            100 * sol["pbe_sol"].N2Fv,
            label=f"Sim. depois da válvula",
        )
        ax.plot(
            sol["pbe_sol"].moc.xi_d * 1e6,
            exp.get_DTG(ID=IDs[1], marco=m, teste=N_esc).freq_v,
            label=f"Exp. depois da válvula",
            dashes=[5, 1],
        )
        ax.text(
            0.702,
            0.855,
            "{3:s} \n N: {4:d} \n H2O: {0:.0f}\% \n ANM: {1:.0f}\% \n Re: {2:.1f} \n".format(
                df["H2O [%]"], df["ANM [%]"], df["Re"], run, N_emul
            )
            + r"$\psi_e$: {0:.1f}\%".format(df["Erro"]),
            fontsize=8,
            horizontalalignment="left",
            verticalalignment="center",
            transform=ax.transAxes,
        )
        ax.set_ylabel("")
        ax.legend().set_visible(False)
        ax.semilogx()
        ax.set_title(f"({enum[i]})", fontsize=11)
        i += 1

    # ax.set_ylabel("Diferença Sim. e Exp. [\%]")
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        ncol=3,
        fontsize=12,
        bbox_to_anchor=(0.5, 1.03),
    )
    fig.supylabel(r"Frequência volumétrica $f_v$ [%]", x=0.05)
    fig.supxlabel(r"Diâmetro da gota [$\mu$m]")

    # y_formatter = ticker.ScalarFormatter(useOffset=False)
    # y_formatter.set_scientific(False)
    if save:
        fig.savefig(
            join(dir, f"{name}_best_worst.pdf"),
            backend="pgf",
            bbox_inches="tight",
            pad_inches=0.05,
        )

    # plt.plot(sol.moc.xi_d*1e6, sol.N2Fv[0], label = 'Extrator 2')
    # ax.yaxis.set_major_formatter(y_formatter)
    plt.close()


def DTGs_plot(data, runs, name, dir, save=False):

    plt_config2(relative_fig_width=0.7)

    # plt.plot(sol.moc.xi_d*1e6, sol.N2Fv[0], label = 'Extrator 2')
    # ax.yaxis.set_major_formatter(y_formatter)

    for i in range(len(data)):

        fig = plt.figure()
        ax = fig.gca()

        df = data.iloc[i]
        run = df["Parâmetro"]
        N_emul = df["N_emul"]
        D43_r = df["Red D43"]
        m = df["Marco"]
        We = df["We"]
        exp = runs[run]["experiments"]  # Pega o experimento
        IDs = exp.compares[data["Compare"].iloc[0]]

        sol = runs[run][df["idict"]]
        N_esc = sol["N_escoam"]
        ax.plot(
            sol["pbe_sol"].moc.xi_d * 1e6,
            exp.get_DTG(ID=IDs[0], marco=m, teste=N_esc).freq_v / 100,
            label="Antes da válvula",
            dashes=[2, 2],
        )
        ax.plot(
            sol["pbe_sol"].moc.xi_d * 1e6,
            sol["pbe_sol"].N2Fv,
            label=f"Sim. depois da válvula",
        )
        ax.plot(
            sol["pbe_sol"].moc.xi_d * 1e6,
            exp.get_DTG(ID=IDs[1], marco=m, teste=N_esc).freq_v / 100,
            label=f"Exp. depois da válvula",
            dashes=[5, 1],
        )
        ax.text(
            0.77,
            0.58,
            "{3:s} \n N: {4:d} \n We: {5:g}\n H2O: {0:.0f}\% \n ANM: {1:.0f}\% \n Re: {2:.1f} \n".format(
                df["H2O [%]"],
                df["ANM [%]"],
                df["Re"],
                run,
                N_emul,
                We,
            )
            + r"$R_{D43}$: "
            + "{0:.2f}".format(D43_r)
            + "\n"
            + r"$\psi_e$: {0:.1f}\%".format(df["Erro"]),
            fontsize=8,
            horizontalalignment="left",
            verticalalignment="center",
            transform=ax.transAxes,
        )
        ax.set_xlabel(r"Diâmetro da gota [$\mu$m]")
        ax.set_ylabel(r"Frequência volumétrica $f_v$ [%]")

        ax.legend()
        ax.semilogx()
        i += 1
        y_formatter = ticker.ScalarFormatter(useOffset=False)
        y_formatter.set_scientific(False)
        # ax.yaxis.set_major_formatter(y_formatter)
        if save:
            fig.savefig(
                join(dir, f"{run}_{N_esc}_{m}.pdf"),
                backend="pgf",
                bbox_inches="tight",
                pad_inches=0.05,
            )
        plt.close()


def doe_mp_plot():
    pass