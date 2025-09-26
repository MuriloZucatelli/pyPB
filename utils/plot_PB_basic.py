from os.path import join, abspath, dirname
from os import getcwd
from sys import path as sph
import pickle
from numpy import abs, pi, set_printoptions
from pbe.setup.helpers import plt_config2
from matplotlib import pyplot as plt
from matplotlib import ticker, rcParams
from matplotlib.ticker import ScalarFormatter
from seaborn import scatterplot, color_palette, move_legend
from pbe.app.dtg_class import Import_flow_DSD2

plt_config2(relative_fig_width=0.7)
set_printoptions(precision=4)

cor = color_palette("Paired", n_colors=5)
rcParams["legend.fontsize"] = 5
y_name_error = r"Erro relativo ANM $\psi$ \% "


def D43Reduction_x_Epsilon(data, name, dir, save=False, hue="MV-01 [%]", dp_name="MV01"):
    fig = plt.figure()
    ax = fig.gca()
    ax.set_xlabel(r"Dissipação de energia cinética média $\overline{\epsilon}$ [m²/s³]")
    ax.set_ylabel(r"Razão $D_{[4,3]}d / D_{[4,3]}a$ " + f"{dp_name}")
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
        x="epsilon [m²/s³]",
        y="D43r",
        style="H2O [%]",
        hue=hue,
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


# def D95_x_Epsilon(data, name, dir, save=False, ylabel=r"Diâmetro $d_{90}$ depois ANM $[\mu m]$"):
#     fig = plt.figure()
#     ax = fig.gca()
#     ax.set_xlabel(r"Dissipação de energia cinética média $\overline{\epsilon}$ [m²/s³]")
#     ax.set_ylabel(ylabel)
#     ax.text(
#         0.3,
#         0.9,
#         r"$L_{diss}$ = 5D",
#         fontsize=10,
#         horizontalalignment="left",
#         verticalalignment="center",
#         transform=ax.transAxes,
#     )

#     scatterplot(
#         data=data,
#         x="epsilon [m²/s³]",
#         y="dv99d MV-01 [um]",
#         style="H2O [%]",
#         hue="ANM [%]",
#         palette=color_palette("bright"),
#         ax=ax,
#     )
#     # ax.hlines(
#     #    y=1.0, xmin=0, xmax=0.9 * ax.get_xlim()[1], linestyles="dashed", linewidth=0.8
#     # )
#     ax.set_xmargin(1 / 100)
#     # ax.semilogx()
#     ax.legend()
#     move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
#     # ax.grid()
#     y_formatter = ticker.ScalarFormatter(useOffset=False)
#     y_formatter.set_scientific(False)
#     if save:
#         fig.savefig(
#             join(dir, f"{name}.pdf"),
#             backend="pgf",
#             bbox_inches="tight",
#             pad_inches=0.05,
#         )
#     return fig


def dXX_x_Epsilon(
    data,
    name,
    dir,
    save=False,
    y="",
    f=99,
    ylabel="",
    base="v",
    hue=None,
    local="depois MV01",
):
    fig = plt.figure()
    ax = fig.gca()
    d = f"d_{{{base+str(f)}}}"
    ylabel = r"Diâmetro ${}$ {} $[\mu m]$".format(d, local)
    ax.set_xlabel(r"Dissipação de energia cinética média $\overline{\epsilon}$ [m²/s³]")
    ax.set_ylabel(ylabel)
    ax.text(
        0.4,
        0.9,
        r"$L_{diss}$ = 5D",
        fontsize=10,
        horizontalalignment="left",
        verticalalignment="center",
        transform=ax.transAxes,
    )

    scatterplot(
        data=data,
        x="epsilon [m²/s³]",
        y=y,
        style="H2O [%]",
        hue=hue,
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


def dXX_x_Epsilon_2x2(
    data,
    name,
    dir,
    save=False,
    y="",
    f=99,
    ylabel="",
    base="v",
    local="depois MV01",
    Ca_crit=None,
):

    fig, axs = plt.subplots(2, 2, figsize=(6, 5), sharey=True, sharex=True)
    axs = axs.flatten()
    h2o_concentrations = data.groupby(by="H2O [%]").groups.keys()
    xlabel = r"Dissipação de energia cinética $\overline{\epsilon}$ [m²/s³]"
    d = f"d_{{{base+str(f)}}}"
    ylabel = r"${}$ {} $[\mu m]$".format(d, local)
    for i, (ax, h2o) in enumerate(zip(axs, h2o_concentrations)):
        hue = "ANM [%]"  # Exemplo de coluna no dataset
        style = "H2O [%]"  # Exemplo de coluna no dataset
        # Filtrando os dados para cada concentração de H2O
        data_h2o = data[data["H2O [%]"] == h2o]
        # Adicionando scatterplot
        scatterplot(
            data=data_h2o,
            x="epsilon [m²/s³]",
            y=y,
            palette="gray",
            ax=ax,
            legend=False,  # Desativar legendas individuais
        )

        # Títulos e ajustes
        ax.text(
            0.4, 0.9, f"$H2O = {h2o}$%", fontsize=10,
            horizontalalignment="left", verticalalignment="center",
            transform=ax.transAxes
        )
        ax.xaxis.set_major_formatter(ScalarFormatter(useOffset=False))
        ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
        ax.grid(True, which="both", linestyle="--", linewidth=0.5)
        # ax.xaxis.set_major_formatter(ScalarFormatter(useOffset=False))
        # ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))

        # Apenas o primeiro subplot terá o rótulo do eixo y
        if i % 2 == 0:
            ax.set_ylabel(ylabel)

        # Apenas a última linha terá o rótulo do eixo x
        if i >= 2:
            ax.set_xlabel(xlabel)

    handles, labels = axs[0].get_legend_handles_labels()
    # fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 1.05), ncol=4)

    fig.tight_layout()

    if save:
        fig.savefig(
            join(dir, f"{name}.pdf"),
            backend="pgf",
            bbox_inches="tight",
            pad_inches=0.05,
        )
    return fig


def dXX_x_Capilaridade_2x2(
    data,
    name,
    dir,
    save=False,
    y="",
    f=99,
    ylabel="",
    base="v",
    local="depois MV01",
    Ca_crit=None,
):

    fig, axs = plt.subplots(2, 2, figsize=(6, 5), sharey=True, sharex=True)
    axs = axs.flatten()
    h2o_concentrations = data.groupby(by="H2O [%]").groups.keys()
    xlabel = r"Número de capilaridade [-]"
    d = f"d_{{{base+str(f)}}}"
    ylabel = r"${}$ {} $[\mu m]$".format(d, local)
    for i, (ax, h2o) in enumerate(zip(axs, h2o_concentrations)):
        hue = "ANM [%]"  # Exemplo de coluna no dataset
        style = "H2O [%]"  # Exemplo de coluna no dataset
        # Filtrando os dados para cada concentração de H2O
        data_h2o = data[data["H2O [%]"] == h2o]
        # Adicionando scatterplot
        scatterplot(
            data=data_h2o,
            x="Ca_1",
            y=y,
            palette="gray",
            ax=ax,
            legend=False,  # Desativar legendas individuais
        )

        # Títulos e ajustes
        ax.text(
            0.4, 0.9, f"$H2O = {h2o}$%", fontsize=10,
            horizontalalignment="left", verticalalignment="center",
            transform=ax.transAxes
        )
        ax.xaxis.set_major_formatter(ScalarFormatter(useOffset=False))
        ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
        ax.grid(True, which="both", linestyle="--", linewidth=0.5)
        # ax.xaxis.set_major_formatter(ScalarFormatter(useOffset=False))
        # ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))

        # Apenas o primeiro subplot terá o rótulo do eixo y
        if i % 2 == 0:
            ax.set_ylabel(ylabel)

        # Apenas a última linha terá o rótulo do eixo x
        if i >= 2:
            ax.set_xlabel(xlabel)

    handles, labels = axs[0].get_legend_handles_labels()
    # fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 1.05), ncol=4)

    fig.tight_layout()

    if save:
        fig.savefig(
            join(dir, f"{name}.pdf"),
            backend="pgf",
            bbox_inches="tight",
            pad_inches=0.05,
        )
    return fig


def d90_d_a_Epsilon(data, name, dir, save=False, hue=None, dp_name=None):
    fig = plt.figure()
    ax = fig.gca()
    ax.set_xlabel(r"Dissipação de energia cinética média $\overline{\epsilon}$ [m²/s³]")
    ax.set_ylabel(r"Diâmetro $d_{90}d$ - $d_{90}a$ " + f"{dp_name}" + r" $[\mu m]$")
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
        x="epsilon [m²/s³]",
        y=f"d90d_a {dp_name} [um]",
        style="H2O [%]",
        hue=hue,
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
    style="Nome",
    hue="MV-01 [%]",
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
    # ylabel = f"Erro relativo {dp_name}" + r"$\psi$ \% "
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
    data, name, dir, x="Red D43", y="Erro", style="Nome", hue="ANM [%]", ylabel=None, save=False, dp_name="MV01",
):
    fig = plt.figure()
    ax = fig.gca()
    ax.set_xlabel(r"Razão $D_{[4,3]}d / D_{[4,3]}a$" + f" {dp_name}")
    ax.set_ylabel(ylabel)
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


def DTG_best_worst(data, runs, name, dir, save=False, nrows=2, ws=0.2, dp_name=None):

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
        nome = df["Nome"]
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
        xx = f"{dp_name} [%]"
        ax.text(
            0.702,
            0.855,
            "{0:s} \n N: {1:d} \n H2O: {2:.0f}\% \n".format(nome, N_emul, df["H2O [%]"])
            + f"{dp_name}: {df[xx]}\% \n"
            + "Re: {0:.1f} \n".format(df["Re"])
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
            join(dir, f"{name}_{dp_name}.pdf"),
            backend="pgf",
            bbox_inches="tight",
            pad_inches=0.05,
        )

    # plt.plot(sol.moc.xi_d*1e6, sol.N2Fv[0], label = 'Extrator 2')
    # ax.yaxis.set_major_formatter(y_formatter)
    plt.close()


def DTGs_plot(data, runs, name, dir, dp_name=None, save=False):

    # plt_config2(relative_fig_width=0.7)
    plt_config2(relative_fig_width=0.7, page_width=180, invert=True, golden_mean=0.8)

    # plt.plot(sol.moc.xi_d*1e6, sol.N2Fv[0], label = 'Extrator 2')
    # ax.yaxis.set_major_formatter(y_formatter)

    for i in range(len(data)):

        fig = plt.figure()
        ax = fig.gca()

        df = data.iloc[i]
        run = df["Parâmetro"]
        nome = df["Nome"]
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
            label=f"Mont. {dp_name}",
            dashes=[2, 2],
        )
        ax.plot(
            sol["pbe_sol"].moc.xi_d * 1e6,
            sol["pbe_sol"].N2Fv,
            label=f"Sim. jus. {dp_name}",
        )
        ax.plot(
            sol["pbe_sol"].moc.xi_d * 1e6,
            exp.get_DTG(ID=IDs[1], marco=m, teste=N_esc).freq_v / 100,
            label=f"Exp. jus. {dp_name}",
            dashes=[5, 1],
        )
        ax.text(
            0.76, #0.77,
            0.6,
            "{3:s} \n N: {4:d} \n We: {5:g}\n H2O: {0:.0f}\% \n {6:s}: {1:.0f}\% \n Re: {2:.1f} \n".format(
                df["H2O [%]"],
                df[f"{dp_name} [%]"],
                df["Re"],
                nome,
                N_emul,
                We,
                dp_name,
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

        ax.legend(loc="upper right", fontsize=9)
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
