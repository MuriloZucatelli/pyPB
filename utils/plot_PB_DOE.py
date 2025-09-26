from os.path import join, abspath, dirname
from os import getcwd
from sys import path as sph
import pickle
from numpy import abs, pi, set_printoptions
from pbe.setup.helpers import plt_config2
from matplotlib import pyplot as plt
from matplotlib import ticker, rcParams
from itertools import cycle
from seaborn import scatterplot, color_palette, move_legend
from pbe.app.dtg_class import Import_flow_DSD2

plt_config2(relative_fig_width=0.7)
set_printoptions(precision=4)

cor = color_palette("Paired", n_colors=5)
rcParams["legend.fontsize"] = 5
y_name_error = r"Erro relativo ANM $\psi$ \% "


def doe_mp_plot(data_DOE, data_base, const, const_name, runs, name, dir, dp_name, save=False):

    plt_config2(relative_fig_width=0.7, page_width=180, invert=True, golden_mean=1.0)
    rcParams["legend.fontsize"] = 8
    linestyles = ["--", "-.", ":"]
    linestyle_cycle = cycle(linestyles)
    # plt.plot(sol.moc.xi_d*1e6, sol.N2Fv[0], label = 'Extrator 2')
    # ax.yaxis.set_major_formatter(y_formatter)

    for i in range(len(data_base)):

        fig = plt.figure()
        ax = fig.gca()

        df = data_base.iloc[i]
        run = df["Parâmetro"]
        N_emul = df["N_emul"]
        D43_r = df["Red D43"]
        m = df["Marco"]
        We = df["We"]
        exp = runs[run]["experiments"]  # Pega o experimento
        IDs = exp.compares[data_base["Compare"].iloc[0]]

        sol_base = runs[run][df["idict"]]
        N_esc = sol_base["N_escoam"]
        ax.plot(
            sol_base["pbe_sol"].moc.xi_d * 1e6,
            exp.get_DTG(ID=IDs[0], marco=m, teste=N_esc).freq_v / 100,
            label="Mont.",
            dashes=[2, 2],
            alpha=0.6,
        )
        ax.plot(
            sol_base["pbe_sol"].moc.xi_d * 1e6,
            exp.get_DTG(ID=IDs[1], marco=m, teste=N_esc).freq_v / 100,
            label=f"Exp. jus.",
            dashes=[5, 1],
        )
        ax.plot(
            sol_base["pbe_sol"].moc.xi_d * 1e6,
            sol_base["pbe_sol"].N2Fv,
            label=f"Sim. jus.",
        )
        ax.text(
            0.76,
            0.6,
            "{3:s} \n N: {4:d} \n We: {5:g}\n H2O: {0:.0f}\% \n {6:s}: {1:.0f}\% \n Re: {2:.1f} \n".format(
                df["H2O [%]"],
                df[f"{dp_name} [%]"],
                df["Re"],
                name,
                N_emul,
                We,
                dp_name,
            )
            + r"$R_{D43}$: "
            + "{0:.2f}".format(D43_r)
            + "\n",#+ r"$\psi_e$: {0:.1f}\%".format(df["Erro"]),
            fontsize=8,
            horizontalalignment="left",
            verticalalignment="center",
            transform=ax.transAxes,
        )
        ax.set_xlabel(r"Diâmetro da gota [$\mu$m]")
        ax.set_ylabel(r"Frequência volumétrica $f_v$ [%]")
        ax.semilogx()
        i += 1
        y_formatter = ticker.ScalarFormatter(useOffset=False)
        y_formatter.set_scientific(False)
        # ax.yaxis.set_major_formatter(y_formatter)

        for j in range(len(data_DOE)):
            df2 = data_DOE.iloc[j]
            run = df2["Parâmetro"]
            assert all(
                df[
                    [
                        "N_emul",
                        "Marco",
                        "Compare",
                        "MV01 [%]",
                        "epsilon",
                        "Re",
                        "We",
                        "Red D43",
                        "H2O [%]",
                    ]
                ]
                == df2[
                    [
                        "N_emul",
                        "Marco",
                        "Compare",
                        "MV01 [%]",
                        "epsilon",
                        "Re",
                        "We",
                        "Red D43",
                        "H2O [%]",
                    ]
                ]
            )
            sol = runs[run][df["idict"]]
            try:
                xDOE = (
                    100
                    * (
                        getattr(sol["pbe_sol"].mp, const)
                        - getattr(sol_base["pbe_sol"].mp, const)
                    )
                    / getattr(sol_base["pbe_sol"].mp, const)
                )
            except Exception:
                if const == "varsigma":
                    xDOE = (
                        100
                        * (
                            sol["pbe_sol"].moc.varsigma.mean()
                            - sol_base["pbe_sol"].moc.varsigma.mean()
                        )
                        / sol_base["pbe_sol"].moc.varsigma.mean()
                    )
            sign = "+" if xDOE >= 0 else "-"
            print(const)
            ax.plot(
                sol["pbe_sol"].moc.xi_d * 1e6,
                sol["pbe_sol"].N2Fv,
                linestyle=next(linestyle_cycle),
                label=f"Sim. jus. {sign}" + "{0:g}\%".format(xDOE) + const_name[const],
            )

        ax.legend(loc="upper left", fontsize=9)
        if save:
            fig.savefig(
                join(dir, f"run_{list(data_DOE['run_id'])}_{dp_name}_{N_esc}_{m}.pdf"),
                backend="pgf",
                bbox_inches="tight",
                pad_inches=0.05,
            )
        plt.close()
