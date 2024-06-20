from os.path import join, abspath, dirname
from os import getcwd
from sys import path as sph
import pickle
from numpy import abs, pi, set_printoptions
from pbe.setup.helpers import plt_config2
from matplotlib import pyplot as plt
from matplotlib import ticker
from seaborn import scatterplot

plt_config2(relative_fig_width=0.7)
set_printoptions(precision=4)


def DTG_exp_x_sim(sol, experiments, i, name=None, dir=None, save=False):

    m = sol["marco"]
    N = sol["N_escoam"]
    exp = sol["exp"]
    IDs = sol["compares"]
    fig = plt.figure()
    ax = fig.gca()
    # plt.plot(sol.moc.xi_d*1e6, sol.N2Fv[0], label = 'Extrator 2')
    ax.plot(
        sol["best_pbe_sol"].moc.xi_d * 1e6,
        experiments.get_DTG(ID=IDs[0], marco=m, teste=N).freq_v / 100,
        label="Antes da válvula",
    )
    ax.plot(
        sol["best_pbe_sol"].moc.xi_d * 1e6,
        sol["best_pbe_sol"].N2Fv,
        label=f"Sim. depois da válvula",
    )
    ax.plot(
        sol["best_pbe_sol"].moc.xi_d * 1e6,
        experiments.get_DTG(ID=IDs[1], marco=m, teste=N).freq_v / 100,
        label=f"Exp. depois da válvula",
    )

    ax.set_xlabel(r"Diametro da gota [$\mu$m]")
    ax.set_ylabel("Frequência volumétrica")

    ax.semilogx()
    ax.text(
        0.75,
        0.6,
        "H2O: {0:.0f}\% \n ANM: {1:.0f}\% \n Re: {2:.1f} \n Erro DTG: {3:.1f}\%".format(
            exp["C_agua [%]"], exp["ANM"], exp["Re"], 100 * sol["opt"]["error"]
        ),
        fontsize=10,
        horizontalalignment="left",
        verticalalignment="center",
        transform=ax.transAxes,
    )
    ax.legend()
    ax.grid()
    y_formatter = ticker.ScalarFormatter(useOffset=False)
    y_formatter.set_scientific(False)
    # ax.yaxis.set_major_formatter(y_formatter)
    if save:
        fig.savefig(
            join(dir, f"{name}_{N}_{m}.pdf"),
            backend="pgf",
            bbox_inches="tight",
            pad_inches=0.05,
        )
    plt.close()
    #return fig


def C_global_x_erro(C_global, Cs, name, dir):

    plt_config2(relative_fig_width=0.5)
    C_global['Error'] *= 100
    for C in Cs:
        fig = plt.figure()
        ax = fig.gca()

        scatterplot(data=C_global, x=C, y="Error", style="method", hue="method", ax=ax)
        ax.set_ylabel("Diferença Sim. e Exp. [\%]")
        ax.semilogx()
        fig.savefig(
            join(dir, f"{name}_{C}.pdf"),
            backend="pgf",
            bbox_inches="tight",
            pad_inches=0.05,
        )


def C_global_x_erro_2x2(C_global, Cs, name, dir):

    plt_config2(relative_fig_width=0.8)
    C_global['Error'] *= 100
    fig, axes = plt.subplots(nrows=2, ncols=2,
                            gridspec_kw={'width_ratios': [10, 10],
                                         'height_ratios': [10, 10],
                                         'hspace': 0.3, #'wspace': 0.0, 
                                         },
                            subplot_kw=(dict(box_aspect=1) if False else None))
    for ax in fig.get_axes():
        #ax.label_outer()
        #ax.margins(0.1)
        pass
    # custom_xlim = (0, 100)
    #plt.margins(x=3, y=3)
    # ax = fig.gca()
    i=0
    for ax in fig.get_axes():
        scatterplot(data=C_global, x=Cs[i], y="Error", style="method", hue="method", ax=ax)
        ax.set_ylabel('')
        ax.legend('')
        ax.semilogx()
        i+=1
    #ax.set_ylabel("Diferença Sim. e Exp. [\%]")
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', ncol=3, fontsize=12)
    fig.supylabel("Diferença Sim. e Exp. [\%]")
    ax.semilogx()
    fig.savefig(
            join(dir, f"{name}_C.pdf"),
            backend="pgf",
            bbox_inches="tight",
            pad_inches=0.05,
        )

"""

plt.plot(sol[0]["pbe_sol"].moc.xi_d * 1e6, sol[0]["pbe_sol"].moc.N[0] / 1e9, label="t0")
plt.plot(sol[0]["pbe_sol"].moc.xi_d * 1e6, sol[0]["pbe_sol"].moc.N[10] / 1e9, label="t1")
plt.plot(sol[0]["pbe_sol"].moc.xi_d * 1e6, sol[0]["pbe_sol"].moc.N[20] / 1e9, label="t2")
plt.semilogx()
plt.legend()


"""
