"""
    Função para pos processar os resultados do modelo do Mitre
    com as constantes dele, ou seja, plot de epsilon e modelos básicos
    Basico porque é o modelo de outros autores sem modificar em nada
    os valores de suas constantes
    Returns:
        plots
        Dataframes
"""

from os.path import join, abspath, dirname
from os import getcwd
from sys import path as sph
import pickle
import gzip
from numpy import abs, pi, set_printoptions, diff
from pandas import DataFrame, concat, ExcelWriter, Series

dir = dirname(__file__)

if __name__ == "__main__":
    sph.append(abspath(join(dir, "..\\..\\..\\..\\2. APP/")))
    sph.append(abspath(join(dir, "..\\..\\..")))

from deffunctions import calc_error, get_dtg, get_properties
from pbe.app.dtg_class import (
    DTGSolution,
    Import_flow_DSD2,
    DTG_experiment,
    get_location,
)
from utils.plot_PB_basic import (
    D43Reduction_x_Epsilon,
    dXX_x_Epsilon,
    dXX_x_Epsilon_2x2,
    dXX_x_Capilaridade_2x2,
    d90_d_a_Epsilon,
    error_x_,
    DTG_best_worst,
    DTGs_plot,
    error_x_D43red,
)
from pbe.setup.helpers import plt_config2

plt_config2(relative_fig_width=0.7)
set_printoptions(precision=4)

SOLUTION_NAME = "MitreOriginal"
MODELO = "Mitre_CCE_CEM_original"
VALVULA = "MV01"
IGNORAR_TESTES = [88]
data1 = "03-12-2024"
data2 = "03-12-2024"
data1 = "29-11-2024"
data2 = "29-11-2024"

s_folder = join(dir, f"..\\solutions\\{SOLUTION_NAME}\\")
plots_out = join(s_folder, "plots")
data_out = join(s_folder, "tabelas")

plots2made = {
    "Erro_X_Teste": True,
    "Erro_X_Teste_otimizado": False,
    "DTGs": False,
    "Convergencia_no_tempo": True,
    "Convergencia_volume_dissipacao": True,
    "Frequencia_quebra_coalescencia": False,
}

solutions = {
    f"sol_Mitre_{VALVULA}_CCE_Original_S1_{data2}.pickle": "CCE S1",
    f"sol_Mitre_{VALVULA}_CCE_Original_S3_{data2}.pickle": "CCE S3",
    f"sol_Mitre_{VALVULA}_CEM_Original_S1_{data2}.pickle": "CEM S1",
}
solutions_otimizado = {
    # f"sol_C&T_murilo_{VALVULA}_{data2}.pickle": "CT opt",
}
solutions_convergencia_tempo = {
    f"sol_Mitre_{VALVULA}_CEM_Original_5ts_{data1}.pickle": "CEM",  # 14
    f"sol_Mitre_{VALVULA}_CEM_Original_10ts_{data1}.pickle": "CEM",
    f"sol_Mitre_{VALVULA}_CEM_Original_20ts_{data1}.pickle": "CEM",
    f"sol_Mitre_{VALVULA}_CEM_Original_40ts_{data1}.pickle": "CEM",
    f"sol_Mitre_{VALVULA}_CEM_Original_50ts_{data1}.pickle": "CEM",
    f"sol_Mitre_{VALVULA}_CEM_Original_70ts_{data1}.pickle": "CEM",
    f"sol_Mitre_{VALVULA}_CEM_Original_100ts_{data1}.pickle": "CEM",
    f"sol_Mitre_{VALVULA}_CEM_Original_150ts_{data1}.pickle": "CEM",  # 21
}

solutions_convergencia_volume = {
    f"sol_Mitre_{VALVULA}_CEM_Original_0.1D_{data1}.pickle": "CEM",  # 4
    f"sol_Mitre_{VALVULA}_CEM_Original_0.5D_{data1}.pickle": "CEM",
    f"sol_Mitre_{VALVULA}_CEM_Original_1D_{data1}.pickle": "CEM",
    f"sol_Mitre_{VALVULA}_CEM_Original_2D_{data1}.pickle": "CEM",
    f"sol_Mitre_{VALVULA}_CEM_Original_2.4D_{data1}.pickle": "CEM",
    f"sol_Mitre_{VALVULA}_CEM_Original_3D_{data1}.pickle": "CEM",
    f"sol_Mitre_{VALVULA}_CEM_Original_4D_{data1}.pickle": "CEM",
    f"sol_Mitre_{VALVULA}_CEM_Original_5D_{data1}.pickle": "CEM",
    f"sol_Mitre_{VALVULA}_CEM_Original_6D_{data1}.pickle": "CEM",
    f"sol_Mitre_{VALVULA}_CEM_Original_7D_{data1}.pickle": "CEM",  # 13
    f"sol_Mitre_{VALVULA}_CEM_Original_8D_{data1}.pickle": "CEM",
    f"sol_Mitre_{VALVULA}_CEM_Original_10D_{data1}.pickle": "CEM",
    f"sol_Mitre_{VALVULA}_CEM_Original_15D_{data1}.pickle": "CEM",
}

# Plot 2:
#        Error, Constants and DTG simulation results
#        1: Error x testes for orignal and Liao, colorized by H2O
#        2: two best and two worst DTG results (2x2 plot)
if plots2made["Erro_X_Teste"]:

    runs = {}
    # runs_name = SOLUTION_NAME
    # folderCTOut = join(plots_out, "C&T_model")

    for solution in solutions:
        with gzip.open(join(dir, s_folder, solution), "rb") as f:
            runs[solution] = pickle.load(f)

    data = get_properties(runs, solutions, IGNORAR_TESTES)
    datared = data.groupby(by=["Nome", "H2O [%]"]).agg(
        {
            "run_id": Series.mode,
            "Erro": ["mean", "max"],
            "Erro_Er": ["mean", "max", "min"],
            "We": "mean",
            "Re": "mean",
        }
    )

    j = 0
    for r in runs:
        sol = runs[r]
        print(f"Getting {r} simulation")
        for i in sol:  # Iterate over each entry
            if not isinstance(i, int):
                continue
            if sol[i]["N_escoam"] in IGNORAR_TESTES:
                continue
            if sol[i].get("opt_flag") is not None:
                print("Optmization case, better use another pos-process")
                break
            dtg_a, d43_a, dtg_d, d43_d = get_dtg(sol["experiments"], sol[i])
        j += 1

    error_x_(
        data,
        f"Erro_x_test_{SOLUTION_NAME}_{VALVULA}",
        plots_out,
        save=True,
        ylabel=f"Erro relativo {VALVULA} " + r"$\psi$ \% ",
        hue=f"{VALVULA} [%]",
    )
    error_x_(
        data,
        f"Erro_x_test_{SOLUTION_NAME}_{VALVULA}_D43_red",
        plots_out,
        save=True,
        ylabel=f"Erro relativo {VALVULA} " + r"$\psi$ \% ",
        size="Red D43",
        hue=f"{VALVULA} [%]",
    )

    error_x_(
        data,
        f"Erro_x_phi_{SOLUTION_NAME}_{VALVULA}",
        plots_out,
        save=True,
        x="H2O [%]",
        hue="Nome",
        style="Nome",
        xlabel=r"Concentração de água [-]",   #$\phi$
        ylabel=f"Erro relativo {VALVULA} " + r"$\psi$ \% ",
    )

    error_x_(
        data,
        f"Erro_x_abertura_{SOLUTION_NAME}_{VALVULA}",
        plots_out,
        save=True,
        x=f"{VALVULA} [%]",
        hue="H2O [%]",
        xlabel=f"Abertura da válvula {VALVULA} %",
        ylabel=f"Erro relativo {VALVULA} " + r"$\psi$ \% ",
    )

    error_x_(
        data,
        f"Erro_x_epsilon_{SOLUTION_NAME}_{VALVULA}",
        plots_out,
        save=True,
        x="epsilon",
        hue="H2O [%]",
        xlabel=r"Dissipação de energia cinética $\overline{\epsilon}$ [m²/s³]",
        ylabel=f"Erro relativo {VALVULA} " + r"$\psi$ \% ",
    )

    error_x_D43red(
        data,
        f"Erro_x_D43Red_{SOLUTION_NAME}_{VALVULA}",
        plots_out,
        save=True,
        hue=f"{VALVULA} [%]",
        ylabel=f"Erro relativo {VALVULA} " + r"$\psi$ \% ",
        dp_name=VALVULA,
    )

    to_plot = []
    x = data.groupby(by="Nome")
    for run in x.groups.keys():
        to_plot = []
        to_plot.append(x.get_group(run)["Erro"].idxmin())
        to_plot.append(x.get_group(run)["Erro"].idxmax())

        DTG_best_worst(
            data.iloc[to_plot],
            runs,
            f"DTG_best_worst_{run}",
            plots_out,
            save=True,
            nrows=1,
            ws=0.02,
            dp_name=VALVULA,
        )

    DTGs_plot(
        data,
        runs,
        f"{SOLUTION_NAME}",
        join(plots_out, f"DTGs_{SOLUTION_NAME}"),
        dp_name=VALVULA,
        save=plots2made["DTGs"],
    )

    with ExcelWriter(join(data_out, f"{MODELO}_{VALVULA}.xlsx")) as writer:
        data.to_excel(writer, sheet_name=f"{SOLUTION_NAME}", index=False)
        datared.to_excel(writer, sheet_name="datared")

#
#
#
# Plot 4:
#       Plot dos testes feitos com diferentes passos de tempo


if plots2made["Convergencia_no_tempo"]:
    #
    # folderCTOut = join(plots_out, "C&T_model")
    #
    runs = {}
    runs_name = "Mitre"
    for solution in solutions_convergencia_tempo:
        with gzip.open(join(dir, s_folder, "Teste convergencia", solution), "rb") as f:
            runs[solution] = pickle.load(f)

    data = get_properties(runs, solutions_convergencia_tempo, IGNORAR_TESTES)

    datared = data.groupby(by="run_id").agg(
        {
            "run_id": Series.mode,
            "Erro": ["mean", "max"],
            "ts": Series.mode,
            "dt": "mean",
        }
    )
    print(datared)
    datared.columns = list(map(" ".join, datared.columns.values))  # Remove multindex

    error_x_(
        datared,
        f"erro_x_timestep_{runs_name}_{VALVULA}",
        plots_out,
        save=True,
        x="ts mode",
        y="Erro mean",
        hue=None,
        xlabel=r"Número de passos de tempo [-]",
        ylabel=f"Erro médio {VALVULA} " + r"$\overline{\psi}$ \% ",
        style=None,
        size=None,
        fw=0.6,
        xtext=0.8,
    )

    error_x_(
        datared,
        f"erro_x_dt_{runs_name}_{VALVULA}",
        plots_out,
        save=True,
        x="dt mean",
        y="Erro mean",
        hue=None,
        xlabel=r"Passo de tempo [s]",
        ylabel=f"Erro médio {VALVULA} " + r"$\overline{\psi}$ \% ",
        style=None,
        size=None,
        fw=0.6,
        xtext=0.8,
    )


#
#
#
#
#

# Plot 5:
#       Plot dos testes feitos com diferentes comprimentos
#       de dissipação


if plots2made["Convergencia_volume_dissipacao"]:
    #
    runs = {}
    runs_name = "Mitre"

    for solution in solutions_convergencia_volume:
        with gzip.open(
            join(dir, s_folder, "Teste convergencia", solution), "rb"
        ) as f:
            runs[solution] = pickle.load(f)

    data = get_properties(runs, solutions_convergencia_volume, IGNORAR_TESTES)

    datared = data.groupby(by="run_id").agg(
        {
            "run_id": Series.mode,
            "ts": Series.mode,
            "Erro": ["mean", "max"],
            "dt": "mean",
            "fdiss": "mean",
        }
    )
    print(datared)
    datared.columns = list(map(" ".join, datared.columns.values))  # Remove multindex

    error_x_(
        datared,
        f"erro_x_fdiss_{runs_name}_{VALVULA}",
        plots_out,
        save=True,
        x="fdiss mean",
        y="Erro mean",
        hue=None,
        xlabel=r"Razão $V_{diss}/D$",
        ylabel=f"Erro médio {VALVULA} " + r"$\overline{\psi}$ \% ",
        style=None,
        size=None,
        fw=0.6,
        xtext=0.6,
        annot="Passos de tempo: {0:.0f} [-]".format(datared["ts mode"].mean()),
    )


#
#
#
#
#
#

# Plot 6:
#       Plot das frequencias de quebra e coalescencia
#       TODO: implementar!


if plots2made["Frequencia_quebra_coalescencia"]:

    "murilo_sol_C&T_liao_17-05-2024.pickle",  # 1
    "murilo_sol_C&T_original_17-05-2024.pickle",  # 2

    xi = []
    gamma = []
    epsilon = []
    D43_r = []
    j = 0
    for r in runs:
        sol = runs[r]
        print(f"Getting {r} simulation")
        for i in sol:  # Iterate over each entry
            if not isinstance(i, int):
                continue
            if sol[i]["N_escoam"] in IGNORAR_TESTES:
                continue
            if sol[i].get("opt_flag") is not None:
                print("Optmization case, better use another pos-process")
                break

            xi.append(sol[i]["pbe_sol"].moc.xi)
            gamma.append(sol[i]["pbe_sol"].moc.gamma)
            beta = sol[i]["pbe_sol"].moc.beta
            Q = sol[i]["pbe_sol"].moc.Q
        j += 1

    data = concat(
        [
            data,
            DataFrame(
                {
                    "Gamma": gamma,
                }
            ),
        ],
        axis=1,
    )
