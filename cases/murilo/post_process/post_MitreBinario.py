"""
    Função para pos processar a simulação sem modelo nenhum
    Returns:
        plots
        Dataframes
"""

from os.path import join, abspath, dirname
from os import getcwd
from sys import path as sph
import pickle
import gzip
from numpy import set_printoptions, diff
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
    D95_x_Epsilon,
    D90_d_a_Epsilon,
    error_x_,
    DTG_best_worst,
    DTGs_plot,
    error_x_D43red,
)
from pbe.setup.helpers import plt_config2

plt_config2(relative_fig_width=0.7)
set_printoptions(precision=4)

SOLUTION_NAME = "binario"

s_folder = join(dir, f"..\\solutions\\Mitre{SOLUTION_NAME}\\")
plots_out = join(s_folder, "plots")
data_out = join(s_folder, "tabelas")
data1 = "12-11-2024"

# data2 = "20-06-2024"
plots2made = {
    "Erro_X_Teste": True,
    "DTG_plot": True,
    "Convergencia_no_tempo": False,
    "Convergencia_volume_dissipacao": False,
    "Frequencia_quebra_coalescencia": False,
}

solutions = {
    f"sol_Mitre_CCE_{SOLUTION_NAME}_S1_{data1}.pickle": "Mitre CCE S1",
    f"sol_Mitre_CEM_{SOLUTION_NAME}_S1_{data1}.pickle": "Mitre CEM S1",
    f"sol_Mitre_CCE_{SOLUTION_NAME}_S3_{data1}.pickle": "Mitre CCE S3",
}

solutions_convergencia_tempo = {
    f"sol_Mitre_CEM_{SOLUTION_NAME}_5ts_{data1}.pickle": "Mitre CEM S1",
    f"sol_Mitre_CEM_{SOLUTION_NAME}_10ts_{data1}.pickle": "Mitre CEM S1",
    f"sol_Mitre_CEM_{SOLUTION_NAME}_20ts_{data1}.pickle": "Mitre CEM S1",
    f"sol_Mitre_CEM_{SOLUTION_NAME}_40ts_{data1}.pickle": "Mitre CEM S1",
    f"sol_Mitre_CEM_{SOLUTION_NAME}_50ts_{data1}.pickle": "Mitre CEM S1",
    f"sol_Mitre_CEM_{SOLUTION_NAME}_70ts_{data1}.pickle": "Mitre CEM S1",
    f"sol_Mitre_CEM_{SOLUTION_NAME}_100ts_{data1}.pickle": "Mitre CEM S1",
    f"sol_Mitre_CEM_{SOLUTION_NAME}_150ts_{data1}.pickle": "Mitre CEM S1",
}

solutions_convergencia_volume = {
    f"sol_Mitre_CEM_{SOLUTION_NAME}_0.1D_{data1}.pickle": "Mitre CEM S1",
    f"sol_Mitre_CEM_{SOLUTION_NAME}_0.5D_{data1}.pickle": "Mitre CEM S1",
    f"sol_Mitre_CEM_{SOLUTION_NAME}_1D_{data1}.pickle": "Mitre CEM S1",
    f"sol_Mitre_CEM_{SOLUTION_NAME}_2D_{data1}.pickle": "Mitre CEM S1",
    f"sol_Mitre_CEM_{SOLUTION_NAME}_2.4D_{data1}.pickle": "Mitre CEM S1",
    f"sol_Mitre_CEM_{SOLUTION_NAME}_3D_{data1}.pickle": "Mitre CEM S1",
    f"sol_Mitre_CEM_{SOLUTION_NAME}_4D_{data1}.pickle": "Mitre CEM S1",
    f"sol_Mitre_CEM_{SOLUTION_NAME}_5D_{data1}.pickle": "Mitre CEM S1",
    f"sol_Mitre_CEM_{SOLUTION_NAME}_6D_{data1}.pickle": "Mitre CEM S1",
    f"sol_Mitre_CEM_{SOLUTION_NAME}_7D_{data1}.pickle": "Mitre CEM S1",
    f"sol_Mitre_CEM_{SOLUTION_NAME}_8D_{data1}.pickle": "Mitre CEM S1",
    f"sol_Mitre_CEM_{SOLUTION_NAME}_10D_{data1}.pickle": "Mitre CEM S1",
    f"sol_Mitre_CEM_{SOLUTION_NAME}_15D_{data1}.pickle": "Mitre CEM S1", # 20
}


#
#
# Plot:
#        Error, Constants and DTG simulation results
#        1: Error x testes, colorized by H2O
#        2: two best and two worst DTG results (2x2 plot)


if plots2made["Erro_X_Teste"]:

    runs = {}
    runs_name = SOLUTION_NAME

    for solution in solutions:
        with gzip.open(join(dir, s_folder, solution), "rb") as f:
            runs[solution] = pickle.load(f)

    data = get_properties(runs, solutions)

    datared = data.groupby(by="run_id").agg(
        {
            "run_id": Series.mode,
            "Erro": ["min", "mean", "max"],
            "ts": Series.mode,
            "dt": "mean",
        }
    )
    print(datared)
    datared.columns = list(map(" ".join, datared.columns.values))  # Remove multindex

    error_x_(
        data,
        f"Erro_x_test_{runs_name}",
        plots_out,
        save=True,
    )
    error_x_(
        data,
        f"Erro_x_test_{runs_name}_D43_red",
        plots_out,
        save=True,
        size="Red D43",
    )

    error_x_(
        data,
        f"Erro_x_phi_{runs_name}",
        plots_out,
        save=True,
        x="H2O [%]",
        hue="Nome",
        style="Nome",
        xlabel=r"Concentração de água $\phi$ [-]",
    )

    error_x_(
        data,
        f"Erro_x_abertura_{runs_name}",
        plots_out,
        save=True,
        x="ANM [%]",
        hue="H2O [%]",
        xlabel=r"Abertura da válvula ANM \%",
    )

    error_x_(
        data,
        f"Erro_x_epsilon_{runs_name}",
        plots_out,
        save=True,
        x="epsilon",
        hue="H2O [%]",
        xlabel=r"Dissipação de energia cinética média $\overline{\epsilon}$ [m²/s³]",
    )

    error_x_D43red(data, f"Erro_x_D43Red_{runs_name}", plots_out, save=True)

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
            ws=0.05,
        )

    if plots2made["DTG_plot"]:
        DTGs_plot(
            data,
            runs,
            f"{runs_name}",
            join(plots_out, "DTGs"),
            save=True,
        )

    datared = data.groupby(by=["Parâmetro", "H2O [%]"]).agg(
        {
            "run_id": Series.mode,
            "Erro": ["mean", "max"],
            "Erro_Er": ["mean", "max", "min"],
            "We": "mean",
            "Re": "mean",
        }
    )

    with ExcelWriter(join(data_out, f"{runs_name}_{data1}.xlsx")) as writer:
        data.to_excel(writer, sheet_name=f"{runs_name}", index=False)
        datared.to_excel(writer, sheet_name="datared")


if plots2made["Convergencia_no_tempo"]:

    runs = {}
    runs_name = SOLUTION_NAME
    for solution in solutions_convergencia_tempo:
        with gzip.open(join(dir, s_folder, "Teste convergencia", solution), "rb") as f:
            runs[solution] = pickle.load(f)

    data = get_properties(runs, solutions_convergencia_tempo)

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
        f"erro_x_timestep_{runs_name}",
        plots_out,
        save=True,
        x="ts mode",
        y="Erro mean",
        hue=None,
        xlabel=r"Número de passos de tempo [-]",
        ylabel=r"Erro médio na ANM $\overline{\psi}$ \% ",
        style=None,
        size=None,
        fw=0.6,
        xtext=0.8,
    )

    error_x_(
        datared,
        f"erro_x_dt_{runs_name}",
        plots_out,
        save=True,
        x="dt mean",
        y="Erro mean",
        hue=None,
        xlabel=r"Passo de tempo [s]",
        ylabel=r"Erro médio na ANM $\overline{\psi}$ \% ",
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

    runs = {}
    runs_name = SOLUTION_NAME
    for solution in solutions_convergencia_volume:
        with gzip.open(
            join(dir, s_folder, "Teste convergencia", solution), "rb"
        ) as f:
            runs[solution] = pickle.load(f)

    data = get_properties(runs, solutions_convergencia_volume)

    datared = data.groupby(by="run_id").agg(
        {"ts": "median", "Erro": "mean", "dt": "mean", "fdiss": "mean"}
    )
    print(datared)

    error_x_(
        datared,
        f"erro_x_fdiss_{runs_name}",
        plots_out,
        save=True,
        x="fdiss",
        y="Erro",
        hue=None,
        xlabel=r"Razão $V_{diss}/D$",
        ylabel=r"Erro médio na ANM $\overline{\psi}$ \% ",
        style=None,
        size=None,
        fw=0.6,
        xtext=0.6,
        annot="Passos de tempo: {0:.0f} [-]".format(datared["ts"].mean()),
    )


#
#
#
#

# Plot 6:
#       Plot das frequencias de quebra e coalescencia
#       TODO: implementar!


if plots2made["Frequencia_quebra_coalescencia"]:

    runs = {}
    for x in range(1, 3):
        with open(join(dir, s_folder, solutions[x]), "rb") as f:
            runs[solutions[x]] = pickle.load(f)

    "murilo_sol_C&T_liao_17-05-2024.pickle",  # 1
    "murilo_sol_C&T_original_17-05-2024.pickle",  # 2

    data = get_properties(runs)

    datared = data.groupby(by="run_id").agg(
        {"ts": "median", "Erro": "mean", "dt": "mean", "fdiss": "mean"}
    )
    print(datared)

    error_x_(
        datared,
        "erro_x_fdiss_Mitre",
        plots_out,
        save=True,
        x="fdiss",
        y="Erro",
        hue=None,
        xlabel=r"Razão $V_{diss}/D$",
        ylabel=r"Erro médio na ANM $\overline{\psi}$ \% ",
        style=None,
        size=None,
        fw=0.6,
        xtext=0.6,
        annot="Passos de tempo: {0:.0f} [-]".format(datared["ts"].mean()),
    )
