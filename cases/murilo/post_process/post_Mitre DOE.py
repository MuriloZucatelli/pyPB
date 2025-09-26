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
from numpy import set_printoptions, diff
from pandas import DataFrame, concat, ExcelWriter, Series
from deffunctions import calc_error, get_dtg, get_properties

dir = dirname(__file__)

if __name__ == "__main__":
    sph.append(abspath(join(dir, "..\\..\\..\\..\\2. APP/"))) # Pasta APP
    sph.append(abspath(join(dir, "..\\..\\..")))  # Pasta pbe

from pbe.app.dtg_class import (
    DTGSolution,
    Import_flow_DSD2,
    DTG_experiment,
    get_location,
)
from utils.plot_PB_basic import (
    error_x_,
)
from utils.plot_PB_DOE import (
    doe_mp_plot,
)

from pbe.setup.helpers import plt_config2

plt_config2(relative_fig_width=0.7)
set_printoptions(precision=4)

SOLUTION_NAME = "MitreDOE2"
MODELO = "Mitre_CEM"
VALVULA = "MV01"
IGNORAR_TESTES = [88]
data = "19-12-2024"

s_folder = join(dir, f"..\\solutions\\{SOLUTION_NAME}\\")
plots_out = join(s_folder, "plots")
data_out = join(s_folder, "tabelas")
DTG_analise = {67: {3}}  # N_escoam 107 marco 3
DTG_analise = {63: {2}}  # N_escoam 104 marco 2
DTG_analise = {81: {3}}  # N_escoam 104 marco 2

plots2made = {
    "Erro_X_Teste": True,
    "Erro_X_Teste_otimizado": False,
    "DTGs": False,
}
solutions = {
    f"sol_{VALVULA}_Mitre_CEM_1_{data}.pickle": "CEM",   # 0
    f"sol_{VALVULA}_Mitre_CEM_2_{data}.pickle": "CEM",   # 1
    f"sol_{VALVULA}_Mitre_CEM_3_{data}.pickle": "CEM",   # 2
    f"sol_{VALVULA}_Mitre_CEM_4_{data}.pickle": "CEM",   # 3
    f"sol_{VALVULA}_Mitre_CEM_5_{data}.pickle": "CEM",   # 4
    f"sol_{VALVULA}_Mitre_CEM_6_{data}.pickle": "CEM",   # 5
    f"sol_{VALVULA}_Mitre_CEM_7_{data}.pickle": "CEM",   # 6
    f"sol_{VALVULA}_Mitre_CEM_8_{data}.pickle": "CEM",   # 7
    f"sol_{VALVULA}_Mitre_CEM_9_{data}.pickle": "CEM",   # 8
    f"sol_{VALVULA}_Mitre_CEM_10_{data}.pickle": "CEM",  # 9
    f"sol_{VALVULA}_Mitre_CEM_11_{data}.pickle": "CEM",  # 10
    f"sol_{VALVULA}_Mitre_CEM_12_{data}.pickle": "CEM",  # 11
    f"sol_{VALVULA}_Mitre_CEM_13_{data}.pickle": "CEM",  # 12
    f"sol_{VALVULA}_Mitre_CEM_14_{data}.pickle": "CEM",  # 13
    f"sol_{VALVULA}_Mitre_CEM_15_{data}.pickle": "CEM",  # 14
    f"sol_{VALVULA}_Mitre_CEM_16_{data}.pickle": "CEM",  # 15
    f"sol_{VALVULA}_Mitre_CEM_17_{data}.pickle": "CEM",  # 16
    f"sol_{VALVULA}_Mitre_CEM_18_{data}.pickle": "CEM",  # 17
    f"sol_{VALVULA}_Mitre_CEM_19_{data}.pickle": "CEM",  # 18
    f"sol_{VALVULA}_Mitre_CEM_20_{data}.pickle": "CEM",  # 19
    f"sol_{VALVULA}_Mitre_CEM_21_{data}.pickle": "CEM",  # 20
    f"sol_{VALVULA}_Mitre_CEM_22_{data}.pickle": "CEM",  # 21
    f"sol_{VALVULA}_Mitre_CEM_23_{data}.pickle": "CEM",  # 22
    f"sol_{VALVULA}_Mitre_CEM_24_{data}.pickle": "CEM",  # 23
    f"sol_{VALVULA}_Mitre_CEM_25_{data}.pickle": "CEM",  # 24
}

# indices da solução que possui variação dos parâmetros
idx_solut_base = 24
idx_solut = {"Cc": [16, 17], "Ce": [18, 19], "Cb": [20, 21], "varsigma": [22, 23]}
const_name = {"Cc": r"$C_c$", "Ce": r"$C_e$", "Cb": r"$C_b$", "varsigma": r"$\varsigma$"}

#
#
#
#
#

# Plot 1:
#       Plot dos testes feitos com diferentes passos de tempo


if plots2made["Erro_X_Teste"]:

    runs = {}
    # runs_name = SOLUTION_NAME
    # folderCTOut = join(plots_out, "C&T_model")

    for solution in solutions:
        with gzip.open(join(dir, s_folder, solution), "rb") as f:
            runs[solution] = pickle.load(f)

    data = get_properties(runs, solutions, IGNORAR_TESTES)
    datared = data.groupby(by=["run_id", "H2O [%]"]).agg(
        {
            # "run_id": Series.mode,
            "Erro": ["min", "mean", "max"],
            "Erro_Er": ["mean", "max", "min"],
            "ts": Series.mode,
            "We": "mean",
            "Re": "mean",
        }
    )

    print(datared)
    datared.columns = list(map(" ".join, datared.columns.values))  # Remove multindex

    error_x_(
        datared,
        f"Erro_x_timestep_{MODELO}_{VALVULA}",
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
        f"Erro_x_dt_{MODELO}_{VALVULA}",
        plots_out,
        save=True,
        x="run_id",
        y="Erro mean",
        hue=None,
        xlabel="i",
        ylabel=f"Erro médio {VALVULA} " + r"$\overline{\psi}$ \% ",
        style=None,
        size=None,
        fw=0.6,
        xtext=0.8,
    )

    if DTG_analise is not None:
        for teste in DTG_analise:
            for m in list(DTG_analise[teste]):
                data_base = data[
                    (data["N_emul"] == teste)
                    & (data["Marco"] == m)
                    & (data["run_id"] == idx_solut_base)
                ]
                for const in idx_solut:
                    idx_solut[const]
                    data[(data["N_emul"] == teste) & (data["Marco"] == m)]
                    data_DOE = data[
                        (data["N_emul"] == teste)
                        & (data["Marco"] == m)
                        & (data["run_id"].isin(idx_solut[const]))
                    ]

                    doe_mp_plot(
                        data_DOE,
                        data_base,
                        const,
                        const_name,
                        runs,
                        f"{MODELO}",
                        join(plots_out, "DTGs"),
                        VALVULA,
                        save=True,
                    )
                # idx_solut

    with ExcelWriter(join(data_out, f"{MODELO}_{VALVULA}.xlsx")) as writer:
        data.to_excel(writer, sheet_name=f"{MODELO}", index=False)
        datared.to_excel(writer, sheet_name="datared")
