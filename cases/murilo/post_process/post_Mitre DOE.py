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
import pickle, gzip
from numpy import set_printoptions, diff
from pandas import DataFrame, concat, ExcelWriter, Series

dir = dirname(__file__)

if __name__ == "__main__":
    sph.append(abspath(join(dir, "..\\..\\..\\..\\2. APP/")))
    sph.append(abspath(join(dir, "..\\..\\..")))
from deffunctions import calc_error, get_dtg
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

s_folder = "..\\solutions\\Mitre_DOE\\"
plots_out = join(dir, "..\\", "plots")
folderModelOut = join(plots_out, "Mitre_model_DOE")
data_out = join(dir, "..\\", "tabelas")
data = "19-06-2024"

plots_made = {1: True}
solutions = [
    f"sol_Mitre_CEM_1_{data}.pickle",
    f"sol_Mitre_CEM_2_{data}.pickle",
    f"sol_Mitre_CEM_3_{data}.pickle",
    f"sol_Mitre_CEM_4_{data}.pickle",
    f"sol_Mitre_CEM_5_{data}.pickle",
    f"sol_Mitre_CEM_6_{data}.pickle",
    f"sol_Mitre_CEM_7_{data}.pickle",
    f"sol_Mitre_CEM_8_{data}.pickle",
    f"sol_Mitre_CEM_9_{data}.pickle",
    f"sol_Mitre_CEM_10_{data}.pickle",
    f"sol_Mitre_CEM_11_{data}.pickle",
    f"sol_Mitre_CEM_12_{data}.pickle",
    f"sol_Mitre_CEM_13_{data}.pickle",
    f"sol_Mitre_CEM_14_{data}.pickle",
    f"sol_Mitre_CEM_15_{data}.pickle",
    f"sol_Mitre_CEM_16_{data}.pickle",
    f"sol_Mitre_CEM_17_{data}.pickle",
    f"sol_Mitre_CEM_18_{data}.pickle",
    f"sol_Mitre_CEM_19_{data}.pickle",
    f"sol_Mitre_CEM_20_{data}.pickle",
    f"sol_Mitre_CEM_21_{data}.pickle",
    f"sol_Mitre_CEM_22_{data}.pickle",
    f"sol_Mitre_CEM_23_{data}.pickle",
    f"sol_Mitre_CEM_24_{data}.pickle",
    f"sol_Mitre_CEM_25_{data}.pickle",
]

#
#
#
#
#

# Plot 1:
#       Plot dos testes feitos com diferentes passos de tempo
# NOTE: Nemul 75 não ta batendo


if plots_made[1]:

    runs = {}
    runs_name = "Mitre_CEM_DOE"
    for x in range(0, 25):
        with gzip.open(join(dir, s_folder, solutions[x]), "rb") as f:
            runs[solutions[x]] = pickle.load(f)

    psi = []
    run = []
    run_id = []
    N_emul = []
    marco = []
    method = []
    model = []
    anm = []
    compare = []
    idict = []
    ch2o = []
    re = []
    we = []
    dt = []
    ts = []
    epsilon = []
    D43_r = []
    j = 0
    for r in runs:
        sol = runs[r]
        print(f"Getting {r} simulation")
        for i in sol:  # Iterate over each entry
            if not isinstance(i, int):
                continue
            if sol[i].get("opt_flag") is not None:
                print("Optmization case, better use another pos-process")
                break
            N_emul.append(int(sol[i]["exp"]["N_emul"]))
            idict.append(i)
            run.append(r)
            marco.append(sol[i]["marco"])
            if sol[i]["compares"] == [2, 3]:
                compare.append("E_ANM")
            anm.append(sol[i]["exp"]["ANM"])
            re.append(sol[i]["exp"]["Re"])
            dtg_a, d43_a, dtg_d, d43_d = get_dtg(sol["experiments"], sol[i])
            D43_r.append(d43_d / d43_a)
            we.append(dtg_a["We"].mean())
            ch2o.append(sol[i]["exp"]["C_agua [%]"])
            # C.append(sol[i]["opt"]["best_fit"])
            epsilon.append(sol[i]["pbe_sol"].cp.epsilon)
            psi.append(calc_error(sol[i], sol["experiments"], "psi"))
            dt.append(diff(sol[i]["pbe_sol"].time).mean())
            ts.append(len(sol[i]["pbe_sol"].time))
            run_id.append(j)
        j += 1

    data = DataFrame(run, columns=["Parâmetro"])
    data = concat(
        [
            data,
            DataFrame(
                {
                    "Erro": psi,
                    "N_emul": N_emul,
                    "Marco": marco,
                    "Compare": compare,
                    "ANM [%]": anm,
                    "epsilon": epsilon,
                    "Re": re,
                    "We": we,
                    "Red D43": D43_r,
                    "H2O [%]": ch2o,
                    "idict": idict,
                    "dt": dt,
                    "ts": ts,
                    "run_id": run_id,
                }
            ),
        ],
        axis=1,
    )
    data["Erro"] *= 100

    datared = data.groupby(by="run_id").agg(
        {
            "run_id": Series.mode,
            "Erro": ["min", "mean", "max"],
            "ts": Series.mode,
            "dt": "mean",
        }
    )
    print(datared)
    datared.columns = list(map(' '.join, datared.columns.values))  # Remove multindex

    error_x_(
        datared,
        f"erro_x_timestep_{runs_name}",
        folderModelOut,
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
        folderModelOut,
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

    

    with ExcelWriter(join(data_out, f"{runs_name}.xlsx")) as writer:
        data.to_excel(writer, sheet_name=f"{runs_name}", index=False)
        datared.to_excel(writer, sheet_name="datared")
