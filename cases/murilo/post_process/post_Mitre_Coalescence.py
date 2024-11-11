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

s_folder = "..\\solutions\\Mitre_Coalescence\\"
plots_out = join(dir, "..\\", "plots")
folderModelOut = join(plots_out, "Mitre_Coalescence")
data_out = join(dir, "..\\", "tabelas")
data1 = "19-06-2024"
data2 = "26-09-2024"
plots_made = {2: True, 4: True, 5: True, 6: False}
# 2: Plot
# 4:  
# 5: 
# 6: 

solutions = [
    f"sol_Mitre_C_CCE_S1_{data2}.pickle",
    f"sol_Mitre_C_CEM_S1_{data2}.pickle",
    f"sol_Mitre_C_CCE_S3_{data2}.pickle",
    # f"sol_Mitre_CEM_5ts_{data1}.pickle",  # 3
    # f"sol_Mitre_CEM_10ts_{data1}.pickle",
    # f"sol_Mitre_CEM_20ts_{data1}.pickle",
    # f"sol_Mitre_CEM_40ts_{data1}.pickle",
    # f"sol_Mitre_CEM_50ts_{data1}.pickle",
    # f"sol_Mitre_CEM_70ts_{data1}.pickle",
    # f"sol_Mitre_CEM_100ts_{data1}.pickle",
    # f"sol_Mitre_CEM_150ts_{data1}.pickle",  # 10
    # f"sol_Mitre_CEM_0.1D_{data1}.pickle",  # 11
    # f"sol_Mitre_CEM_0.5D_{data1}.pickle",
    # f"sol_Mitre_CEM_1D_{data1}.pickle",
    # f"sol_Mitre_CEM_2D_{data1}.pickle",
    # f"sol_Mitre_CEM_2.4D_{data1}.pickle",
    # f"sol_Mitre_CEM_3D_{data1}.pickle",
    # f"sol_Mitre_CEM_4D_{data1}.pickle",
    # f"sol_Mitre_CEM_5D_{data1}.pickle",
    # f"sol_Mitre_CEM_6D_{data1}.pickle",
    # f"sol_Mitre_CEM_7D_{data1}.pickle",
    # f"sol_Mitre_CEM_8D_{data1}.pickle",
    # f"sol_Mitre_CEM_10D_{data1}.pickle",  # 22
    # f"sol_Mitre_CEM_15D_{data1}.pickle",  # 23
]

runs_name = "Mitre_C_100ts_5D"

#
#
# Plot 2:
#        Error, Constants and DTG simulation results
#        1: Error x testes, colorized by H2O
#        2: two best and two worst DTG results (2x2 plot)
if plots_made[2]:
    #
    #
    with gzip.open(join(dir, s_folder, solutions[0]), "rb") as f:
        mitre_fq_S1 = pickle.load(f)
        exp: Import_flow_DSD2 = mitre_fq_S1["experiments"]

    with gzip.open(join(dir, s_folder, solutions[1]), "rb") as f:
        mitre_fqN_S1 = pickle.load(f)

    with gzip.open(join(dir, s_folder, solutions[2]), "rb") as f:
        mitre_fq_S3 = pickle.load(f)

    psi = []
    Er = []
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
    epsilon = []
    D43_r = []
    runs = {
        "Mitre CCE S1": mitre_fq_S1,
        "Mitre CCE S3": mitre_fq_S3,
        "Mitre CEM S1": mitre_fqN_S1,
    }
    models = {"Mitre CCE S1": "Mitre", "Mitre CCE S3": "Mitre", "Mitre CEM S1": "Mitre"}
    
    j = 0
    for r in runs:
        sol = runs[r]
        print(f"Getting {r} simulation")
        for i in sol:  # Iterate over each entry
            if not isinstance(i, int):
                continue
            if i == 0 or i == 1:  # Pular teste 88
                continue
            if sol[i].get("opt_flag") is not None:
                print("Optmization case, better use another pos-process")
                break
            N_emul.append(int(sol[i]["exp"]["N_emul"]))
            idict.append(i)
            model.append(models[r])
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
            epsilon.append(sol[i]["pbe_sol"].cp.epsilon)
            psi.append(calc_error(sol[i], sol["experiments"]))
            Er.append(calc_error(sol[i], sol["experiments"], "Er"))
            run_id.append(j)
        j += 1

    data = DataFrame(run, columns=["Parâmetro"])
    data = concat(
        [
            data,
            DataFrame(
                {
                    "Erro": psi,
                    "Erro_Er": Er,
                    "run_id": run_id,
                    "N_emul": N_emul,
                    "Modelo": model,
                    "Marco": marco,
                    "Compare": compare,
                    "ANM [%]": anm,
                    "epsilon": epsilon,
                    "Re": re,
                    "We": we,
                    "Red D43": D43_r,
                    "H2O [%]": ch2o,
                    "idict": idict,
                }
            ),
        ],
        axis=1,
    )
    data["Erro"] *= 100
    data["Erro_Er"] *= 100

    error_x_(
        data,
        f"Erro_x_test_{runs_name}",
        folderModelOut,
        save=True,
    )
    error_x_(
        data,
        f"Erro_x_test_{runs_name}_D43_red",
        folderModelOut,
        save=True,
        size="Red D43",
    )

    error_x_(
        data,
        f"Erro_x_phi_{runs_name}",
        folderModelOut,
        save=True,
        x="H2O [%]",
        hue="Parâmetro",
        style="Parâmetro",
        xlabel=r"Concentração de água $\phi$ [-]",
    )

    error_x_(
        data,
        f"Erro_x_abertura_{runs_name}",
        folderModelOut,
        save=True,
        x="ANM [%]",
        hue="H2O [%]",
        xlabel=r"Abertura da válvula ANM \%",
    )

    error_x_(
        data,
        f"Erro_x_epsilon_{runs_name}",
        folderModelOut,
        save=True,
        x="epsilon",
        hue="H2O [%]",
        xlabel=r"Dissipação de energia cinética média $\overline{\epsilon}$ [m²/s³]",
    )

    error_x_D43red(data, f"Erro_x_D43Red_{runs_name}", folderModelOut, save=True)

    x = data.groupby(by="Parâmetro")
    for run in x.groups.keys():
        to_plot = []
        to_plot.append(x.get_group(run)["Erro"].idxmin())
        to_plot.append(x.get_group(run)["Erro"].idxmax())

        DTG_best_worst(
            data.iloc[to_plot],
            runs,
            f"DTG_best_worst_{run}",
            folderModelOut,
            save=True,
            nrows=1,
            ws=0.05,
        )

    DTGs_plot(
        data,
        runs,
        f"{runs_name}",
        join(folderModelOut, "DTGs"),
        save=False,
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

    with ExcelWriter(join(data_out, f"{runs_name}_{data2}.xlsx")) as writer:
        data.to_excel(writer, sheet_name=f"{runs_name}", index=False)
        datared.to_excel(writer, sheet_name="datared")


#
#
#
#

# Plot 4:
#       Plot dos testes feitos com diferentes passos de tempo


if plots_made[4]:

    runs = {}
    runs_name = "Mitre_CEM"
    for x in range(3, 11):
        with gzip.open(join(dir, s_folder, 'Teste convergencia', solutions[x]), "rb") as f:
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
            if i == 0 or i == 1:
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
            "Erro": ["mean", "max"],
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


#
#
#
#
#

# Plot 5:
#       Plot dos testes feitos com diferentes comprimentos
#       de dissipação


if plots_made[5]:

    runs = {}
    for x in range(11, 24):
        with gzip.open(join(dir, s_folder, 'Teste convergencia', solutions[x]), "rb") as f:
            runs[solutions[x]] = pickle.load(f)

    psi = []
    run = []
    run_id = []
    C = []
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
    fdiss = []
    epsilon = []
    D43_r = []
    Cs = ["C1", "C2", "C3", "C4"]
    j = 0
    for r in runs:
        sol = runs[r]
        print(f"Getting {r} simulation")
        for i in sol:  # Iterate over each entry
            if not isinstance(i, int):
                continue
            if i == 0 or i == 1:
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
            fdiss.append((sol[i]["pbe_sol"].fdiss))
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
                    "fdiss": fdiss,
                    "run_id": run_id,
                }
            ),
        ],
        axis=1,
    )
    data["Erro"] *= 100

    datared = data.groupby(by="run_id").agg(
        {"ts": "median", "Erro": "mean", "dt": "mean", "fdiss": "mean"}
    )
    print(datared)

    error_x_(
        datared,
        f"erro_x_fdiss_{runs_name}",
        folderModelOut,
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


if plots_made[6]:

    runs = {}
    for x in range(1, 3):
        with open(join(dir, s_folder, solutions[x]), "rb") as f:
            runs[solutions[x]] = pickle.load(f)

    "murilo_sol_C&T_liao_17-05-2024.pickle",  # 1
    "murilo_sol_C&T_original_17-05-2024.pickle",  # 2

    psi = []
    run = []
    run_id = []
    C = []
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
    fdiss = []
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
            if i == 0 or i == 1:
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
            epsilon.append(sol[i]["pbe_sol"].cp.epsilon)
            psi.append(calc_error(sol[i], sol["experiments"], "psi"))
            dt.append(diff(sol[i]["pbe_sol"].time).mean())
            ts.append(len(sol[i]["pbe_sol"].time))
            fdiss.append((sol[i]["pbe_sol"].fdiss))
            xi.append(sol[i]["pbe_sol"].moc.xi)
            gamma.append(sol[i]["pbe_sol"].moc.gamma)
            beta = sol[i]["pbe_sol"].moc.beta
            Q = sol[i]["pbe_sol"].moc.Q
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
                    "fdiss": fdiss,
                    "run_id": run_id,
                }
            ),
        ],
        axis=1,
    )
    data["Erro"] *= 100

    datared = data.groupby(by="run_id").agg(
        {"ts": "median", "Erro": "mean", "dt": "mean", "fdiss": "mean"}
    )
    print(datared)

    error_x_(
        datared,
        "erro_x_fdiss_Mitre",
        folderModelOut,
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
