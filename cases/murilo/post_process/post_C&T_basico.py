"""
    Função para pos processar os resultados do modelo C&T
    com básico, ou seja, plot de epsilon e modelos básicos
    Basico porque é o modelo de outros autores sem modificar em nada
    os valores de suas constantes
    Returns:
        plots
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

s_folder = "..\\solutions\\basico\\"
plots_out = join(dir, "..\\", "plots")
data_out = join(dir, "..\\", "tabelas")
plots_made = {1: True, 2: True, 3: True, 4: True, 5: True, 6: False}

solutions = [
    "murilo_C&T_basico_15-05-2024.pickle",
    "sol_C&T_liao_20-06-2024.pickle",
    "sol_C&T_original_20-06-2024.pickle",
    "murilo_sol_C&T_murilo_22-05-2024.pickle",

    "sol_C&T_liao_0.1D_19-06-2024.pickle",  # 4
    "sol_C&T_liao_0.5D_19-06-2024.pickle",
    "sol_C&T_liao_1D_19-06-2024.pickle",
    "sol_C&T_liao_2D_19-06-2024.pickle",
    "sol_C&T_liao_2.4D_19-06-2024.pickle",
    "sol_C&T_liao_3D_19-06-2024.pickle",
    "sol_C&T_liao_4D_19-06-2024.pickle",
    "sol_C&T_liao_5D_19-06-2024.pickle",
    "sol_C&T_liao_6D_19-06-2024.pickle",
    "sol_C&T_liao_7D_19-06-2024.pickle",  # 13
    "sol_C&T_liao_5ts_19-06-2024.pickle",  # 14
    "sol_C&T_liao_10ts_19-06-2024.pickle",
    "sol_C&T_liao_20ts_19-06-2024.pickle",
    "sol_C&T_liao_40ts_19-06-2024.pickle",
    "sol_C&T_liao_50ts_19-06-2024.pickle",
    "sol_C&T_liao_60ts_19-06-2024.pickle",
    "sol_C&T_liao_100ts_19-06-2024.pickle",
    "sol_C&T_liao_150ts_19-06-2024.pickle",  # 21
]

#
#

# Plot 1:
#        D43 reduction versus epsilon


if plots_made[1]:
    with open(join(dir, s_folder, solutions[0]), "rb") as f:
        ct_basico = pickle.load(f)
    exp_basico: Import_flow_DSD2 = ct_basico["experiments"]
    D43_r = list()
    epsilon = list()
    N_emul = list()
    N_esc = list()
    marco = list()
    anm = list()
    ch2o = list()
    d95 = list()
    d90da = list()
    for i in ct_basico:
        if not isinstance(i, int):
            continue
        if i == 0 or i == 1:
            continue
        sol = ct_basico[i]
        dtg_a, d43_a, dtg_d, d43_d = get_dtg(exp_basico, sol)

        d95.append(
            1e6 * exp_basico.get_Dx(dtg_a, sol["pbe_sol"].moc.xi_d, dff=90, base="num")
        )
        d90a = 1e6 * exp_basico.get_Dx(
            dtg_a, sol["pbe_sol"].moc.xi_d, dff=90, base="vol"
        )
        d90d = 1e6 * exp_basico.get_Dx(
            dtg_d, sol["pbe_sol"].moc.xi_d, dff=90, base="vol"
        )
        d90da.append(d90d - d90a)
        N_emul.append(int(sol["exp"]["N_emul"]))
        N_esc.append(int(sol["N_escoam"]))
        marco.append(sol["marco"])
        anm.append(sol["exp"]["ANM"])
        ch2o.append(sol["exp"]["C_agua [%]"])
        D43_r.append(d43_d / d43_a)
        epsilon.append(sol["pbe_sol"].cp.epsilon)
        # Concentração, reynolds,
    data = DataFrame(
        {
            "Red D43": D43_r,
            "D95": d95,
            "D90d_a": d90da,  # d95 depois - antes
            "epsilon": epsilon,
            "N_emul": N_emul,
            "N_esc": N_esc,
            "marco": marco,
            "H2O [%]": ch2o,
            "ANM [%]": anm,
        }
    )
    D43Reduction_x_Epsilon(data, "Epsilon_x_D43_ANM", plots_out, save=False)
    D95_x_Epsilon(data, "Epsilon_x_D90_ANM", plots_out, save=False)
    D90_d_a_Epsilon(data, "Epsilon_x_D90_d-a_ANM", plots_out, save=False)

    with ExcelWriter(join(data_out, "ct_basico.xlsx")) as writer:
        data.to_excel(writer, sheet_name="ct_basico", index=False)


# Plot 2:
#        Error, Constants and DTG simulation results
#        1: Error x testes for orignal and Liao, colorized by H2O
#        2: two best and two worst DTG results (2x2 plot)
if plots_made[2]:
    #
    folderCTOut = join(plots_out, "C&T_model")
    #
    with gzip.open(join(dir, s_folder, solutions[1]), "rb") as f:
        ct_liao = pickle.load(f)

    with gzip.open(join(dir, s_folder, solutions[2]), "rb") as f:
        ct_original = pickle.load(f)

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
    Cs = ["C1", "C2", "C3", "C4"]
    runs = {"CT Original": ct_original, "CT Liao": ct_liao}
    runs_name = "CT_50timestep"
    models = {"CT Original": "CT", "CT Liao": "CT"}

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
            psi.append(calc_error(sol[i], sol["experiments"], "psi"))
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
                    "N_emul": N_emul,
                    "run_id": run_id,
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
        folderCTOut,
        save=True,
    )
    error_x_(
        data,
        f"Erro_x_test_{runs_name}_D43_red",
        folderCTOut,
        save=True,
        size="Red D43",
    )

    error_x_(
        data,
        f"Erro_x_phi_{runs_name}",
        folderCTOut,
        save=True,
        x="H2O [%]",
        hue="Parâmetro",
        style="Parâmetro",
        xlabel=r"Concentração de água $\phi$ [-]",
    )

    error_x_(
        data,
        f"Erro_x_abertura_{runs_name}",
        folderCTOut,
        save=True,
        x="ANM [%]",
        hue="H2O [%]",
        xlabel=r"Abertura da válvula ANM \%",
    )

    error_x_(
        data,
        f"Erro_x_epsilon_{runs_name}",
        folderCTOut,
        save=True,
        x="epsilon",
        hue="H2O [%]",
        xlabel=r"Dissipação de energia cinética média $\overline{\epsilon}$ [m²/s³]",
    )

    error_x_D43red(data, f"Erro_x_D43Red_{runs_name}", folderCTOut, save=True)

    to_plot = []
    x = data.groupby(by="Parâmetro")
    for run in x.groups.keys():
        to_plot = []
        to_plot.append(x.get_group(run)["Erro"].idxmin())
        to_plot.append(x.get_group(run)["Erro"].idxmax())

        DTG_best_worst(
            data.iloc[to_plot],
            runs,
            f"DTG_best_worst_{run}",
            folderCTOut,
            save=True,
            nrows=1,
            ws=0.02,
        )
    # DTG_best_worst(
    #    data.iloc[to_plot],
    #    runs,
    #    "DTG_best_worst_CT_50timestep",
    #    folderCTOut,
    #    save=False,
    # )

    DTGs_plot(
        data,
        runs,
        f"{runs_name}",
        join(folderCTOut, "DTGs"),
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
    with ExcelWriter(join(data_out, "ct_original_liao.xlsx")) as writer:
        data.to_excel(writer, sheet_name="CT", index=False)
        datared.to_excel(writer, sheet_name="datared")


# Plot 3:
#        Error, Constants and DTG simulation results
#        1: Error x testes for orignal and Liao, colorized by H2O
#        2: two best and two worst DTG results (2x2 plot)
if plots_made[3]:
    with open(join(dir, s_folder, solutions[3]), "rb") as f:
        ct_murilo = pickle.load(f)

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
    D43_r = []
    Cs = ["C1", "C2", "C3", "C4"]
    runs = {"CT opt": ct_murilo}
    models = {"CT opt": "CT"}

    j=0
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
            model.append(models[r])
            run.append(r)
            marco.append(sol[i]["marco"])
            if sol[i]["compares"] == [2, 3]:
                compare.append("E_ANM")
            anm.append(sol[i]["exp"]["ANM"])
            re.append(sol[i]["exp"]["Re"])
            ch2o.append(sol[i]["exp"]["C_agua [%]"])
            # C.append(sol[i]["opt"]["best_fit"])
            dtg_a, d43_a, dtg_d, d43_d = get_dtg(sol["experiments"], sol[i])
            we.append(dtg_a["We"].mean())
            D43_r.append(d43_d / d43_a)
            psi.append(calc_error(sol[i], sol["experiments"], "psi"))
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
                    "N_emul": N_emul,
                    "run_id": run_id,
                    "Modelo": model,
                    "Marco": marco,
                    "Compare": compare,
                    "ANM [%]": anm,
                    "We": we,
                    "Re": re,
                    "H2O [%]": ch2o,
                    "idict": idict,
                    "Red D43": D43_r,
                }
            ),
        ],
        axis=1,
    )
    data["Erro"] *= 100
    data["Erro_Er"] *= 100

    error_x_(
        data,
        "Erro_x_test_CT_optimizado",
        join(plots_out, "C&T_model"),
        save=True,
    )

    error_x_(
        data,
        "Erro_x_phi_CT_optimizado",
        join(plots_out, "C&T_model"),
        save=True,
        x="H2O [%]",
        xlabel=r"Concentração de água $\phi$ [-]",
    )

    to_plot = []
    x = data.groupby(by="Parâmetro")
    for run in x.groups.keys():
        sorted = x.get_group(run)["Erro"].sort_values()
        to_plot.append(sorted.index[0])
        to_plot.append(sorted.index[1])
        to_plot.append(sorted.index[-2])
        to_plot.append(sorted.index[-1])

    # Grafico 2x2 of DTG
    DTG_best_worst(
        data.iloc[to_plot],
        runs,
        "DTG_best_worst_CT_optimizado",
        join(plots_out, "C&T_model"),
        save=True,
    )

    DTGs_plot(
        data,
        runs,
        "DTGs_CT_optimizado",
        join(plots_out, "C&T_model\\DTGs_CT_otimizado"),
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
    name = "ct_otimizado_murilo_22-05-2024.xlsx"
    with ExcelWriter(join(data_out, name)) as writer:
        data.to_excel(writer, sheet_name="CT", index=False)
        datared.to_excel(writer, sheet_name="datared")

#
#
#
#
#
#

# Plot 4:
#       Plot dos testes feitos com diferentes passos de tempo


if plots_made[4]:
    #
    folderCTOut = join(plots_out, "C&T_model")
    #
    runs = {}
    for x in range(14, 22):
        with gzip.open(join(dir, s_folder, solutions[x]), "rb") as f:
            runs[solutions[x]] = pickle.load(f)

    Diff = []
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
            Diff.append(calc_error(sol[i], sol["experiments"], "psi"))
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
                    "Erro": Diff,
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
        {"ts": "median", "Erro": "mean", "dt": "mean"}
    )
    print(datared)

    error_x_(
        datared,
        "erro_x_timestep_CT_liao",
        folderCTOut,
        save=True,
        x="ts",
        y="Erro",
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
        "erro_x_dt_CT_liao",
        folderCTOut,
        save=True,
        x="dt",
        y="Erro",
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
    #
    folderCTOut = join(plots_out, "C&T_model")
    #
    runs = {}
    for x in range(4, 14):
        with gzip.open(join(dir, s_folder, solutions[x]), "rb") as f:
            runs[solutions[x]] = pickle.load(f)

    Diff = []
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
            Diff.append(calc_error(sol[i], sol["experiments"], "psi"))
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
                    "Erro": Diff,
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
        "erro_x_fdiss_CT_liao",
        folderCTOut,
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
    #
    folderCTOut = join(plots_out, "C&T_model")
    #
    runs = {}
    for x in range(1, 3):
        with open(join(dir, s_folder, solutions[x]), "rb") as f:
            runs[solutions[x]] = pickle.load(f)

    "murilo_sol_C&T_liao_17-05-2024.pickle",  # 1
    "murilo_sol_C&T_original_17-05-2024.pickle",  # 2

    Diff = []
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
            Diff.append(calc_error(sol[i], sol["experiments"], "DTG"))
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
                    "Erro": Diff,
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
        "erro_x_fdiss_CT_liao",
        folderCTOut,
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
