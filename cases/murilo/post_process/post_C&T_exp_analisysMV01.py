"""
    Função para pos processar os resultados do modelo C&T
    com básico, ou seja, plot de epsilon e modelos básicos

    Returns:
        plots
"""

from os.path import join, abspath, dirname
from os import getcwd
from sys import path as sph
import pickle
import gzip
from numpy import abs, pi, set_printoptions, linspace, interp, array
from pandas import DataFrame, concat, ExcelWriter

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
from utils.plot_exp_analisys import (
    d_max_kolmogorov_length,
    d_We_Re,
)
from utils.plot_PB_basic import (
    D43Reduction_x_Epsilon,
    D95_x_Epsilon,
    dXX_x_Epsilon,
    dXX_x_Epsilon_2x2,
    d90_d_a_Epsilon,
    error_x_,
    DTG_best_worst,
    DTGs_plot,
    error_x_D43red,
)
from pbe.setup.helpers import plt_config2

plt_config2(relative_fig_width=0.7)
set_printoptions(precision=4)

SOLUTION_NAME = "basico"
VALVULA = "MV01"
IGNORAR_TESTES = [88]
data1 = "20-06-2024"

s_folder = join(dir, f"..\\solutions\\{SOLUTION_NAME}\\")
plots_out = join(s_folder, f"plots\\Exp_analisys_{VALVULA}")
data_out = join(s_folder, "tabelas")

plots2made = {
    "D43r x epsilon": True,
}
solutions = {
    f"sol_C&T_{SOLUTION_NAME}_{VALVULA}_{data1}.pickle": "C&T_basico",
}


def calc_U_max(U, A, D):
    """Calcula U max de acordo com a abertura mínima da válvula
    C é a abertura mínima
        w: Vazao kg/min
        A: Abertura de valvula 0 a 100
    """
    Ab = linspace(0.1, 1, 1000)
    y = linspace(0.0519 / 1000, D, 1000)
    Ab = linspace(0, 1, 1000)  # Abertura de 0 a 1
    y = linspace(0.0 / 1000, D, 1000)  # Altura de 0 a D
    h = interp(A / 100, Ab, y)  # altura interpolada para A
    U_max = U * D / h
    At_min = 0
    return U_max, At_min, h


def calc_L11(Relambda, ReLe):
    return (Relambda**2 / ReLe) * (3 / 2) ** (1 / 2) / 20


# Plot 1:
#        D43 reduction versus epsilon

if plots2made["D43r x epsilon"]:
    runs = {}
    runs_name = SOLUTION_NAME

    for solution in solutions:
        with gzip.open(join(dir, s_folder, solution), "rb") as f:
            runs[solution] = pickle.load(f)

    exp_basico: Import_flow_DSD2 = runs[solution]["experiments"]
    D43_r = list()
    N_esc = list()
    N_emul = list()
    mv01 = []
    mv02 = []
    marco = []
    dP = list()
    w = list()
    we = list()
    re = list()
    ch2o = []
    Umax = list()  # Max velocity in the valve
    U = list()  # Mean velocity in the tube
    ul = list()  # u': velocity flutuation
    h = list()  # altura de abertura de válvula
    k = list()  # Eddy kinetic energy
    eta = list()  # Kolmogorov length scale
    epsilon = list()  # Kolmogorov length scale
    Le = list()  # Integral Length scale
    ReLe = list()  # Reynolds for integral scale
    ReLambda = list()  # Taylor-scale Reynolds number Eq 6.266 do Pope
    ReT = list()  # turbulence reynolds number
    L11_Le = list()  # Ratio of the longitudinal integral lengthscale
    L11 = list()  # L11 longitudinal length scale
    We_pipe = list()  # Boxal Weber
    Re_pipe = list()  # Boxal Reynolds
    dv99d = list()
    dv95d = list()
    dv90d = list()
    dv99d = list()
    dn95d = list()
    dv50d = list()
    d90da = list()

    j = 0
    for r in runs:
        sol = runs[r]
        print(f"Getting {r} simulation")
        for i in sol:
            if not isinstance(i, int):
                continue
            if sol[i]['N_escoam'] in IGNORAR_TESTES:
                continue
            if sol[i].get("opt_flag") is not None:
                print("Optmization case, better use another pos-process")
                break
            dtg_a, d43_a, dtg_d, d43_d = get_dtg(sol["experiments"], sol[i])
            dv99d.append(
                1e6
                * exp_basico.get_Dx(
                    dtg_d, sol[i]["pbe_sol"].moc.xi_d, dff=99, base="vol"
                )
            )
            dv95d.append(
                1e6
                * exp_basico.get_Dx(
                    dtg_d, sol[i]["pbe_sol"].moc.xi_d, dff=95, base="vol"
                )
            )
            dv90d.append(
                1e6
                * exp_basico.get_Dx(
                    dtg_d, sol[i]["pbe_sol"].moc.xi_d, dff=90, base="vol"
                )
            )
            dn95d.append(
                1e6
                * exp_basico.get_Dx(
                    dtg_d, sol[i]["pbe_sol"].moc.xi_d, dff=95, base="num"
                )
            )
            # Arithmetic mean diameter
            dv50d.append(
                1e6
                * exp_basico.get_Dx(
                    dtg_d, sol[i]["pbe_sol"].moc.xi_d, dff=50, base="vol"
                )
            )
            d90a = 1e6 * exp_basico.get_Dx(
                dtg_a, sol[i]["pbe_sol"].moc.xi_d, dff=90, base="vol"
            )
            d90d = 1e6 * exp_basico.get_Dx(
                dtg_d, sol[i]["pbe_sol"].moc.xi_d, dff=90, base="vol"
            )

            D43_r.append(d43_d / d43_a)
            d90da.append(d90d - d90a)
            dP.append(sol[i]["exp"][f"dP-{VALVULA} [Pa]"])
            w.append(sol[i]["exp"]["FT-01[kg/min]"])
            N_emul.append(int(sol[i]["exp"]["N_emul"]))
            N_esc.append(int(sol[i]["N_escoam"]))
            mv01.append(sol[i]["exp"]["ANM"])
            mv02.append(sol[i]["exp"]["Choke"])
            marco.append(sol[i]["marco"])
            we.append(dtg_a["We"].mean())
            re.append(sol[i]["exp"]["Re"])
            ch2o.append(sol[i]["exp"]["C_agua [%]"])

            # Mitre 2014 + Pope 2000
            U_ = sol[i]["pbe_sol"].Q / (pi * sol[i]["pbe_sol"].D ** 2 / 4)
            if VALVULA == "MV02":
                U_max, At_min, h_ = calc_U_max(
                    U_, mv02[-1], sol[i]["pbe_sol"].D
                )
            else:
                U_max, At_min, h_ = calc_U_max(
                    U_, mv01[-1], sol[i]["pbe_sol"].D
                )
            U.append(U_)
            Umax.append(U_max)
            h.append(h_)
            ul.append(0.2 * U_max)
            k.append(1.5 * ul[-1] ** 2)
            nu = sol[i]["pbe_sol"].cp.mu / sol[i]["pbe_sol"].cp.rho
            epsilon.append(sol[i]["pbe_sol"].cp.epsilon)
            eta.append((nu**3 / epsilon[-1]) ** (1 / 4))
            Le.append(k[-1] ** (3 / 2) / epsilon[-1])
            ReLe.append(k[-1] ** (1 / 2) * Le[-1] / nu)
            ReLambda.append(
                (20 * ReLe[-1] / 3) ** (1 / 2)
            )  # Eq 6.64 Pope  # TODO: Ta errado como calcular corretamente o ReLambda
            # ReLambda 10² to 10^5 Hero 2021:
            # ReLambda.append((20 / 3 * k[-1] ** 2 / (epsilon[-1] * nu)) ** (1 / 2))
            L11_Le.append(calc_L11(ReLambda[-1], ReLe[-1]))
            L11.append(L11_Le[-1] * Le[-1])
            ReT.append(ul[-1] * L11[-1] / nu)  # turbulence Reynolds number

            # Boxal
            We_pipe.append(
                sol[i]["pbe_sol"].cp.rho
                * U_**2
                * sol[i]["pbe_sol"].D ** 3
                / sol[i]["pbe_sol"].dp.sigma
            )
            Re_pipe.append(
                sol[i]["pbe_sol"].cp.rho
                * U_
                * sol[i]["pbe_sol"].D
                / sol[i]["pbe_sol"].cp.mu
            )
        j += 1

    data = DataFrame(
        {
            "N_emul": N_emul,
            "N_esc": N_esc,
            "Marco": marco,
            "D43r": D43_r,
            f"dv99d {VALVULA} [um]": dv99d,
            f"dv95d {VALVULA} [um]": dv95d,
            f"dv90d {VALVULA} [um]": dv90d,
            f"dv50d {VALVULA} [um]": dv50d,
            f"d90d_a {VALVULA} [um]": d90da,  # d95 depois - antes
            "dn95d": dn95d,
            "epsilon [m²/s³]": epsilon,
            "We": we,
            "We_pipe": We_pipe,
            "Re": re,
            "Re_pipe": Re_pipe,
            "H2O [%]": ch2o,
            "MV-01 [%]": mv01,
            "MV-02 [%]": mv02,
            f"dP-{VALVULA} [Pa]": dP,
            "FT-01[kg/min]": w,
            "U [m/s]": U,
            "Umax [m/s]": Umax,
            "h [mm]": 1e3 * array(h),
            "ul [m/s]": ul,
            "k [m²/s²]": k,
            "eta [um]": 1e6 * array(eta),
            "Le [mm]": 1e3 * array(Le),
            "ReLe": ReLe,
            "ReLambda": ReLambda,
            "L11_Le": L11_Le,
            "L11 [mm]": 1e3 * array(L11),
            "ReT": ReT,
        }
    )

    d_max_kolmogorov_length(
        data,
        "dmax_x_eta",
        plots_out,
        save=True,
    )
    d_max_kolmogorov_length(
        data,
        "dmax_x_eta_H2O",
        plots_out,
        save=True,
        style="H2O [%]",
    )
    d_We_Re(data, "d_We_Re", plots_out, save=True)

    D43Reduction_x_Epsilon(data, f"Epsilon_x_D43_{VALVULA}", plots_out, save=True)

    dXX_x_Epsilon(
        data,
        f"Epsilon_x_dv90d_{VALVULA}",
        plots_out,
        save=True,
        y=f"dv90d {VALVULA} [um]",
        f=90,
        base="v",
        local=f"depois {VALVULA}",
    )  # v de volumetric n de num

    dXX_x_Epsilon_2x2(
        data,
        f"Epsilon_x_dv90d_{VALVULA}_2x2",
        plots_out,
        save=True,
        y=f"dv90d {VALVULA} [um]",
        f=90,
        base="v",
        local=f"depois {VALVULA}",
    )  # v de volumetric n de num

    dXX_x_Epsilon(
        data,
        f"Epsilon_x_dv95d_{VALVULA}",
        plots_out,
        save=True,
        y=f"dv95d {VALVULA} [um]",
        f=95,
        base="v",
        local="depois MV-01",
    )  # v de volumetric n de num

    dXX_x_Epsilon(
        data,
        f"Epsilon_x_dv99d_{VALVULA}",
        plots_out,
        save=True,
        y=f"dv99d {VALVULA} [um]",
        f=99,
        base="v",
        local="depois {VALVULA}",
    )  # v de volumetric n de num

    d90_d_a_Epsilon(data, f"Epsilon_x_d90_d-a_{VALVULA}", plots_out, save=True)

    with ExcelWriter(join(data_out, f"exp_analisys_{VALVULA}.xlsx")) as writer:
        data.to_excel(writer, sheet_name="exp_analisys", index=False)
        data.transpose().to_excel(writer, sheet_name="exp_analisys_T")
