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
from numpy import abs, pi, set_printoptions, linspace, interp, array
from pandas import DataFrame, concat, ExcelWriter

dir = dirname(__file__)

if __name__ == "__main__":
    sph.append(abspath(join(dir, "..\\..\\..\\..\\2. APP/")))
    sph.append(abspath(join(dir, "..\\..\\..")))
from pbe.app.dtg_class import (
    DTGSolution,
    Import_flow_DSD2,
    DTG_experiment,
    get_location,
)
from utils.plot_exp_analisys import d_max_kolmogorov_length, d_We_Re
from pbe.setup.helpers import plt_config2

plt_config2(relative_fig_width=0.7)
set_printoptions(precision=4)

s_folder = "..\\solutions\\basico\\"
plots_out = join(dir, "..\\", "plots")
data_out = join(dir, "..\\", "tabelas")
plots_made = {1: True, 2: False, 3: False}
solutions = [
    "murilo_C&T_basico_15-05-2024.pickle",
]


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


def get_dtg(exp_basico, sol):
    dtg_a = exp_basico.get_DTG(
        teste=sol["N_escoam"], marco=sol["marco"], ID=sol["compares"][0]
    )
    dtg_d = exp_basico.get_DTG(
        teste=sol["N_escoam"], marco=sol["marco"], ID=sol["compares"][1]
    )
    d43_a = exp_basico.get_D_caracteristico(dtg_a, sol["pbe_sol"].moc.xi_d, dff=[4, 3])
    d43_d = exp_basico.get_D_caracteristico(dtg_d, sol["pbe_sol"].moc.xi_d, dff=[4, 3])
    return dtg_a, d43_a, dtg_d, d43_d


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
    re = list()
    dP = list()
    w = list()
    we = list()
    Umax = list()  # Max velocity in the valve
    U = list()  # Mean velocity in the tube
    ul = list()  # u': velocity flutuation
    h = list()  # altura de abertura de válvula
    k = list()  # Eddy kinetic energy
    eta = list()  # Kolmogorov length scale
    Le = list()  # Integral Length scale
    ReLe = list()  # Reynolds for integral scale
    ReLambda = list()  # Taylor-scale Reynolds number Eq 6.266 do Pope
    ReT = list()  # turbulence reynolds number
    L11_Le = list()  # Ratio of the longitudinal integral lengthscale
    L11 = list()  # L11 longitudinal length scale
    We_pipe = list()  # Boxal Weber
    Re_pipe = list()  # Boxal Reynolds
    d99d = list()
    d50d = list()
    d90da = list()
    for i in ct_basico:
        if not isinstance(i, int):
            continue
        if i == 0 or i == 1:
            continue
        sol = ct_basico[i]

        dtg_a, d43_a, dtg_d, d43_d = get_dtg(exp_basico, sol)

        d99d.append(
            1e6 * exp_basico.get_Dx(dtg_d, sol["pbe_sol"].moc.xi_d, dff=99, base="vol")
        )
        # Arithmetic mean diameter
        d50d.append(
            1e6 * exp_basico.get_Dx(dtg_d, sol["pbe_sol"].moc.xi_d, dff=50, base="vol")
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
        dP.append(sol["exp"]["dPGV-ANM [Pa]"])
        re.append(sol["exp"]["Re"])
        we.append(dtg_a["We"].mean())
        w.append(sol["exp"]["FT-01[kg/min]"])
        # Mitre 2014 + Pope 2000
        U_ = sol["pbe_sol"].Q / (pi * sol["pbe_sol"].D ** 2 / 4)
        U_max, At_min, h_ = calc_U_max(U_, anm[-1], sol["pbe_sol"].D)  # TODO: Melhorar
        U.append(U_)
        Umax.append(U_max)
        h.append(h_)
        ul.append(0.2 * U_max)
        k.append(1.5 * ul[-1] ** 2)
        nu = sol["pbe_sol"].cp.mu / sol["pbe_sol"].cp.rho
        eta.append((nu**3 / epsilon[-1]) ** (1 / 4))
        Le.append(k[-1] ** (3 / 2) / epsilon[-1])
        ReLe.append(k[-1] ** (1 / 2) * Le[-1] / nu) 
        ReLambda.append((20 * ReLe[-1] / 3) ** (1 / 2))  # Eq 6.64 Pope  # TODO: Ta errado como calcular corretamente o ReLambda
        L11_Le.append(calc_L11(ReLambda[-1], ReLe[-1]))
        L11.append(L11_Le[-1] * Le[-1])
        ReT.append(ul[-1] * L11[-1] / nu)  # turbulence Reynolds number

        # Boxal
        We_pipe.append(
            sol["pbe_sol"].cp.rho
            * U_**2
            * sol["pbe_sol"].D ** 3
            / sol["pbe_sol"].dp.sigma
        )
        Re_pipe.append(
            sol["pbe_sol"].cp.rho * U_ * sol["pbe_sol"].D / sol["pbe_sol"].cp.mu
        )

    data = DataFrame(
        {
            "N_emul": N_emul,
            "N_esc": N_esc,
            "Marco": marco,
            "Red D43": D43_r,
            "ANM_d Dv99 [um]": d99d,
            "ANM_d Dv50 [um]": d50d,
            "D90d_a": d90da,  # d95 depois - antes
            "epsilon [m²/s³]": epsilon,
            "We": we,
            "We_pipe": We_pipe,
            "Re": re,
            "Re_pipe": Re_pipe,
            "H2O [%]": ch2o,
            "ANM [%]": anm,
            "dP-ANM [Pa]": dP,
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
        join(plots_out, "Exp_analisys"),
        save=False,
    )
    d_max_kolmogorov_length(
        data,
        "dmax_x_eta_H2O",
        join(plots_out, "Exp_analisys"),
        save=False,
        style="H2O [%]",
    )
    d_We_Re(data, "d_We_Re", join(plots_out, "Exp_analisys"), save=False)

    with ExcelWriter(join(data_out, "exp_analisys.xlsx")) as writer:
        data.to_excel(writer, sheet_name="exp_analisys", index=False)
        data.transpose().to_excel(writer, sheet_name="exp_analisys_T")
