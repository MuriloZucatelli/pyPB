"""
    Função para resolver utilizando os modelos C&T
    com utilização de otimização numérica SHGO
    SOGO2: Single Objective Global Optimization

    Returns:
        C&T_models_{date}.plicke: _description_
"""

from numpy import arange, abs, array, pi, set_printoptions, diff

# from numpy import genfromtxt, zeros_like
# from scipy.optimize import minimize, differential_evolution
from sys import path as sph
from os.path import join, abspath, dirname
from os import getcwd
from pandas import read_excel
from scipy.optimize import minimize, shgo, differential_evolution, dual_annealing
from multiprocessing import Pool
from datetime import date
import pickle
import time


dir = dirname(__file__)
if __name__ == "__main__":
    sph.append(abspath(join(dir, "..\\..\\..\\2. APP/")))
    sph.append(abspath(join(dir, "..\\..")))
from pbe.app.dtg_class import (
    DTGSolution,
    Import_flow_DSD2,
    DTG_experiment,
    get_location,
)

set_printoptions(precision=4)


"""
    Resolve o balanço populacional para o dados do DTG

    Returns:

"""


data = date.today().strftime("%d-%m-%Y")
# Importa dados
pasta = abspath(join(dir, "..\\..\\..", r"6. Compilado\\LP_PB_completo"))
testes_all = {
    88: {2, 3},
    90: {2, 3},
    91: {2, 3},
    92: {2, 3},
    93: {2, 3},
    94: {2, 3, 5, 6},
    95: {1, 2, 3},
    96: {1, 2, 3},
    97: {1, 2},
    98: {1, 2, 3},
    99: {1, 2, 3},
    100: {2, 3},
    101: {2, 3},
    102: {2, 3},
    103: {1, 2, 3},
    104: {2, 3},
    105: {2, 3},
    106: {2, 3},
    107: {2, 3},
    108: {2, 3},
    109: {2, 3},
    110: {2, 3},
    111: {2, 3},
    112: {2, 3},
    113: {2, 3},
    114: {2, 3},
    115: {2, 3},
    116: {2, 3},
    117: {2, 3},
    118: {2, 3},
    119: {2, 3},
    120: {2, 3},
    121: {2, 3},
    122: {2, 3},
    123: {2, 3},
    124: {2, 3},
    125: {2, 3},
    126: {2, 3},
    127: {2, 3},
    128: {2, 3},
}

experiments = Import_flow_DSD2(get_location(pasta), teste=testes_all)


# Seleciona local de avaliação
# experiments.select_DTG(X=["E_ANM", "E_FlowLine"])
experiments.select_DTG(X=["E_ANM"])

# Como obter apenas uma DTG:
# ID vem de:
# experiments.compares['E_ANM'][0] ou [1], 0: antes, [1]: depois
# experiments.get_DTG(teste=88, marco=1, ID=3)


# Importa propriedades
experiments.get_prop(dir)


# Define todos os dados experimentais em um objeto só
experiments.preparaDados()


# Obtem as classe do Bettersizer e define a malha
# Create mesh
x = read_excel(dir + "/classes.xlsx")
# diameter is in micrometer and volume is in mm³
d, v = x["d"].to_numpy(), x["v [mm³]"].to_numpy()
dxi = diff(v)
xi = v[:-1] + diff(v) / 2
# Remover classes zeros
sl = slice(8, -12)
dxi, xi = dxi[sl] / 1e9, xi[sl] / 1e9  # m³
xi_d = (6 / pi * xi) ** (1.0 / 3)  # metro
M = len(xi)


# Reduz as classes experimentais também
experiments.reduce_DTG(M, sl)


# Define a função de cálculo
def PB_solve(M, xi, dxi, C, sol, data, ts=50):
    pbe_solutions = DTGSolution(
        M=M,
        xi=xi,
        dxi=dxi,
        exp=sol["exp"],
        timesteps=ts,
        data=data,
        IDs=sol["compares"],
        marco=sol["marco"],
        model_parameters=array(C),
    )

    return pbe_solutions


def error_function(C, M, xi, dxi, sol, experiments, d43, dsd=None):
    pbe_sol = PB_solve(M, xi, dxi, C, sol, experiments)
    moc_d43 = pbe_sol.moc.d43
    error = abs(moc_d43 - d43) / d43
    if dsd is not None:
        error = sum(abs(pbe_sol.N2Fv - dsd / 100))
        error = sum((pbe_sol.N2Fv - dsd / 100) ** 2) / sum((dsd / 100) ** 2)
    return error


def global_sum_error(C, Ne: int, arguments) -> float:
    erro = [0] * Ne
    for i in range(Ne):
        erro[i] = error_function(C, *arguments[i])
    psi = sum(erro)
    print(
        "constants {0}\t error {1:.2f}%\t ".format(C, 100 * psi),
        # "sim_d43 {0:.2e}\t exper_d43 {1:.2e}".format(moc_d43, d43),
    )
    return psi


def main():

    def prepare_opt(testes, objetivo):
        Ne = sum([len(testes[i]) for i in testes])
        arguments = []
        i = 0
        sol = dict()
        sol["experiments"] = experiments
        sol["testes"] = testes
        sol["objective"] = objetivo
        for N in testes:
            IDs = experiments.compares["E_ANM"]
            for marco in testes[N]:
                dtg = experiments.get_DTG(
                    teste=N,
                    marco=marco,
                    ID=IDs[1],
                )
                d43 = experiments.get_D_caracteristico(dtg, xi_d, dff=[4, 3])
                exp = experiments.dados.loc[
                    (experiments.dados["Marco"] == marco)
                    * (experiments.dados["N_escoam"] == N)
                ].squeeze()
                sol[i] = {
                    "N_escoam": N,
                    "marco": marco,
                    "compares": IDs,
                    "exp": exp,
                    "M": M,
                }
                arguments.append(
                    (M, xi, dxi, sol[i], experiments, d43, dtg["freq_v"]),
                )
                i += 1
        return sol, Ne, arguments

    # Roda a simulação
    def run_opt(testes, bounds, objetivo=None, ts=50):

        n = 2000
        opt = dict()
        sol, Ne, arguments = prepare_opt(testes, objetivo)
        method = "sobol"
        print("Ready to start optmization")

        Copt = shgo(
            global_sum_error,
            bounds=bounds,
            args=(Ne, arguments),
            n=n,
            iters=5,
            sampling_method=method,
            options={"ftol": 0.000001},
            workers=18,
        )
        print(Copt)
        opt["best_fit"] = Copt.x
        if Copt.success:
            opt["local_solutions"] = Copt.xl
        opt["error"] = Copt.fun
        opt["bounds"] = bounds
        opt["n_shgo"] = n
        for i in range(Ne):
            sol[i]["best_pbe_sol"] = PB_solve(
                M, xi, dxi, Copt.x, sol[i], experiments, ts
            )
            sol[i]["opt"] = opt
            sol[i]["opt_flag"] = Copt.success
        return sol

    # Roda a simulação
    def run_opt2(testes, bounds, objetivo=None, ts=50):

        popsize = 400
        opt = dict()
        sol, Ne, arguments = prepare_opt(testes, objetivo)
        #method = "sobol"
        print("Ready to start optmization")

        Copt = differential_evolution(
            global_sum_error,
            bounds=bounds,
            args=(Ne, arguments),
            popsize=popsize,
            #sampling_method=method,
            workers=-1,
        )
        print(Copt)
        opt["best_fit"] = Copt.x
        if Copt.success:
            #opt["local_solutions"] = Copt.xl
            pass
        opt["error"] = Copt.fun
        opt["bounds"] = bounds
        opt["popsize"] = popsize
        for i in range(Ne):
            sol[i]["best_pbe_sol"] = PB_solve(
                M, xi, dxi, Copt.x, sol[i], experiments, ts
            )
            sol[i]["opt"] = opt
            sol[i]["opt_flag"] = Copt.success
        return sol

    testes = {
        91: {2},
        93: {2, 3},
        95: {2},
        96: {2, 3},
        99: {1},
        100: {3},
        101: {2, 3},
        102: {2},
        103: {2},
        104: {2},
        105: {3},
        106: {2},
        107: {3},
        108: {3},
        110: {3},
        111: {3},
        112: {3},
        115: {2, 3},
        116: {2},
        118: {3},
        119: {3},
        120: {2},
        121: {3},
        122: {3},
        123: {2},
        125: {3},
        128: {2},
    }

    bounds = [
        (1e-5 * 0.00481, 1e6 * 0.00481),
        (1e-6 * 0.08, 1e6 * 0.08),
        (1e-5 * 2.8e-6, 1e6 * 2.8e-6),
        (1e-6 * 1.83e9, 1e3 * 1.83e9),
    ]
    #bounds = [(1e-12, 1e12)] * 4

    t1 = time.time()
    args = [
        (
            testes,
            bounds,
            "Teste com erro relativo soma de testes selecionados 29-05-2024",
        )
    ]

    runs = {"murilo_optglobal_difevolution_erro_soma": 0}

    results = [0] * len(runs)
    i = 0
    for arg in args:
        print("Starting task...")
        #results[i] = run_opt(*arg)
        results[i] = run_opt2(*arg)
        i += 1

    for run in runs:
        id = runs[run]
        with open(
            join(dir, "solutions\\optimizations", f"opt_{run}_{data}.pickle"), "wb"
        ) as f:
            pickle.dump(results[id], f)
    t2 = time.time()
    print("time Pool: ", t2 - t1)


if __name__ == "__main__":
    main()
