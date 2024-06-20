"""
    Função para resolver utilizando os modelos C&T
    com constantes originais e do Liao ou outro autor

    Returns:
        murilo_C&T_basico_{date}.plicke: _description_
"""

from numpy import arange, abs, array, pi, set_printoptions, diff, column_stack

# from numpy import genfromtxt, zeros_like
from sys import path as sph
from os.path import join, abspath, dirname
from os import getcwd
from pandas import read_excel
from multiprocessing import Pool
from datetime import date
import pickle
import time
from pymoo.core.problem import (
    ElementwiseEvaluationFunction,
    LoopedElementwiseEvaluation,
    Problem,
    StarmapParallelization,
    ElementwiseProblem,
)
from pymoo.optimize import minimize
from pymoo.algorithms.moo import nsga2
from pymoo.termination.default import DefaultMultiObjectiveTermination
from pymoo.core.callback import Callback
from matplotlib import pyplot as plt

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
experiments.select_DTG(X=["E_ANM", "E_FlowLine"])
experiments.select_DTG(X=["E_ANM"])
# Como obter apenas uma DTG:
# ID vem de:
# experiments.compares['E_ANM'][0] ou [1], 0: antes, [1]: depois
# ID = experiments.compares['E_ANM'][0] ou experiments.compares['E_ANM'][1]
# experiments.get_DTG(teste=90, marco=3, ID=3)


# Importa propriedades
experiments.get_prop(dir)

# Define os dados experimentais em um objeto so
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


# Roda a simulação
def run_sim(testes, objetive=None, C0=None, ts=50):
    i = 0
    sol = dict()
    sol["experiments"] = experiments
    sol["testes"] = testes
    sol["objective"] = objetive
    for N in testes:
        print("Teste número", N, " Marcos: ", list(testes[N]))
        IDs = experiments.compares["E_ANM"]
        for marco in testes[N]:
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
            sol[i]["pbe_sol"] = PB_solve(
                M, xi, dxi, C0, sol[i], experiments, ts=ts  # m³  # m³
            )
            print(sol[i]["pbe_sol"].moc.d43)
            i += 1
    return sol


def error_function(C, M, xi, dxi, sol, experiments, d43, dsd=None):
    pbe_sol = PB_solve(M, xi, dxi, C, sol, experiments)  # mm³ to m³  # mm³ to m³
    # moc_d43 = pbe_sol.moc.d43
    # error = abs(moc_d43 - d43) / d43
    if dsd is not None:
        # error = sum(abs(pbe_sol.N2Fv - dsd / 100))
        error = sum((pbe_sol.N2Fv - dsd / 100) ** 2) / sum((dsd / 100) ** 2)
    print(
        "constants {0}\t error {1:.2f}%\t ".format(C, 100 * error),
        # "sim_d43 {0:.2e}\t exper_d43 {1:.2e}".format(moc_d43, d43),
    )
    return error


class pop_balance_DTG(ElementwiseProblem):

    def __init__(
        self,
        n_var,
        n_obj,
        n_ieq_constr=0,
        n_eq_constr=0,
        xl=None,
        xu=None,
        testes=None,
        objetivo=None,
        vtype=None,
        vars=None,
        elementwise=True,
        elementwise_func=...,
        elementwise_runner=...,
        requires_kwargs=False,
        replace_nan_values_by=None,
        exclude_from_serialization=None,
        callback=None,
        strict=True,
        **kwargs,
    ):

        super().__init__(
            n_var=n_var,
            n_obj=n_obj,
            n_ieq_constr=n_ieq_constr,
            n_eq_constr=n_eq_constr,
            xl=xl,
            xu=xu,
            elementwise_runner=elementwise_runner,
            # vtype,
            # vars,
            # elementwise,
            # elementwise_func,
            # requires_kwargs,
            # replace_nan_values_by,
            # exclude_from_serialization,
            # callback,
            # strict,
            # **kwargs,
        )
        self.n_obj = n_obj
        self.arguments = []
        i = 0
        self.sol = dict()
        self.sol["experiments"] = experiments
        self.sol["testes"] = testes
        self.sol["objective"] = objetivo
        for N in testes:
            # print("Teste número", N, " Marcos: ", list(testes[N]))
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
                self.sol[i] = {
                    "N_escoam": N,
                    "marco": marco,
                    "compares": IDs,
                    "exp": exp,
                    "M": M,
                }
                self.arguments.append(
                    (M, dxi, xi, self.sol[i], experiments, d43, dtg["freq_v"]),
                )
                i += 1

        print("Ready to start optmization")

    def _evaluate(self, x, out, *args, **kwargs):
        listF = []
        for i in range(self.n_obj):
            listF.append(error_function(x, *self.arguments[i]))

        out["F"] = column_stack(listF)


def main():

    # Define os testes a serem feitos
    # testes = testes_all.copy()
    # testes.pop(88)

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

    # O numero de funções objetivos serão a somatória de todos os marcos de todos os testes
    n_obj = sum([len(testes[i]) for i in testes])
    n_var = 4
    # Constantes minimas
    Cl = array([1e-6 * 0.4, 1e-6 * 0.08, 1e-6 * 2.8e-6, 1e-6 * 1.83e9])
    # Constantes máximas
    Cu = array([1e6 * 0.4, 1e6 * 0.08, 1e6 * 2.8e-6, 1e4 * 1.83e9])

    n_proccess = 20
    pool = Pool(n_proccess)
    runner = StarmapParallelization(pool.starmap)

    problem = pop_balance_DTG(
        n_var=n_var,
        n_obj=n_obj,
        xl=Cl,
        xu=Cu,
        testes=testes,
        objetivo="CT_model otimização",
        elementwise_runner=runner,
    )

    algorithm = nsga2.NSGA2(pop_size=100)
    # algorithm = GA(pop_size=100)

    termination = DefaultMultiObjectiveTermination(
        xtol=1e-8,
        cvtol=1e-6,
        ftol=0.0025,
        period=30,
        n_max_gen=1000,
        n_max_evals=100000,
    )

    class MyCallback(Callback):

        def __init__(self) -> None:
            super().__init__()
            self.data["best"] = []

        def notify(self, algorithm):
            self.data["best"].append(algorithm.pop.get("F").min())

    res = minimize(
        problem,
        algorithm,
        termination=termination,
        seed=1,
        verbose=True,
        display=None,
        callback=MyCallback(),
        return_least_infeasible=False,
        save_history=True,
    )
    pool.close()
    with open(
        join(dir, "solutions\\optimizations", f"opt_pymoo_CT_{data}.pickle"), "wb"
    ) as f:
        pickle.dump(res, f)

    val = res.algorithm.callback.data["best"]
    plt.plot(arange(len(val)), val)
    plt.show()


if __name__ == "__main__":
    main()
