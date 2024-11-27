"""
    Função para resolver utilizando os modelos C&T
    com constantes originais e do Liao ou outro autor

    Returns:
        murilo_C&T_basico_{date}.plicke: _description_
"""

from numpy import arange, abs, array, pi, set_printoptions, diff

# from numpy import genfromtxt, zeros_like
from sys import path as sph
from os.path import join, abspath, dirname
from os import getcwd
from pandas import read_excel
from multiprocessing import Pool
from datetime import date
import pickle
import gzip
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


data = date.today().strftime("%d-%m-%Y")
# Importa dados
pasta = abspath(join(dir, "..\\..\\..", r"6. Compilado\\LP_PB2"))
pasta_out = "solutions\\basico"
VALVULA = "MV01"
X_COMPARE = ["E_ANM"]  # ["E_ANM", "E_FlowLine"]

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

#
# CONSTANTES DE MODELO
#
C0_original = [0.4, 0.08, 2.8e-6, 1.83e9]
C0_liao = [0.00481, 0.08, 2.8e-6, 1.83e9]
C_name = ["Cb", "Cepsilon", "Cc", "Ce"]  # ordem das constantes

#
# DEFINIÇÃO DOS MODELOS
#
Cmodel = {
    "breakup": "coulaloglou",
    "coalescence": "coulaloglou",
    "DDSD": "coulaloglou",
    "varsigma": 2.0,
}

#
# DEFINIÇÃO DOS PARÂMETROS DE SIMULAÇÃO
#
CT_liao = (C0_liao, C_name, Cmodel)
CT_original = (C0_original, C_name, Cmodel)

args = [
    (
        testes_all,
        "Calcular o epsilon para cada teste",
        [9.1609e-01, 3.6968e02, 1.3237e-03, 2.6200e9],
        C_name,
        Cmodel,
        10,
        5,
        VALVULA,
    ),
    (testes_all, "Constantes de C&T", *CT_original),
    (testes_all, "Constantes de C&T do Liao", *CT_liao),
    (
        testes_all,
        "Constantes de C&T otimizadas",
        [3.6, 1.4e3, 1.8e-4, 1.4e10],
        C_name,
        Cmodel,
    ),  # 22-05-2024
    (testes_all, "C&T Liao 5 ts", *CT_liao, 5),  # 4
    (testes_all, "C&T Liao 10 ts", *CT_liao, 10),
    (testes_all, "C&T Liao 20 ts", *CT_liao, 20),
    (testes_all, "C&T Liao 40 ts", *CT_liao, 40),
    (testes_all, "C&T Liao 50 ts", *CT_liao, 50),
    (testes_all, "C&T Liao 60 ts", *CT_liao, 60),
    (testes_all, "C&T Liao 100 ts", *CT_liao, 100),
    (testes_all, "C&T Liao 150 ts", *CT_liao, 150),  # 11
    (testes_all, "C&T Liao_0.1D", *CT_liao, 100, 0.1),  # 12
    (testes_all, "C&T Liao_0.5D", *CT_liao, 100, 0.5),
    (testes_all, "C&T Liao_1D", *CT_liao, 100, 1),
    (testes_all, "C&T Liao_2D", *CT_liao, 100, 2),
    (testes_all, "C&T Liao_2.4D", *CT_liao, 100, 2.4),  # 16
    (testes_all, "C&T Liao_3D", *CT_liao, 100, 3),
    (testes_all, "C&T Liao_4D", *CT_liao, 100, 4),
    (testes_all, "C&T Liao_5D", *CT_liao, 100, 5),
    (testes_all, "C&T Liao_6D", *CT_liao, 100, 6),
    (testes_all, "C&T Liao_7D", *CT_liao, 100, 7),  # 21
    (testes_all, "C&T Liao_8D", *CT_liao, 100, 8),  # 22
]

# 22-05-2024: Otimização de Sobol, média dos erros menores q 8%

sols = {
    F"sol_C&T_basico_{VALVULA}": 0,
    #"sol_C&T_original": 1,
    #"sol_C&T_liao": 2,
    #"sol_C&T_murilo": 3,
}
# sols = {
#     "sol_C&T_liao_5ts": 4,
#     "sol_C&T_liao_10ts": 5,
#     "sol_C&T_liao_20ts": 6,
#     "sol_C&T_liao_40ts": 7,
#     "sol_C&T_liao_50ts": 8,
#     "sol_C&T_liao_60ts": 9,
#     "sol_C&T_liao_100ts": 10,
#     "sol_C&T_liao_150ts": 11,
#     "sol_C&T_liao_0.1D": 12,
#     "sol_C&T_liao_0.5D": 13,
#     "sol_C&T_liao_1D": 14,
#     "sol_C&T_liao_2D": 15,
#     "sol_C&T_liao_2.4D": 16,
#     "sol_C&T_liao_3D": 17,
#     "sol_C&T_liao_4D": 18,
#     "sol_C&T_liao_5D": 19,
#     "sol_C&T_liao_6D": 20,
#     "sol_C&T_liao_7D": 21,
#     "sol_C&T_liao_8D": 22,
# }


experiments = Import_flow_DSD2(get_location(pasta), teste=testes_all)


# Seleciona local de avaliação
experiments.select_DTG(X=X_COMPARE)

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
x = read_excel(join(dir, "pb_data\\classes.xlsx"))
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
def PB_solve(M, xi, dxi, mp, sol, data, model, ts, fator=5):
    pbe_solutions = DTGSolution(
        M=M,
        xi=xi,
        dxi=dxi,
        exp=sol["exp"],
        timesteps=ts,
        data=data,
        IDs=sol["compares"],
        marco=sol["marco"],
        model_parameters=mp,
        breakupmodel=model["breakup"],
        coalescencemodel=model["coalescence"],
        DDSDmodel=model["DDSD"],
        varsigma=model["varsigma"],
        fator=fator,
        dp_name=sol["dp_name"]
    )

    return pbe_solutions


# Roda a simulação
def run_sim(testes, objetive=None, C0=None, Cname=None, model=None, ts=100, fator=5, dp_name="MV01"):
    i = 0
    sol = dict()
    sol["experiments"] = experiments
    sol["testes"] = testes
    sol["objective"] = objetive
    sol["model"] = model
    for N in testes:
        print("Teste número", N, " Marcos: ", list(testes[N]))
        IDs = experiments.compares[X_COMPARE[0]]
        for marco in testes[N]:
            exp = experiments.dados.loc[
                (experiments.dados["Marco"] == marco)
                * (experiments.dados["N_escoam"] == N)
            ].squeeze()
            sol[i] = {
                "N_escoam": N,
                "marco": marco,
                "compares": IDs,
                "C0": C0,
                "exp": exp,
                "M": M,
                "mp": Cname,
                "dp_name": dp_name,
            }
            mp = {i: j for i, j in zip(Cname, C0)}
            sol[i]["pbe_sol"] = PB_solve(
                M,
                xi,
                dxi,
                mp,
                sol[i],
                experiments,
                model,
                ts,
                fator=fator,
            )
            print(sol[i]["pbe_sol"].moc.d43)
            i += 1
    return sol


t1 = time.time()

if __name__ == "__main__":
    print("Starting task...")
    # Redução
    args_sol = [args[sols[x]] for x in sols]
    with Pool(10) as pool:
        # perform calculations
        results = pool.starmap(run_sim, args_sol)

    i = 0
    for sol in sols:
        id = sols[sol]
        with gzip.open(join(dir, pasta_out, f"{sol}_{data}.pickle"), "wb") as f:
            pickle.dump(results[i], f)
            i += 1
t2 = time.time()
print("time Pool: ", t2 - t1)
