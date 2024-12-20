"""
    Função para resolver utilizando o modelo do mitre apenas de Coalescencia

    Returns:
        murilo_MitreCoalescence_{date}.plicke:
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
pasta_out = "solutions\\MitreCoalescence"

VALVULA = "MV02"
X_COMPARE = ["E_Choke"]  # ["E_ANM", "E_FlowLine"]

testes_all = {
    88: {2},
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
    128: {3},
    129: {2, 3},
}

# Args são os parâmetros para resolver a PBE
# testes, objetive, C0, Cname, model, ts
# varsigma: number of droplets from breakup
args = [
    (
        testes_all,
        "Modelo de coalescencia do Mitre",
        [1.0 / 1e2],
        ["Cc"],  # ordem das constantes
        {
            "breakup": None,
            "coalescence": "mitre_rigid_interface",
            "DDSD": None,
            "varsigma": 0.0,  # S1: 26.0 ± 0.9. S3: 32.7 ± 16.8
        },
        100,
        5,
        VALVULA,
    ),
    (
        testes_all,
        "Modelos do Mitre com frequencia não constante",
        [0.90 / 1e2, 1 / 1e2],
        ["Cc", "Ce"],  # ordem das constantes
        {
            "breakup": None,
            "coalescence": "mitre_partmobile_interface",
            "DDSD": None,
            "varsigma": 0.0,
        },
        100,
        5,
        VALVULA,
    ),
    (
        testes_all,
        "Modelos do Mitre com frequencia de quebra constante",
        [1.88 / 1e2],  # S3: [1.88 ± 0.06, 1.08 ± 0.07] . 10e-2
        ["Cc"],  # ordem das constantes
        {
            "breakup": None,
            "coalescence": "mitre_rigid_interface",
            "DDSD": None,
            "varsigma": 0.0,
        },
        100,
        5,
        VALVULA,
    ),
]

mitre_CEM = (
    [0.90 / 1e2, 1 / 1e2],
    ["Cc", "Ce"],  # ordem das constantes
    {
        "breakup": None,
        "coalescence": "mitre_partmobile_interface",
        "DDSD": None,
        "varsigma": 0.0,
    },
)
# Teste de convergencia de malha.
args2 = [
    (testes_all, "Mitre_CEM 5ts", *mitre_CEM, 5, 5, VALVULA),  # 3
    (testes_all, "Mitre_CEM 10ts", *mitre_CEM, 10, 5, VALVULA),
    (testes_all, "Mitre_CEM 20ts", *mitre_CEM, 20, 5, VALVULA),
    (testes_all, "Mitre_CEM 30ts", *mitre_CEM, 30, 5, VALVULA),
    (testes_all, "Mitre_CEM 50ts", *mitre_CEM, 50, 5, VALVULA),
    (testes_all, "Mitre_CEM 70ts", *mitre_CEM, 70, 5, VALVULA),
    (testes_all, "Mitre_CEM 100ts", *mitre_CEM, 100, 5, VALVULA),
    (testes_all, "Mitre_CEM 150ts", *mitre_CEM, 150, 5, VALVULA),  # 10
    (testes_all, "Mitre_CEM 0.1D", *mitre_CEM, 100, 0.1, VALVULA),  # 11
    (testes_all, "Mitre_CEM 0.5D", *mitre_CEM, 100, 0.5, VALVULA),
    (testes_all, "Mitre_CEM 1D", *mitre_CEM, 100, 1, VALVULA),
    (testes_all, "Mitre_CEM 2D", *mitre_CEM, 100, 2, VALVULA),
    (testes_all, "Mitre_CEM 2.4D", *mitre_CEM, 100, 2.4, VALVULA),  # 15
    (testes_all, "Mitre_CEM 3D", *mitre_CEM, 100, 3, VALVULA),
    (testes_all, "Mitre_CEM 4D", *mitre_CEM, 100, 4, VALVULA),
    (testes_all, "Mitre_CEM 5D", *mitre_CEM, 100, 5, VALVULA),
    (testes_all, "Mitre_CEM 6D", *mitre_CEM, 100, 6, VALVULA),
    (testes_all, "Mitre_CEM 7D", *mitre_CEM, 100, 7, VALVULA),  # 20
    (testes_all, "Mitre_CEM 8D", *mitre_CEM, 100, 8, VALVULA),  # 21
    (testes_all, "Mitre_CEM 10D", *mitre_CEM, 100, 10, VALVULA),  # 22
    (testes_all, "Mitre_CEM 15D", *mitre_CEM, 100, 15, VALVULA),  # 23
]

args = [*args, *args2]

sols1 = {
    f"sol_Mitre_{VALVULA}_C_CCE_S1": 0,
    f"sol_Mitre_{VALVULA}_C_CEM_S1": 1,
    f"sol_Mitre_{VALVULA}_C_CCE_S3": 2,
}

sols2 = {
    f"sol_Mitre_{VALVULA}_C_CEM_5ts": 3,
    f"sol_Mitre_{VALVULA}_C_CEM_10ts": 4,
    f"sol_Mitre_{VALVULA}_C_CEM_20ts": 5,
    f"sol_Mitre_{VALVULA}_C_CEM_40ts": 6,
    f"sol_Mitre_{VALVULA}_C_CEM_50ts": 7,
    f"sol_Mitre_{VALVULA}_C_CEM_70ts": 8,
    f"sol_Mitre_{VALVULA}_C_CEM_100ts": 9,
    f"sol_Mitre_{VALVULA}_C_CEM_150ts": 10,
    f"sol_Mitre_{VALVULA}_C_CEM_0.1D": 11,
    f"sol_Mitre_{VALVULA}_C_CEM_0.5D": 12,
    f"sol_Mitre_{VALVULA}_C_CEM_1D": 13,
    f"sol_Mitre_{VALVULA}_C_CEM_2D": 14,
    f"sol_Mitre_{VALVULA}_C_CEM_2.4D": 15,
    f"sol_Mitre_{VALVULA}_C_CEM_3D": 16,
    f"sol_Mitre_{VALVULA}_C_CEM_4D": 17,
    f"sol_Mitre_{VALVULA}_C_CEM_5D": 18,
    f"sol_Mitre_{VALVULA}_C_CEM_6D": 19,
    f"sol_Mitre_{VALVULA}_C_CEM_7D": 20,
    f"sol_Mitre_{VALVULA}_C_CEM_8D": 21,
    f"sol_Mitre_{VALVULA}_C_CEM_10D": 22,
    f"sol_Mitre_{VALVULA}_C_CEM_15D": 23}

sols = {**sols1, **sols2}

experiments = Import_flow_DSD2(get_location(pasta), teste=testes_all)


# Seleciona local de avaliação
experiments.select_DTG(X=X_COMPARE)
# Como obter apenas uma DTG:
# experiments.get_DTG(teste=90, marco=3, ID=3)
# ID vem de:
# ID = experiments.compares['E_ANM'][0] ou experiments.compares['E_ANM'][1]
# [0]: antes, [1]: depois


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
        dp_name=sol["dp_name"],
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
