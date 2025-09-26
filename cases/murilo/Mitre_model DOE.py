"""
    Função para resolver utilizando os modelos do mitre

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
pasta_out = "solutions\\MitreDOE2"

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
    129: {2, 3},
}

testes_selec = {
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
testes_selec = {
    122: {3},
}
# Analisar o 107_3 e o 104_2 121_3 122_3
# Args são os parâmetros para resolver a PBE
# testes, objetive, C0, Cname, model, ts
# varsigma: number of droplets from breakup


def nivel_DOE(m, zeta):
    Cc = 0.90 / 1e2 * (1 + m[0] * zeta[0] / 2)
    Ce = 1 / 1e2 * (1 + m[1] * zeta[1] / 2)
    Cb = 0.81 / 1e2 * (1 + m[2] * zeta[2] / 2)
    varsigma = 31.2 * (1 + m[3] * zeta[3] / 2)

    return (
        [Cc, Ce, Cb],
        ["Cc", "Ce", "Cb"],  # ordem das constantes
        {
            "breakup": "mitre_modified",
            "coalescence": "mitre_partmobile_interface",
            "DDSD": "mitre",
            "varsigma": varsigma,
        },
    )


# Testes padronizados
zeta = [52.2 / 100,
        50 / 100,
        75 / 100,
        30.4 / 100]  # Cc, Ce, Cb, Varsigma

args = [
    (testes_selec, "Mitre_CEM_1", *nivel_DOE([-1, -1, -1, -1], zeta), 100, 5, VALVULA),  # 1
    (testes_selec, "Mitre_CEM_2", *nivel_DOE([1, -1, -1, -1], zeta), 100, 5, VALVULA),
    (testes_selec, "Mitre_CEM_3", *nivel_DOE([-1, 1, -1, -1], zeta), 100, 5, VALVULA),
    (testes_selec, "Mitre_CEM_4", *nivel_DOE([1, 1, -1, -1], zeta), 100, 5, VALVULA),  # 4
    (testes_selec, "Mitre_CEM_5", *nivel_DOE([-1, -1, 1, -1], zeta), 100, 5, VALVULA),
    (testes_selec, "Mitre_CEM_6", *nivel_DOE([1, -1, 1, -1], zeta), 100, 5, VALVULA),
    (testes_selec, "Mitre_CEM_7", *nivel_DOE([-1, 1, 1, -1], zeta), 100, 5, VALVULA),
    (testes_selec, "Mitre_CEM_8", *nivel_DOE([1, 1, 1, -1], zeta), 100, 5, VALVULA),
    (testes_selec, "Mitre_CEM_9", *nivel_DOE([-1, -1, -1, 1], zeta), 100, 5, VALVULA),  # 9
    (testes_selec, "Mitre_CEM_10", *nivel_DOE([1, -1, -1, 1], zeta), 100, 5, VALVULA),
    (testes_selec, "Mitre_CEM_11", *nivel_DOE([-1, 1, -1, 1], zeta), 100, 5, VALVULA),
    (testes_selec, "Mitre_CEM_12", *nivel_DOE([1, 1, -1, 1], zeta), 100, 5, VALVULA),
    (testes_selec, "Mitre_CEM_13", *nivel_DOE([-1, -1, 1, 1], zeta), 100, 5, VALVULA),  # 13
    (testes_selec, "Mitre_CEM_14", *nivel_DOE([1, -1, 1, 1], zeta), 100, 5, VALVULA),
    (testes_selec, "Mitre_CEM_15", *nivel_DOE([-1, 1, 1, 1], zeta), 100, 5, VALVULA),
    (testes_selec, "Mitre_CEM_16", *nivel_DOE([1, 1, 1, 1], zeta), 100, 5, VALVULA),  # 16
    (testes_selec, "Mitre_CEM_17", *nivel_DOE([-2, 0, 0, 0], zeta), 100, 5, VALVULA),
    (testes_selec, "Mitre_CEM_18", *nivel_DOE([2, 0, 0, 0], zeta), 100, 5, VALVULA),
    (testes_selec, "Mitre_CEM_19", *nivel_DOE([0, -2, 0, 0], zeta), 100, 5, VALVULA),
    (testes_selec, "Mitre_CEM_20", *nivel_DOE([0, 2, 0, 0], zeta), 100, 5, VALVULA),
    (testes_selec, "Mitre_CEM_21", *nivel_DOE([0, 0, -2, 0], zeta), 100, 5, VALVULA),
    (testes_selec, "Mitre_CEM_22", *nivel_DOE([0, 0, 2, 0], zeta), 100, 5, VALVULA),
    (testes_selec, "Mitre_CEM_23", *nivel_DOE([0, 0, 0, -2], zeta), 100, 5, VALVULA),
    (testes_selec, "Mitre_CEM_24", *nivel_DOE([0, 0, 0, 2], zeta), 100, 5, VALVULA),
    (testes_selec, "Mitre_CEM_25", *nivel_DOE([0, 0, 0, 0], zeta), 100, 5, VALVULA),
]

sols = {
    f"sol_{VALVULA}_Mitre_CEM_1": 0,
    f"sol_{VALVULA}_Mitre_CEM_2": 1,
    f"sol_{VALVULA}_Mitre_CEM_3": 2,
    f"sol_{VALVULA}_Mitre_CEM_4": 3,
    f"sol_{VALVULA}_Mitre_CEM_5": 4,
    f"sol_{VALVULA}_Mitre_CEM_6": 5,
    f"sol_{VALVULA}_Mitre_CEM_7": 6,
    f"sol_{VALVULA}_Mitre_CEM_8": 7,
    f"sol_{VALVULA}_Mitre_CEM_9": 8,
    f"sol_{VALVULA}_Mitre_CEM_10": 9,
    f"sol_{VALVULA}_Mitre_CEM_11": 10,
    f"sol_{VALVULA}_Mitre_CEM_12": 11,
    f"sol_{VALVULA}_Mitre_CEM_13": 12,
    f"sol_{VALVULA}_Mitre_CEM_14": 13,
    f"sol_{VALVULA}_Mitre_CEM_15": 14,
    f"sol_{VALVULA}_Mitre_CEM_16": 15,
    f"sol_{VALVULA}_Mitre_CEM_17": 16,
    f"sol_{VALVULA}_Mitre_CEM_18": 17,
    f"sol_{VALVULA}_Mitre_CEM_19": 18,
    f"sol_{VALVULA}_Mitre_CEM_20": 19,
    f"sol_{VALVULA}_Mitre_CEM_21": 20,
    f"sol_{VALVULA}_Mitre_CEM_22": 21,
    f"sol_{VALVULA}_Mitre_CEM_23": 22,
    f"sol_{VALVULA}_Mitre_CEM_24": 23,
    f"sol_{VALVULA}_Mitre_CEM_25": 24,
}


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
            # print(sol[i]["pbe_sol"].moc.d43)
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
