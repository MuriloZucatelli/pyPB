"""
    Função para pos processar os resultados do modelo C&T
    com utilização de otimização numérica SHGO

    Returns:
        plots
"""

from os.path import join, abspath, dirname
from os import getcwd
from sys import path as sph
import pickle
from numpy import abs, pi, set_printoptions
from pandas import DataFrame, concat

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
from utils.plot_SOGO_comparition import DTG_exp_x_sim, C_global_x_erro, C_global_x_erro_2x2
from pbe.setup.helpers import plt_config2

plt_config2(relative_fig_width=0.7)
set_printoptions(precision=4)

s_folder = "..\\solutions\\optimizations"
plots_out = join(dir, "..\\", "plots")
data_out = join(dir, "..\\", "tabelas")
optmization = [
    "murilo_optglobal_shgo_simplicial_10-05-2024.pickle",
    "murilo_optglobal_sobol_10-05-2024.pickle",
    "murilo_optglobal_shgo_simplicial_11-05-2024.pickle",
]

plots_made = {1: False, 2: True}

with open(join(dir, s_folder, optmization[0]), "rb") as f:
    ct_SHGO_Simplicial = pickle.load(f)
with open(join(dir, s_folder, optmization[1]), "rb") as f:
    ct_SHGO_Sobol = pickle.load(f)
with open(join(dir, s_folder, optmization[2]), "rb") as f:
    ct_SHGO_Simplicial_2 = pickle.load(f)


exp_SHGO_simp: Import_flow_DSD2 = ct_SHGO_Simplicial["experiments"]
exp_SHGO_sobol: Import_flow_DSD2 = ct_SHGO_Sobol["experiments"]
exp_SHGO_simp_2: Import_flow_DSD2 = ct_SHGO_Simplicial_2["experiments"]

# Dados não mutáveis:
exp_SHGO_simp.dados

# Observação, os testes e marcos otimizados podem sofrer alterações
# Não necessariamente serão os mesmos de experiments.marco
# Os testes e marcos rodados estão em ct_teste[i]['teste'] ou
# ct_teste[i]['marco']


# Primeiro tipo de plot
# DTG com melhor solução

if plots_made[1]:
    for i in ct_SHGO_Simplicial:
        if not isinstance(i, int):
            continue
        DTG_exp_x_sim(
            ct_SHGO_Simplicial[i],
            exp_SHGO_simp,
            i,
            "ctSHGO_Simplic",
            join(plots_out, "C&T_model"),
        )


if plots_made[1]:
    for i in ct_SHGO_Sobol:
        if not isinstance(i, int):
            continue
        DTG_exp_x_sim(
            ct_SHGO_Sobol[i],
            exp_SHGO_sobol,
            i,
            "ctSHGO_Sobol",
            join(plots_out, "C&T_model"),
            # save=True,
        )

if plots_made[1]:
    for i in ct_SHGO_Simplicial_2:
        if not isinstance(i, int):
            continue
        DTG_exp_x_sim(
            ct_SHGO_Simplicial_2[i],
            exp_SHGO_simp_2,
            i,
            "ctSHGO_Simplic_2",
            join(plots_out, "C&T_model"),
            # save=True,
        )


# 2° tipo de plot
# Constantes para todos os testes e todos os marcos com distinção
# apenas do otimizador
# Tem C locais sem o erro calculado


Error = []
C = []
N_emul = []
marco = []
opt_flag = []
method = []
model = []
Cs = ["C1", "C2", "C3", "C4"]
# runs = {"ct_SHGO_Simplicial": ct_SHGO_Simplicial, "ct_SHGO_Sobol": ct_SHGO_Sobol}
runs = {
    "ct_SHGO_Simplicial": ct_SHGO_Simplicial,
    "ct_SHGO_Sobol": ct_SHGO_Sobol,
    "ct_SHGO_Simplicial_2": ct_SHGO_Simplicial_2,
}
# models = {"ct_SHGO_Simplicial": 'C&T', "ct_SHGO_Sobol": 'C&T'}
models = {
    "ct_SHGO_Simplicial": "C&T",
    "ct_SHGO_Sobol": "C&T",
    "ct_SHGO_Simplicial_2": "C&T",
}
# methods = {"ct_SHGO_Simplicial": 'SHGO Simplicial', "ct_SHGO_Sobol": 'SHGO Sobol'}
methods = {
    "ct_SHGO_Simplicial": "SHGO Simplicial",
    "ct_SHGO_Sobol": "SHGO Sobol",
    "ct_SHGO_Simplicial_2": "SHGO Simplicial 2",
}

if plots_made[2]:
    for r in runs:
        run = runs[r]
        print(f"Getting {r} optmization")
        for i in run:  # Iterate over each entry
            if not isinstance(i, int):
                continue
            N_emul.append(int(run[i]["exp"]["N_emul"]))
            method.append(methods[r])
            model.append(models[r])
            marco.append(run[i]["marco"])
            opt_flag.append(run[i]["opt_flag"])
            C.append(run[i]["opt"]["best_fit"])
            Error.append(run[i]["opt"]["error"])

    C_global = DataFrame(C, columns=Cs)
    C_global = concat(
        [
            C_global,
            DataFrame(
                {
                    "Error": Error,
                    "N_emul": N_emul,
                    "opt_flag": opt_flag,
                    "method": method,
                    "model": model,
                    "marco": marco,
                }
            ),
        ],
        axis=1,
    )

    #C_global_x_erro(C_global, Cs, "C_x_erro_C&T_model_2", plots_out)
    C_global_x_erro_2x2(C_global, Cs, "C_x_erro_C&T_model_2", plots_out)

    best = C_global.groupby(by='method').get_group('SHGO Sobol')

    print(best[best['Error']< 0.08][['C1', 'C2', 'C2', 'C2', 'Error']].mean())