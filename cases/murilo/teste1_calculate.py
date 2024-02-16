from numpy import arange, abs, array, pi, set_printoptions, diff
#from numpy import genfromtxt, zeros_like
#from scipy.optimize import minimize, differential_evolution
from sys import path as sph
from os.path import join, abspath, dirname
from pandas import read_excel
#import time
#import pickle
#import matplotlib.pyplot as plt

dir = dirname(__file__)
if __name__ == "__main__":
    sph.append(abspath(join(dir, "..\\..\\..\\2. APP/")))
    sph.append(abspath(join(dir, "..\\..")))
from pbe.app.dtg_class import DTGSolution, import_flow_DSD, DTG_experiment, get_location


set_printoptions(precision=4)

"""
    Resolve o balanço populacional para o dados do DTG

    Returns:

"""


# Pasta contendo as distribuições de tamanho de gotas
pasta = r"OneDrive/NEMOG Murilo/00 Circuito DTG/6. Compilado/LP_PB"
experiments = import_flow_DSD(get_location(pasta), teste=[88], marco=[2])
experiments.select_DTG(X=["E_ANM", "E_FlowLine"])
experiments.select_DTG(X=["E_ANM"])
experiments.calc_DP_GV()
experiments.get_prop(dir)
# Como obter apenas uma DTG:
# ID vem de:
# experiments.compares['E_ANM'][0] ou [1], 0: antes, [1]: depois
# experiments.get_DTG(teste=88, ID=2)

# Obtem as classe do Bettersizer e define a malha
# Create mesh
x = read_excel(dir + "/classes.xlsx")
# diameter is in micrometer and volume is in mm³
d, v = x["d"].to_numpy(), x["v [mm³]"].to_numpy()
dxi = diff(v)
xi = v[:-1] + diff(v) / 2
# Remover classes zeros
sl = slice(8, -12)
dxi, xi = dxi[sl], xi[sl]
M = len(xi)


def teste1_solve(C, experiment: DTG_experiment):
    mp = array(C)
    mp[3] *= 1e13
    pbe_solutions = DTGSolution(
        M=20,
        U=experiment.U,
        phi=experiment.phi,
        theta=experiment.theta,
        model_parameters=mp,
    )

    return pbe_solutions


C0 = [0.3372, 0.1098, 1.3237, 0.262]
grids = [20]
alphas = [0.091]
Us = [1.1]
t = arange(0.0, 10.0, 0.001)
experiments: DTG_experiment = []
for i in range(len(grids)):
    experiments.append(DTG_experiment(Us[i], alphas[i], ds[i] * 1e-03))
pbe_sol = []
for e in experiments:
    pbe_sol.append(teste1_solve(C0, experiments))

