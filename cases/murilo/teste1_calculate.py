from numpy import arange, genfromtxt, abs, array, pi, set_printoptions
from scipy.optimize import minimize, differential_evolution
import sys
import os.path as path

dir = path.dirname(__file__)
if __name__ == "__main__":
    sys.path.append(path.abspath(path.join(dir, "..\\..")))
from pbe.app.dtg_class import DTGSolution, import_flow_DSD, DTG_experiment, get_location
import time
import pickle

set_printoptions(precision=4)

"""
    Resolve o balanço populacional para o dados do DTG

    Returns:

"""

# Pasta contendo as distribuições de tamanho de gotas
pasta = r"OneDrive\NEMOG Murilo\00 Circuito DTG\4. Compilado\LP_PB"
experiments = import_flow_DSD(get_location(pasta))


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
