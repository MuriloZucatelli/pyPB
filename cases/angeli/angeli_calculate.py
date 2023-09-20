from numpy import genfromtxt, abs, array, pi
from scipy.optimize import minimize, differential_evolution
from pbe.app.angeli_class import AngeliSolution
import time
import pickle
import os


class angeli_experiment:
    def __init__(self, U, phi, d32, theta=1200.):
        self.phi = phi          # volumetric concentration
        self.theta = theta      # residence time [s]
        self.U = U              # [m/s]
        self.d32 = d32          # [m]


class angeli_numerical:
    """Calculate numerical solution for each model parameters (C) and 
    experimental condition
    """
    def __init__(self, C, experiment: angeli_experiment) -> None:
        v0s = array([0.5, 1.5]) * pi / 6 * experiment.d32**3
        mp = array(C)
        mp[3] *= 1e11
        self.pbe_solutions = [
            AngeliSolution(
                M=50, v0=v0, U=experiment.U, phi=experiment.phi, theta=experiment.theta,
                model_parameters=mp) for v0 in v0s]
        [print(v0) for v0 in v0s]


mm_to_m = 0.001
experiments = []
Us = [1.1, 1.25, 1.3, 1.5, 1.7]             # Fluid flow velocity            [m/s]
alphas = [0.091, 0.05, 0.077, 0.034, 0.05]  # Dispersed phase concentration  [1]
ds = [1.29, 0.97, 1.25, 0.73, 0.69]         # Sauter mean diameter SMD d32   [mm]
for i in range(5):
    experiments.append(angeli_experiment(Us[i], alphas[i], ds[i] * mm_to_m))

# original C&T parameters where multiplying (Nstar**3 * D**2)
# and epsilon = 0.407 * NStar**3 * D**2
# so
s = 0.407
c0 = [0.4 * s**(-1./3.), 0.08 / s**(-2./3.), 2.8 * s**(-1./3.), 1.83 * s]  # CT original constants

simulacoes = []
for i in range(5):
    simulacoes.append(angeli_numerical(c0, experiments[i]))


dir = os.path.dirname(__file__)
with open(os.path.join(dir, 'angeli_results.pickle'), 'wb') as f:
    pickle.dump(simulacoes, f)
