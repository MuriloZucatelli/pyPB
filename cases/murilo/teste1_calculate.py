from numpy import genfromtxt, abs, array, pi, set_printoptions
from scipy.optimize import minimize, differential_evolution
import sys
import os.path as path

dir = path.dirname(__file__)
if __name__ == "__main__":
    sys.path.append(path.abspath(path.join(dir, "..\\..")))
from pbe.app.dtg_class import DTGSolution
import time
import pickle

set_printoptions(precision=4)


class DTG_experiment:
    def __init__(self, d32=0.32e-03, phi=0.117, U=2.72, theta=600.0):
        self.phi = phi
        self.theta = theta
        self.U = U
        self.d32 = d32

    def __repr__(self) -> str:
        phi = "{:.{}f}".format(100 * self.phi, 2)
        d32 = "{:.{}e}".format(1000 * self.d32, 2)
        return (
            f"concentration {phi}%  |  flow velocity {self.U} m/s  |  "
            + f"residence time {self.theta} s  |  SMD diameter {d32} mm"
        )


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
experiments = []
experiments.append(DTG_experiment())
results = []
for e in experiments:
    pbe_sol = teste1_solve(C0, experiments)
    results.append(pbe_sol)