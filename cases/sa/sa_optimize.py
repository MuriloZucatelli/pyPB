from numpy import genfromtxt, abs, array, pi, set_printoptions
from scipy.optimize import minimize, differential_evolution
import sys
import os.path as path

dir = path.dirname(__file__)
if __name__ == "__main__":
    sys.path.append(path.abspath(path.join(dir, "..\\..")))
from pbe.app.simmons_azzopardi_class_copy import SASolution
import time
import pickle

set_printoptions(precision=4)


class sa_experiment:
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


def sa_solve(C, experiment: sa_experiment, v0s):
    mp = array(C)
    mp[3] *= 1e13
    pbe_solutions = [
        SASolution(
            M=20,
            v0=v0,
            U=experiment.U,
            phi=experiment.phi,
            theta=experiment.theta,
            model_parameters=mp,
        )
        for v0 in v0s(experiment.d32)
    ]
    return pbe_solutions


def error_function(C, e: sa_experiment, v0s):
    pbe_sol = sa_solve(C, e, v0s)
    error = sum([abs(s.moc.d32 - e.d32) / e.d32 for s in pbe_sol])
    print(
        "constants {0}\t error {1:.2f}%\t ".format(C, 100 * error),
        "numer_d [{0:.2e}, {1:.2e}]\t exper_d {2:.2e}".format(
            pbe_sol[0].moc.d32, pbe_sol[1].moc.d32, e.d32
        ),
    )
    return error


mm_to_m = 0.001
experiments = []
experiments.append(sa_experiment())
# original C&T parameters where multiplying (Nstar**3 * D**2)
# and epsilon = 0.407 * NStar**3 * D**2
s = 0.407
c0 = [
    0.4 * s ** (-1.0 / 3.0),
    0.08 / s ** (-2.0 / 3.0),
    2.8 * s ** (-1.0 / 3.0),
    1.83 * s,
]  # CT original constants

# c0 = [0.75, 0.22, 5., 4.3]
# c0 optimo
c0 = [0.537, 0.1154, 3.778, 0.745]  # Error 0.02%
c0 = [0.3372, 0.1098, 1.3237, 0.262]  # Error 0.000%  Nelder-Mead
c0 = [0.4369, 0.1165, 2.9701, 0.7528]  # Error 0.01%  Powell

c0 = [0.537, 0.1154, 0.778, 0.745]  # Error 22%
v0s = lambda exp_d32: array([0.5, 1.5]) * pi / 6 * exp_d32**3
results = []
for e in experiments:
    print(e)
    res = dict()
    Copt = minimize(
        lambda c: error_function(c, e, v0s),
        c0,
        method="Powell",
        bounds=[(0, None)] * 4,
        options={
            "disp": True,
            "ftol": 1e-8,  # default: 2.2e-09
            "gtol": 1e-6,  # 1e-8 default is 1e-5
            "maxls": 200,  # Max line search steps (per iteration). Default is 20.
            "maxiter": 1000,
            # "eps": 0.1, # eps: default 1e-8
        },
    )  # default: 15000
    # GLOBAL OPTIMIZATION
    # b = [(0, 1), (0, 1), (0, 10), (0, 10)]
    # Copt = differential_evolution(lambda c: error_function(c, e),
    # bounds=b, maxiter=10, popsize=1)
    res["best_fit"] = Copt.x
    res["setup"] = {"U": e.U, "phi": e.phi, "d32": e.d32, "theta": e.theta}
    res["sim_d32"] = {
        i: sa_solve(Copt.x, e, v0s)[i].moc.d32 for i, _ in enumerate(v0s(e.d32))
    }
    res["v0s"] = {i: v0 for i, v0 in enumerate(v0s(e.d32))}
    res["error"] = error_function(Copt.x, e, v0s)
    results.append(res)

 # "Nelder-Mead" best than L-BFGS-B 
 # Powell works very well but is slow
with open(path.join(dir, "sa_optimization_results.pickle"), "wb") as f:
    pickle.dump(results, f)
