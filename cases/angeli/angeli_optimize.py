from numpy import genfromtxt, abs, array, pi, set_printoptions
from scipy.optimize import minimize, differential_evolution
from pbe.app.angeli_class import AngeliSolution
import time
import os
import pickle
set_printoptions(precision=4)


class angeli_experiment:
    def __init__(self, U, phi, d32, theta=1200.):
        self.phi = phi
        self.theta = theta
        self.U = U
        self.d32 = d32

    def __repr__(self) -> str:
        phi = "{:.{}f}".format(100*self.phi, 2)
        d32 = "{:.{}e}".format(1000*self.d32, 2)
        return f"concentration {phi}%  |  flow velocity {self.U} m/s  |  " + \
               f"residence time {self.theta} s  |  SMD diameter {d32} mm"





def error_function(C, experiment: angeli_experiment):
    v0s = array([0.5, 1.5]) * pi / 6 * experiment.d32**3  # Inicial volume
    mp = array(C)
    mp[3] *= 1e11
    pbe_solutions = [
        AngeliSolution(
            M=30, v0=v0, U=experiment.U, phi=experiment.phi, theta=experiment.theta,
            model_parameters=mp)
        for v0 in v0s]
    error = sum(
        [abs(s.moc.d32 - experiment.d32) / experiment.d32
            for s in pbe_solutions])
    print('constants {0}\t error {1:.2f}%\t '.format(C, 100*error),
          'numer_d [{0:.2e}, {1:.2e}]\t exper_d {2:.2e}'.format(pbe_solutions[0].moc.d32,
                                                                         pbe_solutions[1].moc.d32,
                                                                         experiment.d32))
    return error


mm_to_m = 0.001

experiments = []
Us = [1.1, 1.25, 1.3, 1.5, 1.7]
alphas = [0.091, 0.05, 0.077, 0.034, 0.05]
ds = [1.29, 0.97, 1.25, 0.73, 0.69]
for i in range(5):
    experiments.append(angeli_experiment(Us[i], alphas[i], ds[i] * 1e-03))

# original C&T parameters where multiplying (Nstar**3 * D**2)
# and epsilon = 0.407 * NStar**3 * D**2
# so
s = 0.407
c0 = [0.4 * s**(-1./3.), 0.08 / s**(-2./3.), 2.8 * s**(-1./3.), 1.83 * s]  # CT original constants

results = []
for e in experiments:
    print(e)
    res = dict()

    Copt = minimize(
        lambda c: error_function(c, e), c0,
        method='L-BFGS-B', bounds=[(0, None)] * 4,
        options={'disp': False, 'ftol': 0.01, 'maxiter': 2000})
    res['best_fit'] = Copt.x
    res['setup'] = {'U': e.U, 'phi': e.phi, 'd32': e.d32, 'theta': e.theta}
    results.append(res)

dir = os.path.dirname(__file__)
with open(os.path.join(dir, 'angeli_optimization_results.pickle'), 'wb') as f:
    pickle.dump(results, f)
