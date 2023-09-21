import pickle
import time
from numpy import genfromtxt, abs, array, pi, set_printoptions
from scipy.optimize import minimize, differential_evolution
import sys
import os.path as path
dir = path.dirname(__file__)
if __name__ == '__main__':
    sys.path.append(path.abspath(path.join(dir, '..\..')))
    from pbe.app.angeli_class import AngeliSolution
else:
    from pbe.app.angeli_class import AngeliSolution
set_printoptions(precision=4)
# Tudo ok


with open(path.join(dir, 'angeli_optimization_results.pickle'), 'rb') as f:
    results = pickle.load(f)

print(results)