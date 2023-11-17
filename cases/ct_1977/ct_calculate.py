import sys
import os.path as path

dir = path.dirname(__file__)
if __name__ == "__main__":
    sys.path.append(path.abspath(path.join(dir, "..\\..")))

from pbe.app.ct_class import CTSolution
from numpy import linspace, array
import pickle
import os

"""
This script calculated solutions to PBE problem that are necessary to reproduce
figure 3 from CT publication.
"""
concentrations = [10]  # , 10, 15]
Ns = linspace(3, 6, 10)  # rps
#Ns = array([3])
# Ns = [5.16]
ct_solutions = dict([(
    c, [CTSolution(M=60, Nstar=N, phi=c / 100.0) for N in Ns])
    for c in concentrations])

dir = os.path.dirname(__file__)
with open(os.path.join(dir, 'ct_solutions.pickle'), 'wb') as f:
    pickle.dump(ct_solutions, f)