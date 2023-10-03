from numpy import genfromtxt, abs, array, pi
from scipy.optimize import minimize, differential_evolution
from pbe.app.angeli_class import AngeliSolution
import time
import pickle
import os


class mitre_experiment:
    def __init__(self, U, phi, d32, theta=None):
        self.phi = phi          # volumetric concentration
        self.theta = theta      # residence time [s]
        self.U = U              # [m/s]
        self.d32 = d32          # [m]