# -*- coding: utf-8 -*-
"""Coalescence models

This module contains coalescence models objects for population balance equation

"""


# Importing dependencies
from numpy import arange, sqrt, exp, pi
from pbe.setup.system import Domain, ContinuosPhase, DispersePhase


class coalescenceModels():
    def __init__(self, C,
                 domain: Domain,
                 cp: ContinuosPhase,
                 dp: DispersePhase) -> None:

        self.C = C
        self.domain = domain
        self.cp = cp
        self.dp = dp

    def Qf(self, v1: float = None, v2: float = None) -> float:
        """Coalescence rate

        Args:
            v1 (float): volume of droplet i
            v2 (float): volume of droplet j

        Returns:
            float: _description_
        """
        C3 = self.C[2]
        C4 = self.C[3]
        d_ratio = (v1**(1. / 3) * v2**(1. / 3)) / (v1**(1. / 3) + v2**(1. / 3))

        return C3 * \
            (v1**(2. / 3) + v2**(2. / 3)) * \
            (v1**(2. / 9) + v2**(2. / 9))**0.5 * \
            self.cp.epsilon**(1. / 3) / \
            ((1 + self.dp.phi)) * \
            exp(
                -C4 * self.cp.mu * self.cp.rho * self.cp.epsilon /
                self.dp.sigma**2 /
                (1 + self.dp.phi)**3 *
                d_ratio**4)


    """Frequencia de coalescencia
    """
    def Qf2(self, v1: float = None, v2: float = None) -> float:
        """Coalescence rate

        Args:
            v1 (float): volume of droplet i
            v2 (float): volume of droplet j

        Returns:
            float: _description_
        """
        C3 = self.C[2]
        C4 = self.C[3]
        d_ratio = (v1**(1. / 3) * v2**(1. / 3)) / (v1**(1. / 3) + v2**(1. / 3))
        return C3 * \
            (v1**(2. / 3) + v2**(2. / 3)) * \
            (v1**(2. / 9) + v2**(2. / 9))**0.5 * \
            self.domain.D**(2. / 3) * \
            self.domain.Nstar / (1 + self.dp.phi) * \
            exp(
                -C4 * self.cp.mu * self.cp.rho * self.domain.D**2 /
                self.dp.sigma**2 *
                self.domain.Nstar**3 /
                (1 + self.dp.phi)**3 *
                d_ratio**4)
