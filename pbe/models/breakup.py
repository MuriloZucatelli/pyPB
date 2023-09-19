# -*- coding: utf-8 -*-
"""Breakup models

This module contains breakup models objects for population balance equation

"""


# Importing dependencies
from numpy import arange, sqrt, exp, pi
from pbe.setup.system import Domain, ContinuosPhase, DispersePhase


class breakupModels():
    def __init__(self, C, domain: Domain, cp: ContinuosPhase, dp: DispersePhase) -> None:
        self.C = C
        self.domain = domain
        self.cp = cp
        self.dp = dp

    def gamma(self, v) -> float:
        """Breakup rate

        Args:
            v (float): droplet volume

        Returns:
            float: breakup rate for the parent drop
        """
        C1 = self.C[0]
        C2 = self.C[1]
        return \
            C1 * v**(-2. / 9) * self.cp.epsilon**(1. / 3) / (1 + self.dp.phi) * \
            exp(- C2 * (1 + self.dp.phi)**2 * self.dp.sigma /
                (self.dp.rho * v**(5. / 9) * self.cp.epsilon**(2. / 3)))

    def gamma2(self, v) -> float:
        # Breakup model
        C1 = self.C[0]
        C2 = self.C[1]

        return \
            C1 * v**(-2. / 9) * self.domain.D**(2. / 3) * \
            self.domain.Nstar / (1 + self.dp.phi) * \
            exp(- C2 * (1 + self.dp.phi)**2 * self.dp.sigma /
                (self.dp.rho * v**(5. / 9) *
                 self.domain.D**(4. / 3) * self.domain.Nstar**2))
    
    @staticmethod
    def beta(v1: float, v2: float) -> float:
        """DDSD: dimensionless daughter particle size distribution

        Args:
            v1 (float): volume of droplet 1
            v2 (float): volume of droplet 1

        Returns:
            float: _description_
        """
        return 2.4 / v2 * exp(-4.5 * (2 * v1 - v2)**2 / (v2**2))
