# -*- coding: utf-8 -*-
"""Breakup models

This module contains breakup models objects for population balance equation

"""


# Importing dependencies
from numpy import arange, sqrt, exp, pi


class breakupModels():
    def __init__(self) -> None:
        pass

    def gamma(self, v) -> float:
        """Breakup rate

        Args:
            v (float): droplet volume

        Returns:
            float: breakup rate for the parent drop
        """
        C1 = self.C[0]   # TODO
        C2 = self.C[1]
        return \
            C1 * v**(-2. / 9) * self.epsilon**(1. / 3) / (1 + self.phi) * \
            exp(- C2 * (1 + self.phi)**2 * self.sigma /
                (self.rhod * v**(5. / 9) * self.epsilon**(2. / 3))) # TODO

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
