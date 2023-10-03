# -*- coding: utf-8 -*-

"""
    Module for droplet size distribuition functions
"""


#importing dependences
from numpy import exp, pi, sqrt
from pbe.setup.system import DispersePhase


class analitico:
    """Analitical probability number density functions
    """
    def __init__(self, dp: DispersePhase) -> None:
        self.dp = dp

    def A0(self, v):
        """C&T

        Args:
            v (float): droplet volume

        Returns:
            probability density number distribuition
        """
        return \
            self.dp.phi / (self.dp.v0 * self.dp.sigma0 * sqrt(2 * pi)) * \
            exp(-(v - self.dp.v0)**2 / (2 * self.dp.sigma0**2))
        # Gaussian distribuition
