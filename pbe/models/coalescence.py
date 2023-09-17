# -*- coding: utf-8 -*-
"""Coalescence models

This module contains coalescence models objects for population balance equation

"""


# Importing dependencies
from numpy import arange, sqrt, exp, pi


class coalescenceModels():
    def __init__(self) -> None:
        pass

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
            self.epsilon**(1. / 3) / \
            ((1 + self.phi)) * \
            exp(
                -C4 * self.muc * self.rhoc * self.epsilon /
                self.sigma**2 /
                (1 + self.phi)**3 *
                d_ratio**4)
