# -*- coding: utf-8 -*-
"""Breakup models

This module contains breakup models objects for population balance equation

"""


# Importing dependencies
from numpy import arange, sqrt, exp, pi
import sys
import os.path as path
from scipy.special import gamma

dir = path.dirname(__file__)
if __name__ == "__main__":
    sys.path.append(path.abspath(path.join(dir, "..\\..")))
from pbe.setup.system import Domain, ContinuosPhase, DispersePhase, droplet
from pbe.models import modelParameters


class breakup_models:
    def __init__(
        self,
        name,
        C,
        domain: Domain,
        cp: ContinuosPhase,
        dp: DispersePhase,
        mp: modelParameters = None,
    ) -> None:
        self.C = C
        self.domain = domain
        self.cp = cp
        self.dp = dp
        self.mp = mp

        match name:
            case "coulaloglou":
                self.gamma = self.coulaloglou_gamma_eq36
                self.beta = self.coulaloglou_beta
            case "alopaeus":
                pass
            case "mitre":
                self.gamma = self.mitre_gamma_eq41
            case None:
                self.gamma = None

    def mitre_gamma_eq41(self, v: float) -> float:
        """
           João F. Mitre. Droplet breakage and coalescence models for the flow of water-in-oil
           emulsions through a valve-like element

        Args:
            v (float): droplet volume

        Returns:
            float: _description_
        """
        C1 = self.C[0]
        C2 = self.C[1]
        Cb = self.mp.Cb
        nu = self.cp.mu / self.cp.rho
        d = (6 * v / pi) ** (1 / 3)
        Ca = (
            self.cp.mu
            * d
            * sqrt(self.cp.epsilon / (self.cp.mu / self.cp.rho))
            / (2 * self.dp.sigma)
        )
        # Ca_crit: Modeled by assuming that dcrit is the minimum attainable droplet diameter at
        # the outlet of the test section as most of the data are for breakage dominant conditions
        # Experimental parameters
        CCa = 1.65e-4  # The minimum attainable droplet diameter at the
        C_Re = 3 / 20  # outlet of the test section could be correlated with
        Stk = self.domain.theta / sqrt((self.cp.mu / self.cp.rho) / self.cp.epsilon)
        Re_max = self.cp.Q / (self.domain.l * self.cp.mu / self.cp.rho)
        Ca_crit = CCa * Stk * Re_max ** (-C_Re)
        We_crit = 6

        eta = "CALCULATION ROUTINE"

        if Ca > Ca_crit:
            b = (
                Cb
                * 63.927
                / (We_crit ** (11 / 5))
                * sqrt(self.cp.epsilon / nu)
                * Ca**2.2
                * (d / (2 * eta)) ** (4 / 5)
            )
        else:
            b = 0
        return b

    def coulaloglou_gamma_eq36(self, v) -> float:
        """Coulaloglou Breakup rate. Eq. 4.90 from Guan 2014.
           account for the “damping” effect of droplets on the local turbulent
           intensities at high volume fraction

        Args:
            v (float): droplet volume

        Returns:
            float: breakup rate for the parent drop
        """
        C1 = self.C[0]
        C2 = self.C[1]
        return (
            C1
            * v
            ** (-2.0 / 9)  # should be: d^(2/3), so C1 is different from Liao and Lucas
            * self.cp.epsilon ** (1.0 / 3)
            / (1 + self.dp.phi)
            * exp(
                -C2
                * self.dp.sigma
                * (1 + self.dp.phi) ** 2
                / (self.dp.rho * v ** (5.0 / 9) * self.cp.epsilon ** (2.0 / 3))
            )
        )

    def coulaloglou_gamma_eq36_uniform(self, v) -> float:
        """Coulaloglou Breakup rate Eq. 36
           When the energy is uniformly distributed throughout the vessel: epsilon = K'N*³D²
           To account for the “damping” effect of droplets on the local turbulent intensities at
           high holdup fractions, Coulaloglou and Tavlarides modified the original expression.
        Args:
            v (float): droplet volume

        Returns:
            float: breakup rate for the parent drop
        """
        C1 = self.C[0]
        C2 = self.C[1]

        return (
            C1
            * v ** (-2.0 / 9)
            * self.domain.D ** (2.0 / 3)
            * self.domain.Nstar
            / (1 + self.dp.phi)  # remove / (1 + self.dp.phi)
            * exp(
                -C2
                * (1 + self.dp.phi) ** 2
                * self.dp.sigma
                / (
                    self.dp.rho
                    * v ** (5.0 / 9)
                    * self.domain.D ** (4.0 / 3)
                    * self.domain.Nstar**2
                )
            )
        )

    def coulaloglou_gamma_eq14_uniform(self, v) -> float:
        """Coulaloglou Breakup rate Equation 14
           When the energy is uniformly distributed throughout the vessel: epsilon = K'N*³D²
        Args:
            v (float): droplet volume

        Returns:
            float: breakup rate for the parent drop
        """
        C1 = self.C[0]
        C2 = self.C[1]

        return (
            C1
            * v ** (-2.0 / 9)
            * self.domain.D ** (2.0 / 3)
            * self.domain.Nstar
            * exp(
                -C2
                * self.dp.sigma
                / (
                    self.dp.rho
                    * v ** (5.0 / 9)
                    * self.domain.D ** (4.0 / 3)
                    * self.domain.Nstar**2
                )
            )
        )

    @staticmethod
    def coulaloglou_beta(v1: float, v2: float) -> float:
        """DDSD: dimensionless daughter particle size distribution

        Args:
            v1 (float): volume of droplet 1
            v2 (float): volume of droplet 1

        Returns:
            float: _description_
        """
        return 2.4 / v2 * exp(-4.5 * (2 * v1 - v2) ** 2 / (v2**2))
