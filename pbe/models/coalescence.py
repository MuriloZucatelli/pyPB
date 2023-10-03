# -*- coding: utf-8 -*-
"""Coalescence models

This module contains coalescence models objects for population balance equation

"""


# Importing dependencies
from numpy import arange, sqrt, exp, pi, log
from pbe.setup.system import Domain, ContinuosPhase, DispersePhase, droplet
from pbe.models.modelParameters import modelParameter


class coalescenceModels:
    def __init__(
        self,
        name: str,
        C: list,
        domain: Domain,
        cp: ContinuosPhase,
        dp: DispersePhase,
        mp: modelParameter = None,
        interface: str = "rigid_interface",
    ) -> None:
        """Classe de modelos de coalescencia

        Args:
            C (list): parameters.\n
            name (str). coulaloglou, alopaeus, mitre ou chesters \n
            domain (Domain).\n
            cp (ContinuosPhase).\n
            dp (DispersePhase).\n
            mp (modelParameters, optional): Defaults to None.\n
            interface (str, optional): droplet interface. Defaults to "rigid_interface".
                Options: rigid_interface"  "deformable" or "partially_mobile"\n
        """
        self.C = C
        self.mp = mp
        self.domain = domain
        self.cp = cp
        self.dp = dp
        self.interface = interface

        match name:
            case "coulaloglou":
                self.Q = self.coulaloglou_Q
            case "alopaeus":
                self.Q = None
            case "mitre":
                self.Q = self.mitre_Q
            case "chesters":
                self.Q = self.chesters_Q
            case None:
                self.Q = None

    def coulaloglou_Q(self, v1: float = None, v2: float = None) -> float:
        """Coalescence rate from coulaloglou and Tavlarides

        Args:
            v1 (float): volume of droplet i
            v2 (float): volume of droplet j

        Returns:
            float: _description_
        """
        C3 = self.C[2]
        C4 = self.C[3]
        d_ratio = (v1 ** (1.0 / 3) * v2 ** (1.0 / 3)) / (
            v1 ** (1.0 / 3) + v2 ** (1.0 / 3)
        )

        return (
            C3
            * (v1 ** (2.0 / 3) + v2 ** (2.0 / 3))
            * (v1 ** (2.0 / 9) + v2 ** (2.0 / 9)) ** 0.5
            * self.cp.epsilon ** (1.0 / 3)
            / ((1 + self.dp.phi))
            * exp(
                -C4
                * (self.cp.mu * self.cp.rho * self.cp.epsilon)
                / self.dp.sigma**2
                / (1 + self.dp.phi) ** 3
                * d_ratio**4
            )
        )

    def coulaloglou_Q_vessel(self, v1: float = None, v2: float = None) -> float:
        """Coalescence rate

        Args:
            v1 (float): volume of droplet i
            v2 (float): volume of droplet j

        Returns:
            float: _description_
        """
        C3 = self.C[2]
        C4 = self.C[3]
        d_ratio = (v1 ** (1.0 / 3) * v2 ** (1.0 / 3)) / (
            v1 ** (1.0 / 3) + v2 ** (1.0 / 3)
        )
        return (
            C3
            * (v1 ** (2.0 / 3) + v2 ** (2.0 / 3))
            * (v1 ** (2.0 / 9) + v2 ** (2.0 / 9)) ** 0.5
            * self.domain.D ** (2.0 / 3)
            * self.domain.Nstar
            / (1 + self.dp.phi)
            * exp(
                -C4
                * (self.cp.mu * self.cp.rho * self.domain.D**2)
                / self.dp.sigma**2
                * self.domain.Nstar**3
                / (1 + self.dp.phi) ** 3
                * d_ratio**4
            )
        )

    def chesters_Q(self, v1: float = None, v2: float = None) -> float:
        """Coalescence frequency

        Args:
            v1 (float): volume of droplet i
            v2 (float): volume of droplet j
        Returns:
            float: _description_
        """
        # C3 = self.C[2]
        # C4 = self.C[3]
        Cc = self.mp.Cc
        Ce = self.mp.Ce

        d1 = (6 * v1 / pi) ** (1.0 / 3)
        d2 = (6 * v2 / pi) ** (1.0 / 3)

        # equivalent diameter in the collision of drops i and j
        deq = (0.5 * (d1 ** (-1) + d2 ** (-1))) ** (-1)
        # collision cross-section of particles with particles d and d'
        sec_area = pi / 4 * (d1 + d2) ** 2

        u_d1 = self.cv_chesters_1991(d=d1, epsilon=self.cp.epsilon, nu=self.cp.nu)
        u_d2 = self.cv_chesters_1991(d=d2, epsilon=self.cp.epsilon, nu=self.cp.nu)

        # relative velocity of droplets
        u_r = (u_d1 + u_d2) / 2

        # Collision frequency
        coli_freq = Cc * sec_area * u_r

        A = 0.6e-20  # Joules
        hf = A * deq / (16 * pi * self.dp.sigma)

        # Models for coalescence eficience
        # depends on the type of droplet interface
        if self.interface == "rigid_interface":
            if Ce is None:
                eficiency = exp(-1 / 4 * log(1 / hf))  # TODO revisar
            else:
                hi = 1
                eficiency = exp(-Ce / 4 * log(hi / hf))  # TODO revisar

        elif self.interface == "deformable" or self.interface == "partially_mobile":
            Caeq = (
                self.cp.mu
                * sqrt(self.cp.epsilon / self.cp.nu)
                * deq
                / (2 * self.dp.sigma)
            )
            eficiency = exp(
                -Ce
                * (sqrt(3) / 8)
                * (self.dp.mu / self.cp.mu)
                * Caeq ** (3 / 2)
                * (deq / hf)
            )
        else:
            eficiency = 0  # TODO revisar

        return coli_freq * eficiency

    def mitre_Q(self, v1: float = None, v2: float = None) -> float:
        """Coalescence frequency from Mitre Equation 40.
           modified because all droplets were in turbulent dissipation size range

        Args:
            v1 (float): volume of droplet i
            v2 (float): volume of droplet j
        Returns:
            float: _description_
        """
        # C3 = self.C[2]
        # C4 = self.C[3]
        Cc = self.mp.Cc
        Ce = self.mp.Ce

        d1 = (6 * v1 / pi) ** (1.0 / 3)
        d2 = (6 * v2 / pi) ** (1.0 / 3)

        # equivalent diameter in the collision of drops i and j
        deq = (0.5 * (d1 ** (-1) + d2 ** (-1))) ** (-1)

        # Collision frequency
        coli_freq = Cc * pi / 8 * (d1 + d2) ** 3 * sqrt(self.cp.epsilon / self.cp.nu)

        A = 0.6e-20  # Joules
        hf = A * deq / (16 * pi * self.dp.sigma)

        # Models for coalescence eficience
        # depends on the type of droplet interface
        if self.interface == "rigid_interface":
            if Ce is None:
                eficiency = exp(-1 / 4 * log(1 / hf))  # TODO revisar
            else:
                hi = 1
                eficiency = exp(-Ce / 4 * log(hi / hf))  # TODO revisar

        elif self.interface == "deformable" or self.interface == "partially_mobile":
            Caeq = (
                self.cp.mu
                * sqrt(self.cp.epsilon / self.cp.nu)
                * deq
                / (2 * self.dp.sigma)
            )
            eficiency = exp(
                -Ce
                * (sqrt(3) / 8)
                * (self.dp.mu / self.cp.mu)
                * Caeq ** (3 / 2)
                * (deq / hf)
            )
        else:
            eficiency = 0  # TODO revisar

        return coli_freq * eficiency

    @classmethod
    def cv_chesters_1991(cls, d: float, epsilon: float, nu: float) -> float:
        """Chesters characteristic velocity model 1991

        Args:
            d (float): droplet diameter
            epsilon (float): turbulent dissipation rate
            nu (float): cinematic viscosity

        Returns:
            float: collision
        """
        eta = (nu**3 / epsilon) ** (1 / 4)
        if d < eta:
            return d * sqrt(epsilon / nu)
        else:
            return (epsilon * d) ** (1.0 / 3)

    # @classmethod
    # def ce_chesters_1991(cls, C, deq, dp: DispersePhase)
    # """Chesters coalescence efficiency model 1991
    # """
