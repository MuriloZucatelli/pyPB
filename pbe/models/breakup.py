# -*- coding: utf-8 -*-
"""Breakup models

This module contains breakup models objects for population balance equation

"""


# Importing dependencies
from numpy import arange, sqrt, exp, pi, log, linspace, interp, array, where, zeros_like
import sys
import os.path as path
from scipy.special import gamma

dir = path.dirname(__file__)
if __name__ == "__main__":
    sys.path.append(path.abspath(path.join(dir, "..\\..")))
from pbe.setup.system import Domain, ContinuosPhase, DispersePhase, droplet
from pbe.models.modelParameters import modelParameter


class breakupModels:
    def __init__(
        self,
        name: str,
        mp: modelParameter,
        domain: Domain,
        cp: ContinuosPhase,
        dp: DispersePhase,
        nu=None,
    ) -> None:
        """Classe de modelos de quebra

        Args:
            name (str). coulaloglou, alopaeus, mitre, mb ou cristini\n

            domain (Domain).\n

            cp (ContinuosPhase).\n

            dp (DispersePhase).\n

            mp (modelParameters, optional): Defaults to None.\n
        """
        self.domain = domain
        self.cp = cp
        self.dp = dp
        self.mp = mp
        self.breakupmodelslist = [
            "coulaloglou",
            "alopaeus",
            "mitre_modified",
            "mb",
            "cristini",
        ]
        if name == "coulaloglou":
            self.gamma = self.coulaloglou_gamma_damping
        elif name == "alopaeus":
            self.gamma = None
        elif name == "mitre_modified":
            self.gamma = self.mitre_modified_gamma
        elif name == "mb":
            self.gamma = self.mb_gamma
        elif name == "cristini":
            self.gamma = self.cristini_gamma_1
        elif name is None:
            self.gamma = None
        else:
            raise Exception(
                f"breakup model {name} not found, try one of the brekup models: "
                + "\n"
                + f"{self.breakupmodelslist}"
            )

        # match name:
        #    case "coulaloglou":
        #        self.gamma = self.coulaloglou_gamma_damping
        #    case "alopaeus":
        #       pass
        #  case "mitre":
        #       self.gamma = self.mitre_gamma
        #   case "mb":
        #       self.gamma = self.mb_gamma
        #   case "cristini":
        #        self.gamma = self.cristini_gamma_1
        #    case None:
        #        self.gamma = None

    def mitre_modified_gamma(self, v: float) -> float:
        """
           João F. Mitre. Droplet breakage and coalescence models for
           the flow of water-in-oil emulsions through a valve-like element
           Equation 41
           modified because all droplets were in turbulent dissipation size range
        Eqs. (8), (40), (41), (27) and (28) for the model with constant coalescence
        efficiency and also using Eq. (18) for the coalescence efficiency model
        for drops with deformable and partially mobile interfaces.
        Args:
            v (float): droplet volume

        Returns:
            float: _description_
        """

        Cb = self.get_constant(self.mp.Cb, "Cb")

        d = (6 * array(v) / pi) ** (1 / 3)

        # Capilarity number
        Ca = self.cp.mu * d * sqrt(self.cp.epsilon / self.cp.nu) / (2 * self.dp.sigma)

        # Ca_crit: Modeled by assuming that dcrit is the minimum attainable droplet diameter at
        # the outlet of the test section as most of the data are for breakage dominant conditions
        # Experimental parameters
        CCa = 1.65e-4  # The minimum attainable droplet diameter at the
        C_Re = 3 / 20  # outlet of the test section could be correlated with
        We_crit = 6
        # Eq. 27
        Stk = self.domain.tres / sqrt(self.cp.nu / self.cp.epsilon)
        _, _, h = self.calc_U_max(1, self.domain.A, self.domain.D)
        Re_max = self.domain.Q / (h * self.cp.nu)
        # Re_max = (Q / A_min) * c / nu  # c é a abertura minima da válvula
        Ca_crit = CCa * Stk * Re_max ** (-C_Re)

        eta = (self.cp.nu**3 / self.cp.epsilon) ** (1 / 4)

        i = where(Ca > Ca_crit)
        b = zeros_like(d)

        b[i] = (
            Cb
            * 63.927
            / (We_crit ** (11 / 5))
            * sqrt(self.cp.epsilon / self.cp.nu)
            * Ca[i] ** 2.2
            * (d[i] / (2 * eta)) ** (4 / 5)
        )
        # b[not i] = 0
        return b

    def cristini_gamma_1(self, v: float) -> float:
        """
           Cristini et al. (2003) breakup frequency

        Args:
            v (float): droplet volume

        Returns:
            float: _description_
        """
        # Cristini breakup constant
        Cb_c = self.mp.Cb

        d = (6 * v / pi) ** (1 / 3)
        Ca = self.cp.mu * d * sqrt(self.cp.epsilon / self.cp.nu) / (2 * self.dp.sigma)
        # Ca_crit: Modeled by assuming that dcrit is the minimum attainable droplet diameter at
        # the outlet of the test section as most of the data are for breakage dominant conditions
        # Experimental parameters
        CCa = 1.65e-4  # The minimum attainable droplet diameter at the
        C_Re = 3 / 20  # outlet of the test section could be correlated with
        Stk = self.domain.theta / sqrt((self.cp.mu / self.cp.rho) / self.cp.epsilon)
        Re_max = self.cp.Q / (self.domain.l * self.cp.nu)
        Ca_crit = CCa * Stk * Re_max ** (-C_Re)

        if Ca > Ca_crit:
            b = Cb_c * sqrt(self.cp.epsilon / self.cp.nu) * Ca**3
        else:
            b = 0
        return b

    def mb_gamma(self, v: float) -> float:
        """
           Martinez-Bazán et al. (1999a,b) breakup frequency

        Args:
            v (float): droplet volume

        Returns:
            float: _description_
        """
        Cb = self.mp.Cb
        beta = 8.2  # TODO universalizar esse parâmetro para 'mp class'
        We_crit = 6  # TODO universalizar esse parâmetro para 'mp class'

        d = (6 * v / pi) ** (1.0 / 3)

        We = (
            (self.cp.rho * beta * (self.cp.epsilon * d) ** (2 / 3))
            * d
            / (2 * self.dp.sigma)
        )

        # critical particle diameter from martinez-bazan
        dc_mb = (2 * We_crit * self.dp.sigma / (beta * self.cp.rho)) ** (
            3.0 / 5
        ) * self.cp.epsilon ** (-2.0 / 5)

        if d > dc_mb:
            b = (
                Cb
                * sqrt(
                    beta * (self.cp.epsilon * d) ** (2.0 / 3)
                    - 2 * We_crit * self.dp.sigma / (self.cp.rho * d)
                )
                / d
            )
        else:
            b = 0
        return b

    def coulaloglou_gamma_damping(self, v) -> float:
        """Coulaloglou Breakup rate Equation 36 (1977) or Eq. 4.90 from Guan 2014.
           account for the “damping” effect of droplets on the local turbulent
           intensities at high volume fraction

        Args:
            v (float): droplet volume

        Returns:
            float: breakup rate for the parent drop
        """
        Cb = self.get_constant(self.mp.Cb, "Cb")
        Cepsilon = self.get_constant(self.mp.Cepsilon, "Cepsilon")

        return (
            Cb
            * v
            ** (-2.0 / 9)  # should be: d^(-2/3), so C1 is different from Liao and Lucas
            * self.cp.epsilon ** (1.0 / 3)
            / (1 + self.dp.phi)
            * exp(
                -Cepsilon
                * self.dp.sigma
                * (1 + self.dp.phi) ** 2
                / (self.dp.rho * v ** (5.0 / 9) * self.cp.epsilon ** (2.0 / 3))
            )
        )

    def coulaloglou_gamma_damping_vessel(self, v) -> float:
        """Coulaloglou Breakup rate Equation 36
           When the energy is uniformly distributed throughout the vessel: epsilon = K'N*³D²
           To account for the “damping” effect of droplets on the local turbulent intensities at
           high holdup fractions, Coulaloglou and Tavlarides modified the original expression.
        Args:
            v (float): droplet volume

        Returns:
            float: breakup rate for the parent drop
        """
        if self.mp.Cb is not None:
            Cb = self.mp.Cb
        else:
            raise Exception(f"Parameter Cb not defined in args")
        if self.mp.Cepsilon is not None:
            Cepsilon = self.mp.Cepsilon
        else:
            raise Exception(f"Parameter Cepsilon not defined in args")

        return (
            Cb
            * v ** (-2.0 / 9)
            * self.domain.D ** (2.0 / 3)
            * self.domain.Nstar
            / (1 + self.dp.phi)  # remove / (1 + self.dp.phi)
            * exp(
                -Cepsilon
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

    def coulaloglou_gamma_vessel(self, v) -> float:
        """Coulaloglou Breakup rate Equation 14
           When the energy is uniformly distributed throughout the vessel: epsilon = K'N*³D²
        Args:
            v (float): droplet volume

        Returns:
            float: breakup rate for the parent drop
        """
        if self.mp.Cb is not None:
            Cb = self.mp.Cb
        else:
            raise Exception(f"Parameter Cb not defined in args")
        if self.mp.Cepsilon is not None:
            Cepsilon = self.mp.Cepsilon
        else:
            raise Exception(f"Parameter Cepsilon not defined in args")

        return (
            Cb
            * v ** (-2.0 / 9)
            * self.domain.D ** (2.0 / 3)
            * self.domain.Nstar
            * exp(
                -Cepsilon
                * self.dp.sigma
                / (
                    self.dp.rho
                    * v ** (5.0 / 9)
                    * self.domain.D ** (4.0 / 3)
                    * self.domain.Nstar**2
                )
            )
        )

    @classmethod
    def get_constant(cls, constant, name):
        if constant is not None:
            return constant
        else:
            raise Exception(f"Parameter {name} not defined in args")

    @classmethod
    def calc_U_max(cls, U, A, D):
        """Calcula U max de acordo com a abertura mínima da válvula
        C é a abertura mínima
            w: Vazao kg/min
            A: Abertura de valvula 0 a 100
        """
        Ab = linspace(0.1, 1, 1000)
        y = linspace(0.0519 / 1000, D, 1000)
        Ab = linspace(0, 1, 1000)  # Abertura de 0 a 1
        y = linspace(0.0 / 1000, D, 1000)  # Altura de 0 a D
        h = interp(A / 100, Ab, y)  # altura interpolada para A
        U_max = U * D / h
        At_min = 0
        return U_max, At_min, h


class DDSD:
    def __init__(
        self,
        name,
        mp: modelParameter,
        domain: Domain,
        cp: ContinuosPhase,
        dp: DispersePhase,
        varsigma=None,
    ) -> None:
        self.domain = domain
        self.cp = cp
        self.dp = dp
        self.mp = mp
        self.varsigma = varsigma

        self.DDSDmodelslist = [
            "coulaloglou",
            "alopaeus",
            "mitre",
            "mb",
            "cristini",
        ]

        if name == "coulaloglou":
            self.beta = self.coulaloglou_beta
        elif name == "alopaeus":
            self.beta = None
        elif name == "mitre":
            self.beta = self.mitre_beta
        elif name == "mb":
            self.beta = None
        elif name == "cristini":
            self.beta = None
        elif name is None:
            self.beta = None
        else:
            raise Exception(
                f"Daughter Droplet model {name} not found, try one of the models: "
                + "\n"
                + f"{self.DDSDmodelslist}"
            )

        # match name:
        #    case "coulaloglou":
        #        self.beta = self.coulaloglou_beta
        #    case "alopaeus":
        #        self.beta = None
        #    case "mitre":
        #        self.beta = None
        #    case "mb_gamma":
        #        self.beta = None
        #    case "cristini":
        #        self.beta = None
        #    case None:
        #        self.beta = None

    @staticmethod
    def coulaloglou_beta(v1: float, v2: float, i: int) -> float:
        """DDSD: dimensionless daughter particle size distribution
           from coulaloglou and tavlarides 1977

        Args:
            v1 (float): volume of droplet 1 v daughter
            v2 (float): volume of droplet 2 v' mother

        Returns:
            float: _description_
        """
        return 2.4 / v2 * exp(-4.5 * (2 * v1 - v2) ** 2 / (v2**2))

    def mitre_beta(self, v1: float, v2: float, i: int) -> float:
        """DDSD: dimensionless daughter particle size distribution
           from Mitre
           The only one that gave good results is the
            simple model that assumes that ς equal daughter droplets are
            formed upon breakage
        Args:
            v1 (float): volume of droplet 1 v daughter
            v2 (float): volume of droplet 2 v' mother

        Returns:
            float: _description_
        """
        # delta(v1 - v2/nu)
        # varsigma é o número de gotas geradas na quebra (isso é problemático)

        self.domain.dxi
        self.varsigma
        if abs(v1 - v2 / self.varsigma) < self.domain.dxi[i] / 2:
            P = 1 / self.domain.dxi[i]
        else:
            P = 0
        return P

    # NOTE: Não é satisfeito a integral de P(v,v') para v>v_max/varsigma

    def lehr_beta(v1: float, v2: float, i: int) -> float:
        """Lehr 2002 M-shape DDSD"""
        fbv = v1 / v2  # Vi/Vj
        # num = exp(-9/4 * (log( (2 ** (2/5) * v2 * rho ** (3/5) * epsilon ** (2/5)) / sigma ** (3/5)))** 2)
        # den = exp(-9/4 * (log( (2 ** (2/5) * v2 * rho ** (3/5) * epsilon ** (2/5)) / sigma ** (3/5)) ** 2)
        # beta = 1 / (pi * fbv) * num / den
        pass


def check_DDSD(beta, domain):
    i = 0
    int_P = zeros_like(domain.xi)
    for x1 in domain.xi:
        for x2 in domain.xi:
            P = beta(x1, x2, i)
            int_P[i] += P * domain.dxi[i]
        i += 1
    return int_P
