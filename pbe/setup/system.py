"""Fluid system classes
"""
from numpy import pi
import pandas as pd


class Fluid:
    """Definição basica do fluido"""

    def __init__(self, name: str, rho: float, mu: float) -> None:
        """Inicializa o fluido

        Args:
            name (str): nome
            rho (float): massa especifica
            mu (float): dynamic viscosity
        """
        self.name = name
        self.rho = rho
        self.mu = mu


class Fluid2:
    """Definição basica do fluido"""

    def __init__(self, name: str, rho: pd.DataFrame, mu: pd.DataFrame) -> None:
        """Inicializa o fluido

        Args:
            name (str): nome
            rho (DataFrame): massa especifica
            mu (DataFrame): dynamic viscosity
        """
        self.name = name
        self.rho = rho
        self.mu = mu


class Prop:
    def __init__(self, emulsion: Fluid2, water: Fluid2, oil: Fluid2) -> None:
        """Inicializa o sistema de propriedades"""
        self.emulsion = emulsion
        self.water = water
        self.oil = oil


class ContinuosPhase(Fluid):
    """Fase continua"""

    def __init__(self, name: str, rho: float, mu: float, epsilon: float = None) -> None:
        """Continuos phase properties

        Args:
            name (str): nome
            rho (float): massa especifica
            mu (float): dynamic viscosity
            epsilon (float) : turbulent energy dissipation rate
        """
        super(ContinuosPhase, self).__init__(name, rho, mu)
        self.epsilon = epsilon
        self.nu = mu / rho


class DispersePhase(Fluid):
    """Properties of the dispersed phase"""

    def __init__(
        self,
        name: str,
        phi: float,
        rho: float,
        sigma: float,
        v_max: float = None,
        v0: float = None,
        sigma0: float = None,
        mu: float = None,
    ):
        """Cria o objeto de fase dispersa

        Args:
            name (str)         name
            phi (float)        volume fraction
            rho (float)        density
            sigma (float)      interfacial tension
            vMax (float)       maximum volume for the discretization
            v0 (float)         mean volume for the distributinon function
            sigma0 (float)     standard deviation for the distribution function
        """
        super(DispersePhase, self).__init__(name, rho, mu)
        self.phi = phi
        self.sigma = sigma
        self.v_max = v_max
        self.v0 = v0
        self.sigma0 = sigma0


class Domain:
    def __init__(
        self,
        theta: float = None,
        tres: float = None,
        V: float = None,
        M: float = None,
        xi: float = None,
        dxi: float = None,
        D: float = None,
        Nstar: float = None,
        w: float = None,
        A: float = None,
        Q: float = None,
    ) -> None:
        """domainProperties
         description: is a dictionary with properties of the computational
         domain and discretization parameters
        Args:
            theta (float): mean residence time
            V (float): domain volume
            M (float): number of classes
            D (float, optional): Impeler diameter. Defaults to None.
            Nstar (float, optional): Impeler revolutions per second, sec-1. Defaults to None.
            w (float): Width of channel (Mitre)
            U (float): Possible fluid velocity
            A float: Valve opening %
            Q (float): flow rate [m³/s]
        """
        self.theta = theta
        self.tres = tres
        self.V = V
        self.M = M
        self.xi = xi
        self.dxi = dxi
        self.D = D
        self.Nstar = Nstar
        self.w = w
        self.A = A
        self.Q = Q


class BoundaryInitial:
    def __init__(self, V) -> None:
        """_summary_

        Args:
            V (float): velocity
        """
        self.V = V


class droplet:
    def __init__(self, v: float) -> float:
        self.v = v

    @property
    def v2d(self) -> float:
        """Transform volume of droplet to a spherical type diameter

        Args:
            v (float): volume

        Returns:
            float: droplet spherical mean diameter
        """
        return (6 * self.v / pi) ** (1.0 / 3)


class FLOW:
    def __init__(
        self,
        flow: pd.DataFrame,
        flow_ext: pd.DataFrame,
        flow_mean: pd.DataFrame,
        flow_std: pd.DataFrame,
    ) -> None:
        self.flow = flow
        self.flow_ext = flow_ext
        self.flow_mean = flow_mean
        self.flow_std = flow_std
