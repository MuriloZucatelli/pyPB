"""Fluid system classes
"""
from numpy import pi


class Fluid:
    """ Definição basica do fluido
    """
    def __init__(self, name: str, rho: float, mu: float) -> None:
        """ Inicializa o fluido

        Args:
            name (str): nome
            rho (float): massa especifica
            mu (float): dynamic viscosity
        """
        self.name = name
        self.rho = rho
        self.mu = mu


class ContinuosPhase(Fluid):
    """Fase continua
    """
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
        self.nu = rho/mu


class DispersePhase(Fluid):
    """Properties of the dispersed phase
    """
    def __init__(self, name: str, phi: float, rho: float, sigma: float,
                 v_max: float, v0: float, sigma0: float, mu: float = None):
        """ Cria o objeto de fase dispersa
        
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
    """ domainProperties:
        description: is a dictionary with properties of the computational
        domain and discretization parameters
        Required fields:
        theta       mean residence time
        M           number of classes used
        V           domain volume
        D           Possible impeller diameter [m]
        Nstar       Possible revolutions per second, sec-1
    """
    def __init__(self,
                 theta: float,
                 V: float,
                 M: float,
                 D: float = None,
                 Nstar: float = None) -> None:
        """ domainProperties

        Args:
            theta (float): _description_
            V (float): domain volume
            M (float): number of classes
            D (float, optional): Impeler diameter. Defaults to None.
            Nstar (float, optional): Impeler revolutions. Defaults to None.
        """
        self.theta = theta
        self.V = V
        self.M = M
        self.D = D
        self.Nstar = Nstar
        # TODO: create pipe diameter Dp?


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
