"""Fluid system classes
"""

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
    """
    def __init__(self, theta, V, M):
        self.theta = theta
        self.V = V
        self.M = M
