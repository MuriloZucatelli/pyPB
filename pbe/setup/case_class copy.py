from .moc import MOCSolution
from numpy import arange, sqrt, exp, pi
from pbe.models import breakup


class DomainProperties:
    """    domainProperties:
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


class DispersionProperties:
    """Dictionary with properties of the dispersed phase
        dispersion = dict()
        Args:
            phi (float)        volume fraction
            rho (float)        density
            sigma (float)      interfacial tension
            vMax (float)       maximum volume for the discretization
            v0 (float)         mean volume for the distributinon function
            sigma0 (float)     standard deviation for the distribution function
    """
    def __init__(self, phi, rho, sigma, v_max, v0, sigma0):
        """_summary_

        Args:
            phi (_type_): _description_
            rho (_type_): _description_
            sigma (_type_): _description_
            v_max (_type_): _description_
            v0 (_type_): _description_
            sigma0 (_type_): _description_
        """
        self.phi = phi
        self.rho = rho
        self.sigma = sigma
        self.v_max = v_max
        self.v0 = v0
        self.sigma0 = sigma0


beta = breakup.breakupModels.beta  # Função passada


'''
input:
    dispersion:
        description: is a 
        

    contProperties:
        description: is a dictionary with properties of the continuous phase
        Required fields:
        mu          viscosity
        rho         density
        epsilon     turbulent energy dissipation rate



    model_parameters:
        description: list of four values C1 - C4 representing breakup model
        constants (C1 and C2) and coalescence model constants (C3 and C4)

    time:
        description: discretized time domain
'''


class CaseSolution(MOCSolution):
    def __init__(
            self,
            dispersion,
            contProperties,
            domain,
            model_parameters=None,
            time=arange(0.0, 3600, 0.5)):

        self.contProperties = contProperties
        self.phi = dispersion.phi

        self.muc = contProperties['mu']
        self.rhoc = contProperties['rho']
        self.epsilon = contProperties['epsilon']
        self.rhod = dispersion.rho
        self.sigma = dispersion.sigma

        vmax = dispersion.v_max

        # Feed distribution
        self.v0 = dispersion.v0
        self.sigma0 = dispersion.sigma0

        # Feed
        theta = domain.theta
        M = domain.M  # TODO: Fix variable name, what is M?
        self.Vt = domain.Vt
        if theta is None:
            self.n0 = None
        else:
            self.n0 = self.Vt / theta

        if model_parameters is None:
            self.C = [0.4, 0.08, 2.8, 1.83e13]
        else:
            self.C = model_parameters

        MOCSolution.__init__(
            self, M, time, vmax / M,
            beta=beta, gamma=self.g, Q=self.Qf,
            theta=theta, n0=self.n0, A0=self.A0)

    
    @property
    def pbe_phi(self):
        return self.total_volume / self.Vt
