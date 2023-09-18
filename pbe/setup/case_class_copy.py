from pbe.solvers.moc import MOCSolution
from numpy import arange, sqrt, exp, pi
from pbe.models import breakup, coalescence, DSD
from pbe.setup.system import ContinuosPhase, DispersePhase, Domain




class CaseSolution(MOCSolution):
    def __init__(
            self,
            dispersion: DispersePhase,
            continuo: ContinuosPhase,
            domain: Domain,
            model_parameters=None,
            time=arange(0.0, 3600, 0.5)):
        """Case configuration for solution

        Args:
            dispersion (DispersePhase):
            contProperties (ContinuosPhase):
            domain (Domain):
            model_parameters (dict, optional):
                description: list of four values C1 - C4 representing breakup model 
                constants (C1 and C2) and coalescence model constants (C3 and C4)
            time (_type_, optional): discretized time domain. Defaults to arange(0.0, 3600, 0.5).
        """
        self.contProperties = continuo
        self.phi = dispersion.phi

        self.muc = continuo.mu
        self.rhoc = continuo.rho
        self.epsilon = continuo.epsilon
        self.rhod = dispersion.rho
        self.sigma = dispersion.sigma

        vmax = dispersion.v_max

        # Feed distribution
        self.v0 = dispersion.v0
        self.sigma0 = dispersion.sigma0

        # Feed
        theta = domain.theta
        M = domain.M  # TODO: Fix variable name, what is M?
        self.V = domain.V
        if theta is None:
            self.n0 = None
        else:
            self.n0 = self.V / theta

        if model_parameters is None:
            self.C = [0.4, 0.08, 2.8, 1.83e13]
        else:
            self.C = model_parameters
        
        # Função frequencia de quebra
        beta = breakup.breakupModels.beta

        # Função 
        g = breakup.breakupModels(C=self.C, cp=continuo, dp=dispersion).gamma
        Qf = coalescence.coalescenceModels(C=self.C, cp=continuo, dp=dispersion).Qf
        A0 = DSD.analitico(dp=dispersion).A0
        
        MOCSolution.__init__(
            self, M, time, vmax / M,
            beta=beta, gamma=g, Q=Qf,
            theta=theta, n0=self.n0, A0=A0)

    
    @property
    def pbe_phi(self):
        return self.total_volume / self.V
