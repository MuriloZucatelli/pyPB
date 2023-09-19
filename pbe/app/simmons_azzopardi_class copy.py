from numpy import arange, pi
from pbe.setup.system import Domain, DispersePhase, ContinuosPhase
from pbe.solvers.moc import MOCSolution
from pbe.models import breakup, coalescence, DSD


class SASolution:
    def __init__(self,
                 M=10,        # number of classes
                 U=2.71,      # average fluid velocity
                 phi=0.117,   # [1] holdup
                 vmax=8e-11,    # Maximum droplet volume [m³]
                 v0=4e-11,    # [m³]
                 model_parameters=None,
                 theta=600.):

        self.D = 0.063  # [m] pipe diameter diameter
        self.L = 4.5    # [m] length diameter
        self.phi = phi
        time = arange(0.0, 3600, 0.5)  # Tempo de integração

        self.calc_turbulent_properties(U)  # calculate turbulent properties
        self.n0()
        vmin = None


        # oil
        self.continuousphase = ContinuosPhase(
            name='oil',
            rho=797.,       # [kg/m³]
            mu=1.8e-3,     # [P = kg * m^-1 s^-1]  dynamic viscosity
            epsilon=self.epsilon)


        # water solution
        self.dispersephase = DispersePhase(
            name='water solution',
            phi=self.phi, rho=1166.,  # [kg/m3]
            sigma=1.e-2,  # [P = kg * m^-1 s^-1]
            v_max=6e-11,
            v0=v0,
            sigma0=v0 / 10)

        self.domain = Domain(
            theta=theta, V=pi * self.L * (self.D / 2) ** 2, M=M)

        # Função frequencia de quebra
        beta = breakup.breakupModels.beta

        # Função
        g = breakup.breakupModels(C=self.C,
                                  domain=self.domain,
                                  cp=self.continuousphase,
                                  dp=self.dispersephase).gamma
        Qf = coalescence.coalescenceModels(C=self.C,
                                           domain=self.domain,
                                           cp=self.continuousphase,
                                           dp=self.dispersephase).Qf
        A0 = DSD.analitico(dp=self.dispersephase).A0

        self.moc = MOCSolution(
            M, time, vmax / M, N0=self.N0, xi0=vmin,
            beta=beta, gamma=g, Q=Qf, theta=self.domain.theta,
            n0=self.n0, A0=A0)
        
    def calc_turbulent_properties(self, U):
        self.Re = U * self.D / self.continuousphase.mu * self.continuousphase.rho
        self.TI = 0.16 * self.Re ** (-1. / 8.)
        u_rms = U * self.TI
        k = 3. / 2. * u_rms ** 2
        L_t = 0.038 * self.D
        self.epsilon = 0.09 * k ** (3. / 2.) / L_t


    def N0(self, v):              # total number of drops?
        return 0 * v
    
    def n0(self):
        self.n0 = self.domain.V / self.domain.theta

    @property
    def pbe_phi(self):
        return self.moc.total_volume / self.domain.V
