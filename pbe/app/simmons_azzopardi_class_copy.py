from numpy import arange, pi
from pbe.setup.system import Domain, DispersePhase, ContinuosPhase
from pbe.solvers.moc import MOCSolution
from pbe.models import breakup, coalescence, DSD


class SASolution:
    def __init__(self,
                 M=10,        # number of classes
                 U=2.71,      # average fluid velocity
                 phi=0.117,   # [1] holdup
                 v_max=8e-11,  # Maximum droplet volume [m³]
                 v0=4e-11,    # [m³]
                 model_parameters=None,
                 theta=600.):

        self.D = 0.063  # [m] pipe diameter diameter
        self.L = 4.5    # [m] length diameter
        self.phi = phi
        time = arange(0.0, 3600, 0.5)  # Tempo de integração



        # oil
        self.cp = ContinuosPhase(name='oil',
                                 rho=797.,    # [kg/m³]
                                 mu=1.8e-3)   # [kg*m^-1 s^-1] dyn. viscosity
        self.cp.epsilon = calc_epsilon(U,
                                       self.D,
                                       self.cp)
        # water solution
        self.dp = DispersePhase(name='water solution',
                                phi=self.phi, rho=1166.,  # [kg/m3]
                                sigma=1.e-2,  # [P = kg * m^-1 s^-1]
                                v_max=6e-11,  # DIFERENT FROM INIT
                                v0=v0,
                                sigma0=v0 / 10)

        self.domain = Domain(theta=theta,
                             V=pi * self.L * (self.D / 2) ** 2,
                             M=M)

        self.set_parameters(model_parameters)
        self.nf0()
        vmin = None

        # Função frequencia de quebra
        beta = breakup.breakupModels.beta

        # Função
        g = breakup.breakupModels(C=self.C,
                                  domain=self.domain,
                                  cp=self.cp,
                                  dp=self.dp).gamma
        Qf = coalescence.coalescenceModels(C=self.C,
                                           domain=self.domain,
                                           cp=self.cp,
                                           dp=self.dp).Qf
        A0 = DSD.analitico(dp=self.dp).A0

        self.moc = MOCSolution(M, time, self.dp.v_max / M, 
                               n0=self.n0, xi0=vmin, beta=beta,
                               gamma=g, Q=Qf, theta=self.domain.theta,
                               nf0=self.nf0, A0=A0)

    def set_parameters(self, model_parameters):
        if model_parameters is None:
            self.C = [0.4, 0.08, 2.8, 1.83e13]
        else:
            self.C = model_parameters

    def n0(self, v):              # total number of drops?
        return 0 * v

    def nf0(self):
        self.nf0 = self.domain.V / self.domain.theta

    @property
    def pbe_phi(self):
        return self.moc.total_volume / self.domain.V


# Revisado não testado
def calc_epsilon(U: float, D: float, cp: ContinuosPhase):
    """Calcula as propriedades turbulentas

    Args:
        U (float): average flow velocity
    """
    Re = U * D / cp.mu * cp.rho
    TI = 0.16 * Re ** (-1. / 8.)
    u_rms = U * TI
    k = 3. / 2. * u_rms ** 2
    L_t = 0.038 * D
    epsilon = 0.09 * k ** (3. / 2.) / L_t
    return epsilon
