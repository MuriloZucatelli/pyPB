from numpy import arange, pi
from pbe.setup.system import Domain, DispersePhase, ContinuosPhase
from pbe.solvers.moc import MOCSolution
from pbe.models import breakup, coalescence, DSD

# TODO: implementar uma classe que ja traz todos os parâmetros constantes
#       do caso, pra não precisar, por exemplo, de obter o objeto
#       continuousphase pra depois calcular epsilon. Além de não precisar
#       salvar em self.D ou L parâmetros que não se encaixam nem na classe solução
#       nem em domain ou phase continua ou dispersa


class MitreSolution:
    def __init__(
        self,
        M=10,  # number of classes
        U=1.10,  # average fluid velocity
        dP=None,  # pressure drop
        phi=0.091,  # [1] holdup
        v0=5e-10,  # [m³]
        model_parameters=None,
        theta=600.0,
        omega=0.05,
    ):
        self.phi = phi
        self.omega = omega
        self.Ve = (
            2 * (0.001 * 0.005 * 0.005)
            + omega * 0.005 * 0.005
            + (0.002 * 0.005 * 0.005)
            + (0.003 * 0.005 * 0.005)
        )  # effetive volume in the domain
        self.theta = self.Ve / Q
        time = arange(0.0, 3600, 0.5)  # Tempo de integração

        # oil
        self.cp = ContinuosPhase(
            name="oil", rho=865.0, mu=1.15e-2  # [kg/m³]
        )  # [P = kg * m^-1 s^-1]  dynamic viscosity

        self.cp.epsilon = calc_epsilon(U, self.D, self.cp.mu, self.cp.rho)

        # water solution
        self.dp = DispersePhase(
            name="water solution",
            phi=self.phi,
            rho=1000.0,  # [kg/m3]
            sigma=2.23e-2,  # [P = kg * m^-1 s^-1]
            v_max=v0 * 3,
            v0=v0,
            sigma0=v0 / 10,
        )

        self.domain = Domain(theta=theta, V=pi * self.L * (self.D / 2) ** 2, M=M)

        self.set_parameters(model_parameters)
        self.n0()
        vmin = None

        # Função distribuição de gotas filhas
        beta = breakup.DDSD.mitre_beta

        # Função frequencia de quebra
        g = breakup.breakupModels(
            name="mitre", C=self.C, domain=self.domain, cp=self.cp, dp=self.dp
        ).gamma

        # Função frequencia de Coalescencia
        Qf = coalescence.coalescenceModels(
            name="mitre", C=self.C, domain=self.domain, cp=self.cp, dp=self.dp
        ).Qf
        
        # Função distribuição inicial
        A0 = DSD.analitico(dp=self.dp).A0

        self.moc = MOCSolution(
            M,
            time,
            self.dp.v_max / M,
            N0=self.N0,
            xi0=vmin,
            beta=beta,
            gamma=g,
            Q=Qf,
            theta=self.domain.theta,
            nf0=self.n0,
            A0=A0,
        )

    def set_parameters(self, model_parameters):
        if model_parameters is None:
            self.C = [0.4, 0.08, 2.8, 1.83e13]
        else:
            self.C = model_parameters

    def N0(self, v):  # total number of drops?
        return 0 * v

    def n0(self):
        self.n0 = self.domain.V / self.domain.theta

    @property
    def pbe_phi(self):
        return self.moc.total_volume / self.domain.V


# Function used by class
def calc_epsilon(U: float, D: float, mu: float, rho: float):
    """Calcula as propriedades turbulentas

    Args:
        U (float): average flow velocity
        D (float):
        mu (float):
        rho (float):

    Returns:
        float: epsilon
    """
    Re = U * D / mu * rho
    TI = 0.16 * Re ** (-1.0 / 8.0)
    u_rms = U * TI
    k = 3.0 / 2.0 * u_rms**2
    L_t = 0.038 * D
    epsilon = 0.09 * k ** (3.0 / 2.0) / L_t

    return epsilon
