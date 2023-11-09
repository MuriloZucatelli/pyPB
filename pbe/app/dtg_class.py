import pandas as pd

# from .classes import FLOW
from numpy import arange, pi
from pbe.setup.system import Domain, DispersePhase, ContinuosPhase, FLOW
from pbe.solvers.moc import MOCSolution
from pbe.models import breakup, coalescence
from os.path import join
from pathlib import Path
from os import environ

# TODO: implementar uma classe que ja traz todos os parâmetros constantes
#       do caso, pra não precisar, por exemplo, de obter o objeto
#       continuousphase pra depois calcular epsilon. Além de não precisar
#       salvar em self.D ou L parâmetros que não se encaixam nem na classe solução
#       nem em domain ou phase continua ou dispersa


class DTGSolution:
    def __init__(
        self,
        M=10,  # number of classes
        U=1.10,  # average fluid velocity
        phi=0.1,  # [1] holdup
        v0=5e-10,  # [m³]
        model_parameters=None,
        theta=600.0,
    ):
        self.D = 0.02095  # [m] pipe diameter diameter
        self.L = 2.5  # [m] length diameter
        self.phi = phi
        time = arange(0.0, 10, 0.01)  # Tempo de integração

        # oil
        self.cp = ContinuosPhase(
            name="oil", rho=873.0, mu=1.21e-1  # [kg/m³]
        )  # [P = kg * m^-1 s^-1]  dynamic viscosity

        self.cp.epsilon = calc_epsilon(U, self.D, self.cp.mu, self.cp.rho)

        # water solution
        self.dp = DispersePhase(
            name="water solution",
            phi=self.phi,
            rho=1000.0,  # [kg/m3]
            sigma=4.5e-3,  # [P = kg * m^-1 s^-1]
            # v_max=v0 * 3,
            # v0=v0,
            # sigma0=v0 / 10,
        )

        self.domain = Domain(theta=theta, V=pi * self.L * (self.D / 2) ** 2, M=M)

        self.set_parameters(model_parameters)
        self.nf0()
        vmin = None

        # Função distribuição de gotas filhas
        beta = breakup.DDSD.coulaloglou_beta

        # Função frequencia de quebra
        g = breakup.breakupModels(
            name="coulaloglou", C=self.C, domain=self.domain, cp=self.cp, dp=self.dp
        ).gamma
        # Função frequencia de Coalescencia
        Qf = coalescence.coalescenceModels(
            name="coulaloglou", C=self.C, domain=self.domain, cp=self.cp, dp=self.dp
        ).Q

        # Initial probability density function distribuition
        # A0 = DSD.analitico(dp=self.dp).A0

        self.moc = MOCSolution(
            M,
            time,
            self.dp.v_max / M,
            n0=self.n0,
            xi0=vmin,
            beta=beta,
            gamma=g,
            Q=Qf,
            theta=self.domain.theta,
        )

    def set_parameters(self, model_parameters):
        if model_parameters is None:
            self.C = [0.4, 0.08, 2.8, 1.83e13]
        else:
            self.C = model_parameters

    def n0(self, v):
        return 0 * v

    def nf0(self):
        self.nf0 = self.domain.V / self.domain.theta

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


"""
    Classes to import flow sensor and DSD measures
"""


class import_flow_DSD:
    def __init__(self, folder: str) -> None:
        self.folder = folder
        f1 = join(folder, "flow")
        # Le os dados externos (Tipo,Weber,Concentração,...)
        df1 = pd.read_excel(f1, sheet_name="flow_ext")
        # Le os dados de media e desvio padrao (Tipo,Weber,Concentração,...)
        df2 = pd.read_excel(f1, sheet_name="flow_media")
        df3 = pd.read_excel(f1, sheet_name="flow_std")
        # le os dados do circuito (P,T,V,...)
        df4 = pd.read_excel(f1, sheet_name="flow")

        self.flow = FLOW(flow=df4, flow_ext=df1, flow_mean=df2, flow_std=df3)

        f2 = join(folder, "DTG")
        self.DTG = pd.read_excel(f2, sheet_name="DTG")


class DTG_experiment:
    def __init__(self, d32=0.32e-03, phi=0.117, U=2.72, theta=600.0):
        self.phi = phi
        self.theta = theta
        self.U = U
        self.d32 = d32

    def __repr__(self) -> str:
        phi = "{:.{}f}".format(100 * self.phi, 2)
        d32 = "{:.{}e}".format(1000 * self.d32, 2)
        return (
            f"concentration {phi}%  |  flow velocity {self.U} m/s  |  "
            + f"residence time {self.theta} s  |  SMD diameter {d32} mm"
        )


def get_location(folder, os_envi: str = "USERPROFILE"):
    home = Path.home()
    path = join(home, folder)
    path = join(environ['USERPROFILE'], folder)
    return path
