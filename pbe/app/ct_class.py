from numpy import arange
from pbe.setup.system import Domain, DispersePhase, ContinuosPhase
from pbe.solvers.moc import MOCSolution
from pbe.models import breakup, coalescence, DSD


"""
Case setup based on:

    Colaloglou and Tavlarides (1977)
    "Description of Interaction Processes in Agitated Liquid-Liquid
    Dispersions", J. Chem Eng Vol 32

"""


class CTSolution:
    def __init__(
        self,
        M=10,  # number of classes
        Nstar=4.16,  # [rps] impeller revolutions
        phi=0.15,  # [1] holdup, dispersed phase holdup fraction
        vmax=8e-11,  # Maximum droplet volume [m³]
        v0=4e-11,  # [m³]
        model_parameters=None,
    ) -> None:
        self.vmax = vmax
        self.D = 0.1  # [m]  impeller diameter
        # self.D = 10         # [cm] impeller diameter
        self.Nstar = Nstar
        self.phi = phi
        time = arange(0.0, 3600, 0.5)  # Tempo de integração

        self.domain = Domain(
            theta=600,  # nominal residence time, sec (10 min)
            V=12e-3,  # tank volume [m³]
            M=M,
            D=self.D,
            Nstar=self.Nstar,
        )

        # Feed distribution (DSD)
        self.v0 = vmax / 2  # Average volume
        self.sigma0 = (vmax - self.v0) / 3.3  # Standard deviation
        self.nf0 = self.domain.V / self.domain.theta  # number feed rate
        vmin = None  # 5.23e-7                       # Minimum volume

        self.continuousphase = ContinuosPhase(
            name="water",
            rho=1000,  # [kg/m³]
            # rho=1,         # [g/cm³]
            mu=0.89e-3,  # [P = kg * m^-1 s^-1]
            # mu=0.89e-2     # [P = g * cm^-1 s^-1] # dynamic viscosity
            epsilon=Nstar**3 * self.D**2,
        )
        # contProperties['epsilon'] = 0.407 * Nstar**3 * self.D**2

        # Kerosene-dicholorebenzene
        self.dispersephase = DispersePhase(
            name="Kerosene-dicholorebenzene",
            phi=self.phi,
            rho=972.0,  # [kg/m3]
            # rho=0.972,      # [g/cm3]
            sigma=42.82e-3,  # [P = kg.m^-1.s^-1: N/m]
            # sigma=42.82,     # [P = dynes.cm^-1]
            v_max=self.vmax,  # m³
            v0=self.v0,
            # sigma0=v0 / 10.)
            sigma0=self.sigma0,
        )

        # mm3_to_cm3 = 0.1**3            # Conversão mm³ para cm³
        # vmax = vmax * mm3_to_cm3       # Maximum volume
        # vmax = 0.06 * mm3_to_cm3

        if model_parameters is None:
            # C3=2.8e-6, C4=1.83e9
            self.C = [0.4, 0.08, 2.8, 1.83e13]
        else:
            self.C = model_parameters

        # Função distribuição de gotas filhas
        beta = breakup.DDSD.coulaloglou_beta

        # Função frequencia de quebra
        g = breakup.breakupModels(
            name="coulaloglou",
            C=self.C,
            domain=self.domain,
            cp=self.continuousphase,
            dp=self.dispersephase,
        ).gamma

        # Função frequencia de Coalescencia
        Qf = coalescence.coalescenceModels(
            name="coulaloglou",
            C=self.C,
            domain=self.domain,
            cp=self.continuousphase,
            dp=self.dispersephase,
        ).Q

        # Função distribuição inicial
        A0 = DSD.analitico(dp=self.dispersephase).A0

        self.moc = MOCSolution(
            M,
            time,
            vmax / M,
            n0=self.n0,
            xi0=vmin,
            beta=beta,
            gamma=g,
            Q=Qf,
            theta=self.domain.theta,
            nf0=self.nf0,
            A0=A0,
        )

    def n0(self, v):  # initial density function
        return 0 * v

    @property
    def pbe_phi(self):
        return self.moc.total_volume / self.domain.V
