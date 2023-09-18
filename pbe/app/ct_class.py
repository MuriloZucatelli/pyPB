from pbe.setup.case_class_copy import CaseSolution, Domain, DispersePhase, ContinuosPhase


class CTSolution():
    def __init__(
            self,
            M=10,
            Nstar=4.16,  # [rps] impeller revolutions
            phi=0.15,  # [1] holdup
            v0=4e-11,  # [cm^3]
            model_parameters=None):
        self.D = 0.10  # [m] impeller diameter

        # Water
        cont = ContinuosPhase(
            name='water',
            rho=1000,      # [kg/cm3]
            mu=0.89e-3,    # [P = kg * m^-1 s^-1]
            epsilon=Nstar**3 * self.D**2)
        # contProperties['epsilon'] = 0.407 * Nstar**3 * self.D**2

        # Kerosene-dicholorebenzene
        dispersion = DispersePhase(
            name='Kerosene-dicholorebenzene',
            phi=phi, rho=972.,  # [kg/m3]
            sigma=42.82e-3,  # [P = kg * m^-1 s^-1]
            v_max=6e-11,
            v0=v0,
            sigma0=v0 / 10.)

        domain = Domain(theta=600, V=12e-3, M=M)

        CaseSolution(dispersion, cont, domain,
            model_parameters=model_parameters)
