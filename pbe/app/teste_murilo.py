from numpy import arange, sqrt, exp, pi


class MOCSolution:
    def __init__(self, number_of_classes, t, dxi, N0=None, xi0=None,
                 beta=None, gamma=None, Q=None,
                 theta=None, n0=None, A0=None):
        self.number_of_classes = number_of_classes
        self.A0


class CTSolution(MOCSolution):
    def __init__(
            self,
            M=10,
            phi=0.15,
            v0=4e-11,
            vmax=0.08,
            model_parameters=None):
        self.D = 0.10  # [m] impeller diameter
        self.phi = phi
        self.v0 = v0
        self.sigma0 = (vmax - self.v0) / 3.3
        time = 0
        MOCSolution(M, time, vmax / M, A0=self.A0)

    def A0(self, v):
        return \
            self.phi / (self.v0 * self.sigma0 * sqrt(2 * pi)) * \
            exp(-(v - self.v0)**2 / (2 * self.sigma0**2))  # Gaussian distribuition

print(1)