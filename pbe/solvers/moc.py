from numpy import arange, zeros, pi, zeros_like, dot, array
from numpy import sum as nsum
from scipy.integrate import odeint, quad

"""
Method of classes
"""


class MOCSolution:
    """
    Based on Brooks and Hidy uniform discretisation

    """

    def __init__(
        self,
        number_of_classes,
        t,
        dxi,
        N0=None,
        xi0=None,
        beta=None,
        gamma=None,
        Q=None,
        theta=None,
        n0=None,
        A0=None,
    ):
        """Solution by method of classes

        Args:
            number_of_classes (_type_): number of classes
            t (_type_): time span
            dxi (_type_): grid size
            N0 (_type_, optional): _description_. Defaults to None.
            xi0 (_type_, optional): _description_. Defaults to None.
            beta (_type_, optional): DDSD. Defaults to None.
            gamma (_type_, optional): Breakup frequency. Defaults to None.
            Q (_type_, optional): Coalescence frequency. Defaults to None.
            theta (_type_, optional): mean residence time. Defaults to None.
            n0 (_type_, optional): _description_. Defaults to None.
            A0 (_type_, optional): _description_. Defaults to None.
        """
        self.number_of_classes = number_of_classes
        if xi0 is None:
            self.xi0 = dxi  # define o ξ0, volume minimo
        else:
            self.xi0 = xi0
        self.n0 = n0  # ??
        # inflow and outflow replaced with relaxation to equilibrium
        # process with relaxation time equal to residence time theta
        self.theta = theta  # mean residence time
        # Uniform grid
        self.xi = self.xi0 + dxi * arange(
            self.number_of_classes
        )  # vetor com as classes vk ξ, dxi é o parametro k
        xi_len = len(self.xi)
        if N0 is None:
            N0 = zeros_like(self.xi)
        else:
            N0 = array(
                [
                    quad(N0, self.xi[i] - dxi / 2.0, self.xi[i] + dxi / 2.0)[0]
                    for i in range(number_of_classes)
                ]
            )  # number concentration???  integração para o delta de dirac

        self.nu = 2.0  # Binary breakup
        # Kernels setup avaliando a função beta e gama para cada classe
        if gamma is not None and beta is not None:
            self.gamma = gamma(self.xi)
            self.betadxi = zeros(
                (self.number_of_classes, self.number_of_classes)
            )  # β(ξ,ξ′j)
            for i in range(1, xi_len):
                for j in range(i):
                    self.betadxi[j, i] = beta(self.xi[j], self.xi[i])
                self.betadxi[:, i] = self.betadxi[:, i] / nsum(
                    self.betadxi[:, i]
                )  # normalizando apenas a coluna

        else:
            self.gamma = None
            self.betadxi = None
        # Q: coalescence rate
        if Q is not None:
            self.Q = array(
                [
                    [Q(self.xi[i], self.xi[j]) for j in range(xi_len)]
                    for i in range(xi_len)
                ]
            )
            # self.Q = zeros((self.number_of_classes, self.number_of_classes))
            # for i in range(len(self.xi)):
            #    for j in range(len(self.xi)):
            #        self.Q[i, j] = Q(self.xi[i], self.xi[j])  #
        else:
            self.Q = None

        if A0 is None:
            self.A0 = None
        else:
            self.A0 = array(
                [
                    quad(A0, self.xi[i] - dxi / 2.0, self.xi[i] + dxi / 2.0)[0]
                    for i in range(number_of_classes)
                ]
            )
        # Solve procedure
        self.N = odeint(lambda NN, t: self.RHS(NN, t), N0, t)

    def RHS(self, N, t):
        dNdt = zeros_like(N)

        if self.gamma is not None and self.betadxi is not None:
            # Death breakup term
            dNdt[1:] -= N[1:] * self.gamma[1:]
            # Birth breakup term
            dNdt[:-1] += self.nu * dot(self.betadxi[:-1, 1:], N[1:] * self.gamma[1:])

        if self.Q is not None:
            Cd = zeros_like(dNdt)
            for i in arange(self.number_of_classes // 2):
                ind = slice(i, self.number_of_classes - i - 1)
                Cb = self.Q[i, ind] * N[i] * N[ind]
                # Death coalescence term
                Cd[i] += nsum(Cb)
                # Birth coalescence term
                Cd[(i + 1): (self.number_of_classes - i - 1)] += Cb[1:]
                Cb[0] = 0.5 * Cb[0]
                dNdt[(2 * i + 1):] += Cb

            dNdt -= Cd
        if self.theta is not None:
            dNdt += self.n0 * self.A0 - N / self.theta  #dNdt-= (N - self.N0) / self.theta??

        # print('Time = {0:g}'.format(t))
        return dNdt

    @property
    def number_density(self):
        return self.N / self.xi0  # N/dV, TODO xi0 só é igual a dV quando xi0=dxi

    @property
    def d32(self):
        return (6 / pi * sum(self.N[-1] * self.xi) / sum(self.N[-1])) ** (1.0 / 3)
        # (6/pi * sum(NiVi)/sum(N))^(1/3)

    @property
    def xi_d(self):  # xi: volume, xi_d: diametro
        return (6 / pi * self.xi) ** (1.0 / 3)

    @property
    def total_volume(self):
        return nsum(self.N[-1] * self.xi)  # integral de NiVi, N[-1]?

    @property
    def total_numbers(self):
        return nsum(self.N, axis=1)  # total Droplets per m³
