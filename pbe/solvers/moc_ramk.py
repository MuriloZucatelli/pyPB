from numpy import arange, zeros, pi, zeros_like, dot, array, where, delete, eye
from numpy import sum as nsum
from scipy.integrate import odeint, quad
import matplotlib.pyplot as plt
import numpy as np

"""
Method of classes
"""


class MOCSolution:
    """
    Based on Ramkrishma 1996
    """

    def __init__(
        self,
        M,
        t,
        xi=None,
        dxi=None,
        v=None,
        n0=None,
        N0=None,
        xi0=None,
        beta=None,
        gamma=None,
        Q=None,
        nu=None,
        theta=None,
        nf0=None,
        A0=None,
    ):
        """Solution by method of classes based on Coulaloglou & Tavlarides 1977

        Args:
            M (_type_): number of classes
            t (_type_): time span
            dxi (_type_): grid size
            xi (_type_): grid
            v (_type_): volume classes
            n0 (_type_, optional): number density function.
                initial distribuition number of droplets per m³ of discrete phase. Defaults to None.
            N0 (_type_): Number concentration
                number of droplets by volume of disperse phase
            xi0 (_type_, optional): initial droplet volume. Defaults to None.
            beta (_type_, optional): DDSD. Defaults to None.
            gamma (_type_, optional): Breakup frequency. Defaults to None.
            Q (_type_, optional): Coalescence frequency. Defaults to None.
            theta (_type_, optional): mean residence time. Defaults to None.
            nu(_type, optional): number of droplets formed of a break of a droplet
            nf0 (_type_, optional): number feed rate of drops, sec-1. Defaults to None.
            A0 (_type_, optional): probability density of droplet size in the feed . Defaults to None.
        """
        self.M = M
        self.t = t
        if dxi is not None:
            if xi0 is None and isinstance(self.dxi, float):
                self.xi0 = self.dxi  # define o ξ0, volume minimo
            else:
                self.xi0 = xi0

            # Uniform grid
            if isinstance(self.dxi, float) and xi is None:
                self.xi = self.xi0 + self.dxi * arange(
                    self.M
                )  # vetor com as classes vk ξ, dxi é o parametro k

        self.nf0 = nf0
        # inflow and outflow replaced with relaxation to equilibrium
        # process with relaxation time equal to residence time theta
        self.theta = theta  # mean residence time

        # Non-uniform grid
        if xi is not None:
            v = zeros(self.M + 1)  # Classes de volume v = (x_i-1 + x_i)/2
            dxi = zeros(self.M)
            for i in range(self.M):
                if i == 0:
                    dxi[i] = (xi[0] + xi[1]) / 2 - xi[0]
                    v[0] = max(0, xi[0] - dxi[0])
                    v[1] = (xi[0] + xi[1]) / 2
                elif i == self.M - 1:
                    dxi[i] = xi[-1] - (xi[-1] + xi[-2]) / 2
                    v[i + 1] = xi[i] + dxi[i]
                else:
                    dxi[i] = (xi[i + 1] + xi[i]) / 2 - (xi[i] + xi[i - 1]) / 2
                    v[i + 1] = (xi[i] + xi[i + 1]) / 2
            self.dxi = dxi
            self.xi = xi
            self.v = v
        else:
            pass

        # integration of initial number density function to obtain
        # number concentration N0
        # Number of droplets per m³ of discrete phase (number concentration)
        if n0 is None and N0 is None:
            N0 = zeros_like(self.xi)
        elif N0 is None and isinstance(self.dxi, float):
            N0 = array(
                [
                    quad(n0, self.xi[i] - self.dxi / 2.0, self.xi[i] + self.dxi / 2.0)[
                        0
                    ]
                    for i in range(self.M)
                ]
            )  # initial number concentration

        elif N0 is None and v is not None:
            N0 = array(
                [quad(n0, v[i], v[i + 1])[0] for i in range(self.M)]
            )  # initial number concentration

        # plt.plot(self.xi, N0, label=str(self.M))
        # plt.legend()
        if nu is None:
            self.nu = 2.0  # Binary breakup

        def nik1(v, xi, i, k):
            return (xi[i + 1] - v) / (xi[i + 1] - xi[i]) * beta(v, xi[k])

        def nik2(v, xi, i, k):
            return (v - xi[i - 1]) / (xi[i] - xi[i - 1]) * beta(v, xi[k])

        # Kernels setup avaliando a função beta e gama para cada classe
        if gamma is not None and beta is not None:
            self.gamma = gamma(self.xi)
            self.beta = zeros((self.M, self.M))  # β(ξ,ξ′j)
            self.nik = zeros((self.M, self.M))  # β(ξ,ξ′j)
            for i in range(0, self.M):
                for k in range(i, self.M):
                    if i != k:
                        self.nik[i, k] += quad(
                            lambda v: nik1(v, self.xi, i, k), self.xi[i], self.xi[i + 1]
                        )[0]
                    if i != 0:
                        self.nik[i, k] += quad(
                            lambda v: nik2(v, self.xi, i, k), self.xi[i - 1], self.xi[i]
                        )[0]

        else:
            self.gamma = None
            self.nik = None
        # Q: coalescence rate
        if Q is not None:
            self.D = 1 - 0.5 * eye(self.M)
            Xjk = array([xi + xi[i] for i in range(self.M)])
            cond = dict()
            self.condjbk = dict()
            for i in range(self.M):
                if i == 0:
                    condicao_sup = Xjk <= xi[i + 1]
                    cond[i] = array(where(condicao_sup))
                elif i == self.M - 1:
                    condicao_inf = Xjk >= xi[i - 1]
                    cond[i] = array(where(condicao_inf))
                else:
                    condicao_inf = Xjk >= xi[i - 1]
                    condicao_sup = Xjk <= xi[i + 1]
                    cond[i] = array(where(condicao_inf & condicao_sup))

            # Eliminando os pares onde j é menor do que k
            eta = zeros((self.M, self.M, self.M))
            for i in range(self.M):
                # where j is bigger or equal than k
                jbk = where(cond[i][0, :] - cond[i][1, :] < 0)
                self.condjbk[i] = delete(cond[i], jbk, axis=1)
                j, k = self.condjbk[i][0], self.condjbk[i][1]
                vjk = self.xi[j] + self.xi[k]

                for v, j, k in zip(vjk, j, k):
                    if i < self.M - 1 and i > 0:
                        if v >= self.xi[i] and v <= self.xi[i + 1]:
                            eta[i, j, k] = (self.xi[i + 1] - v) / (
                                self.xi[i + 1] - self.xi[i]
                            )
                        elif v >= self.xi[i - 1] and v <= self.xi[i]:
                            eta[i, j, k] = (v - self.xi[i - 1]) / (
                                self.xi[i] - self.xi[i - 1]
                            )
                    elif i == self.M - 1:
                        if v >= self.xi[i]:
                            pass
                        elif v >= self.xi[i - 1] and v <= self.xi[i]:
                            eta[i, j, k] = (v - self.xi[i - 1]) / (
                                self.xi[i] - self.xi[i - 1]
                            )
                    elif i == 0:
                        if v >= self.xi[i] and v <= self.xi[i + 1]:
                            eta[i, j, k] = (self.xi[i + 1] - v) / (
                                self.xi[i + 1] - self.xi[i]
                            )
                        elif v <= self.xi[i]:
                            pass
            self.eta = eta

            self.Q = array(
                [
                    [Q(self.xi[i], self.xi[j]) for j in range(self.M)]
                    for i in range(self.M)
                ]
            )

        else:
            self.Q = None

        if A0 is None:
            self.A0 = None
        else:
            self.A0 = array(
                [
                    quad(A0, self.xi[i] - dxi / 2.0, self.xi[i] + dxi / 2.0)[0]
                    for i in range(self.M)
                ]
            )
        # Solve procedure
        self.N = odeint(lambda NN, t: self.RHS(NN, t), N0, t)

    def RHS(self, N, t):
        dNdt = zeros_like(N)

        if self.gamma is not None and self.nik is not None:
            # Death breakup term
            dNdt -= N * self.gamma
            # Birth breakup term
            dNdt += self.nu * dot(self.nik, N * self.gamma)

        if self.Q is not None:
            dNdt -= N * dot(self.Q, N)

            Rab = zeros_like(dNdt)
            for i in arange(self.M):
                j, k = self.condjbk[i][0], self.condjbk[i][1]
                Rab[i] = dot(
                    self.D[j, k] * self.eta[i, j, k], self.Q[j, k] * N[j] * N[k]
                )

            dNdt += Rab

        if self.theta is not None:
            dNdt += self.nf0 * self.A0 - N / self.theta

        # print('Time = {0:g}'.format(t))
        return dNdt

    @property
    def number_density(self):
        """Calculate de number density n or f in some authors
           n is the number of droplets of a size range from di to di+1 per
           unit volume of discrete phase

           Units: 1 / (m³ * m³)

        Returns:
            float: number density
        """
        return self.N / self.dxi  # N/dV

    @property
    def d32(self):
        return (6 / pi * sum(self.N[-1] * self.xi) / sum(self.N[-1])) ** (1.0 / 3)
        # (6/pi * sum(NiVi)/sum(N))^(1/3)

    @property
    def xi_d(self):  # xi: volume, xi_d: diametro
        return (6 / pi * self.xi) ** (1.0 / 3)

    @property
    def total_volume(self):
        """Calculate de total volume concentration of discrete phase
            $alpha$ = \\sum $alpha_i = \\sum (N_i v_i)$
        Returns:
            float: total volume concentration
        """
        return nsum(self.N[-1] * self.xi)  # integral de NiVi, N[-1] is last time

    @property
    def initial_total_volume(self):
        """Calculate de total volume concentration of discrete phase
            $alpha$ = \\sum $alpha_i = \\sum (N_i v_i)$
        Returns:
            float: total volume concentration
        """
        return nsum(self.N[0] * self.xi)  # integral de NiVi, N[0] is initial time

    @property
    def total_numbers(self):
        """Calculate the total number of droplets per unit volume of discrete phase
           N
           Units: 1 / m³

        Returns:
            float: number concentration
        """
        # return integral de n de 0 ao infinito
        return nsum(self.N, axis=1)  # total Droplets per m³
