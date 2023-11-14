from numpy import arange, zeros, pi, zeros_like, dot, array
from numpy import sum as nsum
from scipy.integrate import odeint, quad
import matplotlib.pyplot as plt
import numpy as np

"""
Method of classes
"""


class MOCSolution:
    """
    Based on Brocks and Hidy uniform discretisation 1970

    """

    def __init__(
        self,
        M,
        t,
        dxi,
        xi=None,
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
            Nclasses (_type_): number of classes
            t (_type_): time span
            dxi (_type_): grid size
            xi (_type_): grid
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
        self.dxi = dxi
        if xi0 is None and isinstance(self.dxi, float):
            self.xi0 = self.dxi  # define o ξ0, volume minimo
        else:
            self.xi0 = xi0

        self.nf0 = nf0
        # inflow and outflow replaced with relaxation to equilibrium
        # process with relaxation time equal to residence time theta
        self.theta = theta  # mean residence time

        # Uniform grid
        if isinstance(self.dxi, float):
            self.xi = self.xi0 + self.dxi * arange(
                self.M
            )  # vetor com as classes vk ξ, dxi é o parametro k

        # Non-uniform grid
        elif xi is not None:
            self.xi = xi
        else:
            pass

        xi_len = len(self.xi)

        # integration of initial number density function to obtain
        # number concentration N0
        # Number of droplets per m³ of discrete phase (number concentration)
        if n0 is None and N0 is None:
            N0 = zeros_like(self.xi)
        elif N0 is None and dxi is not None:
            N0 = array(
                [
                    quad(n0, self.xi[i] - dxi / 2.0, self.xi[i] + dxi / 2.0)[0]
                    for i in range(M)
                ]
            )  # initial number concentration
        print(sum(N0))
        #plt.plot(self.xi, N0, label=str(Nclasses))
        #plt.legend()
        if nu is None:
            self.nu = 2.0  # Binary breakup

        # Kernels setup avaliando a função beta e gama para cada classe
        if gamma is not None and beta is not None:
            self.gamma = gamma(self.xi)
            self.betadxi = zeros((self.M, self.M))  # β(ξ,ξ′j)
            for i in range(1, xi_len):
                for j in range(i):
                    self.betadxi[j, i] = beta(self.xi[j], self.xi[i])
                self.betadxi[:, i] = self.betadxi[:, i] / nsum(
                    self.betadxi[:, i]
                )  # normalizando apenas a coluna
                # betadxi substitui beta * (v_i+1 - vi)

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
                    for i in range(M)
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
            for i in arange(self.M // 2):
                ind = slice(i, self.M - i - 1)
                Cb = self.Q[i, ind] * N[i] * N[ind]
                # Death coalescence term
                Cd[i] += nsum(Cb)
                # Birth coalescence term
                Cd[(i + 1) : (self.M - i - 1)] += Cb[1:]
                Cb[0] = 0.5 * Cb[0]
                dNdt[(2 * i + 1) :] += Cb

            dNdt -= Cd
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
        # TODO xi0 só é igual a dV quando xi0=dxi
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
        return nsum(self.N, axis=1)  # total Droplets per m³
