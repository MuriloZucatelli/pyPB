from numpy import (
    arange,
    zeros,
    pi,
    zeros_like,
    dot,
    array,
    where,
    delete,
    eye,
    ones,
    diff,
    fill_diagonal,
)

from numpy import sum as nsum
from scipy.integrate import odeint, quad, solve_ivp
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
        varsigma=None,
        theta=None,
        nf0=None,
        A0=None,
    ):
        """

        Args:
            M: number of classes
            t: time span
            dxi: grid size
            xi: grid
            v: volume classes
            n0: number density function.
                initial distribuition number of droplets per m³ of discrete phase. Defaults to None.
            N0: Number concentration
                number of droplets by volume of disperse phase
            xi0: initial droplet volume. Defaults to None.
            beta: DDSD. Defaults to None.
            gamma: Breakup frequency. Defaults to None.
            Q: Coalescence frequency. Defaults to None.
            theta: mean residence time. Defaults to None.
            varsigma: number of droplets formed of a break of a droplet
            nf0: number feed rate of drops, sec-1. Defaults to None.
            A0: probability density of droplet size in the feeds. Defaults to None.
        """
        self.M = M
        self.t = t

        if dxi is not None:
            self.dxi = dxi
        if dxi is not None and xi is None:
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
        # self.theta = theta  # mean residence time

        # Non-uniform grid
        # a = 2 mantem os extremos com dx na mesma medida que os
        # pontos intermediarios, funciona bem pra coalescencia
        # mas pra quebra necessita de um ajuste em dNdT[0]
        # Não encontrado nenhum bug nem o motivo disso acontecer
        # Esse impacto é visto apenas na quebra, mas não foi encontrado
        # nenhum bug ou cálculo errado na implementação
        a = 2
        if a == 1:
            if xi is not None and dxi is None:
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
                self.xi = xi

        if a == 2:
            if xi is not None and dxi is None:
                v = zeros(self.M + 1)  # Classes de volume v = (x_i-1 + x_i)/2
                dxi = zeros(self.M)
                for i in range(self.M):
                    if i == 0:
                        dxi[i] = 2 * ((xi[0] + xi[1]) / 2 - xi[0])
                        v[0] = max(0, xi[0] - dxi[0])
                        v[1] = (xi[0] + xi[1]) / 2
                    elif i == self.M - 1:
                        dxi[i] = 2 * (xi[-1] - (xi[-1] + xi[-2]) / 2)
                        v[i + 1] = xi[i] + dxi[i]
                    else:
                        dxi[i] = (xi[i + 1] + xi[i]) / 2 - (xi[i] + xi[i - 1]) / 2
                        v[i + 1] = (xi[i] + xi[i + 1]) / 2
                self.dxi = dxi
                self.xi = xi
                self.v = v
            elif xi is not None and dxi is not None:
                self.xi = xi
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

        elif callable(N0):
            N0 = array(
                [N0(self.xi[i]) for i in range(M)]
            )  # initial number concentration

        # plt.plot(self.xi, N0, label=str(self.M))
        # plt.legend()
        if varsigma is None:
            self.varsigma = 2.0 * ones(self.M)  # Binary breakup for all droplets
        elif isinstance(varsigma, float) or isinstance(varsigma, int):
            self.varsigma = varsigma * ones(self.M)  # Equal number for all droplets
        else:
            self.varsigma = varsigma

        def nik1(v, xi, i, k):
            # xi[k] > v, xi[k] = v' (Mother droplet)
            return ((xi[i + 1] - v) / (xi[i + 1] - xi[i])) * beta(v, xi[k], i)

        def nik2(v, xi, i, k):
            return ((v - xi[i - 1]) / (xi[i] - xi[i - 1])) * beta(v, xi[k], i)

        # Kernels setup avaliando a função beta e gamma para cada classe
        if gamma is not None and beta is not None:
            self.gamma = gamma(self.xi)
            self.nik = zeros((self.M, self.M))  # β(v,v′)
            for i in range(0, self.M):
                for k in range(i, self.M):
                    if i != k:  # Pq i ñ tem q ser diferente de M
                        nik = lambda v: nik1(v, self.xi, i, k)
                        self.nik[i, k] += quad(nik, self.xi[i], self.xi[i + 1])[0]
                    if i != 0:
                        nik = lambda v: nik2(v, self.xi, i, k)
                        self.nik[i, k] += quad(nik, self.xi[i - 1], self.xi[i])[0]

            self.nik[0, :] *= 2
            mesh = "geometric"
            if diff(self.xi).std() < 0.5 * diff(self.xi).mean():
                mesh = "uniform"
                # NOTE: Necessário para bater com o resultado analítico do
                # Blatz e Tobolsky. Mas faz errar o resultado de Ziff
                # Isso deve ser feito para quando a malha for uniforme
                # Porem quando a malha é geométrica, isso produz um erro
                fill_diagonal(self.nik, 0)
                for i in range(self.M):
                    if nsum(self.nik[:, i]) > 0:
                        self.nik[:, i] = self.nik[:, i] / nsum(self.nik[:, i])

            elif mesh == "geometric":
                pass
            # print(mesh)
        else:
            self.gamma = None
            self.nik = None
        # Q: coalescence rate
        if Q is not None:
            self.D = 1 - 0.5 * eye(self.M)
            Xjk = array([self.xi + self.xi[i] for i in range(self.M)])
            cond = dict()
            self.condjbk = dict()
            for i in range(self.M):
                if i == 0:
                    condicao_sup = Xjk <= self.xi[i + 1]
                    cond[i] = array(where(condicao_sup))
                elif i == self.M - 1:
                    condicao_inf = Xjk >= self.xi[i - 1]
                    cond[i] = array(where(condicao_inf))
                else:
                    condicao_inf = Xjk >= self.xi[i - 1]
                    condicao_sup = Xjk <= self.xi[i + 1]
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
            self.eta = eta  # NOTE: Matriz 3-Dim grande

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
                    quad(A0, self.xi[i] - self.dxi / 2.0, self.xi[i] + self.dxi / 2.0)[
                        0
                    ]
                    for i in range(self.M)
                ]
            )
        # Solve procedure
        # self.N = odeint(lambda NN, t: self.RHS(t, NN), N0, t)
        # NOTE: 1D PB is a initial boundary problem
        self.N = solve_ivp(
            self.RHS, (t[0], t[-1]), N0, t_eval=t, method="Radau"
        )["y"].T
        # Solução no tempo t sol['y'][:,-1]

    def RHS(self, t, N):
        dNdt = zeros_like(N)

        if self.gamma is not None and self.nik is not None:
            # Death breakup term
            dNdt -= N * self.gamma
            # Birth breakup term
            dNdt += dot(self.varsigma * self.nik, N * self.gamma)

        if self.Q is not None:
            dNdt -= N * dot(self.Q, N)

            Rab = zeros_like(dNdt)
            for i in arange(self.M):
                j, k = self.condjbk[i][0], self.condjbk[i][1]
                Rab[i] = dot(self.D[j, k] * self.eta[i, j, k], self.Q[j, k] * N[j] * N[k])

            dNdt += Rab

        # if self.theta is not None:
        #    dNdt += self.nf0 * self.A0 - N / self.theta

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
        """D32 em micras
        O correto, o de cima não bate com o bettersizer
        """
        return sum(self.N[-1] * self.xi_d**3) / sum(self.N[-1] * self.xi_d**2)
        # (6/pi * sum(NiVi)/sum(N))^(1/3)

    @property
    def d43(self):
        """D43 em micras"""
        return sum(self.N[-1] * self.xi_d**4) / sum(self.N[-1] * self.xi_d**3)

    @property
    def xi_d(self):  # xi: volume, xi_d: diametro
        return (6 / pi * self.xi) ** (1.0 / 3)

    @property
    def phase_fraction(self):
        """Calculate de total volume concentration of discrete phase
            $alpha$ = \\sum $alpha_i = \\sum (N_i v_i)$
        Returns:
            float: total volume concentration
        """
        return nsum(self.N[-1] * self.xi)  # integral de NiVi, N[-1] is last time

    @property
    def initial_phase_fraction(self):
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
