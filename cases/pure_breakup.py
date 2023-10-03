from numpy import arange, linspace, array, piecewise
from itertools import cycle
import sys
import os.path as path

dir = path.dirname(__file__)
if __name__ == "__main__":
    sys.path.append(path.abspath(path.join(dir, "..\\")))
from pbe.solvers.moc import MOCSolution
from tests.test_moc import ziff_total_number_solution, ziff_pbe_solution
import matplotlib.pyplot as plt

"""
Case setup based on:

    Ziff, R.M and McGrady, E.D
    "New solutions to the fragmentation equation", J. Phys A, vol. 18, 1985

    We are looking at case 4. from their paper with kernel F(x, y) = x + y.
    This corresponds to choosing:
        beta = 1.0/y,
        Gamma = y^2,
    for our kernels.
"""

grids = [10, 20, 40, 80, 160]  # Número de classes utilizadas na discretização
time = arange(0.0, 10.0, 0.001)  # Tempo e passo
vmax = 1.0  # Volume máximo
v0 = 1.0
pbe_solutions = dict()

for g in grids:
    threshold = vmax / g / 2

    # This is modelling Dirac's delta
    # initial number density function
    def n0_init(x):
        return piecewise(
            x, [x < v0 - threshold, x >= v0 - threshold], [0, g / vmax]
        )  # 0 ou g/vmax   função delta de dirac

    def n0_init2(x):
        return piecewise(
            x, [abs(v0 - x) < threshold, v0 - x >= threshold], [g / vmax, 0]
        )  # 0 ou g/vmax   função delta de dirac generica

    pbe_solutions[g] = MOCSolution(
        g,  # number of classes
        time,
        vmax / g,  # dxi
        n0=n0_init2,
        beta=lambda x, y: 1.0 / y,  # DDSD   era 2.0 / y
        gamma=lambda x: x**2,  # breakage rate
    )

totals = dict((n, pbe_solutions[n].total_numbers) for n in pbe_solutions)
volume = dict((n, pbe_solutions[n].total_volume) for n in pbe_solutions)
print(volume)
v = linspace(0, vmax, 100)
Na = array([ziff_total_number_solution(v, t, vmax) for t in time])

fig = plt.figure()
ax = fig.gca()
linestyles = cycle(["-", "--", ":"])

# Iterando sobre cada resultado para o numero de classes n
for n in sorted(totals):
    ax.loglog(
        time,
        totals[n]
        / totals[n][0],  # Total numbers of droplets / Initial total n drop. t=0
        linestyle=next(linestyles),
        label="MOC with N={0}".format(n),
    )
ax.loglog(time, Na, "-k", linewidth=2, label="Analytical")
ax.legend(loc="lower right", shadow=True)
ax.set_xlabel("t")
ax.set_ylabel("N/N0")
plt.show()


fig = plt.figure()
ax = fig.gca()
markers = cycle(["o", "s", "v", "*", ".", ","])
v = linspace(0, vmax, 10000)

for n in sorted(pbe_solutions):
    ax.semilogy(
        pbe_solutions[n].xi,  # volume pivot
        pbe_solutions[n].number_density[-1],  # last time
        marker=next(markers),
        label="MOC with N={0}".format(n),
    )

ax.semilogy(
    v,
    ziff_pbe_solution(v, time[-1], vmax),  # analytical in last time
    "-k",
    linewidth=2,
    label="Analytical $t=\infty$",
)

ax.legend(loc="lower left", shadow=True)
ax.set_xlabel("Volume")
ax.set_ylabel("Number density function")
plt.show()
