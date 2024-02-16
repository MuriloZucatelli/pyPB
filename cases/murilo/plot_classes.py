from numpy import linspace, array, pi, diff
import sys
import os.path as path
from pandas import read_excel

dir = path.dirname(__file__)
if __name__ == "__main__":
    sys.path.append(path.abspath(path.join(dir, "..\\..")))
from pbe.setup.helpers import plt_config2
import matplotlib.pyplot as plt
import itertools

plt_config2(relative_fig_width=0.7)

# Obtem as classe do Bettersizer e define a malha
# Create mesh
x = read_excel(dir + "/classes.xlsx")
# diameter is in micrometer and volume is in mm³
d, v = x["d"].to_numpy(), x["v [mm³]"].to_numpy()
dxi = diff(v)
xi = v[:-1] + diff(v) / 2
# Remover classes zeros
sl = slice(8, -12)
dxi, xi = dxi[sl], xi[sl]
M = len(xi)

classes = linspace(1, M, M, endpoint=True, dtype=int)  # classes

fig = plt.figure()
ax = fig.gca()
ax.set_xlabel(r"Classe")

# ax.set_ylim([0, 0.5])

marker = itertools.cycle(("s", "v", "o"))
ax.scatter(classes, xi, marker=next(marker), label="xi", color='tab:blue')
ax.set_ylabel(r"Volume médio da classe [mm³]", color='tab:blue')
ax.set_yscale("log")
ax.tick_params(axis='y', labelcolor='tab:blue')

ax2 = ax.twinx()
ax2.scatter(classes, 1e3*(6*xi/pi)**(1/3.0), color='tab:red')
ax2.set_ylabel(r"Diâmetro médio da classe [$\mu$m]", color='tab:red')
ax2.tick_params(axis='y', labelcolor='tab:red')
#ax2.set_yscale("log")
plt.show()

if True:
    ax.grid()
    fig.savefig(
        path.join(dir, "Bettersize_class.pdf"),
        backend="pgf",
        bbox_inches="tight",
        pad_inches=0.05,
    )
