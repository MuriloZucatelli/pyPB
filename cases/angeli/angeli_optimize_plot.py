#! .conda\python.exe
import pickle
import time
from numpy import genfromtxt, abs, array, pi, set_printoptions
from scipy.optimize import minimize, differential_evolution
from matplotlib.pyplot import figure
import itertools
import sys
import os.path as path
dir = path.dirname(__file__)
if __name__ == '__main__':
    sys.path.append(path.abspath(path.join(dir, '..\\..')))
from pbe.app.angeli_class import AngeliSolution
from pbe.setup.helpers import set_plt_params
set_printoptions(precision=4)
set_plt_params(relative_fig_width=0.8)


with open(path.join(dir, 'angeli_optimization_results.pickle'), 'rb') as f:
    results = pickle.load(f)
print(results)


fig = figure()
ax = fig.gca()
ax.set_xlabel(r'$N^*$ [rpm]')
ax.set_ylabel(r'$d_{32}$ [mm]')
ax.set_ylim([0, 0.5])

marker = itertools.cycle(('s', 'v', 'o'))
for c in concentrations:
    ax.plot(
        ct_data[c][:, 0], ct_data[c][:, 1],
        marker=next(marker), linestyle='',
        label=r'C\&T $\phi={0:0.2f}$'.format(c / 100.))

for c in concentrations:
    ax.plot(Ns * 60,
            [cts.moc.d32 * 1000 for cts in ct_solutions[c]],
            label=r'Num. $\phi={0:0.2f}$'.format(c / 100.))
    # plt.plot(
    # Ns * 60, [cts.d32 * 1000 for cts in ct_solutions_g[c]],
    # label=r'Num. gen., $\phi={0:0.2f}$'.format(c / 100.))

handles, labels = ax.get_legend_handles_labels()
first_legend = ax.legend(handles[1:], labels[1:], loc='upper right')
ax.add_artist(first_legend)
second_legend = ax.legend(handles[:1], labels[:1], loc='lower left')
fig.patch.set_alpha(0)
fig.savefig(path.join(dir, "ang-fig.pdf"), bbox_inches='tight')
