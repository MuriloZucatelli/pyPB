from numpy import genfromtxt, linspace
import pickle
import os
import sys
import os.path as path

dir = path.dirname(__file__)
if __name__ == "__main__":
    sys.path.append(path.abspath(path.join(dir, "..\\..")))
from pbe.setup.helpers import set_plt_params
from matplotlib.pyplot import figure
import itertools

"""
This script plot figure 3 from CT publication.
"""
dir = os.path.dirname(__file__)
with open(os.path.join(dir, 'ct_solutions.pickle'), 'rb') as f:
    ct_solutions = pickle.load(f)

concentrations = [10]  # , 10, 15]
Ns = linspace(3, 6, 10)  # rps

ct_data = dict([(
    c, genfromtxt(
        './pbe/app/data/coulaloglou/d32_N_alpha{0}.txt'.format(c)))
    for c in concentrations])

set_plt_params(relative_fig_width=0.7)
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
    ax.plot(
        Ns * 60, [cts.moc.d32 * 1000 for cts in ct_solutions[c]],
        label=r'Num. $\phi={0:0.2f}$'.format(c / 100.),
        marker='d')
    # plt.plot(
    # Ns * 60, [cts.d32 * 1000 for cts in ct_solutions_g[c]],
    # label=r'Num. gen., $\phi={0:0.2f}$'.format(c / 100.))

handles, labels = ax.get_legend_handles_labels()
first_legend = ax.legend(handles[1:], labels[1:], loc='upper right')
ax.add_artist(first_legend)
second_legend = ax.legend(handles[:1], labels[:1], loc='lower left')
fig.patch.set_alpha(0)
fig.savefig(os.path.join(dir, "ct-fig3.pdf"), bbox_inches='tight')
