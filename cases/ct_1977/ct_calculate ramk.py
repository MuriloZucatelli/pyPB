import sys
import os.path as path
from numpy import genfromtxt, linspace, array
import os
import itertools
dir = path.dirname(__file__)
if __name__ == "__main__":
    sys.path.append(path.abspath(path.join(dir, "..\\..")))
from pbe.app.ct_class_ramk import CTSolution
from pbe.setup.helpers import set_plt_params
from matplotlib.pyplot import figure

"""
This script calculated solutions to PBE problem that are necessary to reproduce
figure 3 from CT publication.
"""

concentrations = [5, 10, 15]  # , 10, 15]
Ns = linspace(3, 6, 5)  # rps
#Ns = array([3])
# Ns = [5.16]
ct_solutions = dict([(
    c, [CTSolution(M=40, Nstar=N, phi=c / 100.0) for N in Ns])
    for c in concentrations])


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
fig.savefig(os.path.join(dir, "ct-d32_ramk_t4.pdf"), bbox_inches='tight')
