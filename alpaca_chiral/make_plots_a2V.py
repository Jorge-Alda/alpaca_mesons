import numpy as np
import matplotlib.pyplot as plt
import pathlib
import sys
from time import time

ROOT = pathlib.Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from alpaca_chiral.pseudoscalar_diag import ALPcouplings
from alpaca_chiral.vector_pair import neutral_vector_pair_width

mesons_tex = {
 'omega': r'\omega',
 'phi': r'\phi'
}

mesons = mesons_tex.keys()

def make_plot_cG():
    cG = ALPcouplings(cG=1.0, cuA=np.zeros((3,3)), cdA=np.zeros((3,3)), cuV=np.zeros((3,3)), cdV=np.zeros((3,3)))

    fa = 1000/(16*np.pi**2)

    channels = {}

    ma = np.linspace(0.4, 3.0, 100)

    for i, meson1 in enumerate(mesons):
        for j, meson2 in enumerate(mesons):
            if i <= j:
                channel = f"{mesons_tex[meson1]} {mesons_tex[meson2]}"
                print("Computing channel", channel)
                t0 = time()
                dw = [neutral_vector_pair_width(m, fa, cG, meson1, meson2) for m in ma]
                tf = time()
                if any(dw) and np.max(dw) > 1e-12:
                    channels[channel] = np.array(dw)
                    print('\t' , channel, "computed in", round(tf-t0, 2), "seconds")

    colors = plt.cm.jet(np.linspace(0, 1, len(channels)))

    fig = plt.figure()
    for (channel, color) in zip(channels.keys(), colors):
        plt.semilogy(ma, channels[channel], label=f'${channel}$', color=color)
    total = np.sum([c for c in channels.values()], axis=0)
    plt.semilogy(ma, total, c='k', lw=1.5)
    plt.xlabel(r'$m_a$ [GeV]')
    plt.ylabel(r'$\Gamma(a \to VV)$ [GeV]')
    plt.legend(ncols=3, fontsize=6, borderpad=0.2, labelspacing=0.1, handlelength=1.4)
    plt.ylim(1e-12, 1e-3)
    fig.savefig(pathlib.Path(__file__).parent / 'plots/dw_a2V_cG.pdf', bbox_inches='tight')

if __name__ == '__main__':
    make_plot_cG()
