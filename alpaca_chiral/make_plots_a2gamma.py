import numpy as np
import matplotlib.pyplot as plt
import pathlib
import sys
from time import time

ROOT = pathlib.Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from alpaca_chiral.pseudoscalar_diag import ALPcouplings
from alpaca_chiral.two_photon import two_photon_width, photon_coupling


def make_plot_cG():
    cG = ALPcouplings(cG=1.0, cuA=np.zeros((3,3)), cdA=np.zeros((3,3)), cuV=np.zeros((3,3)), cdV=np.zeros((3,3)))

    fa = 1000/(16*np.pi**2)

    channels = {}

    ma = np.linspace(0.01, 3.0, 150)

    print("Computing channel a -> gamma gamma")
    t0 = time()
    dw = [two_photon_width(m, fa, cG) for m in ma]
    c_eff = [photon_coupling(m, fa, cG) for m in ma]
    tf = time()
    if any(dw) and np.max(dw) > 1e-12:
        channels["a \\to \\gamma \\gamma"] = np.array(dw)
        print('\t' , "a -> gamma gamma", "computed in", round(tf-t0, 2), "seconds")

    colors = ['red']

    fig = plt.figure()
    for (channel, color) in zip(channels.keys(), colors):
        plt.semilogy(ma, channels[channel], label=f'${channel}$', color=color)
    total = np.sum([c for c in channels.values()], axis=0)
    #plt.semilogy(ma, total, c='k', lw=1.5)
    plt.xlabel(r'$m_a$ [GeV]')
    plt.ylabel(r'$\Gamma(a \to \gamma\gamma)$ [GeV]')
    plt.legend(ncols=3, fontsize=6, borderpad=0.2, labelspacing=0.1, handlelength=1.4)
    plt.ylim(1e-12, 1e-3)
    fig.savefig(pathlib.Path(__file__).parent / 'plots/dw_a2gamma_cG.pdf', bbox_inches='tight')

    fig_c = plt.figure()
    plt.plot(ma, np.real(c_eff), label=r'$c_{\gamma\gamma}^{\rm eff}$', color='blue')
    plt.xlabel(r'$m_a$ [GeV]')
    plt.ylabel(r'$c_{\gamma\gamma}^{\rm eff}$')
    plt.ylim(-5, 5)
    fig_c.savefig(pathlib.Path(__file__).parent / 'plots/c_eff_a2gamma_cG.pdf', bbox_inches='tight')

    fig_c_log = plt.figure()
    plt.semilogy(ma, np.abs(c_eff), label=r'$c_{\gamma\gamma}^{\rm eff}$', color='blue')
    plt.xlabel(r'$m_a$ [GeV]')
    plt.ylabel(r'$|c_{\gamma\gamma}^{\rm eff}|$')
    plt.ylim(1e-2, 1e2)
    fig_c_log.savefig(pathlib.Path(__file__).parent / 'plots/c_eff_a2gamma_cG_log.pdf', bbox_inches='tight')

if __name__ == '__main__':
    make_plot_cG()
