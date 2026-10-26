import numpy as np
import matplotlib.pyplot as plt
import pathlib
import sys
from time import time

ROOT = pathlib.Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from alpaca_chiral.pseudoscalar_diag import ALPcouplings
from alpaca_chiral.vector_two_pseudoscalar import vector_two_pseudoscalar_width
from alpaca_chiral.photon_two_pseudoscalar import photon_two_pseudoscalar_width

mesons_tex = {
 'pi+': r'\pi^+',
 'pi-': r'\pi^-',
 'K+': r'K^+',
 'K-': r'K^-',
 #'K0': r'K^0',
 #'K0bar': r'\bar{K}^0'
 'KL': r'K_L',
 'KS': r'K_S',
}

vectors_tex = {
 'omega': r'\omega',
 'phi': r'\phi'
}

pions = ['pi+', 'pi-']
kaons = ['K+', "K-", 'K0', 'K0bar', 'KL', 'KS']

groups = {
    r'\omega \pi \pi': (0, 2, 0),
    r'\omega K K': (0, 0, 2),
    r'\phi \pi \pi': (1, 2, 0),
    r'\phi K K': (1, 0, 2),
    r'\gamma \pi \pi': (2, 2, 0),
    r'\gamma K K': (2, 0, 2),
}

groups_lookup = {v: k for k, v in groups.items()}

mesons = mesons_tex.keys()

def make_plot_cG():
    cG = ALPcouplings(cG=1.0, cuA=np.zeros((3,3)), cdA=np.zeros((3,3)), cuV=np.zeros((3,3)), cdV=np.zeros((3,3)))

    fa = 1000/(16*np.pi**2)

    channels = {}

    ma = np.linspace(0.4, 3.0, 100)

    channels_grouped = {k: np.zeros_like(ma) for k in groups.keys()}

    for i, mesonV in enumerate(vectors_tex.keys()):
        for j, mesonP1 in enumerate(mesons):
            for k, mesonP2 in enumerate(mesons):
                if i <= j <= k:
                    channel = f"{vectors_tex[mesonV]} {mesons_tex[mesonP1]} {mesons_tex[mesonP2]}"
                    print("Computing channel", channel)
                    t0 = time()
                    dw = [vector_two_pseudoscalar_width(m, fa, cG, mesonV, mesonP1, mesonP2) for m in ma]
                    tf = time()
                    if any(dw) and np.max(dw) > 1e-12:
                        non_zero = np.where(np.array(dw) > 1e-12*np.max(dw))[0]
                        n_pions = sum(1 for m in [mesonP1, mesonP2] if m in pions)
                        n_kaons = sum(1 for m in [mesonP1, mesonP2] if m in kaons)
                        particles = groups_lookup[(i, n_pions, n_kaons)]
                        channels_grouped[particles] += np.array(dw)
                        print('\t' , channel, "computed in", round(tf-t0, 2), "seconds. Non-zero points:", len(non_zero), ". Points per second:", len(non_zero)/(tf-t0))

    for j, mesonP1 in enumerate(mesons):
        for k, mesonP2 in enumerate(mesons):
            if j <= k:
                channel = f"\\gamma {mesons_tex[mesonP1]} {mesons_tex[mesonP2]}"
                print("Computing channel", channel)
                t0 = time()
                dw = [photon_two_pseudoscalar_width(m, fa, cG, mesonP1, mesonP2) for m in ma]
                tf = time()
                if any(dw) and np.max(dw) > 1e-12:
                    non_zero = np.where(np.array(dw) > 1e-12*np.max(dw))[0]
                    n_pions = sum(1 for m in [mesonP1, mesonP2] if m in pions)
                    n_kaons = sum(1 for m in [mesonP1, mesonP2] if m in kaons)
                    particles = groups_lookup[(2, n_pions, n_kaons)]
                    channels_grouped[particles] += np.array(dw)
                    print('\t' , channel, "computed in", round(tf-t0, 2), "seconds. Non-zero points:", len(non_zero), ". Points per second:", len(non_zero)/(tf-t0))

    colors = plt.cm.jet(np.linspace(0, 1, len(channels_grouped)))

    fig = plt.figure()
    for (channel, color) in zip(channels_grouped.keys(), colors):
        if np.max(channels_grouped[channel]) < 1e-12:
            continue
        print(channel, "max width:", np.max(channels_grouped[channel]))
        plt.semilogy(ma, channels_grouped[channel], label=f'${channel}$', color=color)
    plt.xlabel(r'$m_a$ [GeV]')
    plt.ylabel(r'$\Gamma(a \to V P P/\gamma P P)$ [GeV]')
    plt.legend(ncols=3, fontsize=6, borderpad=0.2, labelspacing=0.1, handlelength=1.4)
    plt.ylim(1e-12, 1e-3)
    fig.savefig(pathlib.Path(__file__).parent / 'plots/dw_aVPP_cG_WZW.pdf', bbox_inches='tight')

if __name__ == '__main__':
    make_plot_cG()
