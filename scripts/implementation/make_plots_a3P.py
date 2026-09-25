import numpy as np
import matplotlib.pyplot as plt
import alpaca
import pathlib
from time import time

from pseudoscalar_diag import cqA, mixing_shift
from three_meson_width import decay_width

mesons_tex = {
 'pi0': r'\pi^0',
 'pi+': r'\pi^+',
 'pi-': r'\pi^-',
 'eta': r'\eta',
 'eta_prime': r'\eta^\prime',
 'K+': r'K^+',
 'K-': r'K^-',
 #'K0': r'K^0',
 #'K0bar': r'\bar{K}^0'
 'KL': r'K_L',
 'KS': r'K_S',
}

mesons = mesons_tex.keys()

def make_plot_cG():
    cG = alpaca.ALPcouplings({'cG': 1.0}, 5.0, 'VA_below')
    cAhat = lambda ma: -mixing_shift(cG, ma)
    cPhat = lambda ma: -mixing_shift(cG, ma) - cqA(cG)

    fa = 1000/(16*np.pi**2)

    channels = {}

    ma = np.linspace(0.4, 3.0, 100)

    for i, meson1 in enumerate(mesons):
        for j, meson2 in enumerate(mesons):
            for k, meson3 in enumerate(mesons):
                if i <= j <= k:
                    channel = f"{mesons_tex[meson1]} {mesons_tex[meson2]} {mesons_tex[meson3]}"
                    print("Computing channel", channel)
                    t0 = time()
                    dw = [decay_width(m, cAhat(m), cPhat(m), meson1, meson2, meson3, fa, include_vector_contributions=True) for m in ma]
                    tf = time()
                    if any(dw):
                        channels[channel] = dw
                        print('\t' , channel, "computed in", round(tf-t0, 2), "seconds")

    colors = plt.cm.jet(np.linspace(0, 1, len(channels)))

    fig = plt.figure()
    for (channel, color) in zip(channels.keys(), colors):
        plt.semilogy(ma, channels[channel], label=f'${channel}$', color=color)
    plt.xlabel(r'$m_a$ [GeV]')
    plt.ylabel(r'$\Gamma(a \to 3\text{mesons})$')
    plt.legend(ncols=3, fontsize=8)
    plt.ylim(1e-12, 1e-3)
    fig.savefig(pathlib.Path(__file__).parent.parent.parent / 'plots/dw_a3P_cG.pdf', bbox_inches='tight')

def check_KSKS():
    cG = alpaca.ALPcouplings({'cG': 1.0}, 5.0, 'VA_below')
    cAhat = lambda ma: -mixing_shift(cG, ma)
    cPhat = lambda ma: -mixing_shift(cG, ma) - cqA(cG)

    fa = 1000/(16*np.pi**2)
    ma = 2.9

    t0 = time()
    dw = decay_width(ma, cAhat(ma), cPhat(ma), 'pi0', 'KS', 'KS', fa, include_vector_contributions=False)
    tf = time()
    print("Decay width for a -> pi0 KS KS, no vector contribution:", dw, "computed in", f'{tf-t0:.2f}', "seconds")

    t0 = time()
    dw = decay_width(ma, cAhat(ma), cPhat(ma), 'pi0', 'KS', 'KS', fa, include_vector_contributions=True)
    tf = time()
    print("Decay width for a -> pi0 KS KS, with vector contribution:", dw, "computed in", f'{tf-t0:.2f}', "seconds")

    t0 = time()
    dw = decay_width(ma, cAhat(ma), cPhat(ma), 'pi0', 'K0', 'K0bar', fa, include_vector_contributions=False)
    tf = time()
    print("Decay width for a -> pi0 K0 K0bar, no vector contribution:", dw, "computed in", f'{tf-t0:.2f}', "seconds")

    t0 = time()
    dw = decay_width(ma, cAhat(ma), cPhat(ma), 'pi0', 'K0', 'K0bar', fa, include_vector_contributions=True)
    tf = time()
    print("Decay width for a -> pi0 K0 K0bar, with vector contribution:", dw, "computed in", f'{tf-t0:.2f}', "seconds")

if __name__ == '__main__':
    #make_plot_cG()
    check_KSKS()