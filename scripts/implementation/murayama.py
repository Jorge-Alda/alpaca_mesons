import matplotlib
from pathlib import Path
matplotlib.rc_file(Path(__file__).parent / "matplotlibrc")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, LogFormatterMathtext
import numpy as np
if __package__:
    from .only_up_decays import spectrum, two_photon_width, three_pseudoscalar_width, photon_two_pseudoscalar_width, sigma_two_pseudoscalar_width
    from .only_up_inputs import only_up_inputs
else:
    from only_up_decays import spectrum, two_photon_width, three_pseudoscalar_width, photon_two_pseudoscalar_width, sigma_two_pseudoscalar_width
    from only_up_inputs import only_up_inputs


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
 'K_L': r'K_L',
 'K_S': r'K_S',
}

pions = ['pi0', 'pi+', 'pi-']
kaons = ['K+', "K-", 'K0', 'K0bar', 'K_L', 'K_S']
eta = ['eta']
etap = ['eta_prime']

mesons = mesons_tex.keys()

def masses(quartic_coupling=1e-3):
    kappa_u = np.logspace(-3, 0, 500)
    m_pi0 = []
    m_eta = []
    m_etap = []
    m_alp = []
    m_sigma = []
    for k in kappa_u:
        inputs = only_up_inputs(k)
        masses = spectrum(inputs)
        m_pi0.append(masses['pi0'])
        m_eta.append(masses['eta'])
        m_etap.append(masses['eta_prime'])
        m_alp.append(masses['a'])
        f_a, B_0, f_pi = inputs['f_a'], inputs['B_0'], inputs['f_pi']
        up_mass = inputs['hat_S'][0][0, 0]
        sigma_squared = B_0*up_mass*f_pi**2/f_a**2 + 2*quartic_coupling*f_a**2
        m_sigma.append(np.sqrt(sigma_squared))
    plt.figure()
    plt.loglog(kappa_u, m_alp, label='$a$')
    plt.loglog(kappa_u, m_pi0, label=r'$\pi^0$')
    plt.loglog(kappa_u, m_eta, label=r'$\eta$')
    plt.loglog(kappa_u, m_etap, label=r'$\eta^\prime$')
    plt.loglog(kappa_u, m_sigma, label=fr'$\sigma, \lambda = {quartic_coupling:.2e}$')
    plt.legend()
    plt.xlim(1e-3, 1)
    plt.ylim(1e-3, 10)
    plt.xlabel(r'$\kappa_u$')
    plt.ylabel(r'$m$ [GeV]')
    plt.tight_layout(pad=0.5)
    plt.savefig('plots/murayama_spectrum.pdf')
    return m_alp, m_sigma

def alp_mass_axis(axis, kappa_u, alp_masses):
    """Add an ALP mass axis, restricting tick locations to the sampled range (GeV)."""
    alp_masses = np.asarray(alp_masses)
    def forward(values):
        return np.interp(values, kappa_u, alp_masses)
    def inverse(values):
        return np.interp(values, alp_masses, kappa_u)
    mass_axis = axis.secondary_xaxis('top', functions=(forward, inverse))
    minimum_mass, maximum_mass = alp_masses[0], alp_masses[-1]
    decades = np.arange(np.floor(np.log10(minimum_mass)), np.ceil(np.log10(maximum_mass)) + 1)
    major_ticks = np.array([10.**power for power in decades if minimum_mass <= 10.**power <= maximum_mass])
    minor_ticks = np.array([multiplier * 10.**power for power in decades for multiplier in range(2, 10)
                            if minimum_mass <= multiplier * 10.**power <= maximum_mass])
    mass_axis.xaxis.set_major_locator(FixedLocator(major_ticks))
    mass_axis.xaxis.set_minor_locator(FixedLocator(minor_ticks))
    mass_axis.xaxis.set_major_formatter(LogFormatterMathtext())
    mass_axis.set_xlabel(r'$m_a$ [GeV]')
    return mass_axis

def gammagamma(m_alp):
    kappa_u = np.logspace(-3, 0, 500)
    Gamma = []
    for k in kappa_u:
        inputs = only_up_inputs(k)
        fa = inputs['f_a']
        G = two_photon_width(inputs)
        Gamma.append(G)
    plt.figure()
    plt.loglog(kappa_u, Gamma, label=r'$a\to \gamma \gamma$ [WZW only]')
    plt.xlabel(r'$\kappa_u$')
    plt.ylabel(r'$\Gamma$ [GeV]')
    plt.xlim(1e-3, 1)
    alp_mass_axis(plt.gca(), kappa_u, m_alp)
    plt.legend()
    plt.tight_layout(pad=0.1)
    plt.savefig('plots/murayama_photons.pdf')

def three_mesons(m_alp, quartic_coupling=1e-3):
    kappa_u = np.logspace(-3, 0, 500)

    groups = {
        r'3\pi': (3, 0, 0, 0),
        r'\eta \pi \pi': (2, 0, 1, 0),
        r'\eta^\prime \pi \pi': (2, 0, 0, 1),
        r'\eta \eta \pi': (1, 0, 2, 0),
        r'\eta \eta^\prime \pi': (1, 0, 1, 1),
        r'\eta^\prime \eta^\prime \pi': (1, 0, 0, 2),
        r'3\eta': (0, 0, 3, 0),
        r'\eta \eta \eta^\prime': (0, 0, 2, 1),
        r'\eta \eta^\prime \eta^\prime': (0, 0, 1, 2),
        r'3\eta^\prime': (0, 0, 0, 3),
        r'K K \pi': (1, 2, 0, 0),
        r'K K \eta': (0, 2, 1, 0),
        r'K K \eta^\prime': (0, 2, 0, 1)
    }

    groups_lookup = {v: k for k, v in groups.items()}

    channels_grouped = {k: np.zeros_like(kappa_u) for k in groups.keys()}

    for i, meson1 in enumerate(mesons):
        for j, meson2 in enumerate(mesons):
            for k, meson3 in enumerate(mesons):
                if i <= j <= k:
                    channel = f"{mesons_tex[meson1]} {mesons_tex[meson2]} {mesons_tex[meson3]}"
                    print("Computing channel", channel)

                    dw = []
                    for k in kappa_u:
                        inputs = only_up_inputs(k)
                        fa = inputs['f_a']
                        f_a, B_0, f_pi = inputs['f_a'], inputs['B_0'], inputs['f_pi']
                        up_mass = inputs['hat_S'][0][0, 0]
                        sigma_squared = B_0*up_mass*f_pi**2/f_a**2 + 2*quartic_coupling*f_a**2
                        dw.append(three_pseudoscalar_width(inputs, [meson1, meson2, meson3], sigma_squared) )
                    if any(dw) and np.max(dw) > 1e-12:
                        n_pions = sum(1 for m in [meson1, meson2, meson3] if m in pions)
                        n_kaons = sum(1 for m in [meson1, meson2, meson3] if m in kaons)
                        n_eta = sum(1 for m in [meson1, meson2, meson3] if m == 'eta')
                        n_etap = sum(1 for m in [meson1, meson2, meson3] if m == 'eta_prime')
                        particles = groups_lookup[(n_pions, n_kaons, n_eta, n_etap)]
                        channels_grouped[particles] += np.array(dw)

    colors = plt.cm.jet(np.linspace(0, 1, len(channels_grouped)))

    plt.figure()
    for (channel, color) in zip(channels_grouped.keys(), colors):
        plt.loglog(kappa_u, channels_grouped[channel], label=f'${channel}$', color=color)
    plt.title(rf'LO-$\chi$PT only, $\lambda={quartic_coupling:.2e}$')

    plt.xlabel(r'$\kappa_u$')
    plt.ylabel(r'$\Gamma(a\to 3P)$ [GeV]')
    plt.xlim(1e-3, 1)
    plt.legend(ncols=3, fontsize=6, borderpad=0.2, labelspacing=0.1, handlelength=1.4)
    alp_mass_axis(plt.gca(), kappa_u, m_alp)
    plt.tight_layout(pad=0.1)
    plt.savefig('plots/murayama_3P.pdf')

def gammaPP(m_alp):
    kappa_u = np.logspace(-3, 0, 500)
    Gammagammapipi = []
    GammagammaKK = []
    for k in kappa_u:
        inputs = only_up_inputs(k)
        fa = inputs['f_a']
        Gammagammapipi.append(photon_two_pseudoscalar_width(inputs, ['pi+', 'pi-']))
        GammagammaKK.append(photon_two_pseudoscalar_width(inputs, ['K+', 'K-']))

    plt.figure()
    plt.loglog(kappa_u, Gammagammapipi, label=r'$\gamma\pi^+\pi^-$')
    plt.loglog(kappa_u, GammagammaKK, label=r'$\gamma K^+ K^-$')
    plt.title(r'WZW only')
    plt.xlabel(r'$\kappa_u$')
    plt.ylabel(r'$\Gamma(a\to\gamma PP)$ [GeV]')
    plt.xlim(1e-3, 1)
    plt.legend(fontsize=8)
    alp_mass_axis(plt.gca(), kappa_u, m_alp)
    plt.tight_layout(pad=0.1)
    plt.savefig('plots/murayama_gammaPP.pdf')

def sigma_PP(m_sigma, quartic_coupling=1e-3):
    kappa_u = np.logspace(-3, 0, 500)
    dw_channels = {}
    for i, meson1 in enumerate(mesons):
        for j, meson2 in enumerate(mesons):
            if j > i:
                continue
            dw = np.array([sigma_two_pseudoscalar_width(k, quartic_coupling, [meson1, meson2]) for k in kappa_u])
            if any(dw): #and np.max(dw) > 1e-12:
                dw_channels[mesons_tex[meson1] + ' ' + mesons_tex[meson2]] = dw

    colors = plt.cm.jet(np.linspace(0, 1, len(dw_channels)))
    plt.figure()
    for (channel, color) in zip(dw_channels.keys(), colors):
        plt.loglog(kappa_u, dw_channels[channel], label=f'${channel}$', color=color)
    plt.title(rf'LO-$\chi$PT only, $\lambda={quartic_coupling:.2e}$')
    plt.xlabel(r'$\kappa_u$')
    plt.ylabel(r'$\Gamma(\sigma\to 2P)$ [GeV]')
    plt.xlim(1e-3, 1)
    plt.legend(ncols=3, fontsize=6, borderpad=0.2, labelspacing=0.1, handlelength=1.4)
    # Find each crossing of round mass values without inverting the full curve.
    mass_axis = plt.gca().twiny()
    mass_axis.set_xscale('log')
    mass_axis.set_xlim(1e-3, 1)
    masses = np.asarray(m_sigma)
    exponents = range(int(np.floor(np.log10(np.min(masses)))),
                      int(np.ceil(np.log10(np.max(masses)))) + 1)
    round_masses = [coefficient * 10.**exponent for exponent in exponents
                    for coefficient in (1, 2, 5)]
    crossings = []
    for mass in round_masses:
        for index in range(len(masses) - 1):
            if (masses[index] - mass) * (masses[index + 1] - mass) <= 0 and masses[index] != masses[index + 1]:
                fraction = (mass - masses[index]) / (masses[index + 1] - masses[index])
                log_position = np.log10(kappa_u[index]) + fraction * (np.log10(kappa_u[index + 1]) - np.log10(kappa_u[index]))
                crossings.append((log_position, mass))
    crossings.sort()
    selected = []
    for position, mass in crossings:
        if not selected or position - selected[-1][0] >= 0.16:
            selected.append((position, mass))
    mass_axis.set_xticks([10.**position for position, mass in selected])
    mass_axis.set_xticklabels([f'{mass:g}' for position, mass in selected])
    mass_axis.set_xlabel(r'$m_\sigma$ [GeV]')
    plt.tight_layout(pad=0.1)
    plt.savefig('plots/murayama_sigma2P.pdf')

if __name__ == '__main__':
    m_alp, m_sigma = masses()
    gammagamma(m_alp)
    three_mesons(m_alp)
    gammaPP(m_alp)
    sigma_PP(m_sigma)
