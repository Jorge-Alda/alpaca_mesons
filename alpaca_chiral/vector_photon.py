"""ALP decays to omega gamma and phi gamma from the eleven odd operators.

The central provisional single-trace fit supplies the vector masses, mixing,
and photon residues. The calculation is at tree level and leading order in
the electromagnetic coupling, with narrow external vectors. It includes
direct VJP and photon-converted VVP terms. Masses, fa and widths are in GeV.
"""
import numpy as np

from .pseudoscalar_diag import ALPcouplings, _no_mixing_couplings
from .vector_inputs import _read_vector_inputs

masses_PDG = {
    'omega': 782.66e-3,
    'phi': 1019.461e-3,
}

def vector_photon_amplitude(ma: float, fa: float, couplings: ALPcouplings,
                            vector: str) -> complex:
    """Return G in GeV^-1 for a -> V gamma, with V = omega or phi.

    M = G epsilon_{mu nu rho sigma} eps_V^mu eps_gamma^nu q_V^rho
    q_gamma^sigma, with epsilon^0123 = +1 and metric (+---). The input
    ALPcouplings is transformed to Alda:2026nrd's no-mixing basis using the
    existing pseudoscalar machinery. No additional ALP--meson pole or local
    anomalous a V gamma term is included. The analytic amplitude is also
    returned below threshold; no fit uncertainty is propagated.
    """
    fitted = _read_vector_inputs()
    axial, pseudoscalar, _ = _no_mixing_couplings(couplings, ma)
    quark_masses = fitted['quark_mass_matrix']
    mass_anticommutator = (pseudoscalar@quark_masses
                          + quark_masses@pseudoscalar)
    charges = np.diag([2/3, -1/3, -1/3])
    photon_response = sum(
        fitted['photon_vector_residues'][name]
        * fitted['vector_generators'][name]/mass**2
        for name, mass in fitted['vector_masses'].items())
    vector_mass = fitted['vector_masses'][vector]
    generator = fitted['vector_generators'][vector]
    odd = fitted['odd_coefficients']
    direct_axial_trace = np.trace(generator@charges@axial)
    direct_mass_trace = np.trace(generator@charges@mass_anticommutator)
    converted_axial_trace = np.trace(generator@photon_response@axial)
    converted_mass_trace = np.trace(
        generator@photon_response@mass_anticommutator)
    # c4, c7 and d4 vanish for the neutral real-photon contractions.
    direct = -8/(fitted['operator_mass']*fa)*(
        direct_axial_trace*(
            odd['c1']*(ma**2-vector_mass**2)
            + odd['c2']*(ma**2+vector_mass**2)
            - odd['c5']*(ma**2-vector_mass**2)
            - 2*odd['c6']*vector_mass**2)
        + 4*odd['c3']*direct_mass_trace)
    # Both VVP vector assignments are included in this expression.
    converted = 8*np.sqrt(2)/fa*(
        converted_axial_trace*((odd['d1']-odd['d3'])*ma**2
                               + odd['d3']*vector_mass**2)
        + 4*odd['d2']*converted_mass_trace)
    electric_charge = np.sqrt(4*np.pi*fitted['alpha_em'])
    return complex(electric_charge/vector_mass*(direct+converted))


def vector_photon_width(ma: float, fa: float, couplings: ALPcouplings,
                        vector: str) -> float:
    """Return Gamma(a -> V gamma) in GeV, with V = omega or phi.

    The width is |G|^2 (ma^2-M_V^2)^3/(32 pi ma^3) above the fitted vector
    mass and zero at and below threshold. There is no identical-particle
    factor or initial spin average. Daughter widths are not integrated.
    """
    vector_mass = masses_PDG[vector]
    if ma <= vector_mass:
        return 0.0
    amplitude = vector_photon_amplitude(ma, fa, couplings, vector)
    return float(abs(amplitude)**2*(ma**2-vector_mass**2)**3
                 / (32*np.pi*ma**3))
