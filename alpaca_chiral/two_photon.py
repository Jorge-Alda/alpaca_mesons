"""ALP two-real-photon coefficient and width in the Bauer:2020jbp convention.

L=alpha_em c_gamma a F Ftilde/(4 pi fa). The no-mixing anomaly term
and original-eleven-operator VJP/VVP vector conversions interfere coherently.
The resonance inputs are the central provisional single-trace fit. Real
photons use stable propagators at q^2=0 and leading electromagnetic order.
Loops and resummed photon-vector mixing are omitted.
"""
import numpy as np

from .pseudoscalar_diag import ALPcouplings, _no_mixing_couplings
from .vector_inputs import _read_vector_inputs


def photon_coupling(ma: float, fa: float, couplings: ALPcouplings) -> complex:
    """Return dimensionless c_gamma_eff including anomaly and vector exchange.

    couplings.cgamma is the independent coefficient in the input quark basis;
    its default is zero. Masses and fa are in GeV. The WZW/meson mixing term
    is carried by the transformed photon coefficient and appears only once.
    """
    fitted = _read_vector_inputs()
    axial, pseudoscalar, _ = _no_mixing_couplings(couplings, ma)
    charges = np.diag([2/3, -1/3, -1/3])
    photon_local = couplings['cgamma']-6*np.trace(charges@charges@pseudoscalar)
    derivative_source = 2*axial/fa
    quark_masses = fitted['quark_mass_matrix']
    mass_source = -2*(pseudoscalar@quark_masses+quark_masses@pseudoscalar)/fa
    odd = fitted['odd_coefficients']
    response = {name: residue/fitted['vector_masses'][name]**2
                for name, residue in fitted['photon_vector_residues'].items()}
    one_vector, two_vector = 0j, 0j
    for first, first_response in response.items():
        generator = fitted['vector_generators'][first]
        ordered = 4/fitted['operator_mass']*(
            ma**2/2*(odd['c1']+odd['c2']-odd['c5'])*np.trace(generator@charges@derivative_source)
            -2*odd['c3']*np.trace(generator@charges@mass_source))
        one_vector -= 4*np.sqrt(2)*first_response*ordered
        for second, second_response in response.items():
            product = generator@fitted['vector_generators'][second]
            ordered_pair = 4*ma**2*(odd['d1']-odd['d3'])*np.trace(product@derivative_source)-16*odd['d2']*np.trace(product@mass_source)
            two_vector += 2*first_response*second_response*ordered_pair
    return complex(photon_local-4*np.pi**2*fa*(one_vector+two_vector))


def two_photon_width(ma: float, fa: float, couplings: ALPcouplings) -> float:
    """Return Gamma(a -> gamma gamma) in GeV, alpha_em^2 |c|^2 ma^3/(64 pi^3 fa^2)."""
    alpha_em = _read_vector_inputs()['alpha_em']
    return float(alpha_em**2*abs(photon_coupling(ma, fa, couplings))**2*ma**3/(64*np.pi**3*fa**2))
