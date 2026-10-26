"""Leading WZW scalar for a photon and two physical pseudoscalars.

All masses and decay constants use GeV. Couplings must already be in the
no-mixing basis of Alda:2026nrd; this module performs no basis conversion.
The calculation and sign are recorded in
memory/alp/photon-two-pseudoscalar/wzw-alp-photon-pp-amplitude.md.

Run the example from the repository root with::

    .venv/bin/python -m scripts.implementation.alp_photon_pp_wzw
"""
import math

import numpy as np

from alpaca_chiral.photon_two_pseudoscalar import _canonical_meson_name
from alpaca_chiral.pseudoscalar_diag import ALPcouplings
from alpaca_chiral.three_pseudoscalar import F_PI, _pseudoscalars
from alpaca_chiral.vector_inputs import _read_vector_inputs


def photon_two_pseudoscalar_wzw_amplitude(
        ma: float, fa: float, couplings: ALPcouplings,
        meson_1: str, meson_2: str, *, orientation: int = 1) -> complex:
    """Return the no-mixing WZW scalar coefficient in GeV^-3.

    The invariant amplitude is epsilon(eps_gamma,p1,p2,k) times the
    returned coefficient, using metric (+---) and epsilon^0123=+1.
    This local leading-order coefficient has no Dalitz dependence. ``ma``
    only imposes the physical threshold; the coefficient is zero at and
    below it, as in the public resonance amplitude.

    ``couplings.cuA`` and ``couplings.cdA`` contain the already transformed
    no-mixing axial coefficients in the existing ALPcouplings convention.
    The Alda u,d,s axial matrix is one half of these entries, including the
    ds and sd coefficients. Neither translate nor a mixing solver is called.
    Gluonic, photon and vector couplings do not enter this direct vertex.

    Meson names are pi0, pi+, pi-, eta, eta_prime, K+, K-, K0, K0bar,
    KS and KL (plus the aliases of the public photon PP amplitude). Neutral
    mesons use the fitted physical generators, and KS/KL use the same
    complex kaon projections as the library, including epsilon_K.

    ``orientation=+1`` selects the explicit common WZW orientation of the
    derivation record. It is a sign convention, independent of the epsilon
    tensor convention, and must be matched when combining odd amplitudes.
    The retained approximation is tree level, first order in 1/fa, and
    leading anomalous chiral order. Loops and rescattering are omitted.
    """
    first_name, second_name = map(_canonical_meson_name, (meson_1, meson_2))
    particles = _pseudoscalars()
    first_mass, first = particles[first_name]
    second_mass, second = particles[second_name]
    if ma <= first_mass + second_mass:
        return 0j
    up_axial, down_axial = couplings['cuA'], couplings['cdA']
    axial = .5*np.array([
        [up_axial[0, 0], 0, 0],
        [0, down_axial[0, 0], down_axial[0, 1]],
        [0, down_axial[1, 0], down_axial[1, 1]],
    ], dtype=complex)
    charge = np.diag([2/3, -1/3, -1/3])
    flavour = np.trace(charge @ (
        axial @ first @ second + first @ second @ axial
        + second @ axial @ first - axial @ second @ first
        - second @ first @ axial - first @ axial @ second))
    electric_charge = math.sqrt(4*math.pi*_read_vector_inputs()['alpha_em'])
    return complex(orientation*3*electric_charge*flavour
                   / (12*math.pi**2*fa*F_PI**2))


if __name__ == '__main__':
    example_couplings = ALPcouplings(
        cG=0., cuA=np.diag([2., 0., 0.]), cuV=np.zeros((3, 3)),
        cdA=np.zeros((3, 3)), cdV=np.zeros((3, 3)))
    print(photon_two_pseudoscalar_wzw_amplitude(
        1., 1000., example_couplings, 'pi+', 'pi-'))
