"""Original-eleven-operator tree ALP amplitudes to neutral vector pairs.

Masses and fa are in GeV. The central provisional single-trace fit fixes
odd coefficients and neutral mixing. The vectors are stable external states;
no spectral integration, loops or uncertainty propagation is included.
"""
import numpy as np

from .pseudoscalar_diag import ALPcouplings, _no_mixing_couplings
from .vector_inputs import _read_vector_inputs

masses_PDG = {
    'omega': 782.66e-3,
    'phi': 1019.461e-3,
}

def neutral_vector_pair_amplitude(ma: float, fa: float, couplings: ALPcouplings,
                                  vector_1: str, vector_2: str) -> complex:
    """Return G in GeV^-1 for the selected rho0, omega or phi pair.

    M=G epsilon(eps1,eps2,q1,q2), epsilon^0123=+1, metric (+---).
    The commuting neutral-vector contractions include VVP and both VJP
    photon-conversion diagrams. Charged vector channels are outside this
    implementation. No second meson-pole contribution is added.
    """
    fitted = _read_vector_inputs()
    axial, pseudoscalar, _ = _no_mixing_couplings(couplings, ma)
    derivative_source = 2*axial/fa
    quark_masses = fitted['quark_mass_matrix']
    mass_source = -2*(pseudoscalar@quark_masses+quark_masses@pseudoscalar)/fa
    odd = fitted['odd_coefficients']
    first, second = vector_1, vector_2
    first_mass, second_mass = (fitted['vector_masses'][name] for name in (first, second))
    generators = fitted['vector_generators'][first]@fitted['vector_generators'][second]
    derivative_trace = np.trace(generators@derivative_source)
    mass_trace = np.trace(generators@mass_source)
    momentum_product = (ma**2-first_mass**2-second_mass**2)/2
    two_vector = 4*(odd['d1']*ma**2-2*odd['d3']*momentum_product)*derivative_trace-16*odd['d2']*mass_trace

    def vector_photon(vector_squared, photon_squared):
        """Return the ordered VJP contraction before photon conversion."""
        return 4/fitted['operator_mass']*(derivative_trace*(
            odd['c1']*(photon_squared+momentum_product)
            +odd['c2']*(vector_squared+momentum_product)
            -odd['c5']*momentum_product-odd['c6']*vector_squared)-2*odd['c3']*mass_trace)

    conversion = -2*np.sqrt(2)*4*np.pi*fitted['alpha_em']*(
        fitted['photon_vector_residues'][second]*vector_photon(first_mass**2, second_mass**2)
        +fitted['photon_vector_residues'][first]*vector_photon(second_mass**2, first_mass**2))
    return complex((two_vector+conversion)/(first_mass*second_mass))


def neutral_vector_pair_width(ma: float, fa: float, couplings: ALPcouplings,
                              vector_1: str, vector_2: str) -> float:
    """Return the selected neutral-vector partial width in GeV.

    vector_1 and vector_2 are rho0, omega or phi. Includes the identical-
    particle factor 1/2 and stable-vector two-body thresholds.
    """
    first_mass, second_mass = (masses_PDG[name]
                               for name in (vector_1, vector_2))
    if ma <= first_mass+second_mass:
        return 0.0
    amplitude = neutral_vector_pair_amplitude(ma, fa, couplings, vector_1, vector_2)
    phase_space = (ma**2-(first_mass+second_mass)**2)*(ma**2-(first_mass-second_mass)**2)
    symmetry = 2 if vector_1 == vector_2 else 1
    return float(abs(amplitude)**2*phase_space**1.5/(32*np.pi*ma**3*symmetry))
