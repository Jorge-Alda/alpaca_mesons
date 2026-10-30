"""KKpi odd-coefficient form factors and normalized electromagnetic integrals.

GomezDumm:2010ra Eqs. (30),(35), Appendix A; canonical epsilon convention.
All masses/decay constants are in GeV, invariants in GeV^2. The current is
H^mu=i epsilon^{mu nu rho sigma} p1_nu p2_rho p3_sigma F, with F in GeV^-3.
The resonance part of the charged-current F3 is implemented without OPE
relations. To obtain an electromagnetic observable, supply its explicit
isospin map and anomaly term. No such map or experimental Dalitz projection
is silently inferred from the BaBar isovector K*(892)K table.
"""

from dataclasses import dataclass
from typing import Callable

import numpy as np
from numpy.polynomial.legendre import leggauss


ODD_COEFFICIENTS = tuple([f"c{i}" for i in range(1, 8)]
                         + [f"d{i}" for i in range(1, 5)]
                         + [f"g{i}" for i in range(1, 6)])
GEV_MINUS_TWO_TO_NB = 389379.3656


@dataclass(frozen=True)
class KKPiInputs:
    """Exact-isospin inputs; every inverse propagator is supplied by caller.

    Propagators return 1/(M^2-x-Sigma(x)) in GeV^-2 on the upper rim.
    Their self energies, masses and subtraction prescriptions are external.
    """

    kaon_mass: float
    pion_mass: float
    pion_decay_constant: float
    vector_source_coupling: float
    vector_pion_coupling: float
    vector_mass_normalization: float
    rho_propagator: Callable
    omega_propagator: Callable
    kstar_propagator: Callable


def resonance_scalar_basis(q_squared, pair_s, pair_t, inputs, mode,
                           component="full"):
    """Return the 16 canonical odd-coefficient gradients of weak-current F3.

    Modes are ``charged`` (K+ K- pi-) and ``neutral`` (K- K0 pi0), in
    GomezDumm:2010ra's invariant ordering. ``kstar`` keeps precisely the
    K*(892) denominators of the source formulas. This algebraic selection
    is not the experimental Dalitz-fit projection. The WZW constant and
    excited-vector amplitudes must be supplied separately. Couplings use
    the canonical convention, c=-c_reference etc.; linear scalar gradients
    retain the reference formula's signs under the tensor conversion.
    """
    if mode not in ("charged", "neutral") or component not in ("full", "kstar"):
        raise ValueError("Choose charged/neutral mode and full/kstar component")
    q_squared, pair_s, pair_t = np.broadcast_arrays(q_squared, pair_s, pair_t)
    kaon_squared = inputs.kaon_mass**2
    pion_squared = inputs.pion_mass**2
    pair_u = q_squared-pair_s-pair_t+2*kaon_squared+pion_squared
    shape = q_squared.shape+(len(ODD_COEFFICIENTS),)

    def direct(pair, first_mass, second_mass, third_mass):
        result = np.zeros(shape, complex)
        result[..., 0] = q_squared-pair+third_mass
        result[..., 1] = -q_squared+pair+third_mass
        result[..., 2] = 8*third_mass
        result[..., 3] = 8*(first_mass-second_mass)
        result[..., 4] = q_squared+pair-third_mass
        result[..., 5] = -2*pair
        return result

    def double(pair, third_mass):
        result = np.zeros(shape, complex)
        result[..., 7] = third_mass
        result[..., 8] = 8*third_mass
        result[..., 9] = q_squared+pair-third_mass
        return result

    kstar_t = np.asarray(inputs.kstar_propagator(pair_t))[..., None]
    direct_terms = direct(pair_t, kaon_squared, pion_squared, kaon_squared)*kstar_t
    double_terms = double(pair_t, kaon_squared)*kstar_t
    local_terms = np.zeros(shape, complex)
    if mode == "neutral":
        kstar_u = np.asarray(inputs.kstar_propagator(pair_u))[..., None]
        direct_terms -= direct(pair_u, kaon_squared, pion_squared, kaon_squared)*kstar_u
        double_terms -= double(pair_u, kaon_squared)*kstar_u
        local_terms[..., 11] = pair_t-pair_u
        local_terms[..., 12] = 2*(pair_t-pair_u)
        local_terms[..., 13] = -(pair_t-pair_u)
        prefactor = 2*np.sqrt(2)*inputs.vector_pion_coupling
    else:
        if component == "full":
            omega_s = np.asarray(inputs.omega_propagator(pair_s))[..., None]
            direct_terms += direct(pair_s, kaon_squared, kaon_squared, pion_squared)*omega_s
            double_terms += double(pair_s, pion_squared)*omega_s
        local_terms[..., 11] = pair_s+pair_t-3*kaon_squared-pion_squared
        local_terms[..., 12] = 2*(pair_s+pair_t-q_squared-kaon_squared)
        local_terms[..., 13] = -local_terms[..., 11]
        local_terms[..., 14] = 2*(kaon_squared+pion_squared)
        local_terms[..., 15] = 2*kaon_squared
        prefactor = -4*inputs.vector_pion_coupling
    if component == "kstar":
        local_terms[...] = 0
    rho = np.asarray(inputs.rho_propagator(q_squared))[..., None]
    scalar = (direct_terms-2*inputs.vector_source_coupling/inputs.vector_pion_coupling
              *rho*local_terms-np.sqrt(2)*inputs.vector_source_coupling
              *inputs.vector_mass_normalization*rho*double_terms)
    return prefactor/(inputs.vector_mass_normalization*inputs.pion_decay_constant**3)*scalar


def triangle(first, second, third):
    """Kallen function of squared masses/invariants, in GeV^4."""
    return (first-second-third)**2-4*second*third


def electromagnetic_cross_section_matrix(energy, masses, scalar_basis,
                                         alpha, pair_order=32, angle_order=32):
    """Integrate EM scalar basis into a real symmetric cross-section matrix.

    ``scalar_basis(Q2,s12,s13)`` returns (...,n) complex scalar coefficients
    in GeV^-3 including any fixed amplitude as a separate column. For real
    parameters v the Born cross section is v.T @ result @ v, in nb. The
    caller must supply an electromagnetic current, including its physical
    charge/isospin factors, and every coherent excited-resonance amplitude.
    There is one distinct charge final state; charge conjugates and K_S
    projections require explicit factors in the supplied scalar.
    """
    first_mass, second_mass, third_mass = masses
    if energy <= sum(masses):
        # The basis dimension is obtained without probing an unphysical pole.
        return None
    q_squared = energy**2
    pair_lower = (first_mass+second_mass)**2
    pair_upper = (energy-third_mass)**2
    nodes, weights = leggauss(pair_order)
    angles, angular_weights = leggauss(angle_order)
    pair_s = (pair_upper+pair_lower)/2+(pair_upper-pair_lower)/2*nodes
    pair_weights = (pair_upper-pair_lower)/2*weights
    outer_triangle = np.maximum(triangle(q_squared, pair_s, third_mass**2), 0)
    inner_triangle = np.maximum(triangle(pair_s, first_mass**2, second_mass**2), 0)
    pair_span = np.sqrt(outer_triangle*inner_triangle)/pair_s
    pair_t = (first_mass**2+third_mass**2
              +(pair_s+first_mass**2-second_mass**2)
              *(q_squared-pair_s-third_mass**2)/(2*pair_s))[:, None]
    pair_t = pair_t+pair_span[:, None]*angles[None, :]/2
    scalar = np.asarray(scalar_basis(q_squared, pair_s[:, None], pair_t), complex)
    if scalar.shape[:2] != pair_t.shape or scalar.ndim != 3:
        raise ValueError("Electromagnetic scalar basis must have shape (ns,nz,ncoeff)")
    # |epsilon p1 p2 p3|^2 = lambda(Q2,s,m3²) lambda(s,m1²,m2²)
    #                           * (1-z²)/(16s).
    current_squared = (outer_triangle*inner_triangle/(16*pair_s))[:, None]*(1-angles**2)
    dalitz_measure = pair_weights[:, None]*angular_weights[None, :]*pair_span[:, None]/2
    weight = (GEV_MINUS_TWO_TO_NB*alpha**2/(48*np.pi*q_squared**3)
              *dalitz_measure*current_squared)
    return np.einsum("ij,ija,ijb->ab", weight, scalar.conj(), scalar).real


def energy_bin_matrices(edges, masses, scalar_basis, alpha, energy_order=8,
                        pair_order=32, angle_order=32):
    """Uniform dE bin averages; use explicit experimental weights if available.

    Uniform bin averaging is an approximation for published bin estimates.
    No ISR luminosity or acceptance response is available in Table 7.
    """
    nodes, weights = leggauss(energy_order)
    matrices = []
    for lower, upper in np.asarray(edges, float):
        energies = (upper+lower)/2+(upper-lower)/2*nodes
        values = [electromagnetic_cross_section_matrix(energy, masses, scalar_basis,
                  alpha, pair_order, angle_order) for energy in energies]
        template = next((value for value in values if value is not None), None)
        if template is None:
            raise ValueError("A closed bin needs an explicit scalar-basis dimension")
        matrices.append(sum(weight*(np.zeros_like(template) if value is None else value)
                            for weight, value in zip(weights, values))/2)
    return np.asarray(matrices)


def predict_from_matrices(coefficients, matrices):
    """Evaluate the coherent quadratic prediction for real coefficients."""
    coefficients = np.asarray(coefficients, float)
    return np.einsum("a,...ab,b->...", coefficients, matrices, coefficients)


def babar_table7_prediction(*args, **kwargs):
    """Reject a direct inclusive-to-Dalitz-component identification."""
    raise NotImplementedError(
        "BaBar:2007ceh Table 7 requires its K*(892)K Dalitz-component "
        "projection, charge/isospin normalization and excited-vector model. "
        "The inclusive KKpi F3 integral is a different observable; its "
        "VPPP coefficient dependence cannot be assigned to Table 7.")
