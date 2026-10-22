"""Selected strong a -> (omega, phi) PP widths, in GeV.

No-mixing c_q^A, c_q^P and c_q^{VS} are supplied directly. The original
seven VJP/four VVP provisional fit supplies d1,d2,d3, vector generators,
internal vector masses and total widths. G_V is a required independent input
in GeV: it was not determined by that fit. Electromagnetic conversion,
rescattering, loops and the additional vector-mass-spurion interaction are
omitted, as in memory/alp/three-meson/alp-vpp-amplitude.md.

Internal vectors use constant-width denominators r^2-M_R^2+i M_R Gamma_R,
with shared fitted total widths; charged real components share their pole.
Internal pseudoscalars use the real LO masses of Alda:2026nrd. Physical masses
are used for external kinematics. K_S/K_L use the repository epsilon_K and
common neutral-kaon mass. Every identical denominator is combined before
quadrature, with relative roundoff cancellation made exact.

Example::

    import numpy as np
    from scripts.implementation.alp_vpp_width import alp_vpp_width
    width = alp_vpp_width(1.5, np.diag([1., 0., 0.]), np.zeros((3, 3)),
                          ('omega', 'pi+', 'pi-'), 1000., g_v=0.065)

The example value of G_V is a caller choice, not a fit determination.
"""
from dataclasses import dataclass
import json
import math

import numpy as np
from scipy.integrate import quad

from .alda_neutral_inputs import SNAPSHOT, alda_2026nrd_neutral_inputs
from .alp_vv_width import read_preliminary_fit_parameters
from .three_meson_width import (
    EPSILON_K, _MESONS, _canonical_name, _physical_neutral_kaon,
    _kaon_mixing_coefficients,
)
from .vector_meson_exchange import _aggregate_flavor_coefficients


def _anticommutator(first, second):
    return first @ second + second @ first


def _commutator(first, second):
    return first @ second - second @ first


def _couplings(values):
    """Accept a diagonal three-vector or full dimensionless source matrix."""
    matrix = np.asarray(values, dtype=complex)
    if matrix.shape == (3,):
        matrix = np.diag(matrix)
    if matrix.shape != (3, 3):
        raise ValueError('quark couplings must have shape (3,) or (3, 3)')
    if not np.all(np.isfinite(matrix)):
        raise ValueError('quark couplings must contain finite entries')
    if not np.allclose(matrix, matrix.conj().T, rtol=1e-12, atol=1e-12):
        raise ValueError('quark couplings must be Hermitian')
    # A neutral ALP cannot have charged u-d or u-s source entries.
    if np.any(matrix[0, 1:]) or np.any(matrix[1:, 0]):
        raise ValueError('neutral ALP sources may contain diagonal and ds/sd entries only')
    return matrix


def _positive_number(value, name):
    """Convert a positive finite real input in GeV."""
    if isinstance(value, (bool, complex)) or not np.isscalar(value):
        raise ValueError(name + ' must be a positive finite real number')
    value = float(value)
    if not math.isfinite(value) or value <= 0:
        raise ValueError(name + ' must be a positive finite real number')
    return value


def _real_components(row, column):
    """Hermitian generators with Tr(T_i T_j)=delta_ij."""
    transition = np.zeros((3, 3), dtype=complex)
    transition[row, column] = 1.
    return ((transition + transition.T) / math.sqrt(2),
            -1j * (transition - transition.T) / math.sqrt(2))


def _meson_inputs(epsilon_k):
    """Return Alda generators (trace norm two) and LO masses in GeV."""
    neutral = alda_2026nrd_neutral_inputs(mass_scheme='tree')
    generators = {name: meson.generator for name, meson in _MESONS.items()}
    tree_masses = {}
    for name, key in (('pi0', 'pi0'), ('eta', 'eta'), ('eta_prime', 'etap')):
        generators[name] = math.sqrt(2)*neutral[key]['TP']
        tree_masses[name] = neutral[key]['tree_mass']
    quark_masses = json.loads(SNAPSHOT.read_text())['B0mq_gev2']
    for names, row, column in ((('pi+', 'pi-'), 0, 1),
                               (('K+', 'K-'), 0, 2),
                               (('K0', 'K0bar', 'KS', 'KL'), 1, 2)):
        for name in names:
            tree_masses[name] = math.sqrt(quark_masses[row] + quark_masses[column])
    for name in ('KS', 'KL'):
        generators[name] = _physical_neutral_kaon(name, epsilon_k).generator
    return generators, tree_masses


def _flavor_projection(name, epsilon_k):
    """Return outgoing flavor-bra coefficients, including complex kaon mixing."""
    if name not in ('KS', 'KL'):
        return ((name, 1.+0j),)
    p_mixing, q_mixing = _kaon_mixing_coefficients(epsilon_k)
    return (('K0', p_mixing.conjugate()),
            ('K0bar', (1 if name == 'KS' else -1)*q_mixing.conjugate()))


@dataclass(frozen=True)
class VPPAmplitude:
    """Prepared invariant coefficients for one ordered V,P1,P2 state.

    ``terms`` stores (coefficient index, invariant index, pole mass, width,
    constant numerator, linear numerator). Invariant indices 0,1,2 denote
    s,t,u. Coefficient indices denote F_odd,F_1,F_2. All units are GeV.
    """
    alp_mass: float
    masses: tuple
    final_state: tuple
    terms: tuple

    def __call__(self, s, t):
        """Return the coherent F_odd, F_1, F_2 coefficients at (s,t)."""
        invariants = (s, t, self.alp_mass**2 + sum(m**2 for m in self.masses) - s - t)
        result = np.zeros(3, dtype=complex)
        for coefficient, invariant, mass, width, constant, linear in self.terms:
            virtuality = invariants[invariant]
            result[coefficient] += (constant + linear * virtuality) / (
                virtuality - mass**2 + 1j * mass * width)
        return tuple(result)

    def squared(self, s, t):
        """Return the polarization-summed squared amplitude (no spin average)."""
        vector_mass, first_mass, second_mass = self.masses
        u = self.alp_mass**2 + sum(m**2 for m in self.masses) - s - t
        vector_first = (u - vector_mass**2 - first_mass**2) / 2
        vector_second = (t - vector_mass**2 - second_mass**2) / 2
        first_second = (s - first_mass**2 - second_mass**2) / 2
        gram = np.array([[vector_mass**2, vector_first, vector_second],
                         [vector_first, first_mass**2, first_second],
                         [vector_second, first_second, second_mass**2]])
        odd, first, second = self(s, t)
        transverse_first = vector_first**2 / vector_mass**2 - first_mass**2
        transverse_second = vector_second**2 / vector_mass**2 - second_mass**2
        transverse_cross = vector_first * vector_second / vector_mass**2 - first_second
        return max(0., float(max(0., np.linalg.det(gram)) * abs(odd)**2
                             + transverse_first * abs(first)**2
                             + transverse_second * abs(second)**2
                             + 2 * transverse_cross * (first * second.conjugate()).real))


def build_vpp_amplitude(alp_mass, c_q_A, c_q_P, final_state, f_a, *, g_v,
                        c_q_VS=None, fit_parameters=None, epsilon_k=EPSILON_K):
    """Prepare only the requested ('omega' or 'phi', P1, P2) amplitude.

    G_V and f_a are in GeV. c_q_VS is restricted to ds/sd by the retained
    |Delta S|=1 even-source calculation; diagonal entries commute with chi.
    K_L/K_S aliases follow three_meson_width. A zero amplitude has no terms.
    The amplitude uses f_pi=0.093 GeV from the stored Alda:2026nrd input.
    """
    alp_mass = _positive_number(alp_mass, 'alp_mass')
    f_a = _positive_number(f_a, 'f_a')
    g_v = _positive_number(g_v, 'g_v')
    axial, pseudoscalar = _couplings(c_q_A), _couplings(c_q_P)
    supplied = None if c_q_VS is None else _couplings(c_q_VS)
    vector, first_name, second_name = final_state
    if vector not in ('omega', 'phi'):
        raise ValueError('the derived external vector is omega or phi')
    first_name, second_name = map(_canonical_name, (first_name, second_name))
    fitted = fit_parameters if fit_parameters is not None else read_preliminary_fit_parameters()
    physical_vector_mass = {'omega': 0.78266, 'phi': 1.019461}[vector]
    masses = (physical_vector_mass, _MESONS[first_name].mass, _MESONS[second_name].mass)
    state = (vector, first_name, second_name)
    if alp_mass <= sum(masses) or _MESONS[first_name].charge + _MESONS[second_name].charge:
        return VPPAmplitude(alp_mass, masses, state, ())
    vector_source = np.zeros((3, 3), complex)
    if supplied is not None:
        vector_source[1, 2], vector_source[2, 1] = supplied[1, 2], supplied[2, 1]
    generators, tree_masses = _meson_inputs(epsilon_k)
    external = fitted.vector_generators[vector]
    chi = fitted.quark_mass_matrix
    pion_constant = json.loads(SNAPSHOT.read_text())['fpi_gev']
    odd = fitted.odd_coefficients
    vector_groups = [(fitted.vector_masses[name], fitted.shared_total_widths[name],
                      (fitted.vector_generators[name],)) for name in ('rho0', 'omega', 'phi')]
    for row, column, width_name in ((0, 1, 'rhoplus'), (0, 2, 'Kstarplus'), (1, 2, 'Kstar0')):
        mass = fitted.even_parameters['m_rho_GeV' if row == 0 and column == 1 else 'm_Kstar_GeV']
        vector_groups.append((mass, fitted.shared_total_widths[width_name], _real_components(row, column)))
    contributions = {}

    def add(coefficient, invariant, mass, width, constant, linear=0j):
        values = contributions.setdefault((coefficient, invariant, mass, width), ([], []))
        values[0].append(projection*constant)
        values[1].append(projection*linear)

    # Project complete flavor amplitudes; physical generators already contain
    # conjugated mixing weights and must not be conjugated a second time.
    for first_flavor, first_weight in _flavor_projection(first_name, epsilon_k):
        for second_flavor, second_weight in _flavor_projection(second_name, epsilon_k):
            projection = first_weight*second_weight
            first = generators[first_flavor].conj().T
            second = generators[second_flavor].conj().T
            for mass, width, components in vector_groups:
                for internal in components:
                    anticommutator = _anticommutator(external, internal)
                    for invariant, source, mass_source, source_mass_squared, trace, factor in (
                        (0, 2*axial/f_a, -2*_anticommutator(pseudoscalar, chi)/f_a,
                         alp_mass**2, np.trace(internal @ _commutator(first, second)),
                         -g_v/(math.sqrt(2)*pion_constant**2*physical_vector_mass)),
                        (1, -second/pion_constant,
                         _anticommutator(second, chi)/pion_constant,
                         tree_masses[second_flavor]**2, np.trace(internal @ _commutator(first, axial)),
                         math.sqrt(2)*g_v/(pion_constant*f_a*physical_vector_mass)),
                        (2, -first/pion_constant,
                         _anticommutator(first, chi)/pion_constant,
                         tree_masses[first_flavor]**2, np.trace(internal @ _commutator(second, axial)),
                         -math.sqrt(2)*g_v/(pion_constant*f_a*physical_vector_mass))):
                        derivative_trace = np.trace(anticommutator @ source) / 2
                        mass_trace = np.trace(anticommutator @ mass_source) / 2
                        constant = 4*((odd['d1']*source_mass_squared
                                       - odd['d3']*(source_mass_squared-physical_vector_mass**2))*derivative_trace
                                      - 4*odd['d2']*mass_trace)
                        add(0, invariant, mass, width, factor*trace*constant,
                            factor*trace*4*odd['d3']*derivative_trace)
            pseudoscalar_groups = [(tree_masses[name], (generators[name],))
                                   for name in ('pi0', 'eta', 'eta_prime')]
            for row, column, name in ((0, 1, 'pi+'), (0, 2, 'K+'), (1, 2, 'K0')):
                pseudoscalar_groups.append((tree_masses[name], tuple(
                    math.sqrt(2)*generator for generator in _real_components(row, column))))
            for mass, components in pseudoscalar_groups:
                for internal in components:
                    for coefficient, invariant, spectator, emitted in ((1, 2, second, first), (2, 1, first, second)):
                        flavor = (np.trace(_anticommutator(spectator, internal) @ _commutator(vector_source, chi))
                                  * np.trace(external @ _commutator(internal, emitted)))
                        add(coefficient, invariant, mass, 0.,
                            -1j*math.sqrt(2)*g_v*physical_vector_mass/(8*pion_constant**2*f_a)*flavor)
    terms = []
    for key, (constants, linears) in contributions.items():
        constant, linear = map(_aggregate_flavor_coefficients, (constants, linears))
        if constant != 0j or linear != 0j:
            terms.append((*key, constant, linear))
    return VPPAmplitude(alp_mass, masses, state, tuple(terms))


def alp_vpp_width(alp_mass, c_q_A, c_q_P, final_state, f_a, *, g_v,
                   c_q_VS=None, fit_parameters=None, epsilon_k=EPSILON_K,
                   epsabs=0., epsrel=1e-5, limit=150):
    """Return one partial width in GeV; forbidden states return 0 before quad.

    Uses nested adaptive quadrature in s and t. Internal vector pole masses
    inside the integration interval are supplied as subdivision points.
    IntegrationWarning is propagated by SciPy if convergence is not achieved.
    """
    amplitude = build_vpp_amplitude(alp_mass, c_q_A, c_q_P, final_state, f_a,
                                    g_v=g_v, c_q_VS=c_q_VS, fit_parameters=fit_parameters,
                                    epsilon_k=epsilon_k)
    if not amplitude.terms:
        return 0.
    alp_mass = amplitude.alp_mass
    vector_mass, first_mass, second_mass = amplitude.masses
    total_squared = alp_mass**2 + sum(mass**2 for mass in amplitude.masses)

    def kallen(first, second, third):
        return first**2 + second**2 + third**2 - 2*(first*second + first*third + second*third)

    def integrate_t(s):
        center = vector_mass**2 + second_mass**2 + (
            (alp_mass**2-vector_mass**2-s)*(s+second_mass**2-first_mass**2)/(2*s))
        half_range = math.sqrt(max(0., kallen(alp_mass**2, vector_mass**2, s))
                               * max(0., kallen(s, first_mass**2, second_mass**2))) / (2*s)
        lower, upper = center-half_range, center+half_range
        points = set()
        for _, invariant, mass, _, _, _ in amplitude.terms:
            pole = mass**2 if invariant == 1 else total_squared-s-mass**2
            if invariant in (1, 2) and lower < pole < upper:
                points.add(pole)
        return quad(lambda t: amplitude.squared(s, t), lower, upper,
                    points=sorted(points), epsabs=epsabs, epsrel=epsrel, limit=limit)[0]

    lower, upper = (first_mass+second_mass)**2, (alp_mass-vector_mass)**2
    points = sorted({mass**2 for _, invariant, mass, _, _, _ in amplitude.terms
                     if invariant == 0 and lower < mass**2 < upper})
    integral = quad(integrate_t, lower, upper, points=points,
                    epsabs=epsabs, epsrel=epsrel, limit=limit)[0]
    identical_factor = 2 if amplitude.final_state[1] == amplitude.final_state[2] else 1
    return integral / (256*math.pi**3*alp_mass**3*identical_factor)
