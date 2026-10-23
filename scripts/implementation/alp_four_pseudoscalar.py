"""Axial/pseudoscalar a -> 4P vector exchange in the no-mixing basis.

Masses and widths are in GeV. Source matrices are dimensionless. The original
preliminary odd fit and conditional G_V fit supply the mesonic parameters.
The retained tree amplitudes, constant vector widths and phase-space measure
are documented in memory/alp/four-pseudoscalar/. All crossed diagrams interfere.

Example::

    amplitude = build_four_pseudoscalar_amplitude(
        1.5, [1., 0., 0.], [0., 0., 0.], ('pi+', 'pi-', 'pi+', 'pi-'), 1000.)
    width = alp_four_pseudoscalar_width(
        1.5, [1., 0., 0.], [0., 0., 0.], ('pi+', 'pi-', 'pi+', 'pi-'), 1000.)

The width uses deterministic adaptive five-dimensional cubature. No narrow
width approximation or extra ALP--meson rotation is made.
"""
from collections import Counter
from dataclasses import dataclass
from itertools import combinations, permutations, product
import json
import math
from pathlib import Path
import warnings

import numpy as np
from scipy.integrate import cubature

from .alda_neutral_inputs import SNAPSHOT
from .alp_vv_width import read_preliminary_fit_parameters
from .alp_vpp_width import (
    _anticommutator, _commutator, _couplings, _meson_inputs,
    _flavor_projection, _real_components,
)
from .three_meson_width import EPSILON_K, _MESONS, _canonical_name
from .vector_meson_exchange import DEFAULT_VECTOR_POLES, _aggregate_flavor_coefficients

PAIR_LABELS = tuple(combinations(range(4), 2))
PAIR_INDEX = {pair: index for index, pair in enumerate(PAIR_LABELS)}
DEFAULT_GV_FIT_PATH = Path(__file__).resolve().parents[2] / 'data/vector_sector/gv_conditional_fit.json'


def read_four_pseudoscalar_g_v(fit_path=None):
    """Return the positive conditional fitted G_V, in GeV."""
    return float(json.loads(Path(fit_path or DEFAULT_GV_FIT_PATH).read_text())['parameters']['G_V_GeV'])


def _parity(labels):
    """Permutation parity relative to the external momentum ordering."""
    return (-1) ** sum(labels[i] > labels[j] for i in range(4) for j in range(i+1, 4))


def _kallen(first, second, third):
    """Kallen polynomial for squared masses, returning GeV^4."""
    return (first-second-third)**2 - 4*second*third


@dataclass(frozen=True)
class FourPseudoscalarAmplitude:
    """Prepared scalar F in M=epsilon(p1,p2,p3,p4) F.

    Six pair invariants use the order (12,13,14,23,24,34). ``terms`` holds
    two invariant indices, their complex squared poles, and the constant and
    two linear numerator coefficients. Indices 6+i mean (p_a-p_i)^2.
    Flavor sums and equal-pole cancellations are completed during preparation.
    ``scalar`` accepts arrays with final dimension six and returns GeV^-5.
    """
    alp_mass: float
    masses: tuple
    final_state: tuple
    terms: tuple

    @property
    def symmetry_factor(self):
        """Product of identical-species factorials in the labeled width."""
        return math.prod(math.factorial(count) for count in Counter(self.final_state).values())

    def scalar(self, pair_invariants):
        """Return coherent F for one point or a batch of pair invariants."""
        pairs = np.asarray(pair_invariants, dtype=float)
        triples = []
        mass_squared = np.square(self.masses)
        for omitted in range(4):
            retained = [index for index in range(4) if index != omitted]
            triples.append(sum(pairs[..., PAIR_INDEX[pair]] for pair in combinations(retained, 2))
                           - sum(mass_squared[index] for index in retained))
        invariants = np.concatenate((pairs, np.stack(triples, axis=-1)), axis=-1)
        result = np.zeros(pairs.shape[:-1], dtype=complex)
        # Evaluate each distinct propagator once for the whole batch.
        reciprocals = {}
        for first, second, first_pole, second_pole, constant, first_linear, second_linear in self.terms:
            for index, pole in ((first, first_pole), (second, second_pole)):
                if (index, pole) not in reciprocals:
                    reciprocals[index, pole] = 1 / (invariants[..., index] - pole)
            result += ((constant + first_linear*invariants[..., first]
                        + second_linear*invariants[..., second])
                       * reciprocals[first, first_pole] * reciprocals[second, second_pole])
        return result

    __call__ = scalar

    def squared(self, pair_invariants):
        """Return |M|^2 in GeV^-2 from the four-momentum Gram determinant."""
        pairs = np.asarray(pair_invariants, dtype=float)
        gram = np.zeros(pairs.shape[:-1] + (4, 4))
        for index, mass in enumerate(self.masses):
            gram[..., index, index] = mass**2
        for index, (first, second) in enumerate(PAIR_LABELS):
            gram[..., first, second] = gram[..., second, first] = (
                pairs[..., index]-self.masses[first]**2-self.masses[second]**2)/2
        return np.maximum(0., -np.linalg.det(gram)) * np.abs(self.scalar(pairs))**2

    def amplitude(self, momenta):
        """Return M for four contravariant momenta (E,px,py,pz).

        Uses epsilon_0123=-1 for metric (+---), epsilon^0123=+1.
        A common overall epsilon sign does not affect the partial width.
        """
        momenta = np.asarray(momenta, dtype=float)
        pairs = []
        for first, second in PAIR_LABELS:
            total = momenta[..., first, :] + momenta[..., second, :]
            pairs.append(total[..., 0]**2-np.sum(total[..., 1:]**2, axis=-1))
        return -np.linalg.det(momenta) * self.scalar(np.stack(pairs, axis=-1))


def build_four_pseudoscalar_amplitude(alp_mass, c_q_A, c_q_P, final_state, f_a, *,
                                      fit_parameters=None, g_v=None, epsilon_k=EPSILON_K):
    """Prepare a selected ordered four-pseudoscalar amplitude.

    Physical external masses and Alda:2026nrd generators are used; the strong
    odd numerator uses the emitted meson's LO mass. Charged poles follow the
    PDG input of the a -> 3P implementation; neutral poles follow the fit.
    Zero terms encode thresholds, charge restrictions, and vanishing grouped
    flavor coefficients. K_S/K_L are projected coherently before aggregation.
    """
    state = tuple(map(_canonical_name, final_state))
    if len(state) != 4:
        raise ValueError('final_state must contain four pseudoscalars')
    masses = tuple(_MESONS[name].mass for name in state)
    axial, pseudoscalar = _couplings(c_q_A), _couplings(c_q_P)
    empty = FourPseudoscalarAmplitude(float(alp_mass), masses, state, ())
    if alp_mass <= sum(masses) or sum(_MESONS[name].charge for name in state):
        return empty
    if not np.any(axial) and not np.any(pseudoscalar):
        return empty
    if all(name in ('pi0', 'eta', 'eta_prime') for name in state):
        return empty
    fitted = fit_parameters if fit_parameters is not None else read_preliminary_fit_parameters()
    g_v = read_four_pseudoscalar_g_v() if g_v is None else float(g_v)
    pion_constant = json.loads(SNAPSHOT.read_text())['fpi_gev']
    generators, tree_masses = _meson_inputs(epsilon_k)
    vector_components = []
    for name in ('rho0', 'omega', 'phi'):
        pole = DEFAULT_VECTOR_POLES[name]
        vector_components.append((fitted.vector_generators[name], fitted.vector_masses[name]**2
                                  - 1j*fitted.vector_masses[name]*pole.width))
    for row, column, name in ((0, 1, 'rho+'), (0, 2, 'K*+'), (1, 2, 'K*0')):
        pole = DEFAULT_VECTOR_POLES[name]
        vector_components.extend((generator, pole.mass**2-1j*pole.mass*pole.width)
                                 for generator in _real_components(row, column))
    odd = fitted.odd_coefficients
    chi = fitted.quark_mass_matrix
    contributions = {}

    def add(first, second, first_pole, second_pole, factor, derivative, mass_trace, source_mass_squared):
        if factor == 0 or (derivative == 0 and mass_trace == 0):
            return
        constant = factor*(4*(odd['d1']-odd['d3'])*source_mass_squared*derivative
                           - 16*odd['d2']*mass_trace)
        linear = factor*4*odd['d3']*derivative
        if (second, second_pole.real, second_pole.imag) < (first, first_pole.real, first_pole.imag):
            first, second, first_pole, second_pole = second, first, second_pole, first_pole
        values = contributions.setdefault((first, second, first_pole, second_pole), ([], []))
        values[0].append(constant)
        values[1].append(linear)

    for projection in product(*(_flavor_projection(name, epsilon_k) for name in state)):
        names, weights = zip(*projection)
        weight = math.prod(weights)
        outgoing = [generators[name].conj().T for name in names]
        pair_traces = {(first, second): [np.trace(vector @ _commutator(outgoing[first], outgoing[second]))
                                         for vector, _ in vector_components]
                       for first, second in PAIR_LABELS}
        axial_traces = [[np.trace(vector @ _commutator(axial, meson))
                         for vector, _ in vector_components] for meson in outgoing]
        for first_vector, (first_generator, first_pole) in enumerate(vector_components):
            for second_vector, (second_generator, second_pole) in enumerate(vector_components):
                anticommutator = _anticommutator(first_generator, second_generator)
                derivative = np.trace(anticommutator @ axial)/f_a
                mass_trace = -np.trace(anticommutator @ _anticommutator(pseudoscalar, chi))/f_a
                for first_pair, second_pair, sign in (((0, 1), (2, 3), 1),
                                                     ((0, 2), (1, 3), -1),
                                                     ((0, 3), (1, 2), 1)):
                    factor = (-weight*sign*g_v**2/(2*pion_constant**4)
                              *pair_traces[first_pair][first_vector]
                              *pair_traces[second_pair][second_vector])
                    add(PAIR_INDEX[first_pair], PAIR_INDEX[second_pair], first_pole, second_pole,
                        factor, derivative, mass_trace, alp_mass**2)
                for first, second in permutations(range(4), 2):
                    third, fourth = [index for index in range(4) if index not in (first, second)]
                    factor = (weight*_parity((first, second, third, fourth))*g_v**2
                              /(pion_constant**3*f_a)*axial_traces[first][first_vector]
                              *pair_traces[third, fourth][second_vector])
                    if factor == 0:
                        continue
                    derivative = -np.trace(anticommutator @ outgoing[second])/(2*pion_constant)
                    mass_trace = np.trace(anticommutator @ _anticommutator(outgoing[second], chi))/(2*pion_constant)
                    add(6+first, PAIR_INDEX[third, fourth], first_pole, second_pole,
                        factor, derivative, mass_trace, tree_masses[names[second]]**2)
    terms = []
    for key, (constants, linears) in contributions.items():
        constant = _aggregate_flavor_coefficients(constants)
        linear = _aggregate_flavor_coefficients(linears)
        if constant != 0 or linear != 0:
            terms.append((*key, constant, linear, linear))
    return FourPseudoscalarAmplitude(float(alp_mass), masses, state, tuple(terms))


def four_pseudoscalar_integrand(amplitude):
    """Return vectorized width integrand on the five-dimensional unit cube.

    Coordinates are s12, s34, cos(theta1), cos(theta2), phi mapped to [0,1].
    The exact epsilon squared and phase-space Jacobian are evaluated directly;
    no Gram determinants, flavor matrices or four-vectors enter cubature.
    """
    parent_squared = amplitude.alp_mass**2
    first, second, third, fourth = np.square(amplitude.masses)
    lower_pair = (amplitude.masses[0]+amplitude.masses[1])**2
    upper_pair = (amplitude.alp_mass-amplitude.masses[2]-amplitude.masses[3])**2
    lower_second_pair = (amplitude.masses[2]+amplitude.masses[3])**2
    normalization = 1/(32768*math.pi**6*amplitude.alp_mass**3*amplitude.symmetry_factor)

    def integrand(coordinates):
        coordinates = np.asarray(coordinates)
        pair = lower_pair+(upper_pair-lower_pair)*coordinates[..., 0]
        upper_second_pair = (amplitude.alp_mass-np.sqrt(pair))**2
        other_pair = lower_second_pair+(upper_second_pair-lower_second_pair)*coordinates[..., 1]
        first_cosine, second_cosine = 2*coordinates[..., 2]-1, 2*coordinates[..., 3]-1
        azimuth = 2*math.pi*coordinates[..., 4]
        parent_kallen = np.maximum(0., _kallen(parent_squared, pair, other_pair))
        first_kallen = np.maximum(0., _kallen(pair, first, second))
        second_kallen = np.maximum(0., _kallen(other_pair, third, fourth))
        first_sine = np.sqrt(np.maximum(0., 1-first_cosine**2))
        second_sine = np.sqrt(np.maximum(0., 1-second_cosine**2))
        first_fraction = (pair+first-second)/(2*pair)
        third_fraction = (other_pair+third-fourth)/(2*other_pair)
        pair_dot = (parent_squared-pair-other_pair)/2
        longitudinal_dot = pair_dot*first_cosine*second_cosine
        transverse_dot = np.sqrt(pair*other_pair)*first_sine*second_sine*np.cos(azimuth)
        first_third_dot = (first_fraction*third_fraction*pair_dot
            + first_fraction*np.sqrt(second_kallen)*np.sqrt(parent_kallen)*second_cosine/(4*other_pair)
            + third_fraction*np.sqrt(first_kallen)*np.sqrt(parent_kallen)*first_cosine/(4*pair)
            + np.sqrt(first_kallen*second_kallen)*(longitudinal_dot-transverse_dot)/(4*pair*other_pair))
        first_other_dot = first_fraction*pair_dot+np.sqrt(first_kallen*parent_kallen)*first_cosine/(4*pair)
        third_pair_dot = third_fraction*pair_dot+np.sqrt(second_kallen*parent_kallen)*second_cosine/(4*other_pair)
        invariants = np.stack((pair, first+third+2*first_third_dot,
            first+fourth+2*(first_other_dot-first_third_dot),
            second+third+2*(third_pair_dot-first_third_dot),
            second+fourth+2*(pair_dot-first_other_dot-third_pair_dot+first_third_dot), other_pair), axis=-1)
        epsilon_squared = (parent_kallen*first_kallen*second_kallen/(64*pair*other_pair)
                           *first_sine**2*second_sine**2*np.sin(azimuth)**2)
        jacobian = ((upper_pair-lower_pair)*(upper_second_pair-lower_second_pair)*8*math.pi)
        measure = np.sqrt(parent_kallen*first_kallen*second_kallen)/(pair*other_pair)
        return normalization*jacobian*measure*epsilon_squared*np.abs(amplitude.scalar(invariants))**2
    return integrand


def alp_four_pseudoscalar_width(alp_mass, c_q_A, c_q_P, final_state, f_a, *,
                                fit_parameters=None, g_v=None, epsilon_k=EPSILON_K,
                                rtol=1e-3, atol=0., max_subdivisions=10000,
                                return_info=False):
    """Return the partial width in GeV, optionally with SciPy cubature info.

    Forbidden channels return zero before invoking the integrator. ``rtol``
    and ``atol`` control deterministic Genz--Malik adaptive cubature. A failed
    convergence emits RuntimeWarning; inspect ``return_info=True`` for the
    numerical error and status. The zero-channel info is None.
    """
    amplitude = build_four_pseudoscalar_amplitude(
        alp_mass, c_q_A, c_q_P, final_state, f_a,
        fit_parameters=fit_parameters, g_v=g_v, epsilon_k=epsilon_k)
    if not amplitude.terms:
        return (0., None) if return_info else 0.
    result = cubature(four_pseudoscalar_integrand(amplitude), np.zeros(5), np.ones(5),
                      rule='genz-malik', rtol=rtol, atol=atol,
                      max_subdivisions=max_subdivisions)
    if result.status != 'converged':
        warnings.warn(f'four-pseudoscalar cubature: {result.status}; estimated absolute error {result.error}',
                      RuntimeWarning, stacklevel=2)
    width = float(result.estimate)
    return (width, result) if return_info else width
