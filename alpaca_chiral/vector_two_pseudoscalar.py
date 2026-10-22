"""Strong tree ALP decays to omega/phi and two pseudoscalars.

Masses, fa, G_V and widths are in GeV. External particles use PDG masses;
internal vectors use fitted masses and constant PDG widths. Pseudoscalar
poles and generators use the library's LO Alda:2026nrd inputs. The narrow
external vector approximation omits rescattering, loops, electromagnetic
conversion and the additional vector mass spurion interaction.
"""
from dataclasses import dataclass
import math

import numpy as np
from scipy.integrate import quad

from .pseudoscalar_diag import ALPcouplings, _no_mixing_couplings, cqV_no_mixing, fit_results, mass_sol
from .three_pseudoscalar import F_PI, _pseudoscalars
from .vector_inputs import _read_vector_inputs


def _anticommutator(first, second):
    """Return the flavour anticommutator."""
    return first @ second + second @ first


def _commutator(first, second):
    """Return the flavour commutator."""
    return first @ second - second @ first


def _sum_flavour_coefficients(coefficients):
    """Cancel roundoff residuals before evaluating a common propagator."""
    total = sum(coefficients)
    scale = sum(abs(coefficient) for coefficient in coefficients)
    return total if abs(total) > 64*np.finfo(float).eps*scale else 0j


def _real_components(row, column):
    """Hermitian generators with Tr(T_i T_j)=delta_ij."""
    transition = np.zeros((3, 3), dtype=complex)
    transition[row, column] = 1.
    return ((transition + transition.T) / math.sqrt(2),
            -1j * (transition - transition.T) / math.sqrt(2))


def _meson_inputs():
    """Return outgoing RChT generators and LO pole masses in GeV."""
    generators = {name: generator / math.sqrt(2)
                  for name, (_, generator) in _pseudoscalars().items()}
    neutral_squared, _ = mass_sol(fit_results['B0mu'], fit_results['B0md'],
                                  fit_results['B0ms'], fit_results['m02'])
    tree_masses = dict(zip(('pi0', 'eta', 'eta_prime'), np.sqrt(neutral_squared)))
    quark_masses = np.diag(fit_results['B0mq'])
    for names, row, column in ((('pi+', 'pi-'), 0, 1),
                               (('K+', 'K-'), 0, 2),
                               (('K0', 'K0bar', 'KS', 'KL'), 1, 2)):
        for name in names:
            tree_masses[name] = math.sqrt(quark_masses[row] + quark_masses[column])
    return generators, tree_masses


@dataclass(frozen=True)
class _VectorTwoPseudoscalarAmplitude:
    """Prepared invariant coefficients for one ordered V,P1,P2 state.

    ``terms`` stores (coefficient index, invariant index, pole mass, width,
    constant numerator, linear numerator). Invariant indices 0,1,2 denote
    s,t,u. Coefficient indices denote F_odd,F_1,F_2. All units are GeV.
    """
    ma: float
    masses: tuple
    final_state: tuple
    terms: tuple

    def __call__(self, s, t):
        """Return the coherent F_odd, F_1, F_2 coefficients at (s,t)."""
        invariants = (s, t, self.ma**2 + sum(m**2 for m in self.masses) - s - t)
        result = np.zeros(3, dtype=complex)
        for coefficient, invariant, mass, width, constant, linear in self.terms:
            virtuality = invariants[invariant]
            result[coefficient] += (constant + linear * virtuality) / (
                virtuality - mass**2 + 1j * mass * width)
        return tuple(result)

    def squared(self, s, t):
        """Return the polarization-summed squared amplitude (no spin average)."""
        vector_mass, first_mass, second_mass = self.masses
        u = self.ma**2 + sum(m**2 for m in self.masses) - s - t
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


def vector_two_pseudoscalar_amplitude(ma: float, fa: float, couplings: ALPcouplings,
                                      vector: str, meson_1: str, meson_2: str):
    """Return callable (F_odd,F_1,F_2) for a -> V P1 P2.

    V is omega or phi. P names follow three_pseudoscalar_amplitude.
    Invariants are s=(p1+p2)^2 and t=(qV+p2)^2. The amplitude is
    F_odd epsilon(epsV,qV,p1,p2)+F_1 epsV.p1+F_2 epsV.p2;
    the returned object also supplies its polarization sum as squared(s,t).
    ALPcouplings is rotated using pseudoscalar_diag.py, as in the other
    decay functions. G_V is taken from the central vector parameter file.
    Only ds/sd vector sources enter the retained |Delta S|=1 even amplitude.
    No additional ALP meson pole is added to the no-mixing sources.
    """
    axial, pseudoscalar, _ = _no_mixing_couplings(couplings, ma)
    cqVS = np.asarray(cqV_no_mixing(couplings, ma), dtype=complex)
    first_name, second_name = meson_1, meson_2
    fitted = _read_vector_inputs()
    g_v = fitted['even_parameters']['G_V_GeV']
    physical_vector_mass = {'omega': 0.78266, 'phi': 1.019461}[vector]
    particles = _pseudoscalars()
    masses = (physical_vector_mass, particles[first_name][0], particles[second_name][0])
    state = (vector, first_name, second_name)
    charges = {'pi+': 1, 'K+': 1, 'pi-': -1, 'K-': -1}
    if ma <= sum(masses) or charges.get(first_name, 0)+charges.get(second_name, 0):
        return _VectorTwoPseudoscalarAmplitude(ma, masses, state, ())
    vector_source = np.zeros((3, 3), complex)
    vector_source[1, 2], vector_source[2, 1] = cqVS[1, 2], cqVS[2, 1]
    generators, tree_masses = _meson_inputs()
    first, second = generators[first_name], generators[second_name]
    external = fitted['vector_generators'][vector]
    chi = fitted['quark_mass_matrix']
    operator_mass = fitted['operator_mass']
    odd = fitted['odd_coefficients']
    vector_groups = [(fitted['vector_masses'][name], {'rho0': .1474, 'omega': .00868, 'phi': .004249}[name],
                      (fitted['vector_generators'][name],)) for name in ('rho0', 'omega', 'phi')]
    for row, column, width_name in ((0, 1, 'rhoplus'), (0, 2, 'Kstarplus'), (1, 2, 'Kstar0')):
        mass = fitted['even_parameters']['m_rho_GeV' if row == 0 and column == 1 else 'm_Kstar_GeV']
        vector_groups.append((mass, {'rhoplus': .1491, 'Kstarplus': .0514, 'Kstar0': .0473}[width_name], _real_components(row, column)))
    contributions = {}

    def add(coefficient, invariant, mass, width, constant, linear=0j):
        """Collect numerator coefficients with the same pole and invariant."""
        values = contributions.setdefault((coefficient, invariant, mass, width), ([], []))
        values[0].append(constant)
        values[1].append(linear)

    for mass, width, components in vector_groups:
        for internal in components:
            anticommutator = _anticommutator(external, internal)
            for invariant, source, mass_source, source_mass_squared, trace, factor in (
                (0, 2*axial/fa, -2*_anticommutator(pseudoscalar, chi)/fa,
                 ma**2, np.trace(internal @ _commutator(first, second)),
                 -math.sqrt(2)*g_v/(F_PI**2*operator_mass)),
                (1, -math.sqrt(2)*second/F_PI,
                 -math.sqrt(2)*_anticommutator(second, chi)/F_PI,
                 tree_masses[second_name]**2, np.trace(internal @ _commutator(first, axial)),
                 2*g_v/(F_PI*fa*operator_mass)),
                (2, -math.sqrt(2)*first/F_PI,
                 -math.sqrt(2)*_anticommutator(first, chi)/F_PI,
                 tree_masses[first_name]**2, np.trace(internal @ _commutator(second, axial)),
                 -2*g_v/(F_PI*fa*operator_mass))):
                derivative_trace = np.trace(anticommutator @ source) / 2
                mass_trace = np.trace(anticommutator @ mass_source) / 2
                constant = 4*((odd['d1']*source_mass_squared
                               - odd['d3']*(source_mass_squared-operator_mass**2))*derivative_trace
                              - 4*odd['d2']*mass_trace)
                add(0, invariant, mass, width, factor*trace*constant,
                    factor*trace*4*odd['d3']*derivative_trace)
    pseudoscalar_groups = [(tree_masses[name], (generators[name],))
                           for name in ('pi0', 'eta', 'eta_prime')]
    for row, column, name in ((0, 1, 'pi+'), (0, 2, 'K+'), (1, 2, 'K0')):
        pseudoscalar_groups.append((tree_masses[name], _real_components(row, column)))
    for mass, components in pseudoscalar_groups:
        for internal in components:
            for coefficient, invariant, spectator, emitted in ((1, 2, second, first), (2, 1, first, second)):
                flavor = (np.trace(_anticommutator(spectator, internal) @ _commutator(vector_source, chi))
                          * np.trace(external @ _commutator(internal, emitted)))
                add(coefficient, invariant, mass, 0.,
                    -1j*g_v*operator_mass/(math.sqrt(2)*F_PI**2*fa)*flavor)
    terms = []
    for key, (constants, linears) in contributions.items():
        constant, linear = map(_sum_flavour_coefficients, (constants, linears))
        if constant != 0j or linear != 0j:
            terms.append((*key, constant, linear))
    return _VectorTwoPseudoscalarAmplitude(ma, masses, state, tuple(terms))


def vector_two_pseudoscalar_width(ma: float, fa: float, couplings: ALPcouplings,
                                  vector: str, meson_1: str, meson_2: str) -> float:
    """Return Gamma(a -> V P1 P2) in GeV for one selected channel.

    Uses nested Dalitz quadrature with relative tolerance 1e-5, internal pole
    subdivisions and the identical-pseudoscalar factor 1/2!. Threshold and
    charge-forbidden channels return zero before integration. SciPy reports
    quadrature convergence warnings to the caller.
    """
    amplitude = vector_two_pseudoscalar_amplitude(ma, fa, couplings, vector, meson_1, meson_2)
    epsabs, epsrel, limit = 0., 1e-5, 150
    if not amplitude.terms:
        return 0.
    ma = amplitude.ma
    vector_mass, first_mass, second_mass = amplitude.masses
    total_squared = ma**2 + sum(mass**2 for mass in amplitude.masses)

    def kallen(first, second, third):
        """Return the Kallen polynomial for squared invariant arguments."""
        return first**2 + second**2 + third**2 - 2*(first*second + first*third + second*third)

    def integrate_t(s):
        """Integrate the polarization sum over t at fixed s."""
        center = vector_mass**2 + second_mass**2 + (
            (ma**2-vector_mass**2-s)*(s+second_mass**2-first_mass**2)/(2*s))
        half_range = math.sqrt(max(0., kallen(ma**2, vector_mass**2, s))
                               * max(0., kallen(s, first_mass**2, second_mass**2))) / (2*s)
        lower, upper = center-half_range, center+half_range
        points = set()
        for _, invariant, mass, _, _, _ in amplitude.terms:
            pole = mass**2 if invariant == 1 else total_squared-s-mass**2
            if invariant in (1, 2) and lower < pole < upper:
                points.add(pole)
        return quad(lambda t: amplitude.squared(s, t), lower, upper,
                    points=sorted(points), epsabs=epsabs, epsrel=epsrel, limit=limit)[0]

    lower, upper = (first_mass+second_mass)**2, (ma-vector_mass)**2
    points = sorted({mass**2 for _, invariant, mass, _, _, _ in amplitude.terms
                     if invariant == 0 and lower < mass**2 < upper})
    integral = quad(integrate_t, lower, upper, points=points,
                    epsabs=epsabs, epsrel=epsrel, limit=limit)[0]
    identical_factor = 2 if amplitude.final_state[1] == amplitude.final_state[2] else 1
    return integral / (256*math.pi**3*ma**3*identical_factor)
