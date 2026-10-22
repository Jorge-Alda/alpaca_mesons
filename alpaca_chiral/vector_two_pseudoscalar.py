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
    """Return outgoing Alda generators (trace norm two) and LO masses in GeV."""
    generators = {name: generator
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
        return _prepare_amplitude(ma, masses, state, ())
    vector_source = np.zeros((3, 3), complex)
    vector_source[1, 2], vector_source[2, 1] = cqVS[1, 2], cqVS[2, 1]
    generators, tree_masses = _meson_inputs()
    first, second = generators[first_name], generators[second_name]
    external = fitted['vector_generators'][vector]
    chi = fitted['quark_mass_matrix']
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
                 -g_v/(math.sqrt(2)*F_PI**2*physical_vector_mass)),
                (1, -second/F_PI,
                 _anticommutator(second, chi)/F_PI,
                 tree_masses[second_name]**2, np.trace(internal @ _commutator(first, axial)),
                 math.sqrt(2)*g_v/(F_PI*fa*physical_vector_mass)),
                (2, -first/F_PI,
                 _anticommutator(first, chi)/F_PI,
                 tree_masses[first_name]**2, np.trace(internal @ _commutator(second, axial)),
                 -math.sqrt(2)*g_v/(F_PI*fa*physical_vector_mass))):
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
                    -1j*math.sqrt(2)*g_v*physical_vector_mass/(8*F_PI**2*fa)*flavor)
    terms = []
    for key, (constants, linears) in contributions.items():
        constant, linear = map(_sum_flavour_coefficients, (constants, linears))
        if constant != 0j or linear != 0j:
            terms.append((*key, constant, linear))
    return _prepare_amplitude(ma, masses, state, tuple(terms))

@dataclass(frozen=True)
class _VectorTwoPseudoscalarAmplitude:
    """Scalar-prepared coefficients and polarization sum for one VPP state.

    The original terms remain available for identical pole subdivisions.
    Prepared terms contain coefficient index, complex pole and numerator.
    Masses and widths use GeV throughout.
    """
    ma: float
    masses: tuple
    final_state: tuple
    terms: tuple
    prepared_terms: tuple
    squared_masses: tuple
    total_squared: float

    def t_limits(self, s):
        """Return the physical t interval in GeV squared at fixed s."""
        vector_squared, first_squared, second_squared = self.squared_masses
        parent_squared = self.ma*self.ma
        center = vector_squared+second_squared + (
            (parent_squared-vector_squared-s)*(s+second_squared-first_squared)/(2*s))
        first_kallen = (parent_squared-vector_squared-s)**2-4*vector_squared*s
        second_kallen = (s-first_squared-second_squared)**2-4*first_squared*second_squared
        half_range = math.sqrt(max(0., first_kallen)*max(0., second_kallen))/(2*s)
        return center-half_range, center+half_range

    def coefficients_at_fixed_s(self, s):
        """Prepare F_odd,F_1,F_2 evaluation with all s poles computed once."""
        fixed = [0j, 0j, 0j]
        for coefficient, pole, constant, linear in self.prepared_terms[0]:
            fixed[coefficient] += (constant+linear*s)/(s-pole)
        fixed_odd, fixed_first, fixed_second = fixed
        crossed_offset = self.total_squared-s
        t_terms, u_terms = self.prepared_terms[1:]

        def coefficients(t):
            """Evaluate the three coherent coefficients with scalar arithmetic."""
            odd, first, second = fixed_odd, fixed_first, fixed_second
            for coefficient, pole, constant, linear in t_terms:
                contribution = (constant+linear*t)/(t-pole)
                if coefficient == 0:
                    odd = odd+contribution
                elif coefficient == 1:
                    first = first+contribution
                else:
                    second = second+contribution
            u = crossed_offset-t
            for coefficient, pole, constant, linear in u_terms:
                contribution = (constant+linear*u)/(u-pole)
                if coefficient == 0:
                    odd = odd+contribution
                elif coefficient == 1:
                    first = first+contribution
                else:
                    second = second+contribution
            return odd, first, second

        return coefficients

    def __call__(self, s, t):
        """Return coherent F_odd,F_1,F_2 with the reference normalization."""
        return self.coefficients_at_fixed_s(s)(t)

    def at_fixed_s(self, s, lower=None, upper=None):
        """Prepare the scalar polarization sum for repeated inner evaluations.

        The exact Gram determinant is s*(t-t_minus)*(t_plus-t)/4.
        Its t dependence retains the angular zero at either boundary.
        """
        if lower is None:
            lower, upper = self.t_limits(s)
        coefficients = self.coefficients_at_fixed_s(s)
        vector_squared, first_squared, second_squared = self.squared_masses
        inverse_vector_squared = 1/vector_squared
        first_second = (s-first_squared-second_squared)/2
        first_offset = (self.total_squared-s-vector_squared-first_squared)/2
        second_offset = (vector_squared+second_squared)/2
        odd_factor = s/4

        def squared(t):
            """Evaluate the polarization sum without NumPy arrays or determinant."""
            vector_first = first_offset-t/2
            vector_second = t/2-second_offset
            odd, first, second = coefficients(t)
            gram = max(0., odd_factor*(t-lower)*(upper-t))
            transverse_first = vector_first*vector_first*inverse_vector_squared-first_squared
            transverse_second = vector_second*vector_second*inverse_vector_squared-second_squared
            transverse_cross = vector_first*vector_second*inverse_vector_squared-first_second
            return max(0., float(gram*abs(odd)**2 + transverse_first*abs(first)**2
                                 + transverse_second*abs(second)**2
                                 + 2*transverse_cross*(first*second.conjugate()).real))

        return squared

    def squared(self, s, t):
        """Return the pointwise polarization sum for direct callers."""
        return self.at_fixed_s(s)(t)


def _prepare_amplitude(ma, masses, final_state, terms):
    """Prepare complex poles and squared GeV masses for scalar evaluation."""
    groups = [[], [], []]
    for coefficient, invariant, mass, width, constant, linear in terms:
        groups[invariant].append((coefficient, complex(mass*mass, -mass*width),
                                  complex(constant), complex(linear)))
    squared_masses = tuple(mass*mass for mass in masses)
    return _VectorTwoPseudoscalarAmplitude(
        ma, masses, final_state, terms, tuple(tuple(group) for group in groups),
        squared_masses, ma*ma+sum(squared_masses))


def vector_two_pseudoscalar_width(
        ma: float, fa: float, couplings: ALPcouplings,
        vector: str, meson_1: str, meson_2: str, *,
        order: int = 16, epsrel: float = 1e-5) -> float:
    """Return Gamma(a -> VPP) in GeV using Gauss-Legendre quadrature.

    Retains the complete polarization sum and constant-width propagators.
    Pole subdivisions, order doubling, adaptive fallback and ``epsrel``
    follow ``three_pseudoscalar_width``. Closed and forbidden
    channels return zero before constructing any quadrature rule.
    """
    from ._three_body_quadrature import integrate_three_body

    amplitude = vector_two_pseudoscalar_amplitude(
        ma, fa, couplings, vector, meson_1, meson_2)
    if not amplitude.terms:
        return 0.
    vector_mass, first_mass, second_mass = amplitude.masses

    def limits(s):
        """Array-compatible physical t boundaries in GeV squared."""
        center = vector_mass**2+second_mass**2+(ma**2-vector_mass**2-s)*(s+second_mass**2-first_mass**2)/(2*s)
        product = ((ma-vector_mass)**2-s)*((ma+vector_mass)**2-s)
        product *= (s-(first_mass+second_mass)**2)*(s-(first_mass-second_mass)**2)
        half_range = np.sqrt(np.maximum(0., product))/(2*s)
        return center-half_range, center+half_range

    def fixed_squared(s):
        """Prepare the exact polarization sum with NumPy array evaluation."""
        lower, upper = limits(s)
        coefficients = amplitude.coefficients_at_fixed_s(s)
        vector_squared, first_squared, second_squared = amplitude.squared_masses
        first_second = (s-first_squared-second_squared)/2

        def squared(t):
            """Evaluate all inner nodes together, retaining interference."""
            vector_first = (amplitude.total_squared-s-t-vector_squared-first_squared)/2
            vector_second = (t-vector_squared-second_squared)/2
            odd, first, second = coefficients(t)
            gram = np.maximum(0., s/4*(t-lower)*(upper-t))
            return np.maximum(0., gram*np.abs(odd)**2
                +(vector_first**2/vector_squared-first_squared)*np.abs(first)**2
                +(vector_second**2/vector_squared-second_squared)*np.abs(second)**2
                +2*(vector_first*vector_second/vector_squared-first_second)
                *np.real(first*np.conjugate(second)))
        return squared

    poles = tuple(tuple(pole for _, pole, _, _ in group)
                  for group in amplitude.prepared_terms)
    integral = integrate_three_body(
        fixed_squared, limits, (first_mass+second_mass)**2,
        (ma-vector_mass)**2, amplitude.total_squared, poles,
        order=order, epsrel=epsrel)
    symmetry = 2 if amplitude.final_state[1] == amplitude.final_state[2] else 1
    return integral/(256*math.pi**3*ma**3*symmetry)
