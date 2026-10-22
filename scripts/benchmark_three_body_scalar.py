"""Validate and time approved production scalar preparation for a -> 3P and a -> VPP.

All masses, decay constants and widths are GeV. Couplings below are assumed
illustrative inputs. The reference and optimized widths retain their original
adaptive tolerances and finite-width models. No approximation is introduced.

Tests: .venv/bin/python -m unittest scripts.benchmark_three_body_scalar
Benchmark: .venv/bin/python -m scripts.benchmark_three_body_scalar
Timings are medians of three calls and include amplitude preparation in widths.
"""
from collections import Counter
from dataclasses import dataclass
from itertools import permutations
from math import factorial

from scipy.integrate import quad

import json
import math
import statistics
import time
import unittest

import numpy as np
import sympy as sp

from alpaca_chiral import (
    ALPcouplings, three_pseudoscalar_amplitude, three_pseudoscalar_width,
    vector_two_pseudoscalar_amplitude, vector_two_pseudoscalar_width,
)
from alpaca_chiral.three_pseudoscalar import F_PI, _pseudoscalars, _vector_poles
from alpaca_chiral.pseudoscalar_diag import _no_mixing_couplings, fit_results
from alpaca_chiral.vector_two_pseudoscalar import _VectorTwoPseudoscalarAmplitude


def _commutator(first, second):
    """Return the flavour commutator used by the archived 3P evaluator."""
    return first @ second - second @ first


def reference_three_pseudoscalar_amplitude(ma: float, fa: float, couplings: ALPcouplings,
                                   meson_1: str, meson_2: str, meson_3: str):
    """Return callable M(s12,s23), including G_V exchange and matching.

    The ordered invariants are s12=(p1+p2)^2 and s23=(p2+p3)^2.
    Supported particles: pi+, pi-, pi0, K+, K-, K0, K0bar, KS, KL,
    eta, eta_prime. Outgoing generators are Hermitian conjugated.
    All six labelled permutations retain identical-particle multiplicities.
    """
    particles = _pseudoscalars()
    masses, generators = zip(*(particles[name] for name in (meson_1, meson_2, meson_3)))
    axial, pseudoscalar, _ = _no_mixing_couplings(couplings, ma)
    quark_masses = np.asarray(fit_results['B0mq'], dtype=complex)
    mass_source = pseudoscalar @ quark_masses + quark_masses @ pseudoscalar
    derivative_traces = [0j, 0j, 0j]
    mass_trace = 0j
    for first, second, differentiated in permutations(range(3)):
        derivative_traces[differentiated] += np.trace(
            axial @ _commutator(generators[first],
                                _commutator(generators[second], generators[differentiated])))
        mass_trace += np.trace(generators[first] @ generators[second]
                              @ generators[differentiated] @ mass_source)
    exchange_coefficients = {}
    for spectator in range(3):
        first, second = (index for index in range(3) if index != spectator)
        for mass, width, generator in _vector_poles().values():
            flavour = np.trace(generator.conj().T @ _commutator(generators[first], generators[second]))
            flavour *= np.trace(generator @ _commutator(generators[spectator], axial))
            if flavour != 0:
                key = spectator, first, second, mass, width
                exchange_coefficients.setdefault(key, []).append(flavour)
    # Cancel flavour coefficients sharing the same pole and kinematics first.
    exchanges = []
    for key, coefficients in exchange_coefficients.items():
        flavour = sum(coefficients)
        coefficient_scale = sum(abs(coefficient) for coefficient in coefficients)
        if abs(flavour) > 64*np.finfo(float).eps*coefficient_scale:
            exchanges.append((*key, flavour))

    def amplitude(s12, s23):
        """Evaluate the invariant amplitude for the requested ordered channel."""
        invariants = {(0, 1): s12, (1, 2): s23,
                      (0, 2): ma**2+sum(mass**2 for mass in masses)-s12-s23}
        result = -mass_trace/(6*fa*F_PI)
        for differentiated in range(3):
            pair = tuple(index for index in range(3) if index != differentiated)
            result += derivative_traces[differentiated]*(
                ma**2+masses[differentiated]**2-invariants[pair])/(12*fa*F_PI)
        for spectator, first, second, mass, width, flavour in exchanges:
            pair_invariant = invariants[(first, second)]
            first_invariant = invariants[tuple(sorted((first, spectator)))]
            second_invariant = invariants[tuple(sorted((second, spectator)))]
            first_squared, second_squared, spectator_squared = (
                masses[index]**2 for index in (first, second, spectator))
            pole_numerator = pair_invariant/4*(
                (pair_invariant+second_squared-first_squared)*(first_invariant-first_squared-spectator_squared)
                -(pair_invariant+first_squared-second_squared)*(second_invariant-second_squared-spectator_squared))
            matching = ((first_invariant-first_squared-spectator_squared)*(ma**2+second_squared-first_invariant)
                        -(second_invariant-second_squared-spectator_squared)*(ma**2+first_squared-second_invariant))/4
            result += flavour/(2*F_PI*fa*mass**2)*(
                pole_numerator/(mass**2-pair_invariant-1j*mass*width)+matching)
        return complex(result)
    return amplitude



def reference_three_pseudoscalar_width(ma: float, fa: float, couplings: ALPcouplings,
                              meson_1: str, meson_2: str, meson_3: str) -> float:
    """Return Gamma(a -> P1 P2 P3) in GeV, including vector interference.

    Identical particles contribute 1/prod(n!). Numerical integration uses
    relative tolerance 1e-5. Thresholds use observed pseudoscalar masses.
    """
    names = meson_1, meson_2, meson_3
    particles = _pseudoscalars()
    first_mass, second_mass, third_mass = (particles[name][0] for name in names)
    if ma <= first_mass+second_mass+third_mass:
        return 0.0
    amplitude = reference_three_pseudoscalar_amplitude(ma, fa, couplings, *names)

    def integrate_second_invariant(s12):
        """Integrate over the second Dalitz invariant at fixed s12."""
        center = second_mass**2+third_mass**2+(
            (s12+second_mass**2-first_mass**2)*(ma**2-s12-third_mass**2))/(2*s12)
        product = ((s12-(first_mass+second_mass)**2)*(s12-(first_mass-second_mass)**2)
                   *((ma-third_mass)**2-s12)*((ma+third_mass)**2-s12))
        half_width = np.sqrt(max(0., product))/(2*s12)
        return quad(lambda s23: abs(amplitude(s12, s23))**2,
                    center-half_width, center+half_width, epsabs=0., epsrel=1e-5, limit=200)[0]

    integral = quad(integrate_second_invariant, (first_mass+second_mass)**2,
                    (ma-third_mass)**2, epsabs=0., epsrel=1e-5, limit=200)[0]
    symmetry = np.prod([factorial(number) for number in Counter(names).values()])
    return float(integral/(256*np.pi**3*ma**3*symmetry))


@dataclass(frozen=True)
class _ReferenceVectorTwoPseudoscalarAmplitude:
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


def reference_vector_two_pseudoscalar_amplitude(*arguments):
    """Use production flavour preparation with the former pointwise evaluator."""
    amplitude = vector_two_pseudoscalar_amplitude(*arguments)
    return _ReferenceVectorTwoPseudoscalarAmplitude(
        amplitude.ma, amplitude.masses, amplitude.final_state, amplitude.terms)



def reference_vector_two_pseudoscalar_width(ma: float, fa: float, couplings: ALPcouplings,
                                  vector: str, meson_1: str, meson_2: str) -> float:
    """Return Gamma(a -> V P1 P2) in GeV for one selected channel.

    Uses nested Dalitz quadrature with relative tolerance 1e-5, internal pole
    subdivisions and the identical-pseudoscalar factor 1/2!. Threshold and
    charge-forbidden channels return zero before integration. SciPy reports
    quadrature convergence warnings to the caller.
    """
    amplitude = reference_vector_two_pseudoscalar_amplitude(ma, fa, couplings, vector, meson_1, meson_2)
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

def coupling_inputs():
    """Return diagonal and complex flavour-changing axial/vector sources."""
    zero = np.zeros((3, 3))
    diagonal = ALPcouplings(.1, np.diag([.8, 0., 0.]), zero,
                            np.diag([-.3, .15, 0.]), zero)
    changing = zero.astype(complex)
    changing[0, 1], changing[1, 0] = .2+.1j, .2-.1j
    mixed = ALPcouplings(.1, np.diag([.8, 0., 0.]), zero, changing, changing)
    return diagonal, mixed


def cases():
    """Return family, mass, coupling and ordered final-state benchmark inputs.

    Cover below-pole and accessible-vector-pole regions, repeated particles,
    coherent neutral kaons, and both odd and even VPP contributions.
    """
    diagonal, changing = coupling_inputs()
    return [
        ('3P', .6, diagonal, ('pi0', 'pi+', 'pi-')),
        ('3P', 1.5, diagonal, ('pi0', 'pi+', 'pi-')),
        ('3P', 1.7, diagonal, ('pi0', 'K+', 'K-')),
        ('3P', 1.3, diagonal, ('pi0', 'pi0', 'eta')),
        ('3P', 1.5, diagonal, ('pi0', 'KS', 'KL')),
        ('3P', 1.3, changing, ('pi+', 'pi-', 'KS')),
        ('VPP', 1.1, diagonal, ('omega', 'pi+', 'pi-')),
        ('VPP', 1.8, diagonal, ('omega', 'pi+', 'pi-')),
        ('VPP', 2.2, diagonal, ('phi', 'K+', 'K-')),
        ('VPP', 1.8, changing, ('omega', 'K+', 'pi-')),
        ('VPP', 1.8, changing, ('omega', 'pi0', 'KS')),
    ]


def functions(family):
    """Return local slow references and approved production factories."""
    if family == '3P':
        return (reference_three_pseudoscalar_amplitude, three_pseudoscalar_amplitude,
                reference_three_pseudoscalar_width, three_pseudoscalar_width)
    return (reference_vector_two_pseudoscalar_amplitude, vector_two_pseudoscalar_amplitude,
            reference_vector_two_pseudoscalar_width, vector_two_pseudoscalar_width)


def t_limits(parent, masses, s):
    """Return second-invariant boundaries for ordered V/P spectator, P1, P2."""
    spectator, first, second = masses
    center = spectator*spectator+second*second + (
        (parent*parent-spectator*spectator-s)*(s+second*second-first*first)/(2*s))
    product = ((parent-spectator)**2-s)*((parent+spectator)**2-s)
    product *= (s-(first+second)**2)*(s-(first-second)**2)
    half_range = math.sqrt(max(0., product))/(2*s)
    return center-half_range, center+half_range


def grouped_points(family, mass, state):
    """Return interior, near-boundary and accessible-pole Dalitz coordinates."""
    particles = _pseudoscalars()
    if family == '3P':
        first, second, third = (particles[name][0] for name in state)
        masses = third, first, second
    else:
        masses = ({'omega': .78266, 'phi': 1.019461}[state[0]],
                  particles[state[1]][0], particles[state[2]][0])
    lower_s, upper_s = (masses[1]+masses[2])**2, (mass-masses[0])**2
    values = [lower_s+(upper_s-lower_s)*fraction for fraction in (.0001, .1, .4, .8, .9999)]
    values += [pole_mass**2 for pole_mass, _, _ in _vector_poles().values()
               if lower_s < pole_mass**2 < upper_s]
    points = []
    for s in sorted(set(values)):
        lower, upper = t_limits(mass, masses, s)
        # For 3P this helper gives s23 with the same spectator choice.
        points.append((s, tuple(lower+(upper-lower)*fraction
                                for fraction in (.0001, .2, .5, .9, .9999))))
    return points


class ThreeBodyScalarTests(unittest.TestCase):
    """Symbolic identities, pointwise agreement, complete widths and zeros."""

    def test_symbolic_three_pseudoscalar(self):
        parent, first, second, spectator, pair, first_cross, pole = sp.symbols(
            'parent_squared first_squared second_squared spectator_squared pair first_cross pole')
        second_cross = parent+first+second+spectator-pair-first_cross
        first_dot = first_cross-first-spectator
        second_dot = second_cross-second-spectator
        pole_numerator = pair*((pair+second-first)*first_dot
                               -(pair+first-second)*second_dot)/4
        matching = (first_dot*(parent+second-first_cross)
                    -second_dot*(parent+first-second_cross))/4
        reduced = (pair*(first_cross-second_cross)+(second-first)*(parent-spectator))/4
        self.assertEqual(sp.expand(matching-reduced), 0)
        self.assertEqual(sp.expand(pole_numerator-pair*matching), 0)
        self.assertEqual(sp.cancel(pole_numerator/(pole-pair)+matching
                                   -reduced*pole/(pole-pair)), 0)

    def test_symbolic_gram(self):
        parent, vector, first, second, s, t = sp.symbols(
            'parent_squared vector_squared first_squared second_squared s t', nonzero=True)
        u = parent+vector+first+second-s-t
        vector_first = (u-vector-first)/2
        vector_second = (t-vector-second)/2
        first_second = (s-first-second)/2
        gram = sp.Matrix([[vector, vector_first, vector_second],
                          [vector_first, first, first_second],
                          [vector_second, first_second, second]]).det()
        center = vector+second+(parent-vector-s)*(s+second-first)/(2*s)
        first_kallen = (parent-vector-s)**2-4*vector*s
        second_kallen = (s-first-second)**2-4*first*second
        factorized = s*(first_kallen*second_kallen/(4*s*s)-(t-center)**2)/4
        self.assertEqual(sp.factor(gram-factorized), 0)

    def test_pointwise(self):
        for family, mass, couplings, state in cases():
            reference_factory, optimized_factory, _, _ = functions(family)
            arguments = mass, 100., couplings, *state
            reference, optimized = reference_factory(*arguments), optimized_factory(*arguments)
            for s, values in grouped_points(family, mass, state):
                fixed_squared = optimized.at_fixed_s(s) if family == 'VPP' else None
                for t in values:
                    with self.subTest(family=family, mass=mass, state=state, s=s, t=t):
                        np.testing.assert_allclose(optimized(s, t), reference(s, t), rtol=2e-11, atol=1e-15)
                        if family == 'VPP':
                            np.testing.assert_allclose(optimized.squared(s, t), reference.squared(s, t),
                                                       rtol=2e-8, atol=1e-21)
                            np.testing.assert_allclose(fixed_squared(t), reference.squared(s, t),
                                                       rtol=2e-8, atol=1e-21)

    def test_vpp_boundary_zeros(self):
        diagonal = coupling_inputs()[0]
        for mass, state in ((1.1, ('omega', 'pi+', 'pi-')),
                            (1.8, ('omega', 'pi+', 'pi-')),
                            (2.2, ('phi', 'K+', 'K-'))):
            reference = reference_vector_two_pseudoscalar_amplitude(mass, 100., diagonal, *state)
            optimized = vector_two_pseudoscalar_amplitude(mass, 100., diagonal, *state)
            self.assertTrue(all(term[0] == 0 for term in reference.terms))
            for s, _ in grouped_points('VPP', mass, state):
                lower, upper = optimized.t_limits(s)
                interior = reference.squared(s, (lower+upper)/2)
                for t in (lower, upper):
                    self.assertEqual(optimized.squared(s, t), 0.)
                    self.assertLessEqual(reference.squared(s, t), 2e-9*interior)

    def test_vpp_identical_bose_symmetric_terms(self):
        # The retained physical identical PP channels vanish. Synthetic symmetric
        # odd and even coefficients also exercise nonzero identical-mass algebra.
        mass, daughter = 1.8, .1349768
        terms = ((0, 1, .775, .15, 1j, .2), (0, 2, .775, .15, -1j, -.2),
                 (1, 1, .5, 0., .3j, 0.), (2, 2, .5, 0., .3j, 0.))
        reference = _ReferenceVectorTwoPseudoscalarAmplitude(mass, (.78266, daughter, daughter),
                                                   ('omega', 'pi0', 'pi0'), terms)
        squared_masses = tuple(value*value for value in reference.masses)
        groups = tuple(tuple((index, complex(pole*pole, -pole*width), complex(constant), complex(linear))
                             for index, invariant, pole, width, constant, linear in terms if invariant == group)
                       for group in range(3))
        optimized = _VectorTwoPseudoscalarAmplitude(
            mass, reference.masses, reference.final_state, terms, groups, squared_masses,
            mass*mass+sum(squared_masses))
        for s, values in grouped_points('VPP', mass, reference.final_state):
            for t in values:
                np.testing.assert_allclose(optimized.squared(s, t), reference.squared(s, t),
                                           rtol=2e-8, atol=1e-18)

    def test_widths(self):
        for family, mass, couplings, state in cases():
            _, _, reference_width, optimized_width = functions(family)
            with self.subTest(family=family, mass=mass, state=state):
                reference = reference_width(mass, 100., couplings, *state)
                optimized = optimized_width(mass, 100., couplings, *state)
                if reference == 0.:
                    self.assertEqual(optimized, 0.)
                else:
                    # The official width now uses independently converged quadrature
                    # with relative tolerance 1e-5 rather than the same adaptive rule.
                    np.testing.assert_allclose(optimized, reference, rtol=3e-5, atol=0.)

    def test_thresholds_and_zeros(self):
        diagonal, changing = coupling_inputs()
        zero = np.zeros((3, 3))
        absent = ALPcouplings(0., zero, zero, zero, zero)
        inputs = [('3P', .3, diagonal, ('pi0', 'pi+', 'pi-')),
                  ('3P', 1.3, absent, ('pi0', 'pi+', 'pi-')),
                  ('VPP', 1., diagonal, ('omega', 'pi+', 'pi-')),
                  ('VPP', 1.8, diagonal, ('omega', 'pi0', 'pi0')),
                  ('VPP', 1.8, changing, ('omega', 'pi+', 'pi0')),
                  ('VPP', 1.8, absent, ('omega', 'pi+', 'pi-'))]
        for family, mass, couplings, state in inputs:
            _, _, reference, optimized = functions(family)
            self.assertEqual(reference(mass, 100., couplings, *state), 0.)
            self.assertEqual(optimized(mass, 100., couplings, *state), 0.)


def timed(function, repeats=3):
    """Return the value and median elapsed seconds of repeated evaluations."""
    durations = []
    for _ in range(repeats):
        start = time.perf_counter()
        value = function()
        durations.append(time.perf_counter()-start)
    return value, statistics.median(durations)


def benchmark():
    """Print per-channel integrand and complete-width timings as JSON."""
    records = []
    for family, mass, couplings, state in cases():
        reference_factory, optimized_factory, reference_width, optimized_width = functions(family)
        arguments = mass, 100., couplings, *state
        reference, optimized = reference_factory(*arguments), optimized_factory(*arguments)
        points = grouped_points(family, mass, state)

        def evaluate_reference():
            total = 0.
            for _ in range(50):
                for s, values in points:
                    for t in values:
                        total += reference.squared(s, t) if family == 'VPP' else abs(reference(s, t))**2
            return total

        def evaluate_optimized():
            total = 0.
            for _ in range(50):
                for s, values in points:
                    if family == 'VPP':
                        squared = optimized.at_fixed_s(s)
                        for t in values:
                            total += squared(t)
                    else:
                        for t in values:
                            total += abs(optimized(s, t))**2
            return total

        _, reference_integrand_seconds = timed(evaluate_reference)
        _, optimized_integrand_seconds = timed(evaluate_optimized)
        reference_value, reference_seconds = timed(lambda: reference_width(*arguments))
        optimized_value, optimized_seconds = timed(lambda: optimized_width(*arguments))
        records.append(dict(family=family, ma_GeV=mass, final_state=state,
                            width_GeV=reference_value,
                            relative_width_difference=(abs(optimized_value/reference_value-1)
                                                       if reference_value else 0.),
                            original_width_seconds=reference_seconds,
                            optimized_width_seconds=optimized_seconds,
                            width_speedup=reference_seconds/optimized_seconds,
                            integrand_speedup=reference_integrand_seconds/optimized_integrand_seconds))
    print(json.dumps(records, indent=2))


if __name__ == '__main__':
    result = unittest.main(exit=False, argv=['benchmark_three_body_scalar'])
    if not result.result.wasSuccessful():
        raise SystemExit(1)
    benchmark()
