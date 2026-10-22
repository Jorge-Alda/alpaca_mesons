"""Verify and time photon PP scalar preparation; all masses and widths are GeV.

Run: .venv/bin/python -m scripts.benchmark_photon_pp_scalar
Run tests: .venv/bin/python -m unittest scripts.benchmark_photon_pp_scalar
The benchmark uses the retained odd model and the original adaptive tolerances.
Timings are medians of three runs; amplitude construction is included in widths.
"""
import json
import math
import statistics
import time
import unittest

import numpy as np
from scipy.integrate import quad

from alpaca_chiral import (
    ALPcouplings, photon_two_pseudoscalar_amplitude,
    photon_two_pseudoscalar_width,
)
from alpaca_chiral.photon_two_pseudoscalar import _photon_two_pseudoscalar_t_limits


class ReferencePhotonPPAmplitude:
    """Slow scalar evaluation from unprepared terms, retained only for tests.

    This reproduces the former evaluator independently of prepared poles,
    invariant groups, mass combinations and fixed-s evaluation.
    """

    def __init__(self, amplitude):
        self.ma = amplitude.ma
        self.masses = amplitude.masses
        self.final_state = amplitude.final_state
        self.terms = amplitude.terms

    def __call__(self, s, t):
        """Evaluate the former invariant sum in GeV^-3."""
        if not self.terms:
            return 0j
        invariants = (s, t, self.ma**2 + sum(m*m for m in self.masses)-s-t)
        return sum(((constant+linear*invariants[index]) /
                    (invariants[index]-mass**2+1j*mass*width)
                    for index, mass, width, constant, linear in self.terms), 0j)

    def squared(self, s, t):
        """Evaluate the former photon polarization sum."""
        if not self.terms:
            return 0.
        lower, upper = _photon_two_pseudoscalar_t_limits(self.ma, *self.masses, s)
        return max(0., s*(t-lower)*(upper-t)/4)*abs(self(s, t))**2


def reference_width(*arguments):
    """Former nested quadrature using the slow test-only scalar evaluator.

    Amplitude preparation is included, as in the public width timing. The
    comparison tests scalar arithmetic and integration, not flavour preparation.
    """
    amplitude = ReferencePhotonPPAmplitude(photon_two_pseudoscalar_amplitude(*arguments))
    if not amplitude.terms:
        return 0.
    ma, (first_mass, second_mass) = amplitude.ma, amplitude.masses
    total_squared = ma**2+first_mass**2+second_mass**2

    def integrate_t(s):
        """Integrate the former pointwise evaluator at fixed s."""
        lower, upper = _photon_two_pseudoscalar_t_limits(ma, first_mass, second_mass, s)
        points = set()
        for index, mass, _, _, _ in amplitude.terms:
            pole = mass**2 if index == 1 else total_squared-s-mass**2
            if index in (1, 2) and lower < pole < upper:
                points.add(pole)
        return quad(lambda t: amplitude.squared(s, t), lower, upper,
                    points=sorted(points), epsabs=0., epsrel=1e-5, limit=150)[0]

    lower, upper = (first_mass+second_mass)**2, ma**2
    points = sorted({mass**2 for index, mass, _, _, _ in amplitude.terms
                     if index == 0 and lower < mass**2 < upper})
    integral = quad(integrate_t, lower, upper, points=points,
                    epsabs=0., epsrel=1e-5, limit=150)[0]
    symmetry = 2 if amplitude.final_state[0] == amplitude.final_state[1] else 1
    return integral/(256*math.pi**3*ma**3*symmetry)


def coupling_inputs():
    """Deterministic diagonal and Hermitian ds axial sources."""
    diagonal = ALPcouplings(0.1, np.diag([0.8, 0., 0.]), np.zeros((3, 3)),
                            np.diag([-0.3, 0.15, 0.]), np.zeros((3, 3)))
    changing = ALPcouplings(0., np.zeros((3, 3)), np.zeros((3, 3)),
                            np.array([[0., .2+.1j, 0.], [.2-.1j, 0., 0.],
                                      [0., 0., 0.]]), np.zeros((3, 3)))
    return diagonal, changing


def cases():
    """Cover charged/neutral pairs, unequal masses, and coherent kaon states."""
    diagonal, changing = coupling_inputs()
    return [(0.6, diagonal, 'pi+', 'pi-'),
            (1.5, diagonal, 'pi+', 'pi-'),
            (1.5, diagonal, 'K+', 'K-'),
            (1.5, diagonal, 'KL', 'KS'),
            (1.5, changing, 'K+', 'pi-')]


def dalitz_points(amplitude):
    """Deterministic interior and endpoint coordinates, including pole regions."""
    first_mass, second_mass = amplitude.masses
    low_s, high_s = (first_mass+second_mass)**2, amplitude.ma**2
    values = [low_s+(high_s-low_s)*fraction for fraction in (.001, .1, .4, .8, .999)]
    values += [mass**2 for index, mass, _, _, _ in amplitude.terms
               if index == 0 and low_s < mass**2 < high_s]
    for s in values:
        low_t, high_t = _photon_two_pseudoscalar_t_limits(
            amplitude.ma, first_mass, second_mass, s)
        for fraction in (0., .001, .2, .5, .9, .999, 1.):
            yield s, low_t+(high_t-low_t)*fraction


class PhotonPPScalarTests(unittest.TestCase):
    """Compare independent reference evaluations and complete quadratures."""

    def test_pointwise(self):
        for mass, couplings, first, second in cases():
            arguments = mass, 100., couplings, first, second
            optimized = photon_two_pseudoscalar_amplitude(*arguments)
            original = ReferencePhotonPPAmplitude(optimized)
            self.assertTrue(original.terms)
            for s, t in dalitz_points(original):
                with self.subTest(mass=mass, pair=(first, second), s=s, t=t):
                    np.testing.assert_allclose(optimized(s, t), original(s, t), rtol=2e-13, atol=1e-15)
                    np.testing.assert_allclose(optimized.at_fixed_s(s)(t), original(s, t), rtol=2e-13, atol=1e-15)
                    np.testing.assert_allclose(optimized.squared(s, t), original.squared(s, t), rtol=2e-12, atol=1e-18)

    def test_widths(self):
        for mass, couplings, first, second in cases():
            arguments = mass, 100., couplings, first, second
            with self.subTest(mass=mass, pair=(first, second)):
                original = reference_width(*arguments)
                optimized = photon_two_pseudoscalar_width(*arguments)
                self.assertGreater(original, 0.)
                np.testing.assert_allclose(optimized, original, rtol=2e-10, atol=0.)

    def test_exact_zeros(self):
        couplings = coupling_inputs()[0]
        for mass, first, second in ((.2, 'pi+', 'pi-'), (1.5, 'pi+', 'pi0'),
                                    (1.5, 'KS', 'KS'), (1.5, 'pi0', 'eta')):
            arguments = mass, 100., couplings, first, second
            self.assertEqual(photon_two_pseudoscalar_width(*arguments), 0.)
            amplitude = photon_two_pseudoscalar_amplitude(*arguments)
            self.assertEqual(amplitude(object(), object()), 0j)
            self.assertEqual(amplitude.squared(object(), object()), 0.)


def timed(function, repeats=3):
    """Return the result and median elapsed seconds over repeated calls."""
    durations = []
    for _ in range(repeats):
        start = time.perf_counter()
        result = function()
        durations.append(time.perf_counter()-start)
    return result, statistics.median(durations)


def benchmark():
    """Print machine-readable timings and relative agreement for each channel."""
    records = []
    for mass, couplings, first, second in cases():
        arguments = mass, 100., couplings, first, second
        optimized = photon_two_pseudoscalar_amplitude(*arguments)
        original = ReferencePhotonPPAmplitude(optimized)
        points = list(dalitz_points(original))
        # Reproduce a nested integrator's reuse of each fixed s over many t calls.
        grouped_points = [(s, [t for other_s, t in points if other_s == s])
                          for s in dict.fromkeys(s for s, _ in points)]

        def original_integrands():
            return sum(original.squared(s, t) for _ in range(100)
                       for s, ts in grouped_points for t in ts)

        def optimized_integrands():
            total = 0.
            for _ in range(100):
                for s, ts in grouped_points:
                    lower, upper = optimized.t_limits(s)
                    scalar = optimized.at_fixed_s(s)
                    for t in ts:
                        total += max(0., s/4*(t-lower)*(upper-t))*abs(scalar(t))**2
            return total

        _, reference_integrand_time = timed(original_integrands)
        _, optimized_integrand_time = timed(optimized_integrands)
        reference_width_value, reference_time = timed(lambda: reference_width(*arguments))
        optimized_width, optimized_time = timed(lambda: photon_two_pseudoscalar_width(*arguments))
        records.append(dict(ma_GeV=mass, pair=[first, second], width_GeV=reference_width_value,
                            relative_width_difference=abs(optimized_width/reference_width_value-1),
                            original_width_seconds=reference_time,
                            optimized_width_seconds=optimized_time,
                            width_speedup=reference_time/optimized_time,
                            integrand_speedup=reference_integrand_time/optimized_integrand_time))
    print(json.dumps(records, indent=2))


if __name__ == '__main__':
    test_result = unittest.main(exit=False, argv=['benchmark_photon_pp_scalar'])
    if not test_result.result.wasSuccessful():
        raise SystemExit(1)
    benchmark()
