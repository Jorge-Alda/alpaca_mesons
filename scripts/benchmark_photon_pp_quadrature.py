"""Test and time hybrid photon PP quadrature; masses and widths are in GeV.

Tests: .venv/bin/python -m unittest scripts.benchmark_photon_pp_quadrature
Timing: .venv/bin/python -m scripts.benchmark_photon_pp_quadrature
A local copy of the former adaptive integration supplies the timing reference.
"""
import argparse
import json
import math
import statistics
import time
import unittest
from unittest.mock import patch

import numpy as np
from scipy.integrate import quad

from alpaca_chiral import ALPcouplings
from alpaca_chiral.photon_two_pseudoscalar import (
    photon_two_pseudoscalar_amplitude,
    photon_two_pseudoscalar_width,
)
from scripts.benchmark_photon_pp_scalar import coupling_inputs


def cases():
    """Include threshold, broad/narrow poles and unequal-mass crossed poles."""
    diagonal, changing = coupling_inputs()
    return [(.280, diagonal, 'pi+', 'pi-'),
            (.6, diagonal, 'pi+', 'pi-'),
            (1.05, diagonal, 'pi+', 'pi-'),
            (1.5, diagonal, 'pi+', 'pi-'),
            (1.05, diagonal, 'K+', 'K-'),
            (1.5, diagonal, 'K+', 'K-'),
            (1.5, diagonal, 'KL', 'KS'),
            (1.5, changing, 'K+', 'pi-'),
            (2., changing, 'K+', 'pi-')]


def adaptive_reference_width(
        ma: float, fa: float, couplings: ALPcouplings,
        meson_1: str, meson_2: str) -> float:
    """Return the former production adaptive width in GeV for comparisons.

    Uses the original width's tolerances, pole subdivisions and symmetry
    factor. Bounds and the s-channel scalar are computed once per outer
    quadrature call; complex poles and mass combinations once per width.
    No phase-space or propagator approximation is introduced.
    """
    amplitude = photon_two_pseudoscalar_amplitude(
        ma, fa, couplings, meson_1, meson_2)
    if not amplitude.terms:
        return 0.
    s_poles, t_poles, u_poles = (
        tuple(sorted({pole.real for pole, _, _ in group}))
        for group in amplitude.invariant_terms)

    def integrate_t(s):
        """Integrate the prepared polarization sum over t at fixed s."""
        lower, upper = amplitude.t_limits(s)
        scalar = amplitude.at_fixed_s(s)
        quarter_s = s/4
        points = sorted({pole for pole in t_poles if lower < pole < upper}
                        | {amplitude.total_squared-s-pole for pole in u_poles
                           if lower < amplitude.total_squared-s-pole < upper})

        def squared(t):
            """Polarization sum with fixed-s bounds already evaluated."""
            return max(0., quarter_s*(t-lower)*(upper-t))*abs(scalar(t))**2

        return quad(squared, lower, upper, points=points,
                    epsabs=0., epsrel=1e-5, limit=150)[0]

    lower, upper = amplitude.threshold_squared, ma*ma
    points = [pole for pole in s_poles if lower < pole < upper]
    integral = quad(integrate_t, lower, upper, points=points,
                    epsabs=0., epsrel=1e-5, limit=150)[0]
    symmetry = 2 if amplitude.final_state[0] == amplitude.final_state[1] else 1
    return integral/(256*math.pi**3*ma**3*symmetry)


def tightened_adaptive_width(*arguments):
    """Independent scalar quadrature at 2e-7 relative tolerance in both axes."""
    amplitude = photon_two_pseudoscalar_amplitude(*arguments)
    if not amplitude.terms:
        return 0.

    def integrate_t(s):
        lower, upper = amplitude.t_limits(s)
        poles = set()
        for invariant, mass, _, _, _ in amplitude.terms:
            point = mass**2 if invariant == 1 else amplitude.total_squared-s-mass**2
            if invariant in (1, 2) and lower < point < upper:
                poles.add(point)
        return quad(lambda t: amplitude.squared(s, t), lower, upper,
                    points=sorted(poles), epsabs=0., epsrel=2e-7, limit=200)[0]

    lower, upper = amplitude.threshold_squared, amplitude.ma**2
    poles = sorted({mass**2 for invariant, mass, _, _, _ in amplitude.terms
                    if invariant == 0 and lower < mass**2 < upper})
    integral = quad(integrate_t, lower, upper, points=poles,
                    epsabs=0., epsrel=2e-7, limit=200)[0]
    symmetry = 2 if amplitude.final_state[0] == amplitude.final_state[1] else 1
    return integral/(256*math.pi**3*amplitude.ma**3*symmetry)


class PhotonPPQuadratureTests(unittest.TestCase):
    """Compare complete coherent widths, order convergence and exact zeros."""

    def test_width_agreement_and_order_convergence(self):
        for mass, couplings, first, second in cases():
            arguments = mass, 100., couplings, first, second
            with self.subTest(mass=mass, pair=(first, second)):
                reference = adaptive_reference_width(*arguments)
                initial = photon_two_pseudoscalar_width(*arguments)
                refined = photon_two_pseudoscalar_width(
                    *arguments, order=32, epsrel=2e-6)
                tightened = tightened_adaptive_width(*arguments)
                self.assertGreater(reference, 0.)
                np.testing.assert_allclose(initial, tightened, rtol=1e-5, atol=0.)
                np.testing.assert_allclose(initial, reference, rtol=3e-5, atol=0.)
                np.testing.assert_allclose(refined, reference, rtol=3e-5, atol=0.)
                np.testing.assert_allclose(initial, refined, rtol=1e-5, atol=0.)

    def test_exact_zeros(self):
        diagonal = coupling_inputs()[0]
        for mass, first, second in ((.2, 'pi+', 'pi-'), (1.5, 'pi+', 'pi0'),
                                    (1.5, 'KS', 'KS'), (1.5, 'pi0', 'eta')):
            arguments = mass, 100., diagonal, first, second
            with self.subTest(mass=mass, pair=(first, second)):
                self.assertEqual(adaptive_reference_width(*arguments), 0.)
                self.assertEqual(photon_two_pseudoscalar_width(*arguments), 0.)

    def test_adaptive_fallback(self):
        arguments = .6, 100., coupling_inputs()[0], 'pi+', 'pi-'
        with patch('alpaca_chiral.photon_two_pseudoscalar.quad', wraps=quad) as adaptive:
            result = photon_two_pseudoscalar_width(
                *arguments, order=1, epsrel=1e-9)
        self.assertGreater(adaptive.call_count, 1)
        np.testing.assert_allclose(result, tightened_adaptive_width(*arguments),
                                   rtol=1e-7, atol=0.)

    def test_meson_exchange(self):
        couplings = coupling_inputs()[1]
        first = photon_two_pseudoscalar_width(1.5, 100., couplings, 'K+', 'pi-')
        second = photon_two_pseudoscalar_width(1.5, 100., couplings, 'pi-', 'K+')
        np.testing.assert_allclose(first, second, rtol=1e-5, atol=0.)


def timed_width(function, arguments, repeats):
    """Return width and median seconds, including amplitude preparation."""
    durations = []
    for _ in range(repeats):
        start = time.perf_counter()
        width = function(*arguments)
        durations.append(time.perf_counter()-start)
    return width, statistics.median(durations)


def benchmark(repeats=3):
    """Return deterministic case data and machine-dependent wall-clock times."""
    records = []
    for mass, couplings, first, second in cases():
        arguments = mass, 100., couplings, first, second
        # Warm imports, fitted inputs and caches before both measurements.
        photon_two_pseudoscalar_width(*arguments)
        reference, adaptive_seconds = timed_width(adaptive_reference_width, arguments, repeats)
        optimized, optimized_seconds = timed_width(photon_two_pseudoscalar_width,
                                                   arguments, repeats)
        records.append(dict(ma_GeV=mass, pair=[first, second], width_GeV=reference,
                            relative_difference=abs(optimized/reference-1),
                            adaptive_seconds=adaptive_seconds,
                            gauss_legendre_seconds=optimized_seconds,
                            speedup=adaptive_seconds/optimized_seconds))
    return dict(repeats=repeats, initial_order=16, epsrel=1e-5, results=records)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repeats', type=int, default=3)
    arguments = parser.parse_args()
    print(json.dumps(benchmark(arguments.repeats), indent=2))
