"""Tests for numerical scalar Omnes functions.

Invocation:
    .venv/bin/python -m unittest scripts.implementation.test_omnes
"""

import unittest

import numpy as np

from scripts.implementation.omnes import OmnesFunction, PairOmnesFSI


class OmnesFunctionTests(unittest.TestCase):
    def setUp(self) -> None:
        self.threshold = 0.4
        self.endpoint = 2.0
        self.delta = 0.7
        self.constant = OmnesFunction(
            phase=lambda s: self.delta,
            threshold=self.threshold,
            phase_input_max=self.endpoint,
            subtraction_point=0.0,
            tail="constant",
            asymptotic_phase=self.delta,
        )

    def exact_constant_phase(self, s: float) -> complex:
        denominator = complex(self.threshold - s, -0.0 if s <= self.threshold else -1e-300)
        return ((self.threshold / denominator) ** (self.delta / np.pi))

    def test_subtraction_normalization(self) -> None:
        self.assertEqual(self.constant(0.0), 1.0 + 0.0j)

    def test_constant_phase_analytic_benchmark_below_cut(self) -> None:
        for s in (-0.3, 0.1, 0.25):
            np.testing.assert_allclose(
                self.constant(s), self.exact_constant_phase(s), rtol=2e-11, atol=1e-13
            )
            self.assertAlmostEqual(self.constant(s).imag, 0.0, places=14)

    def test_constant_phase_analytic_benchmark_on_cut(self) -> None:
        for s in (0.7, 1.3, 2.5):
            expected_magnitude = (self.threshold / abs(self.threshold - s)) ** (
                self.delta / np.pi
            )
            value = self.constant(s)
            self.assertAlmostEqual(abs(value), expected_magnitude, places=11)
            self.assertAlmostEqual(np.angle(value), self.delta, places=11)

    def test_nonconstant_phase_has_input_phase_on_cut(self) -> None:
        phase = lambda s: 0.2 + 0.3 * (s - self.threshold) / (
            self.endpoint - self.threshold
        )
        omega = OmnesFunction(
            phase=phase,
            threshold=self.threshold,
            phase_input_max=self.endpoint,
            tail="constant",
            asymptotic_phase=phase(self.endpoint),
        )
        for s in (0.6, 1.1, 1.8):
            self.assertAlmostEqual(np.angle(omega(s)), phase(s), places=10)

    def test_truncated_tail_and_endpoint_validation(self) -> None:
        omega = OmnesFunction(
            phase=lambda s: 0.4 * (self.endpoint - s),
            threshold=self.threshold,
            phase_input_max=self.endpoint,
            tail="truncate",
        )
        self.assertAlmostEqual(np.angle(omega(0.8)), 0.48, places=10)
        self.assertAlmostEqual(omega(2.5).imag, 0.0, places=13)
        with self.assertRaisesRegex(ValueError, "one-sided"):
            omega(self.threshold)
        with self.assertRaisesRegex(ValueError, "one-sided"):
            omega(self.endpoint)

    def test_invalid_tail_inputs_are_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "asymptotic_phase"):
            OmnesFunction(lambda s: 0.2, 0.4, 2.0)
        with self.assertRaisesRegex(ValueError, "continuously"):
            OmnesFunction(lambda s: 0.2, 0.4, 2.0, asymptotic_phase=0.3)
        with self.assertRaisesRegex(ValueError, "subtraction_point"):
            OmnesFunction(
                lambda s: 0.2,
                0.4,
                2.0,
                subtraction_point=0.5,
                asymptotic_phase=0.2,
            )

    def test_pair_dispatch_and_direct_fsi_signature(self) -> None:
        other = OmnesFunction(
            phase=lambda s: 0.3,
            threshold=0.42,
            phase_input_max=self.endpoint,
            asymptotic_phase=0.3,
        )
        dispatch = PairOmnesFSI(
            {("K+", "pi0"): self.constant, ("K0", "pi+"): other}
        )
        self.assertEqual(dispatch(0.2, ("pi0", "K+")), self.constant(0.2))
        self.assertEqual(self.constant(0.2, ("K+", "pi0")), self.constant(0.2))
        with self.assertRaisesRegex(KeyError, "no Omnes"):
            dispatch(0.2, ("K-", "pi0"))


if __name__ == "__main__":
    unittest.main()
