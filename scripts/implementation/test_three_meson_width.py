"""Tests for the numerical ALP three-meson widths.

Invocation:
    .venv/bin/python -m unittest scripts.implementation.test_three_meson_width
"""

import unittest

import numpy as np

from scripts.implementation import decay_width


class ThreeMesonWidthTests(unittest.TestCase):
    def setUp(self) -> None:
        self.identity = np.eye(3)
        self.zero = np.zeros((3, 3))

    def test_forbidden_channels_short_circuit_before_coupling_validation(self) -> None:
        malformed = np.zeros((2, 2))
        self.assertEqual(
            decay_width(2.0, malformed, malformed, "pi+", "pi0", "pi0", -1.0),
            0.0,
        )
        self.assertEqual(
            decay_width(2.0, self.zero, self.zero, "K+", "pi-", "pi0", -1.0),
            0.0,
        )
        self.assertEqual(
            decay_width(0.3, malformed, malformed, "pi+", "pi-", "pi0", -1.0),
            0.0,
        )

    def test_invalid_couplings_are_rejected_for_allowed_channel(self) -> None:
        with self.assertRaisesRegex(ValueError, "shape"):
            decay_width(0.6, np.eye(2), self.identity, "pi+", "pi-", "pi0", 1.0)
        nonhermitian = np.array([[0, 1, 0], [0, 0, 0], [0, 0, 0]])
        with self.assertRaisesRegex(ValueError, "Hermitian"):
            decay_width(0.6, nonhermitian, self.identity, "pi+", "pi-", "pi0", 1.0)

    def test_width_is_zero_when_both_couplings_vanish(self) -> None:
        self.assertEqual(
            decay_width(0.7, self.zero, self.zero, "pi0", "pi0", "pi0", 10.0),
            0.0,
        )

    def test_width_scales_quadratically_with_coupling_and_inverse_fa(self) -> None:
        base = decay_width(
            0.7, self.zero, self.identity, "pi0", "pi0", "pi0", 10.0
        )
        doubled_coupling = decay_width(
            0.7, self.zero, 2.0 * self.identity, "pi0", "pi0", "pi0", 10.0
        )
        doubled_fa = decay_width(
            0.7, self.zero, self.identity, "pi0", "pi0", "pi0", 20.0
        )
        self.assertGreater(base, 0.0)
        self.assertAlmostEqual(doubled_coupling / base, 4.0, places=11)
        self.assertAlmostEqual(doubled_fa / base, 0.25, places=11)

    def test_two_identical_result_is_independent_of_input_order(self) -> None:
        coupling = np.diag([1.0, -0.3, 0.2])
        widths = [
            decay_width(1.4, coupling, coupling, *channel, 100.0)
            for channel in (
                ("pi0", "pi0", "eta"),
                ("pi0", "eta", "pi0"),
                ("eta", "pi0", "pi0"),
            )
        ]
        np.testing.assert_allclose(widths, widths[0], rtol=2e-12, atol=0.0)

    def test_three_distinct_result_is_permutation_invariant(self) -> None:
        c_a = np.array(
            [[0.4, 0.2j, 0.1], [-0.2j, -0.3, 0.05j], [0.1, -0.05j, 0.2]]
        )
        c_p = np.array(
            [[0.3, 0.1, -0.2j], [0.1, 0.2, 0.04], [0.2j, 0.04, -0.1]]
        )
        channels = (
            ("pi+", "pi-", "pi0"),
            ("pi0", "pi+", "pi-"),
            ("pi-", "pi0", "pi+"),
        )
        widths = [decay_width(0.8, c_a, c_p, *channel, 100.0) for channel in channels]
        self.assertGreater(widths[0], 0.0)
        np.testing.assert_allclose(widths, widths[0], rtol=2e-8, atol=0.0)

    def test_allowed_strange_channel_is_computed(self) -> None:
        coupling = np.diag([0.8, -0.4, 0.1])
        width = decay_width(
            1.5, coupling, coupling, "K+", "K-", "pi0", 100.0
        )
        self.assertGreaterEqual(width, 0.0)
        self.assertTrue(np.isfinite(width))

    def test_unit_strangeness_requires_a_nonzero_ds_coupling(self) -> None:
        channel = ("K+", "pi-", "pi0")
        self.assertEqual(
            decay_width(1.0, self.identity, self.identity, *channel, 100.0),
            0.0,
        )
        ds_coupling = self.zero.astype(complex)
        ds_coupling[1, 2] = 0.4 + 0.1j
        ds_coupling[2, 1] = 0.4 - 0.1j
        width = decay_width(1.0, ds_coupling, self.zero, *channel, 100.0)
        self.assertGreater(width, 0.0)


if __name__ == "__main__":
    unittest.main()
