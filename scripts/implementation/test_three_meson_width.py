"""Tests for the numerical ALP three-meson widths.

Invocation:
    .venv/bin/python -m unittest scripts.implementation.test_three_meson_width
"""

import unittest
from unittest.mock import patch

import numpy as np
from scipy.integrate import quad_vec

from scripts.implementation import (
    EPSILON_K,
    build_amplitude,
    build_vector_exchange_amplitude,
    build_vector_matching_amplitude,
    decay_width,
)


class ThreeMesonWidthTests(unittest.TestCase):
    def setUp(self) -> None:
        self.identity = np.eye(3)
        self.zero = np.zeros((3, 3))

    @staticmethod
    def _manual_vector_width_components(ma, c_a, c_p, channel, fa):
        """Integrate coherent, ChPT-only, and vector-only squared amplitudes."""
        chipt = build_amplitude(ma, c_a, c_p, *channel, fa)
        pole = build_vector_exchange_amplitude(ma, c_a, c_p, *channel, fa)
        matching = build_vector_matching_amplitude(ma, c_a, c_p, *channel, fa)
        m1, m2, m3 = chipt.masses

        def kallen(x, y, z):
            return x * x + y * y + z * z - 2.0 * (x * y + x * z + y * z)

        def integrate_s23(s12):
            center = m2**2 + m3**2 + (
                (s12 + m2**2 - m1**2) * (ma**2 - s12 - m3**2)
            ) / (2.0 * s12)
            half_width = np.sqrt(
                max(
                    0.0,
                    kallen(s12, m1**2, m2**2)
                    * kallen(ma**2, s12, m3**2),
                )
            ) / (2.0 * s12)

            def components(s23):
                chipt_value = chipt(s12, s23)
                vector_value = pole(s12, s23) + matching(s12, s23)
                return np.array(
                    [
                        abs(chipt_value + vector_value) ** 2,
                        abs(chipt_value) ** 2,
                        abs(vector_value) ** 2,
                    ]
                )

            return quad_vec(
                components,
                center - half_width,
                center + half_width,
                epsabs=1e-25,
                epsrel=2e-9,
            )[0]

        integrals = quad_vec(
            integrate_s23,
            (m1 + m2) ** 2,
            (ma - m3) ** 2,
            epsabs=1e-25,
            epsrel=2e-8,
            limit=200,
        )[0]
        return integrals / (256.0 * np.pi**3 * ma**3)

    def test_default_path_is_backward_compatible_and_skips_vector_options(self) -> None:
        args = (
            0.8,
            np.diag([1.0, 0.0, 0.0]),
            self.zero,
            "pi0",
            "pi+",
            "pi-",
            1000.0,
        )
        expected = decay_width(*args)
        with patch(
            "scripts.implementation.vector_meson_exchange.build_vector_exchange_amplitude",
            side_effect=AssertionError("vector builder called"),
        ), patch(
            "scripts.implementation.vector_meson_exchange.build_vector_matching_amplitude",
            side_effect=AssertionError("matching builder called"),
        ):
            actual = decay_width(
                *args,
                include_vector_contributions=False,
                g_v="unused",
                pole_parameters="unused",
                vector_f_pi=-1.0,
            )
        self.assertEqual(actual, expected)

    def test_vector_flag_requires_a_bool(self) -> None:
        with self.assertRaisesRegex(TypeError, "must be a bool"):
            decay_width(
                0.8,
                self.identity,
                self.zero,
                "pi0",
                "pi+",
                "pi-",
                100.0,
                include_vector_contributions=1,
            )

    def test_enabled_vectors_match_coherent_integral_and_include_interference(self) -> None:
        ma, fa = 0.8, 1000.0
        c_a = np.diag([1.0, 0.0, 0.0])
        channel = ("pi0", "pi+", "pi-")
        coherent, chipt_only, vector_only = self._manual_vector_width_components(
            ma, c_a, self.zero, channel, fa
        )
        width = decay_width(
            ma,
            c_a,
            self.zero,
            *channel,
            fa,
            include_vector_contributions=True,
        )
        np.testing.assert_allclose(width, coherent, rtol=3e-8, atol=0.0)
        self.assertGreater(
            abs(width - (chipt_only + vector_only)),
            1e-3 * width,
        )

    def test_repeated_vector_width_has_symmetry_factor_and_permutation_invariance(self) -> None:
        c_a = np.array(
            [[0.4, 0.2j, 0.1], [-0.2j, -0.3, 0.05j], [0.1, -0.05j, 0.2]]
        )
        channel = ("pi0", "pi0", "K0")
        labeled_coherent = self._manual_vector_width_components(
            0.9, c_a, self.zero, channel, 100.0
        )[0]
        ordered = decay_width(
            0.9,
            c_a,
            self.zero,
            *channel,
            100.0,
            include_vector_contributions=True,
        )
        permuted = decay_width(
            0.9,
            c_a,
            self.zero,
            "pi_0",
            "K0",
            "pi0",
            100.0,
            include_vector_contributions=True,
        )
        np.testing.assert_allclose(ordered, labeled_coherent / 2.0, rtol=3e-8, atol=0.0)
        np.testing.assert_allclose(permuted, ordered, rtol=3e-8, atol=0.0)

    def test_vector_options_are_forwarded_and_validated_only_when_enabled(self) -> None:
        args = (
            0.8,
            np.diag([1.0, 0.0, 0.0]),
            self.zero,
            "pi0",
            "pi+",
            "pi-",
            1000.0,
        )
        chipt = decay_width(*args)
        zero_vector = decay_width(
            *args, include_vector_contributions=True, g_v=0.0
        )
        np.testing.assert_allclose(zero_vector, chipt, rtol=3e-9, atol=0.0)

        with self.assertRaisesRegex(TypeError, "real scalars"):
            decay_width(
                *args, include_vector_contributions=True, g_v="invalid"
            )
        with self.assertRaisesRegex(TypeError, "mapping"):
            decay_width(
                *args,
                include_vector_contributions=True,
                pole_parameters="invalid",
            )
        with self.assertRaisesRegex(ValueError, "f_pi"):
            decay_width(
                *args,
                include_vector_contributions=True,
                vector_f_pi=-1.0,
            )

    def test_all_identical_vector_contribution_is_zero(self) -> None:
        c_a = np.array(
            [[0.4, 0.2j, 0.1], [-0.2j, -0.3, 0.05j], [0.1, -0.05j, 0.2]]
        )
        c_p = np.diag([1.0, -0.3, 0.2])
        args = (0.7, c_a, c_p, "pi0", "pi0", "pi0", 100.0)
        chipt = decay_width(*args)
        with_vectors = decay_width(
            *args, include_vector_contributions=True
        )
        np.testing.assert_allclose(with_vectors, chipt, rtol=3e-9, atol=0.0)

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
        self.assertEqual(
            decay_width(
                0.3,
                malformed,
                malformed,
                "pi+",
                "pi-",
                "pi0",
                -1.0,
                include_vector_contributions=True,
                g_v="not inspected",
                pole_parameters="not inspected",
                vector_f_pi=-1.0,
            ),
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

    def test_two_identical_amplitude_keeps_input_mandelstam_order(self) -> None:
        coupling = np.diag([1.0, -0.3, 0.2])
        ordered = build_amplitude(
            1.4, coupling, coupling, "pi0", "pi0", "eta", 100.0
        )
        permuted = build_amplitude(
            1.4, coupling, coupling, "pi0", "eta", "pi0", 100.0
        )
        permuted_again = build_amplitude(
            1.4, coupling, coupling, "eta", "pi0", "pi0", 100.0
        )
        s12 = 0.3
        s23 = 0.5
        invariant_sum = ordered.alp_mass**2 + sum(mass**2 for mass in ordered.masses)
        s13 = invariant_sum - s12 - s23
        self.assertAlmostEqual(ordered(s12, s23), permuted(s13, s23))
        self.assertAlmostEqual(ordered(s12, s23), permuted_again(s13, s12))

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

    def test_three_distinct_amplitude_keeps_input_mandelstam_order(self) -> None:
        c_a = np.array(
            [[0.4, 0.2j, 0.1], [-0.2j, -0.3, 0.05j], [0.1, -0.05j, 0.2]]
        )
        c_p = np.array(
            [[0.3, 0.1, -0.2j], [0.1, 0.2, 0.04], [0.2j, 0.04, -0.1]]
        )
        ordered = build_amplitude(0.8, c_a, c_p, "pi+", "pi-", "pi0", 100.0)
        cyclic = build_amplitude(0.8, c_a, c_p, "pi0", "pi+", "pi-", 100.0)
        s12 = 0.12
        s23 = 0.18
        invariant_sum = ordered.alp_mass**2 + sum(mass**2 for mass in ordered.masses)
        s13 = invariant_sum - s12 - s23
        self.assertAlmostEqual(ordered(s12, s23), cyclic(s13, s12))

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

    def test_charged_mode_projects_equally_to_ks_and_kl(self) -> None:
        coupling = np.diag([0.8, -0.3, 0.15])
        flavor = decay_width(
            1.6, coupling, self.zero, "pi+", "K-", "K0", 100.0
        )
        short = decay_width(
            1.6, coupling, self.zero, "pi+", "K-", "k_s", 100.0
        )
        long = decay_width(
            1.6, coupling, self.zero, "pi+", "K-", "k_l", 100.0
        )
        self.assertGreater(flavor, 0.0)
        normalization = np.sqrt(2.0 * (1.0 + abs(EPSILON_K) ** 2))
        p_mixing = (1.0 + EPSILON_K) / normalization
        np.testing.assert_allclose(
            short, abs(p_mixing) ** 2 * flavor, rtol=2e-12, atol=0.0
        )
        np.testing.assert_allclose(
            long, abs(p_mixing) ** 2 * flavor, rtol=2e-12, atol=0.0
        )
        cp_limit = decay_width(
            1.6,
            coupling,
            self.zero,
            "pi+",
            "K-",
            "KS",
            100.0,
            epsilon_k=0.0,
        )
        np.testing.assert_allclose(cp_limit, 0.5 * flavor, rtol=2e-12, atol=0.0)

    def test_two_neutral_kaons_include_interference_and_symmetry(self) -> None:
        coupling = np.diag([0.8, -0.3, 0.15])
        flavor = decay_width(
            1.6, coupling, self.zero, "pi0", "K0", "K0bar", 100.0
        )
        short_short = decay_width(
            1.6, coupling, self.zero, "pi0", "KS", "KS", 100.0
        )
        long_long = decay_width(
            1.6, coupling, self.zero, "pi0", "KL", "KL", 100.0
        )
        short_long = decay_width(
            1.6, coupling, self.zero, "pi0", "KS", "KL", 100.0
        )
        normalization = np.sqrt(2.0 * (1.0 + abs(EPSILON_K) ** 2))
        p_mixing = (1.0 + EPSILON_K) / normalization
        q_mixing = (1.0 - EPSILON_K) / normalization
        expected = 2.0 * abs(p_mixing * q_mixing) ** 2 * flavor
        np.testing.assert_allclose(short_short, expected, rtol=2e-12, atol=0.0)
        np.testing.assert_allclose(long_long, expected, rtol=2e-12, atol=0.0)
        self.assertEqual(short_long, 0.0)


if __name__ == "__main__":
    unittest.main()
