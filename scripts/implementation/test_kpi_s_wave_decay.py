"""Tests for the charge-resolved K-pi S-wave implementation.

Invocation:
    .venv/bin/python -m unittest scripts.implementation.test_kpi_s_wave_decay
"""

import unittest

import numpy as np

from scripts.implementation.kpi_s_wave_decay import (
    _kernels_from_inputs,
    build_s_wave_model,
    kpi_thresholds,
    s_wave_decay_width,
    s_wave_decay_width_fast,
    supported_channels,
)
from scripts.implementation.omnes import OmnesFunction, PairOmnesFSI, TabulatedOmnes
from scripts.implementation.three_meson_width import B0_MQ, EPSILON_K, F_PI, _MESONS


class KPiSWaveDecayTests(unittest.TestCase):
    def setUp(self) -> None:
        self.c_a = np.diag([0.8, -0.3, 0.15])
        self.c_p = np.diag([0.2, 0.4, -0.1])
        self.zero = np.zeros((3, 3))

    def test_unity_fsi_restores_contact_amplitude(self) -> None:
        model = build_s_wave_model(
            1.6, self.c_a, self.c_p, ("pi0", "K-", "K+"), 100.0
        )
        s = 0.85
        t = 0.92
        unity = lambda invariant, pair: 1.0 + 0.0j
        self.assertAlmostEqual(model.amplitude(s, t, unity), model.contact_amplitude(s, t))

        g_k, g_pi, g_kbar = (_MESONS[name].generator for name in model.channel)
        m_k, m_pi, m_kbar = model.masses
        u = model.alp_mass**2 + m_k**2 + m_pi**2 + m_kbar**2 - s - t
        comm = lambda left, right: left @ right - right @ left
        tensor = lambda first, second, diff: np.trace(
            self.c_a
            @ (comm(first, comm(second, diff)) + comm(second, comm(first, diff)))
        )
        p_matrix = self.c_p @ B0_MQ + B0_MQ @ self.c_p
        ordered = sum(
            (a @ b @ c for a, b, c in __import__("itertools").permutations((g_k, g_pi, g_kbar))),
            np.zeros((3, 3), dtype=complex),
        )
        direct = (
            tensor(g_k, g_pi, g_kbar) * (model.alp_mass**2 + m_kbar**2 - s)
            + tensor(g_k, g_kbar, g_pi) * (model.alp_mass**2 + m_pi**2 - u)
            + tensor(g_pi, g_kbar, g_k) * (model.alp_mass**2 + m_k**2 - t)
        ) / (12.0 * model.f_a * F_PI) - np.trace(ordered @ p_matrix) / (
            6.0 * model.f_a * F_PI
        )
        np.testing.assert_allclose(model.contact_amplitude(s, t), direct, rtol=2e-13)

    def test_isospin_limit_matches_analytic_charged_kernel(self) -> None:
        sqrt2 = np.sqrt(2.0)
        k_plus = np.zeros((3, 3), dtype=complex)
        k_plus[0, 2] = sqrt2
        k_minus = k_plus.T
        pi_zero = np.diag([1.0, -1.0, 0.0]).astype(complex)
        ma, fa, m_k, m_pi = 1.7, 120.0, 0.495, 0.135
        mhat, ms = 0.008, 0.24
        b0mq = np.diag([mhat, mhat, ms])
        kernels = _kernels_from_inputs(
            ma,
            fa,
            ("K+", "pi0", "K-"),
            (m_k, m_pi, m_k),
            (k_plus, pi_zero, k_minus),
            self.c_a,
            self.c_p,
            b0mq,
        )
        invariant = 0.91
        delta_a = self.c_a[0, 0] - self.c_a[2, 2]
        expected = (
            delta_a * (3.0 * invariant - ma**2 - 3.0 * m_k**2)
            - 2.0
            * (2.0 * mhat * self.c_p[0, 0] + ms * self.c_p[2, 2])
        ) / (6.0 * fa * F_PI)
        np.testing.assert_allclose(kernels[0](invariant), expected, rtol=2e-13)
        np.testing.assert_allclose(kernels[1](invariant), expected, rtol=2e-13)

    def test_all_charge_channels_are_finite_and_distinct(self) -> None:
        values = []
        for channel in supported_channels():
            model = build_s_wave_model(1.6, self.c_a, self.c_p, channel, 100.0)
            values.append(model.contact_amplitude(0.8, 0.9))
            self.assertTrue(np.isfinite(values[-1]))
        self.assertGreater(len({complex(round(v.real, 14), round(v.imag, 14)) for v in values}), 1)

    def test_physical_thresholds_are_split(self) -> None:
        thresholds = kpi_thresholds()
        self.assertAlmostEqual(thresholds[("K+", "pi0")], 0.6286538)
        self.assertAlmostEqual(thresholds[("K0", "pi+")], 0.63718139)
        self.assertNotEqual(thresholds[("K+", "pi0")], thresholds[("K0", "pi+")])

    def test_threshold_and_validation(self) -> None:
        self.assertEqual(
            s_wave_decay_width(1.0, self.c_a, self.c_p, ("K+", "K-", "pi0"), 100.0),
            0.0,
        )
        with self.assertRaisesRegex(ValueError, "diagonal"):
            off_diagonal = self.c_a.copy()
            off_diagonal[0, 1] = off_diagonal[1, 0] = 0.1
            build_s_wave_model(1.5, off_diagonal, self.c_p, ("K+", "K-", "pi0"), 100.0)
        with self.assertRaisesRegex(ValueError, "unsupported"):
            build_s_wave_model(1.5, self.c_a, self.c_p, ("K+", "pi-", "pi0"), 100.0)

    def test_width_scaling_and_nontrivial_fsi(self) -> None:
        channel = ("K+", "K-", "pi0")
        base = s_wave_decay_width(1.5, self.c_a, self.c_p, channel, 100.0)
        doubled = s_wave_decay_width(1.5, 2 * self.c_a, 2 * self.c_p, channel, 100.0)
        dressed = s_wave_decay_width(
            1.5, self.c_a, self.c_p, channel, 100.0, lambda s, pair: 2.0
        )
        self.assertGreater(base, 0.0)
        self.assertAlmostEqual(doubled / base, 4.0, places=9)
        self.assertAlmostEqual(dressed / base, 4.0, places=9)

    def test_fast_width_matches_direct_contact_integration(self) -> None:
        channel = ("K+", "K-", "pi0")
        direct = s_wave_decay_width(1.5, self.c_a, self.c_p, channel, 100.0)
        fast = s_wave_decay_width_fast(
            1.5,
            self.c_a,
            self.c_p,
            channel,
            100.0,
            relative_tolerance=1e-8,
        )
        np.testing.assert_allclose(fast, direct, rtol=1e-8, atol=0.0)

    def test_fast_width_with_tabulated_omnes(self) -> None:
        threshold = 0.4
        endpoint = 1.2
        phase = lambda s: 0.4 * (s - threshold) / (endpoint - threshold)
        direct_omega = OmnesFunction(
            phase,
            threshold,
            endpoint,
            asymptotic_phase=phase(endpoint),
        )
        tabulated = TabulatedOmnes(
            direct_omega,
            minimum=(0.493677 + 0.1349768) ** 2,
            maximum=(1.2 - 0.493677) ** 2,
            relative_tolerance=1e-6,
        )
        fsi = PairOmnesFSI(
            {("K+", "pi0"): tabulated, ("K-", "pi0"): tabulated}
        )
        channel = ("K+", "K-", "pi0")
        direct = s_wave_decay_width(1.2, self.c_a, self.c_p, channel, 100.0, fsi)
        fast = s_wave_decay_width_fast(
            1.2,
            self.c_a,
            self.c_p,
            channel,
            100.0,
            fsi,
            relative_tolerance=1e-6,
        )
        np.testing.assert_allclose(fast, direct, rtol=1e-6, atol=0.0)

    def test_physical_neutral_kaon_projection_for_charged_mode(self) -> None:
        flavor = s_wave_decay_width(
            1.6, self.c_a, self.c_p, ("K0", "pi+", "K-"), 100.0
        )
        short = s_wave_decay_width(
            1.6, self.c_a, self.c_p, ("K-", "pi+", "k_s"), 100.0
        )
        long = s_wave_decay_width(
            1.6, self.c_a, self.c_p, ("k_l", "K-", "pi+"), 100.0
        )
        normalization = np.sqrt(2.0 * (1.0 + abs(EPSILON_K) ** 2))
        p_mixing = (1.0 + EPSILON_K) / normalization
        np.testing.assert_allclose(
            short, abs(p_mixing) ** 2 * flavor, rtol=2e-12, atol=0.0
        )
        np.testing.assert_allclose(
            long, abs(p_mixing) ** 2 * flavor, rtol=2e-12, atol=0.0
        )

    def test_two_physical_neutral_kaons_keep_multiplicity(self) -> None:
        flavor = s_wave_decay_width(
            1.6, self.c_a, self.c_p, ("K0", "pi0", "K0bar"), 100.0
        )
        short_short = s_wave_decay_width(
            1.6, self.c_a, self.c_p, ("KS", "KS", "pi0"), 100.0
        )
        long_long = s_wave_decay_width(
            1.6, self.c_a, self.c_p, ("pi0", "KL", "KL"), 100.0
        )
        short_long = s_wave_decay_width(
            1.6, self.c_a, self.c_p, ("KL", "pi0", "KS"), 100.0
        )
        normalization = np.sqrt(2.0 * (1.0 + abs(EPSILON_K) ** 2))
        p_mixing = (1.0 + EPSILON_K) / normalization
        q_mixing = (1.0 - EPSILON_K) / normalization
        expected = 2.0 * abs(p_mixing * q_mixing) ** 2 * flavor
        np.testing.assert_allclose(short_short, expected, rtol=2e-12, atol=0.0)
        np.testing.assert_allclose(long_long, expected, rtol=2e-12, atol=0.0)
        np.testing.assert_allclose(short_long, 0.0, atol=1e-38)


if __name__ == "__main__":
    unittest.main()
