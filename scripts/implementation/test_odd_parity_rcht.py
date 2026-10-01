"""Tests for the non-HLS odd-parity RChT numerical kernel.

Invocation:
    .venv/bin/python -m unittest scripts.implementation.test_odd_parity_rcht
"""

import unittest
from dataclasses import FrozenInstanceError, replace

import numpy as np

from scripts.implementation.odd_parity_rcht import (
    ALPFlavorSpurions,
    OPEParameters,
    OPEPrior,
    OPEResult,
    OddParityCouplings,
    VJPCouplings,
    VVPCouplings,
    anomaly_seed,
    fixed_ope_prior,
    fixed_ope_prior_scales,
    ideal_vector_generators,
    odd_parity_flavor_kernels,
    ope_linear_system,
    ope_prior_chi2,
    ope_residuals,
    photon_anomaly_seed,
    photon_vector_transform,
    transform_anomalous_bilinear,
    vjp_source_coefficients,
    vvp_source_coefficients,
    weighted_vjp_coefficients,
    weighted_vvp_coefficients,
)


class FlavorKernelTests(unittest.TestCase):
    def setUp(self):
        self.m_a_inverse_C = np.diag([0.2, -0.1, 0.3])
        self.cap = np.array([[0.4, 0.1j, 0.0], [-0.1j, -0.2, 0.0], [0.0, 0.0, 0.1]]) # Why a charge violating coupling matrix????
        self.mass = np.diag([0.01, 0.02, 0.4])
        self.spurions = ALPFlavorSpurions(
            m_a_inverse_C=self.m_a_inverse_C,
            c_q_AP=self.cap,
            b0_mq=self.mass,
        )
        self.ta = np.array([[0.0, 1.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
        self.tb = np.array([[0.3, 0.0, 0.0], [0.2j, -0.1, 0.0], [0.0, 0.0, 0.5]])

    def test_direct_formulas_and_nonhermitian_sources(self):
        np.testing.assert_array_equal(
            self.spurions.m_a_inverse_C, self.m_a_inverse_C
        )
        result = odd_parity_flavor_kernels(self.ta, self.tb, self.spurions)
        anticommutator = self.ta @ self.tb + self.tb @ self.ta
        c_q_AP_plus_m_a_inverse_C = self.cap + self.m_a_inverse_C
        np.testing.assert_allclose(
            result.derivative_trace,
            np.trace(anticommutator @ self.m_a_inverse_C),
        )
        np.testing.assert_allclose(
            result.mass_trace,
            np.trace(
                anticommutator
                @ (
                    c_q_AP_plus_m_a_inverse_C @ self.mass
                    + self.mass @ c_q_AP_plus_m_a_inverse_C
                )
            ),
        )

        vjp = vjp_source_coefficients(result, 0.3, 10.0).coefficients
        np.testing.assert_allclose(vjp, [0.12 * result.derivative_trace, 0.12 * result.derivative_trace,
                                        -0.24 * result.mass_trace, 0.0, 0.12 * result.derivative_trace,
                                        0.12 * result.derivative_trace, 0.12 * result.derivative_trace])
        vvp = vvp_source_coefficients(result, 10.0).coefficients
        np.testing.assert_allclose(vvp, [-0.2 * result.derivative_trace, 0.4 * result.mass_trace,
                                        -0.2 * result.derivative_trace, -0.2 * result.derivative_trace])

    def test_unitary_flavor_covariance(self):
        rng = np.random.default_rng(9206)
        q, _ = np.linalg.qr(rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)))
        rotate = lambda matrix: q @ matrix @ q.conj().T
        original = odd_parity_flavor_kernels(self.ta, self.tb, self.spurions)
        transformed = odd_parity_flavor_kernels(
            rotate(self.ta), rotate(self.tb),
            ALPFlavorSpurions(
                rotate(self.m_a_inverse_C), rotate(self.cap), rotate(self.mass)
            ),
        )
        np.testing.assert_allclose(
            [transformed.derivative_trace, transformed.mass_trace],
            [original.derivative_trace, original.mass_trace], atol=2e-16,
        )

    def test_d2_survives_without_m_a_inverse_C_and_scaling(self):
        spurions = ALPFlavorSpurions(np.zeros((3, 3)), np.eye(3), self.mass)
        kernels = odd_parity_flavor_kernels(np.eye(3), np.eye(3), spurions)
        self.assertEqual(kernels.derivative_trace, 0.0)
        self.assertNotEqual(kernels.mass_trace, 0.0)
        base = vvp_source_coefficients(kernels, 7.0).coefficients
        np.testing.assert_allclose(vvp_source_coefficients(kernels, 14.0).coefficients, base / 2)
        weighted = weighted_vvp_coefficients(kernels, VVPCouplings(1, 2, 3, 4), 7.0).coefficients
        np.testing.assert_allclose(weighted, base * [1, 2, 3, 4])

    def test_weighted_vjp_and_immutability(self):
        kernels = odd_parity_flavor_kernels(self.ta, self.tb, self.spurions)
        couplings = VJPCouplings(1, 2, 3, 4, 5, 6, 7)
        source = vjp_source_coefficients(kernels, 0.2, 5.0).coefficients
        weighted = weighted_vjp_coefficients(kernels, couplings, 0.2, 5.0, 0.8).coefficients
        np.testing.assert_allclose(weighted, source * np.arange(1, 8) / 0.8)
        with self.assertRaises(ValueError):
            weighted[0] = 0.0
        with self.assertRaises(FrozenInstanceError):
            couplings.c1 = 2.0

    def test_flavor_validation(self):
        with self.assertRaisesRegex(ValueError, "Hermitian"):
            ALPFlavorSpurions(self.ta, self.cap, self.mass)
        with self.assertRaisesRegex(ValueError, "shape"):
            odd_parity_flavor_kernels(np.eye(2), self.tb, self.spurions)
        bad = self.ta.copy()
        bad[0, 0] = np.nan
        with self.assertRaisesRegex(ValueError, "finite"):
            odd_parity_flavor_kernels(bad, self.tb, self.spurions)
        with self.assertRaisesRegex(ValueError, "positive"):
            vvp_source_coefficients(odd_parity_flavor_kernels(self.ta, self.tb, self.spurions), 0.0)
        with self.assertRaisesRegex(ValueError, "finite and real"):
            VJPCouplings(1j, 2, 3, 4, 5, 6, 7)


class AnomalyTests(unittest.TestCase):
    def test_identical_distinct_photon_and_reparametrization(self):
        q = np.diag([2 / 3, -1 / 3, -1 / 3])
        ca = np.diag([0.2, -0.1, 0.4])
        m_a_inverse_C = np.diag([0.05, 0.02, -0.03])
        trace = np.trace((q @ q + q @ q) @ (ca + m_a_inverse_C))
        same = anomaly_seed(q, q, ca, m_a_inverse_C, 1.3, 3, True)
        distinct = anomaly_seed(q, q, ca, m_a_inverse_C, 1.3, 3, False)
        np.testing.assert_allclose(same, 1.3 + 3 * trace)
        np.testing.assert_allclose(distinct, 1.3 + 6 * trace)
        np.testing.assert_allclose(
            same, photon_anomaly_seed(q, ca, m_a_inverse_C, 1.3, 3)
        )

        kappa = np.diag([0.03, -0.02, 0.01])
        shifted_c = 1.3 - 3 * np.trace((q @ q + q @ q) @ kappa)
        np.testing.assert_allclose(
            same,
            anomaly_seed(q, q, ca + kappa, m_a_inverse_C, shifted_c, 3, True),
        )

    def test_anomaly_validation(self):
        eye = np.eye(3)
        with self.assertRaisesRegex(TypeError, "explicit boolean"):
            anomaly_seed(eye, eye, eye, eye, 0, 3, 1)
        with self.assertRaisesRegex(ValueError, "positive integer"):
            anomaly_seed(eye, eye, eye, eye, 0, 3.0, True)
        with self.assertRaisesRegex(ValueError, "Hermitian"):
            anomaly_seed(eye, eye, np.triu(np.ones((3, 3))), eye, 0, 3, True)


class OPETests(unittest.TestCase):
    def setUp(self):
        self.parameters = OPEParameters(3, 0.77, 0.15, 0.092)

    def test_exact_zero_point_and_null_directions(self):
        a, b = ope_linear_system(self.parameters)
        point = np.zeros(11)
        point[4] = b[2]
        point[5] = 0.0
        point[1] = point[0] + point[4]
        point[2] = -point[0] / 4
        point[7] = b[3]
        point[9] = b[4]
        couplings = OddParityCouplings(VJPCouplings(*point[:7]), VVPCouplings(*point[7:]))
        np.testing.assert_allclose(ope_residuals(couplings, self.parameters).residuals, 0.0, atol=2e-17)
        null_directions = (
            np.array([1, 1, -0.25, 0, 0, 0, 0, 0, 0, 0, 0]),
            np.array([0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0]),
            np.array([0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0]),
            np.array([0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0]),
            np.array([0, 0, 0, 0, 0, 0, 0, 8, -1, 0, 0]),
            np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]),
        )
        self.assertEqual(np.linalg.matrix_rank(np.stack(null_directions)), 6)
        for direction in null_directions:
            np.testing.assert_allclose(a @ (point + 1.7 * direction), b, atol=2e-17)

    def test_fixed_prior_numerical_reference_and_scaling(self):
        parameters = OPEParameters(3, 0.770, np.sqrt(2.0) * 0.092, 0.092)
        expected_scales = np.array([
            0.01987533408744907,
            0.01987533408744907,
            0.01987533408744907,
            0.29134790486234543,
            0.22884790486234544,
        ])
        expected_sigmas = np.array([
            0.006625111362483023,
            0.006625111362483023,
            0.006625111362483023,
            0.09711596828744847,
            0.07628263495411514,
        ])
        np.testing.assert_allclose(fixed_ope_prior_scales(parameters), expected_scales, rtol=2e-15)
        baseline = fixed_ope_prior(parameters)
        np.testing.assert_allclose(np.sqrt(np.diag(baseline.covariance)), expected_sigmas, rtol=2e-15)
        tight = fixed_ope_prior(parameters, truncation_fraction=0.2)
        broad = fixed_ope_prior(parameters, truncation_fraction=0.5)
        np.testing.assert_allclose(np.sqrt(np.diag(tight.covariance)), 0.2 * expected_scales)
        np.testing.assert_allclose(np.sqrt(np.diag(broad.covariance)), 0.5 * expected_scales)
        self.assertTrue(np.all(np.linalg.eigvalsh(fixed_ope_prior(parameters, rho=0.5).covariance) > 0.0))

    def test_fixed_prior_correlation_spd_interval(self):
        for rho in (-0.25, 1.0, -0.3, 1.2):
            with self.assertRaisesRegex(ValueError, "-1/4 < rho < 1"):
                fixed_ope_prior(self.parameters, rho=rho)
        with self.assertRaisesRegex(ValueError, "positive"):
            fixed_ope_prior(self.parameters, truncation_fraction=0.0)

    def test_prior_manual_chi2_and_rejections(self):
        covariance = np.diag([1.0, 2.0, 3.0, 4.0, 5.0])
        residuals = np.arange(1.0, 6.0)
        result = OPEResult(tuple("abcde"), residuals)
        self.assertAlmostEqual(ope_prior_chi2(result, OPEPrior(covariance)), float(np.sum(residuals**2 / np.diag(covariance))))
        with self.assertRaisesRegex(ValueError, "symmetric"):
            OPEPrior(np.triu(np.ones((5, 5))))
        with self.assertRaisesRegex(ValueError, "positive definite"):
            OPEPrior(np.diag([1, 1, 1, 1, 0]))
        with self.assertRaisesRegex(ValueError, "positive integer"):
            OPEParameters(True, 1, 1, 1)
        with self.assertRaisesRegex(ValueError, "positive"):
            OPEParameters(3, -1, 1, 1)


class PhotonVectorTests(unittest.TestCase):
    def test_congruence_and_quadratic_form_invariance(self):
        rng = np.random.default_rng(451)
        rotation, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        transform = photon_vector_transform([0.1, -0.03, 0.2], rotation)
        raw = rng.normal(size=(4, 4)); lam = raw + raw.T
        physical = transform_anomalous_bilinear(lam, transform)
        np.testing.assert_allclose(physical, transform.s.T @ lam @ transform.s)
        x = rng.normal(size=4)
        np.testing.assert_allclose(x @ physical @ x, (transform.s @ x) @ lam @ (transform.s @ x))
        np.testing.assert_allclose(transform.k_sqrt @ transform.k_sqrt,
                                   np.eye(3) - np.outer(transform.epsilon, transform.epsilon) / transform.z)

    def test_zero_mixing_and_block_behavior(self):
        rotation = np.array([[0.0, -1.0], [1.0, 0.0]])
        transform = photon_vector_transform(np.zeros(2), rotation)
        np.testing.assert_allclose(transform.s, np.block([[np.ones((1, 1)), np.zeros((1, 2))],
                                                          [np.zeros((2, 1)), rotation]]))
        lam = np.diag([2.0, 3.0, 4.0])
        result = transform_anomalous_bilinear(lam, transform)
        self.assertEqual(result[0, 0], 2.0)
        np.testing.assert_allclose(result[1:, 1:], rotation.T @ lam[1:, 1:] @ rotation)

    def test_transform_validation(self):
        with self.assertRaisesRegex(ValueError, "real"):
            photon_vector_transform([1j], [[1]])
        with self.assertRaisesRegex(ValueError, "orthogonal"):
            photon_vector_transform([0.1, 0.2], np.ones((2, 2)))
        transform = photon_vector_transform([0.1], [[1]])
        with self.assertRaisesRegex(ValueError, "must match"):
            replace(transform, z=transform.z + 1.0)
        with self.assertRaisesRegex(ValueError, "symmetric"):
            transform_anomalous_bilinear([[1, 2], [0, 1]], transform)
        with self.assertRaisesRegex(ValueError, "shape"):
            transform_anomalous_bilinear(np.eye(3), transform)


class IdealGeneratorTests(unittest.TestCase):
    def test_normalization_conjugation_and_fresh_copies(self):
        first, second = ideal_vector_generators(), ideal_vector_generators()
        expected = {"rho+", "rho-", "rho0", "K*+", "K*-", "K*0", "K*0bar", "omega", "phi"}
        self.assertEqual(set(first), expected)
        for generator in first.values():
            np.testing.assert_allclose(np.trace(generator @ generator.conj().T), 1.0)
        for particle, antiparticle in (("rho+", "rho-"), ("K*+", "K*-"), ("K*0", "K*0bar")):
            np.testing.assert_allclose(first[particle].conj().T, first[antiparticle])
        first["rho+"][0, 1] = 7.0
        self.assertEqual(second["rho+"][0, 1], 1.0)


if __name__ == "__main__":
    unittest.main()
