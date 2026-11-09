"""Tests for the correlated phi-to-kaon residual transformation."""
import unittest
import numpy as np

from scripts.diagnose_joint_phi_kaon_residuals import channel_residuals, diagnostic


class PhiKaonResidualTests(unittest.TestCase):
    def test_correlated_residuals_preserve_norm_and_covariance(self):
        covariance = np.array([[4., -1.2], [-1.2, 1.]])
        differences = np.array([1.5, -0.2])
        pulls = np.linalg.solve(np.linalg.cholesky(covariance), differences)
        restored, rotated, rotation = channel_residuals(pulls, covariance)
        np.testing.assert_allclose(restored, differences, atol=1e-14)
        np.testing.assert_allclose(rotation.T @ rotation, np.eye(2), atol=1e-14)
        np.testing.assert_allclose(rotated @ rotated, differences @ np.linalg.solve(covariance, differences), atol=1e-14)
        np.testing.assert_allclose(np.sign(rotated[0]), 1.)

    def test_channel_swap_changes_only_channel_order(self):
        covariance = np.array([[4., -1.2], [-1.2, 1.]])
        differences = np.array([1.5, -0.2])
        pulls = np.linalg.solve(np.linalg.cholesky(covariance), differences)
        _, rotated, _ = channel_residuals(pulls, covariance)
        swapped_covariance = covariance[::-1, ::-1]
        swapped_pulls = np.linalg.solve(np.linalg.cholesky(swapped_covariance), differences[::-1])
        _, swapped_rotated, _ = channel_residuals(swapped_pulls, swapped_covariance)
        np.testing.assert_allclose(swapped_rotated, rotated[::-1], atol=1e-14)

    def test_recorded_endpoint_reproduces_block_chi2(self):
        result = diagnostic()
        self.assertAlmostEqual(result['rotated_chi2'], result['stored_block_chi2'], places=10)


if __name__ == '__main__':
    unittest.main()
