"""Tests of the common-width Gaussian diagnostic."""
import unittest

import numpy as np
from scipy.optimize import minimize_scalar

from scripts.diagnose_joint_kstar_width import common_width_fit


class CommonWidthTest(unittest.TestCase):
    def test_correlated_measurements_match_independent_minimization(self):
        observed = np.array([0.0514, 0.0473])
        covariance = np.array([[0.0008**2, -0.2 * 0.0008 * 0.0005],
                               [-0.2 * 0.0008 * 0.0005, 0.0005**2]])
        fitted, minimum = common_width_fit(observed, covariance)
        def chi2(width):
            residual = observed - width
            return residual @ np.linalg.solve(covariance, residual)
        numerical = minimize_scalar(chi2, bracket=(0.04, 0.05, 0.06), tol=1e-13)
        self.assertAlmostEqual(fitted, numerical.x, places=9)
        self.assertAlmostEqual(minimum, numerical.fun, places=11)

    def test_equal_measurements_have_zero_minimum(self):
        fitted, minimum = common_width_fit([0.05, 0.05], np.diag([0.001**2, 0.002**2]))
        self.assertAlmostEqual(fitted, 0.05)
        self.assertEqual(minimum, 0.0)


if __name__ == '__main__':
    unittest.main()
