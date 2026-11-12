"""Tests for scripts/derive_alp_lo_arrow_matrix_poles.py (arrow-matrix spectrum and isolated-pole approximations).

Run: .venv/bin/python -m unittest discover -s scripts -p test_derive_alp_lo_arrow_matrix_poles.py
"""

import unittest

import mpmath

from derive_alp_lo_arrow_matrix_poles import (
    alp_fixed_point_map,
    alp_pole_iteration,
    arrow_determinant,
    arrow_matrix,
    derive,
    interlacing_holds,
    isolated_meson_poles,
    sorted_eigenvalues,
)


class ArrowMatrixPoleTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mpmath.mp.dps = 60
        cls.result = derive()

    def test_all_checks_pass(self):
        for label, passed in self.result["checks"].items():
            self.assertTrue(passed, label)

    def test_two_level_form_is_exact_with_one_coupled_meson(self):
        # only meson 0 has tr(lambda_n C) != 0: Mt_n(s) = M and the frozen roots are the exact poles
        # (tr_lambda_C: list of tr(lambda_n C) over the mesons, GeV^2)
        unmixed = [mpmath.mpf(1), mpmath.mpf(3), mpmath.mpf(5)]
        tr_lambda_C = [mpmath.mpf(2), mpmath.mpf(0), mpmath.mpf(0)]
        epsilon, alp_entry = mpmath.mpf("1.5"), mpmath.mpf("1.2")
        lower, upper, frozen, _ = isolated_meson_poles(0, unmixed, tr_lambda_C, epsilon, alp_entry)
        poles = sorted_eigenvalues(arrow_matrix(unmixed, tr_lambda_C, epsilon, alp_entry))
        self.assertEqual(frozen, alp_entry)
        self.assertLess(abs(lower - poles[0]), mpmath.mpf("1e-50"))
        self.assertTrue(any(abs(upper - pole) < mpmath.mpf("1e-50") for pole in poles))

    def test_determinant_formula_is_product_over_poles(self):
        unmixed = [mpmath.mpf(value) for value in ("0.2", "0.9", "1.7")]
        tr_lambda_C = [mpmath.mpf(value) for value in ("0.4", "-0.8", "0.3")]
        epsilon, alp_entry = mpmath.mpf("0.7"), mpmath.mpf("1.1")
        poles = sorted_eigenvalues(arrow_matrix(unmixed, tr_lambda_C, epsilon, alp_entry))
        point = mpmath.mpf("0.37")
        self.assertLess(abs(arrow_determinant(point, unmixed, tr_lambda_C, epsilon, alp_entry)
                            - mpmath.fprod(point - pole for pole in poles)), mpmath.mpf("1e-50"))

    def test_interlacing_rejects_wrong_order(self):
        self.assertTrue(interlacing_holds([1, 2], [0, 1.5, 3]))
        self.assertFalse(interlacing_holds([1, 2], [1.5, 1.8, 3]))

    def test_iteration_fixed_point_is_alp_pole(self):
        unmixed = [mpmath.mpf(value) for value in (1, 2)]
        tr_lambda_C = [mpmath.mpf(value) for value in ("0.5", "0.5")]
        epsilon, alp_entry = mpmath.mpf(1), mpmath.mpf(6)
        iterate, _, converged = alp_pole_iteration(unmixed, tr_lambda_C, epsilon, alp_entry)
        self.assertTrue(converged)
        self.assertLess(abs(alp_fixed_point_map(iterate, unmixed, tr_lambda_C, epsilon, alp_entry) - iterate),
                        mpmath.mpf("1e-40"))
        self.assertLess(abs(iterate - sorted_eigenvalues(arrow_matrix(unmixed, tr_lambda_C, epsilon, alp_entry))[-1]),
                        mpmath.mpf("1e-40"))


if __name__ == "__main__":
    unittest.main()
