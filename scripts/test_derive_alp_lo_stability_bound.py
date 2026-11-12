"""Tests for the analytical LO bound m_a0^2 >= (f_pi/f_a)^2 B (derive_alp_lo_stability_bound.py)."""

import unittest

import numpy as np
import sympy as sp

from derive_alp_lo_stability_bound import (
    axial_source_independence,
    bound_coefficient,
    bound_coefficient_without_topological_term,
    charged_entries_check,
    completed_square_remark,
    eigenvalue_sign_checks,
    flavour_diagonal_derivation,
    flavour_violating_derivation,
    illustration,
    illustrative_values_check,
    light_eigenvalue_check,
    lo_inputs,
    smallest_mass_eigenvalue,
    sufficient_bound_and_limits,
    topological_part_is_rank_one_square,
)


class StructureTest(unittest.TestCase):
    def test_mass_matrix_independent_of_axial_source(self):
        self.assertTrue(axial_source_independence())

    def test_topological_part_is_rank_one_square(self):
        self.assertTrue(topological_part_is_rank_one_square())


class ClosedFormTest(unittest.TestCase):
    def test_flavour_diagonal_derivation(self):
        for label, passed in flavour_diagonal_derivation()["checks"].items():
            self.assertTrue(passed, label)

    def test_flavour_violating_entries(self):
        for label, passed in flavour_violating_derivation()["checks"].items():
            self.assertTrue(passed, label)

    def test_charge_violating_entries_all_generators(self):
        self.assertTrue(charged_entries_check(seed=11, sample_count=2))

    def test_sufficient_bound_and_limits(self):
        for label, passed in sufficient_bound_and_limits()["checks"].items():
            self.assertTrue(passed, label)

    def test_reproduces_previous_illustrative_values(self):
        for label, passed in illustrative_values_check().items():
            self.assertTrue(passed, label)


class EigenvalueTest(unittest.TestCase):
    def test_sign_change_at_bound(self):
        for threshold, above, below in eigenvalue_sign_checks(seed=3, sample_count=6):
            self.assertGreater(above, 0, threshold)
            self.assertLess(below, 0, threshold)

    def test_sign_change_at_vanishing_mixing_point(self):
        """c_q = m_0^2 c_G/(N_f chi_q): b = 0, bound = -M^2_aa^top < 0, and m_a0^2 = 0 is stable."""
        chi = (sp.Rational(3, 200), sp.Rational(1, 40), sp.Rational(47, 100))
        singlet, gluon = sp.Rational(4, 5), sp.Rational(1, 2)
        pseudoscalar = sp.diag(*[singlet * gluon / (3 * value) for value in chi])
        inputs = lo_inputs(sp.zeros(3), pseudoscalar, chi, gluon, sp.Rational(23, 250), sp.Rational(1, 5), singlet, 0)
        epsilon = inputs["pion_decay_constant"] / inputs["alp_decay_constant"]
        threshold = epsilon**2 * bound_coefficient(pseudoscalar, chi, singlet, gluon)
        self.assertEqual(sp.simplify(threshold + 2 * epsilon**2 * singlet * gluon**2 / 3), 0)
        self.assertGreater(smallest_mass_eigenvalue(inputs), 0)
        for offset, sign in ((sp.Rational(1, 10**6), 1), (-sp.Rational(1, 10**6), -1)):
            value = smallest_mass_eigenvalue(dict(inputs, bare_alp_mass_squared=threshold + offset))
            self.assertEqual(np.sign(value), sign)

    def test_sufficient_bound_is_stable(self):
        """m_a0^2 = epsilon^2 B_no_top (1 + 1e-9) gives M^2 > 0 at random points."""
        generator = np.random.default_rng(17)
        for _ in range(5):
            chi = tuple(sp.Rational(int(generator.integers(1, 60)), 100) for _ in range(3))
            pseudoscalar = sp.diag(*[sp.Rational(int(generator.integers(-9, 10)), 4) for _ in range(3)])
            inputs = lo_inputs(sp.zeros(3), pseudoscalar, chi, sp.Rational(int(generator.integers(-6, 7)), 4),
                               sp.Rational(23, 250), sp.Rational(1, 5), sp.Rational(4, 5), 0)
            epsilon = inputs["pion_decay_constant"] / inputs["alp_decay_constant"]
            value = epsilon**2 * bound_coefficient_without_topological_term(pseudoscalar, chi)
            shifted = value * (1 + sp.Rational(1, 10**9)) + sp.Rational(1, 10**12)
            self.assertGreater(smallest_mass_eigenvalue(dict(inputs, bare_alp_mass_squared=shifted)), 0)

    def test_light_eigenvalue_leading_order(self):
        ratios = light_eigenvalue_check()
        self.assertLess(abs(ratios[-1] - 1), 1e-4)
        self.assertLess(abs(ratios[-1] - 1), abs(ratios[0] - 1))


class RemarkAndIllustrationTest(unittest.TestCase):
    def test_completed_square_remark(self):
        self.assertTrue(completed_square_remark())

    def test_illustration_scales_as_inverse_decay_constant_squared(self):
        _, singlet, rows = illustration()
        self.assertGreater(singlet, 0)
        for row in rows:
            self.assertLessEqual(row["B"], row["B_no_top"] + 1e-15)
            (_, first), (_, second), (_, third) = row["bounds"]
            self.assertNotEqual(row["B"], 0)
            self.assertAlmostEqual(first / second, 100.0, places=9)
            self.assertAlmostEqual(second / third, 100.0, places=9)


if __name__ == "__main__":
    unittest.main()
