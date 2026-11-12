"""Tests for the frame (in)dependence of the LO bare-ALP-mass bound
(derive_alp_lo_stability_bound_invariance.py). Runtime about two minutes.

Run: .venv/bin/python -m unittest discover -s scripts -p test_derive_alp_lo_stability_bound_invariance.py
"""

import unittest

import numpy as np
import sympy as sp

from derive_alp_lo_stability_bound_invariance import (
    RATIONAL_PARAMETERS,
    dressed_mass_real_trace_linear_convention,
    dressed_mass_second_order_linear_convention,
    eigenvalue_sign_checks,
    general_frame_reduction,
    invariant_form_checks,
    nonlinear_vacuum_energy_checks,
    previous_bracket_in_frame,
    rational_frame_checks,
    rational_frames,
    schur_complement_shift_invariance,
    singular_value_formula_numeric_check,
    symbolic_redefinition,
)


def assert_all(test, checks):
    for label, passed in checks.items():
        with test.subTest(label=label):
            test.assertTrue(passed, label)


class CongruenceAndFrameStructureTest(unittest.TestCase):
    def test_schur_complement_invariant_under_shift(self):
        self.assertTrue(schur_complement_shift_invariance())

    def test_mass_matrix_in_any_frame(self):
        assert_all(self, general_frame_reduction())


class SymbolicRedefinitionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = symbolic_redefinition()

    def test_all_checks(self):
        assert_all(self, self.result["checks"])

    def test_previous_bracket_changes_by_nonzero_amount(self):
        self.assertNotEqual(sp.expand(self.result["bracket difference"]), 0)


class DressedSourceTest(unittest.TestCase):
    def test_second_order_dressed_mass_closed_form(self):
        self.assertTrue(dressed_mass_second_order_linear_convention())

    def test_real_trace_of_second_order_dressed_mass(self):
        assert_all(self, dressed_mass_real_trace_linear_convention())

    def test_invariant_form(self):
        assert_all(self, invariant_form_checks())

    def test_singular_value_formula_against_mpmath(self):
        formula, finite_difference = singular_value_formula_numeric_check()
        self.assertAlmostEqual(formula, finite_difference, delta=1e-9 * abs(formula))


class RationalFramesTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = rational_frame_checks()

    def test_all_checks(self):
        assert_all(self, self.result["checks"])

    def test_schur_bound_identical_in_all_frames(self):
        values = {row["Schur bound / eps^2"] for row in self.result["rows"].values()}
        self.assertEqual(values, {self.result["invariant bound"]})

    def test_previous_bracket_is_frame_dependent(self):
        frames = rational_frames()
        values = {sp.nsimplify(previous_bracket_in_frame(frame, RATIONAL_PARAMETERS["anomaly_mass"]))
                  for frame in frames.values()}
        self.assertEqual(len(values), len(frames))


class NumericalTest(unittest.TestCase):
    def test_eigenvalue_signs_at_invariant_threshold(self):
        for label, row in eigenvalue_sign_checks().items():
            with self.subTest(frame=label):
                self.assertGreater(row["above"], 0)
                self.assertLess(row["below"], 0)
                if not np.isclose(row["previous-bracket threshold"], row["threshold"], rtol=1e-12, atol=0):
                    self.assertEqual(np.sign(row["at previous +"]), np.sign(row["at previous -"]))

    def test_vacuum_energy_frame_independent_and_curvature(self):
        result = nonlinear_vacuum_energy_checks()
        for ratio, linear, rotated in result["finite x"]:
            self.assertAlmostEqual(linear, rotated, delta=1e-12)
        for value in result["second derivative"].values():
            self.assertAlmostEqual(value, result["expected second derivative"],
                                   delta=1e-5 * abs(result["expected second derivative"]))


if __name__ == "__main__":
    unittest.main()
