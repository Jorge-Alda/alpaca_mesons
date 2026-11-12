"""Tests for the LO finite ALP-meson mixing written with the operator M_a of Alda:2026nrd.

Run: .venv/bin/python -m unittest discover -s scripts -p test_derive_alp_finite_mixing_operator_inverse.py
"""

import unittest

import sympy as sp

from derive_alp_finite_mixing_operator_inverse import (
    SEEDS,
    closed_form_inverse_of_invariant_matrix,
    invariant_matrix,
    lo_pole_operator,
    partially_symbolic_checks,
    rational_checks,
    rational_inputs,
    step_M2_aa_can,
    step_canonical_source,
    step_component_completeness,
    step_dependence_on_couplings,
    step_explicit_inverse,
    step_flavor_diagonal_quartic,
    step_one_alp_amplitude,
    step_pole_equation_matrix_form,
    step_alp_inverse_propagator,
    step_weak_mixing_limit,
    symbolic_inputs,
)


class FiniteMixingOperatorInverseTest(unittest.TestCase):
    def assert_checks(self, checks, context=""):
        for label, passed in checks.items():
            if not label.startswith("_"):
                self.assertTrue(bool(passed), f"{context}: {label}")

    @classmethod
    def setUpClass(cls):
        cls.symbolic = symbolic_inputs()

    def test_completeness(self):
        self.assert_checks(step_component_completeness())

    def test_step1_matrix_form_of_pole_equation(self):
        self.assert_checks(step_pole_equation_matrix_form(self.symbolic))

    def test_step2_canonical_source(self):
        self.assert_checks(step_canonical_source(self.symbolic))

    def test_step3_explicit_inverse(self):
        self.assert_checks(step_explicit_inverse(self.symbolic))

    def test_step5_M2_aa_can(self):
        self.assert_checks(step_M2_aa_can(self.symbolic))

    def test_step6_alp_inverse_propagator_symbolic(self):
        self.assert_checks(step_alp_inverse_propagator(self.symbolic))

    def test_steps_4_and_7_partially_symbolic(self):
        for step, checks in partially_symbolic_checks().items():
            self.assert_checks(checks, step)

    def test_coupling_dependence_and_shift(self):
        self.assert_checks(step_dependence_on_couplings(self.symbolic))

    def test_weak_mixing_limit(self):
        self.assert_checks(step_weak_mixing_limit(self.symbolic))

    def test_flavor_diagonal_quartic(self):
        self.assert_checks(step_flavor_diagonal_quartic())

    def test_one_alp_amplitude_identity(self):
        self.assert_checks(step_one_alp_amplitude())

    def test_rational_seeds(self):
        for seed in SEEDS:
            for step, checks in rational_checks(seed).items():
                self.assert_checks(checks, f"seed {seed}, step {step}")

    def test_control_wrong_sign_of_shift_is_rejected(self):
        """Control: M_a[M_a^{-1}(calC)] = calC holds exactly and fails for -M_a^{-1}(calC), so the zero test can fail."""
        inputs = rational_inputs(1)
        squared_momentum = sp.Symbol("s", real=True)
        inverse = closed_form_inverse_of_invariant_matrix(inputs, squared_momentum)
        residual = lo_pole_operator(inverse, inputs["chi_values"], inputs["anomaly_mass_squared"], squared_momentum) - invariant_matrix(inputs)
        self.assertEqual(residual.applyfunc(sp.cancel), sp.zeros(3))
        wrong = (lo_pole_operator(-inverse, inputs["chi_values"], inputs["anomaly_mass_squared"], squared_momentum)
                 - invariant_matrix(inputs)).applyfunc(sp.cancel)
        self.assertNotEqual(wrong, sp.zeros(3))


if __name__ == "__main__":
    unittest.main()
