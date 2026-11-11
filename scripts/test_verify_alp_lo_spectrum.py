"""Tests for the solver checks on LO, single-pole and field-change data (tasks 1.2-1.4)."""

import unittest

import sympy as sp

from derive_alp_inverse_propagator_congruence import neutral_test_quadratic_form
from solve_alp_meson_spectrum import solve_spectrum_exact, solve_spectrum_numerical
from derive_alp_lo_neutral_blocks import illustrative_lo_inputs
from verify_alp_lo_spectrum import (
    LO_SOURCE_STRUCTURES,
    derive,
    field_change_checks_exact,
    lo_blocks_numerical,
    neutral_kaon_pole_checks,
    neutral_versus_full_checks,
    single_pole_checks_exact,
    single_pole_checks_numerical,
)


class VerifyLoSpectrumTest(unittest.TestCase):
    def test_all_checks_pass(self):
        for group, checks in derive().items():
            for label, value in checks.items():
                if isinstance(value, bool):
                    self.assertTrue(value, f"{group}: {label}")

    def test_single_pole_record_on_every_exact_pole(self):
        data = neutral_test_quadratic_form()
        entries = single_pole_checks_exact(solve_spectrum_exact(data["kinetic"], data["mass_squared"]))
        self.assertEqual(len(entries), 4)

    def test_alp_like_pole_of_lo_blocks(self):
        for source_structure in LO_SOURCE_STRUCTURES:
            kinetic, mass_squared = lo_blocks_numerical(source_structure)
            self.assertEqual(kinetic.shape, (6, 6))
            spectrum = solve_spectrum_numerical(kinetic, mass_squared)
            checks = single_pole_checks_numerical(spectrum)
            self.assertAlmostEqual(checks["LSZ factor N^(-1/2)"], checks["ALP component v_a"], places=10)

    def test_neutral_system_agrees_with_full_system(self):
        for source_structure in LO_SOURCE_STRUCTURES:
            for label, passed in neutral_versus_full_checks(source_structure).items():
                self.assertTrue(passed, f"{source_structure}: {label}")

    def test_full_system_check_detects_charge_violation(self):
        checks = neutral_versus_full_checks("charge_violating")
        self.assertFalse(checks["neutral restriction of the 10x10 D^-1(s) = 6x6 D^-1(s)"])

    def test_kaon_pole_multiplicity(self):
        expected = {"flavor_diagonal": 2, "cp_conserving": 1, "neutral": 1}
        for source_structure, count in expected.items():
            spectrum = solve_spectrum_numerical(*lo_blocks_numerical(source_structure))
            checks = neutral_kaon_pole_checks(spectrum, illustrative_lo_inputs(source_structure))
            self.assertEqual(checks["number of poles at m_K0^2"], count, source_structure)
            self.assertTrue(checks["its pole vectors have zero ALP and (P_3, P_8, P_0) components"], source_structure)

    def test_degenerate_pole_under_field_change(self):
        conversion = sp.Matrix([[1, sp.Rational(1, 2), 0], [0, 1, sp.Rational(-1, 3)], [sp.Rational(1, 4), 0, 1]])
        kinetic = conversion.T * conversion
        mass_squared = conversion.T * sp.diag(sp.Rational(1, 5), sp.Rational(1, 5), 2) * conversion
        for label, passed in field_change_checks_exact(kinetic, mass_squared, seed=11).items():
            self.assertTrue(passed, label)


if __name__ == "__main__":
    unittest.main()
