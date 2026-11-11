"""Tests for the LO kinetic and mass matrices of the neutral mesons and the ALP."""

import unittest

import sympy as sp

from derive_alp_lo_neutral_blocks import (
    ALL_GENERATOR_INDICES,
    NEUTRAL_FIELD_LABELS,
    blocks_from_explicit_lagrangian,
    closed_form_blocks,
    derive,
    field_labels,
    illustrative_lo_inputs,
)

ALP = 5


class LoNeutralBlocksTest(unittest.TestCase):
    def test_closed_forms_match_lagrangian_expansions(self):
        for source_structure in ("flavor_diagonal", "cp_conserving", "neutral", "charge_violating"):
            result = derive(source_structure)
            self.assertEqual(result["kinetic"].shape, (6, 6))
            self.assertEqual(result["full_kinetic"].shape, (10, 10))
            for label, passed in result["checks"].items():
                self.assertTrue(passed, f"{source_structure}: {label}")

    def test_default_basis_is_neutral(self):
        self.assertEqual(field_labels(), NEUTRAL_FIELD_LABELS)
        self.assertEqual(NEUTRAL_FIELD_LABELS, ("P_3", "P_8", "P_0", "P_6", "P_7", "a"))

    def test_flavor_diagonal_sources_mix_only_p3_p8_p0(self):
        kinetic, mass_squared = closed_form_blocks(**illustrative_lo_inputs("flavor_diagonal"))
        for position in (3, 4):
            self.assertEqual(kinetic[position, ALP], 0)
            self.assertEqual(mass_squared[position, ALP], 0)
        for position in (0, 1, 2):
            self.assertNotEqual(kinetic[position, ALP], 0)
            self.assertNotEqual(mass_squared[position, ALP], 0)

    def test_real_ds_entries_decouple_p7(self):
        kinetic, mass_squared = closed_form_blocks(**illustrative_lo_inputs("cp_conserving"))
        for position in (0, 1, 2, 3, 5):
            self.assertEqual(kinetic[4, position], 0)
            self.assertEqual(mass_squared[4, position], 0)
        kinetic, mass_squared = closed_form_blocks(**illustrative_lo_inputs("neutral"))
        self.assertNotEqual(kinetic[4, ALP], 0)
        self.assertNotEqual(mass_squared[4, ALP], 0)

    def test_comparison_detects_a_wrong_entry(self):
        # Guard against a vacuous comparison: a different c_G changes M^2_ia and M^2_aa.
        inputs = illustrative_lo_inputs("neutral")
        _, mass_squared = closed_form_blocks(**inputs, generator_indices=ALL_GENERATOR_INDICES)
        _, mass_wrong = blocks_from_explicit_lagrangian(**{**inputs, "gluon_coupling": sp.Rational(1, 3)})
        self.assertNotEqual(sp.simplify(mass_squared[9, 9] - mass_wrong[9, 9]), 0)
        self.assertNotEqual(sp.simplify(mass_squared[0, 9] - mass_wrong[0, 9]), 0)

    def test_kinetic_alp_entry_has_the_mixing_square(self):
        kinetic, _ = closed_form_blocks(**illustrative_lo_inputs("neutral"))
        mixing = kinetic[:ALP, ALP]
        self.assertEqual(sp.simplify(kinetic[ALP, ALP] - 1 - (mixing.T * mixing)[0]), 0)

    def test_charge_violating_kinetic_alp_entry_needs_charged_mesons(self):
        # With charged K_ia, the 6x6 K_aa (full trace) exceeds 1 + sum over neutral K_ia^2.
        kinetic, _ = closed_form_blocks(**illustrative_lo_inputs("charge_violating"))
        mixing = kinetic[:ALP, ALP]
        self.assertNotEqual(sp.simplify(kinetic[ALP, ALP] - 1 - (mixing.T * mixing)[0]), 0)


if __name__ == "__main__":
    unittest.main()
