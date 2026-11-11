"""Tests for the reduction of the LO neutral system (P_3, P_8, P_0, P_6, P_7, a)."""

import unittest

import numpy as np
import sympy as sp

from derive_alp_lo_neutral_blocks import closed_form_blocks, illustrative_lo_inputs
from derive_alp_neutral_block_reduction import derive


class NeutralBlockReductionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = derive()

    def test_all_checks_pass(self):
        for part in ("a", "b", "c"):
            for label, passed in self.result[part]["checks"].items():
                self.assertTrue(passed, f"({part}) {label}")

    def test_kaon_pole_factor_of_determinant(self):
        # det D(s) has the factor (s - m_K0^2) for generic neutral rational sources.
        inputs = illustrative_lo_inputs("neutral")
        kinetic, mass_squared = closed_form_blocks(**inputs)
        squared_momentum = sp.Symbol("s")
        kaon_mass_squared = (inputs["quark_spurion"][1, 1] + inputs["quark_spurion"][2, 2]) / 2
        determinant = sp.expand((squared_momentum * kinetic - mass_squared).det(method="berkowitz"))
        self.assertEqual(sp.simplify(determinant.subs(squared_momentum, kaon_mass_squared)), 0)
        self.assertNotEqual(sp.simplify(sp.diff(determinant, squared_momentum).subs(squared_momentum, kaon_mass_squared)), 0)

    def test_no_constant_orthogonal_kaon_rotation_decouples_generic_case(self):
        # For complex, non-aligned d-s entries no angle makes P_7' decouple.
        kinetic, mass_squared = closed_form_blocks(**illustrative_lo_inputs("neutral"))
        kinetic_vector = np.array([float(kinetic[3, 5]), float(kinetic[4, 5])])
        mass_vector = np.array([float(mass_squared[3, 5]), float(mass_squared[4, 5])])
        self.assertGreater(abs(kinetic_vector[0] * mass_vector[1] - kinetic_vector[1] * mass_vector[0]), 1e-6)


if __name__ == "__main__":
    unittest.main()
