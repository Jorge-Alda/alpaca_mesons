"""Tests for the pole-alignment shift as a constant-generator chiral redefinition."""

import unittest

import sympy as sp

from derive_alp_pole_alignment_chiral_generator import derive, general_generator_linear_map, source_data


class PoleAlignmentChiralGeneratorTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = derive()

    def test_all_checks(self):
        for label, passed in self.result["checks"].items():
            self.assertTrue(passed, label)

    def test_any_redefinition_has_shift_type_linear_part(self):
        self.assertTrue(general_generator_linear_map())

    def test_singlet_trace(self):
        data = source_data()
        w_0 = sp.Symbol("w_0", real=True)
        expected = -data["f_a"] / data["f_pi"] * sp.sqrt(6) * w_0
        self.assertEqual(sp.simplify(self.result["trace of G_R - G_L"] - expected), 0)


if __name__ == "__main__":
    unittest.main()
