"""Tests for scripts/derive_alp_k0_doublet_exact_pole.py (K^0 doublet exact pole at diagonal chi_0).

Run: .venv/bin/python -m unittest discover -s scripts -p test_derive_alp_k0_doublet_exact_pole.py
"""

import unittest

import sympy as sp

from derive_alp_k0_doublet_exact_pole import (
    LAMBDA_6,
    LAMBDA_7,
    derive,
    doublet_generators,
    strange_down_part,
)


class KaonDoubletTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = derive()

    def test_all_checks_pass(self):
        for label, passed in self.result["checks"].items():
            self.assertTrue(passed, label)

    def test_real_entry_couples_only_lambda_6(self):
        coupled, decoupled = doublet_generators(sp.Rational(-3, 7))
        self.assertEqual(coupled, -LAMBDA_6)
        self.assertEqual(decoupled, LAMBDA_7)

    def test_strange_down_part(self):
        matrix = sp.Matrix([[1, 2, 3], [2, 4, 5 + sp.I], [3, 5 - sp.I, 6]])
        part = strange_down_part(matrix)
        self.assertEqual(part, sp.Matrix([[0, 0, 0], [0, 0, 5 + sp.I], [0, 5 - sp.I, 0]]))


if __name__ == "__main__":
    unittest.main()
