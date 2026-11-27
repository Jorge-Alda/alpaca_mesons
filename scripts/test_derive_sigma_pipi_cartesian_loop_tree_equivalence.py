"""Unit tests for derive_sigma_pipi_cartesian_loop_tree_equivalence.py (exact symbolic identities only).

Run: PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_cartesian_loop_tree_equivalence
"""
import unittest

import sympy as sp

import derive_sigma_pipi_cartesian_loop_tree as cartesian_tree
import derive_sigma_pipi_cartesian_loop_tree_equivalence as equivalence


class MapTest(unittest.TestCase):
    def test_map_series(self):
        for order in (1, 2):
            self.assertTrue(equivalence.map_check(order))

    def test_first_order_relation(self):
        relation = equivalence.cartesian_coefficients_from_polar(1)
        for (r, q), (real, imaginary) in relation.items():
            self.assertEqual(real, equivalence.polar_symbol('S', r, q))
            self.assertEqual(imaginary, equivalence.polar_symbol('P', r, q))


class TreeEquivalenceTest(unittest.TestCase):
    def test_on_shell_tree(self):
        self.assertTrue(equivalence.tree_equivalence().is_zero())

    def test_letters_are_distinct(self):
        """Without the coefficient relation the Cartesian and polar trees are different letters."""
        difference = cartesian_tree.tree_amplitude() - equivalence.polar_tree_input()
        self.assertFalse(difference.is_zero())

    def test_charged_legs_have_no_alp_component(self):
        legs = cartesian_tree.external_legs()
        self.assertEqual(legs['pi+'].alp_field, 0)
        self.assertEqual(legs['pi-'].alp_field, 0)


if __name__ == '__main__':
    unittest.main()
