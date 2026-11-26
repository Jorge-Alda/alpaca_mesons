"""Exact check of the radial--ALP bubbles 5.1--5.4: compact products against the Lagrangian expansion.

Run:  PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_radial_loop_class_five_bubbles
"""
import unittest

import sympy as sp

import derive_radial_loop_class_five_bubbles as calculation
import derive_radial_loop_generator_products as products
import expand_radial_loop_lagrangian as lagrangian


class LagrangianExpansionAgreement(unittest.TestCase):
    def test_compact_products_equal_lagrangian_expansion(self):
        for label, compact in calculation.diagram_products():
            with self.subTest(label=label):
                _, blocks, lines = products.DIAGRAM_DATA[label]
                self.assertFalse(compact.atoms(sp.Float))
                self.assertEqual(lagrangian.compact_difference(blocks, lines, compact), {})


if __name__ == '__main__':
    unittest.main()
