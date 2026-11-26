"""Exact checks of the quartic vertex with two internal ALP legs and of diagrams 3.1--3.3.

Run:  PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_radial_loop_three_one
"""
import unittest

import sympy as sp

import derive_radial_loop_generator_products as products
import derive_radial_loop_three_one as calculation
import derive_radial_loop_trace_contractions as common
import expand_radial_loop_lagrangian as lagrangian
from group_radial_loop_trace_permutations import polynomial_difference, trace_polynomial


class QuarticVertexChecks(unittest.TestCase):
    def test_quartic_vertex_equals_lagrangian_expansion(self):
        legs = (lagrangian.INTERNAL_ALP, lagrangian.INTERNAL_ALP,
                (str(calculation.lambda_x), calculation.xi_x, 0, 'x'),
                (str(calculation.lambda_y), calculation.xi_y, 0, 'y'))
        expected = lagrangian.vertex_polynomial(legs, common.ma**2, calculation.GENERAL_DOTS)
        self.assertEqual(polynomial_difference(expected, trace_polynomial(calculation.quartic_vertex)), {})

    def test_kinetic_part_is_single_commutator_trace(self):
        """Only the pairing with lambda_a once in each commutator survives, since [lambda_a, lambda_a] = 0."""
        kinetic = sp.expand(calculation.quartic_vertex).coeff(calculation.k_xy)
        expected = -common.tr(common.Commutator(common.la, calculation.lambda_x)
                              * common.Commutator(common.la, calculation.lambda_y)) / (12 * common.fpi**2)
        self.assertEqual(polynomial_difference(trace_polynomial(kinetic), trace_polynomial(expected)), {})


class ClassThreeChecks(unittest.TestCase):
    def test_compact_products_equal_lagrangian_expansion(self):
        for label, compact in calculation.diagram_products():
            with self.subTest(label=label):
                _, blocks, lines = products.DIAGRAM_DATA[label]
                self.assertFalse(compact.atoms(sp.Float))
                self.assertEqual(lagrangian.compact_difference(blocks, lines, compact), {})


if __name__ == '__main__':
    unittest.main()
