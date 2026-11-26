"""Exact checks of the radial vertices and of diagrams 6.1 and 6.7, including the closed forms
displayed in the records murayama-radial-loop-six-one.md and murayama-radial-loop-six-seven.md.

Run:  PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_radial_loop_six_seven
"""
import unittest

import sympy as sp

import derive_radial_loop_generator_products as products
import derive_radial_loop_six_seven as calculation
import derive_radial_loop_trace_contractions as common
import expand_radial_loop_lagrangian as lagrangian
from group_radial_loop_trace_permutations import PermutationSum, polynomial_difference, trace_polynomial

B0, fpi, fa, lam, ms, s = common.B0, common.fpi, common.fa, common.lam, common.ms, common.s
tr, S, P = common.tr, common.S, common.P
xi_a, xi_0, k_a0 = common.xi_a, common.xi_0, common.k_a0
la, l0, lp, lm = common.la, common.l0, common.lp, common.lm
cubic_bracket = lam - B0 * fpi**2 * tr(S[3, 0]) / fa**4
charged_trace = tr(S[1, 0] * PermutationSum(lp, lm))


def exact_difference(first, second):
    return polynomial_difference(trace_polynomial(first), trace_polynomial(second))


class RadialVertexChecks(unittest.TestCase):
    def test_radial_cubic_equals_lagrangian_expansion(self):
        expected = lagrangian.vertex_polynomial((lagrangian.RADIAL,) * 3, 0)
        self.assertEqual(polynomial_difference(expected, trace_polynomial(calculation.radial_cubic)), {})

    def test_general_leg_vertices_equal_lagrangian_expansion(self):
        legs = ((calculation.leg_x[0].name, calculation.leg_x[1], 0, 'x'),
                (calculation.leg_y[0].name, calculation.leg_y[1], 0, 'y'))
        for compact, radial_count in ((calculation.radial_vertex, 1), (calculation.radial_pair_vertex, 2)):
            expected = lagrangian.vertex_polynomial((lagrangian.RADIAL,) * radial_count + legs, 0,
                                                    calculation.GENERAL_DOTS)
            self.assertEqual(polynomial_difference(expected, trace_polynomial(compact)), {})

    def test_radial_cubic_closed_form(self):
        self.assertEqual(exact_difference(calculation.radial_cubic, -6 * fa * cubic_bracket), {})


class DiagramChecks(unittest.TestCase):
    def test_compact_products_equal_lagrangian_expansion(self):
        for label, compact in calculation.diagram_products():
            with self.subTest(label=label):
                _, blocks, lines = products.DIAGRAM_DATA[label]
                self.assertFalse(compact.atoms(sp.Float))
                self.assertEqual(lagrangian.compact_difference(blocks, lines, compact), {})

    def test_six_one_record_display(self):
        bracket = (2 * B0 * fa * fpi * (xi_a * tr(P[2, 1] * l0) + xi_0 * tr(P[2, 1] * la))
                   - B0 * fa**2 * tr(S[2, 0] * PermutationSum(la, l0)) + 4 * B0 * fpi**2 * xi_a * xi_0 * tr(S[2, 2])
                   - 2 * fa**2 * xi_a * xi_0 * k_a0)
        display = 3 * B0 * cubic_bracket / (fa**4 * (ms**2 - s)) * charged_trace * bracket
        self.assertEqual(exact_difference(products.diagram_total('6.1'), display), {})

    def test_six_seven_record_display(self):
        bracket = (2 * B0 * fa * fpi * (xi_a * tr(P[1, 1] * l0) + xi_0 * tr(P[1, 1] * la))
                   - B0 * fa**2 * tr(S[1, 0] * PermutationSum(la, l0)) + 4 * B0 * fpi**2 * xi_a * xi_0 * tr(S[1, 2])
                   - 4 * fa**2 * xi_a * xi_0 * k_a0)
        display = -9 * B0 * cubic_bracket**2 / (fa**2 * (ms**2 - s)**2) * charged_trace * bracket
        self.assertEqual(exact_difference(products.diagram_total('6.7'), display), {})

    def test_six_one_up_only_sources(self):
        """With the up-only sources no S_2s or P_2s is nonzero; only the ALP kinetic term remains."""
        zero_sources = {S[2, t]: 0 for t in range(7)} | {P[2, t]: 0 for t in range(7)}
        reduced = products.diagram_total('6.1').xreplace(zero_sources)
        display = -6 * B0 * cubic_bracket / (fa**2 * (ms**2 - s)) * charged_trace * xi_a * xi_0 * k_a0
        self.assertEqual(exact_difference(reduced, display), {})


if __name__ == '__main__':
    unittest.main()
