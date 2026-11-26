"""Unit tests for derive_sigma_pipi_loop_tree.py (exact symbolic identities only).

Run: PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_loop_tree
"""
import unittest

import sympy as sp

import sigma_pipi_loop_algebra as algebra
import derive_sigma_pipi_loop_tree as tree_script
from sigma_pipi_loop_algebra import B0, EPSILON, F_PI, M_PI, M_SIGMA


class KinematicsTest(unittest.TestCase):
    def test_dot_products(self):
        dots = tree_script.EXTERNAL_DOTS
        self.assertEqual(sp.simplify(dots[frozenset(('p1', 'p2'))] - (M_SIGMA ** 2 - 2 * M_PI ** 2) / 2), 0)
        self.assertEqual(sp.simplify(dots[frozenset(('ps', 'p1'))] + M_SIGMA ** 2 / 2), 0)
        self.assertEqual(sp.simplify(dots[frozenset(('ps', 'p2'))] + M_SIGMA ** 2 / 2), 0)

    def test_momentum_conservation(self):
        dot = algebra.make_dot(tree_script.external_basis_dot)
        total = {'ps': 1, 'p1': 1, 'p2': 1}
        self.assertEqual(sp.simplify(dot(total, total)), 0)


class ChargeRuleTest(unittest.TestCase):
    def test_charge_rule_proof(self):
        for length in (1, 2, 3, 4):
            derivation, trace = tree_script.charge_rule_proof(length)
            self.assertEqual(derivation, {})
            self.assertEqual(trace, {})

    def test_pole_operator_preserves_neutral_matrices(self):
        difference, trace_part = tree_script.pole_operator_commutes_with_charge()
        self.assertEqual(difference, {})
        self.assertEqual(trace_part, {})

    def test_charged_words_dropped(self):
        value = algebra.TraceSum.single(('S10', 'l+')) + algebra.TraceSum.single(('S10', 'l+', 'l-'))
        filtered = algebra.apply_charge_rule(value)
        self.assertEqual(list(filtered.terms), [(('S10', 'l+', 'l-'),)])


class TreeTest(unittest.TestCase):
    def test_tree_from_vertex_rule(self):
        tree = tree_script.tree_amplitude()
        expected = (algebra.TraceSum.single(('S10', 'l+', 'l-'), -B0 * EPSILON / (2 * F_PI))
                    + algebra.TraceSum.single(('S10', 'l-', 'l+'), -B0 * EPSILON / (2 * F_PI)))
        self.assertTrue((tree - expected).is_zero())
        self.assertTrue((tree_script.tree_amplitude_compact().to_trace_sum() - tree).is_zero())
        algebra.assert_no_float(tree)

    def test_up_only_agreement(self):
        ours, theirs = tree_script.up_only_tree_comparison()
        self.assertEqual(sp.simplify(ours - theirs), 0)
        self.assertEqual(sp.simplify(ours + B0 * EPSILON * sp.Symbol('m_u', positive=True) / F_PI), 0)

    def test_assumptions_are_linked(self):
        for text, source in tree_script.ASSUMPTIONS:
            self.assertTrue(source)
        sources = ' '.join(source for _, source in tree_script.ASSUMPTIONS)
        self.assertIn('radial-decay-polar-cartesian-general-sources.md', sources)
        self.assertIn('alp-lo-charged-meson-decoupling.md', sources)


if __name__ == '__main__':
    unittest.main()
