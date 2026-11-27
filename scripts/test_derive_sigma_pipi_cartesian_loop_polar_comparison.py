"""Unit tests for derive_sigma_pipi_cartesian_loop_polar_comparison.py (exact symbolic identities; no numerical matrices).

Run: PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_cartesian_loop_polar_comparison
(about ten minutes: both estimates are rebuilt from their scripts).
"""
import pathlib
import unittest

import sympy as sp

import derive_sigma_pipi_cartesian_loop_polar_comparison as comparison
import derive_sigma_pipi_cartesian_loop_tree_equivalence as tree_equivalence
from sigma_pipi_loop_algebra import B0, F_PI, M_MESON
from sigma_pipi_cartesian_loop_algebra import ALP_MASS_COEFFICIENT, UP_QUARK_MASS
from test_derive_sigma_pipi_cartesian_loop_traces import FORBIDDEN_COMPONENTS, FORBIDDEN_LOOP_VALUES

RECORD = (pathlib.Path(__file__).resolve().parent.parent / 'memory' / 'alp' / 'radial-mode' / 'finite-mixing'
          / 'two-pion-loop' / 'cartesian-estimate' / 'cartesian-sigma-pipi-loop-polar-comparison.md')
STATEMENTS = ('single graph', 'same sign', 'could cancel', 'depends on sources', 'exactly zero')


def polar(kind, radial, alp):
    return tree_equivalence.polar_symbol(kind, radial, alp)


class RelationTest(unittest.TestCase):
    def test_map(self):
        self.assertTrue(tree_equivalence.map_check(comparison.RELATION_ORDER))

    def test_relation_through_second_order(self):
        """Taylor expansion of 1 + (s + i p)/f_a = (1 + sigma/f_a) e^{i a/f_a}, done by hand:
        x = X + Y^2/2 + ..., y = Y - X Y + ..., so S^C_02 = S_02 + S_10/2, S^C_11 = S_11 - S_01, S^C_20 = S_20."""
        relation = comparison.coefficient_relation()
        expected = {(0, 0): polar('S', 0, 0), (1, 0): polar('S', 1, 0), (0, 1): polar('S', 0, 1),
                    (2, 0): polar('S', 2, 0), (1, 1): polar('S', 1, 1) - polar('S', 0, 1),
                    (0, 2): polar('S', 0, 2) + polar('S', 1, 0) / 2}
        for key, value in expected.items():
            self.assertEqual(sp.expand(relation[key] - value), 0, key)

    def test_parameters(self):
        self.assertTrue(all(comparison.parameter_checks().values()), comparison.parameter_checks())

    def test_up_only_sources(self):
        self.assertTrue(comparison.up_only_correspondence_holds())
        values, constant_is_letter = comparison.up_only_correspondence()
        self.assertTrue(constant_is_letter)
        self.assertEqual(len(values), sum(n + 1 for n in range(1, comparison.UP_ONLY_ORDER + 1)))


class ClassTest(unittest.TestCase):
    def test_classes(self):
        both, cartesian_only, polar_only = comparison.all_classes()
        self.assertEqual((cartesian_only, polar_only), ([], []))
        cartesian = sum(len(comparison.cartesian_classes()[key]) for key in both)
        polar_count = sum(len(comparison.polar_classes()[key]) for key in both)
        self.assertEqual((cartesian, polar_count), (37, 37))

    def test_d4_terms_coincide(self):
        """C3, C15, Cpi3 (Cartesian) and 3, 15, pi2 (polar): the excluded terms agree routing by routing."""
        for data in comparison.d4_comparison().values():
            self.assertTrue(data['same routing by routing'])
            self.assertTrue(data['same set of routing values'])

    def test_class_powers(self):
        """Regression: in r < 1 / polar region 1 every class has the same leading power in both."""
        both, _, _ = comparison.all_classes()
        for key in both:
            leading = comparison.class_leading(key)
            self.assertEqual(leading['cartesian'][0], leading['polar'][0], key)
        self.assertTrue(comparison.region_changes()[2])

    def test_totals(self):
        for region_name, indices in comparison.REGION_PAIRS:
            for index in indices:
                data = comparison.total_comparison(region_name, index)
                self.assertEqual((data['cartesian'][0], data['polar'][0]), (4, 4))
                self.assertFalse(data['difference'].is_zero())

    def test_routing_dependent_graphs(self):
        cartesian, polar_names = comparison.routing_dependent_graphs()
        self.assertEqual(cartesian, ['C12', 'C13', 'C15'])
        self.assertEqual(len(polar_names), 10)

    def test_statements(self):
        both, _, _ = comparison.all_classes()
        for key in both:
            statements = comparison.class_statements(key)
            self.assertIn(statements['cartesian'][0], STATEMENTS)
            self.assertIn(statements['polar'][0], STATEMENTS)


class UpOnlyTest(unittest.TestCase):
    def test_values_read_from_the_scripts(self):
        totals = comparison.up_only_totals()
        power, value, graphs = totals['cartesian']
        self.assertEqual(power, 2)
        self.assertEqual(sp.simplify(value + B0 ** 2 * UP_QUARK_MASS ** 2
                                     / (8 * sp.pi ** 2 * F_PI ** 2 * M_MESON ** 2)), 0)
        self.assertEqual(graphs, ['Csigma9'])
        power, value, graphs = totals['polar']
        self.assertEqual(power, 4)
        self.assertEqual(sp.simplify(value + 21 * B0 * UP_QUARK_MASS / (64 * sp.pi ** 2 * F_PI ** 2)), 0)
        self.assertEqual(sorted(graphs), ['13', '2', 'sigma4'])


class FurtherTest(unittest.TestCase):
    def test_vacuum_condition(self):
        (power, coefficient), d4_up_only = comparison.vacuum_condition()
        self.assertEqual(power, 4)
        self.assertEqual(sp.simplify(coefficient + ALP_MASS_COEFFICIENT / (32 * sp.pi ** 2 * F_PI ** 2)), 0)
        self.assertEqual(d4_up_only['C15'], 2)


class RecordTest(unittest.TestCase):
    def test_record_equals_script_output(self):
        text = RECORD.read_text()
        for block in comparison.markdown_blocks():
            self.assertIn(block, text)

    def test_record_content(self):
        text = RECORD.read_text()
        for token in FORBIDDEN_COMPONENTS + FORBIDDEN_LOOP_VALUES:
            self.assertNotIn(token, text, token)
        for statement in ('No loop integral', 'authoritative', 'not expected', '-21/(64\\pi^2)',
                          'after both derivations were complete'):
            self.assertIn(statement, text, statement)
        self.assertNotIn('action', text.lower().replace('interaction', ''))
        for shorthand in ('T_0', '\\mathcal N', '\\alpha', '\\beta'):
            self.assertNotIn(shorthand, text, shorthand)


if __name__ == '__main__':
    unittest.main()
