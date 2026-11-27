"""Unit tests for derive_sigma_pipi_cartesian_loop_estimate.py (exact symbolic identities; no numerical matrices).

Run: PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_cartesian_loop_estimate
"""
import pathlib
import unittest

import sympy as sp

import derive_sigma_pipi_cartesian_loop_estimate as estimate
from sigma_pipi_loop_algebra import EPSILON
from test_derive_sigma_pipi_cartesian_loop_traces import FORBIDDEN_COMPONENTS, FORBIDDEN_LOOP_VALUES

RECORD = (pathlib.Path(__file__).resolve().parent.parent / 'memory' / 'alp' / 'radial-mode' / 'finite-mixing'
          / 'two-pion-loop' / 'cartesian-estimate' / 'cartesian-sigma-pipi-loop-estimate.md')


class AssemblyTest(unittest.TestCase):
    def test_sum_of_parts(self):
        self.assertTrue(estimate.sum_of_parts_check())

    def test_pion_residues_equal(self):
        self.assertTrue(estimate.pion_residues_equal())

    def test_no_float(self):
        self.assertTrue(estimate.no_float())

    def test_tree(self):
        tree = estimate.estimate()['tree']
        expected = estimate.TREE_TRACE.scale(-estimate.B0 * EPSILON / (2 * estimate.F_PI))
        self.assertTrue((tree - expected).is_zero())


class DimensionTest(unittest.TestCase):
    def test_maximum_is_homogeneous(self):
        self.assertTrue(estimate.maximum_is_homogeneous())

    def test_dimensions(self):
        checks = estimate.dimension_checks()
        failed = [name for name, value in checks.items() if not value]
        self.assertEqual(failed, [])

    def test_wrong_dimension_fails(self):
        """Negative control: M_1^est does not scale as c^0 or c^2."""
        total = estimate.estimate()['total']
        self.assertFalse(estimate.has_dimension(total, 0))
        self.assertFalse(estimate.has_dimension(total, 2))


class RatioTest(unittest.TestCase):
    def test_generic_ratio(self):
        value = estimate.ratio()
        self.assertIsNotNone(value)
        self.assertTrue(value.has(sp.Abs))

    def test_vanishing_tree(self):
        tree, value, total = estimate.vanishing_tree_case()
        self.assertTrue(tree.is_zero())
        self.assertIsNone(value)
        self.assertFalse(total.is_zero())
        self.assertEqual(len(total.terms), 57)


class SummaryTest(unittest.TestCase):
    def test_summary(self):
        data = estimate.summary()
        self.assertEqual(data['kept 1PI graphs'], 26)
        self.assertEqual(data['kept s self-energies'], 6)
        self.assertEqual(data['kept pi self-energies (each charged leg)'], 5)
        self.assertEqual(set(data['source letters in M_1^est']),
                         {'SC01', 'SC02', 'SC10', 'SC11', 'SC12', 'SC20', 'SC21', 'SC30', 'PC01', 'PC11', 'PC20'})
        self.assertNotIn('SC00', data['source letters in M_1^est'])
        self.assertEqual(data['graphs with S^C_00'], ['C12', 'C13', 'C15'])
        self.assertEqual(data['products with S^C_00 in M_1^est'], 0)
        table, vanish, absent = estimate.source_by_scale('SC00')
        self.assertTrue(vanish)
        self.assertTrue(absent)
        self.assertEqual({label for counts in table.values() for label in counts}, {'C12', 'C13', 'C15'})
        self.assertEqual(len(table), 2)
        self.assertEqual(data['products in M_1^est'],
                         data['products in the 1PI sum'] + data['products in the external-leg term']
                         - len(set(estimate.estimate()['one particle irreducible sum'].terms)
                               & set(estimate.estimate()['external legs'].terms)))


class RecordTest(unittest.TestCase):
    def test_record_equals_script_output(self):
        text = RECORD.read_text()
        self.assertIn(estimate.markdown_summary(), text)

    def test_record_content(self):
        text = RECORD.read_text()
        for token in FORBIDDEN_COMPONENTS + FORBIDDEN_LOOP_VALUES:
            self.assertNotIn(token, text, token)
        self.assertIn('No loop integral', text)
        self.assertIn('-21/(64\\pi^2)', text)
        self.assertIn('authoritative', text)
        self.assertNotIn('action', text.lower().replace('interaction', ''))


if __name__ == '__main__':
    unittest.main()
