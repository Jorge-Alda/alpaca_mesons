"""Unit tests for derive_sigma_pipi_loop_estimate.py (exact symbolic identities only).

Run: PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_loop_estimate
"""
import unittest

import sympy as sp

import sigma_pipi_loop_algebra as algebra
import derive_sigma_pipi_loop_estimate as estimate


ESTIMATE_TERMS = 144


class EstimateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = estimate.assemble()

    def test_dimensions(self):
        for name, value in estimate.dimension_checks(self.result).items():
            self.assertTrue(value, name)

    def test_vanishing_tree(self):
        checks = estimate.vanishing_tree_checks(self.result)
        self.assertTrue(checks['M_0 vanishes'])
        self.assertTrue(checks['leg terms vanish'])
        self.assertIsNone(checks['R'])
        self.assertTrue(checks['M_1^est nonzero'])

    def test_charged_residues_equal(self):
        residues = self.result['residues']
        self.assertTrue((residues['pi+'] - residues['pi-']).is_zero())

    def test_exact(self):
        algebra.assert_no_float(self.result['estimate'])
        self.assertEqual(len(self.result['contributions']), 26)
        self.assertEqual(len(self.result['estimate'].terms), ESTIMATE_TERMS)
        self.assertEqual(estimate.sources_in_estimate(self.result),
                         ['P01', 'P11', 'P20', 'S01', 'S02', 'S10', 'S11', 'S12', 'S20', 'S21', 'S30'])

    def test_larger_mass_rescaling(self):
        """max(c x, c y, ...) = c max(x, y, ...) is applied exactly, and the
        conditions of Piecewise coefficients are homogeneous."""
        m_p = algebra.M_MESON
        for larger in (sp.Max(algebra.M_SIGMA, algebra.M_A), sp.Max(m_p, algebra.M_A),
                       sp.Max(m_p, algebra.M_A, algebra.M_SIGMA)):
            self.assertEqual(sp.simplify(estimate.rescale_value(larger ** -2) - estimate.SCALE ** -2 * larger ** -2), 0)
            self.assertEqual(estimate.rescale_value(algebra.M_SIGMA / larger), algebra.M_SIGMA / larger)
        condition = sp.LessThan(algebra.M_SIGMA ** 2 - algebra.M_PI ** 2, sp.Max(m_p, algebra.M_A) ** 2)
        branch = sp.Piecewise((algebra.M_A ** 2, condition), (algebra.M_SIGMA ** 2, True))
        self.assertTrue(estimate.graphs.expression_is_zero(estimate.rescale_value(branch)
                                                           - estimate.SCALE ** 2 * branch))
        with self.assertRaises(AssertionError):
            estimate.rescale_value(sp.Piecewise((1, sp.LessThan(algebra.M_SIGMA, algebra.M_A ** 2)), (0, True)))

    def test_table(self):
        rows = estimate.graph_table(self.result)
        self.assertEqual(len(rows), 26 + 6 + 5)
        dependent = {row['label'] for row in rows if row['routing']['numerical'] or row['routing']['mass-dependent']}
        self.assertEqual(len(dependent), 10)
        self.assertEqual(sum(row['routing']['mass-dependent'] for row in rows), 80)

    def test_contributions_are_routing_maxima(self):
        for entry in self.result['entries']:
            if entry['group'] == 'one_particle_irreducible':
                self.assertTrue((self.result['contributions'][entry['label']] - entry['value']).is_zero())

    def test_estimate_is_sum_of_parts(self):
        total = self.result['leg_term']
        for value in self.result['contributions'].values():
            total = total + value
        self.assertTrue((total - self.result['estimate']).is_zero())
        legs = (self.result['residues']['sigma'] + self.result['residues']['pi+']
                + self.result['residues']['pi-']).scale(sp.Rational(1, 2))
        self.assertTrue((self.result['tree'] * legs - self.result['leg_term']).is_zero())

    def test_ratio_defined(self):
        self.assertIsNotNone(estimate.ratio(self.result))


if __name__ == '__main__':
    unittest.main()
