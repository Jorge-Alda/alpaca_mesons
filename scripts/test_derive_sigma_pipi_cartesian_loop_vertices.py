"""Unit tests for derive_sigma_pipi_cartesian_loop_vertices.py (exact symbolic identities only; no numerical matrices).

Run: PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_cartesian_loop_vertices
"""
import unittest

import sympy as sp

import sigma_pipi_cartesian_loop_algebra as cartesian
import derive_sigma_pipi_cartesian_loop_vertices as vertices
from sigma_pipi_loop_algebra import B0, EPSILON, F_PI, M_SIGMA, QUARTIC, TraceSum, assert_no_float, symbolic_dot
from sigma_pipi_cartesian_loop_algebra import UP_QUARK_MASS


class ConstructionTest(unittest.TestCase):
    def test_direct_expansion_equals_closed_rules(self):
        results = vertices.compare_constructions()
        self.assertEqual(len(results), 15)
        for key, (direct_ok, tree_ok) in results.items():
            self.assertTrue(direct_ok, key)
            self.assertTrue(tree_ok, key)

    def test_no_float(self):
        for n, radial in vertices.PATTERNS:
            assert_no_float(vertices.closed_vertex(n, radial, eliminate_lambda=True))


class SymmetryAndDimensionTest(unittest.TestCase):
    def test_relabelling(self):
        self.assertTrue(all(vertices.relabelling_symmetry().values()))

    def test_relabelling_detects_asymmetry(self):
        value = TraceSum.single(('SC00', 'x1', 'x2', 'x3'))
        swapped = vertices.relabel(value, {1: 2, 2: 1, 3: 3})
        self.assertFalse((swapped - value).is_zero())

    def test_dimensions(self):
        self.assertTrue(all(vertices.dimension_check().values()))


class RadialPotentialTest(unittest.TestCase):
    def test_table(self):
        f_a = F_PI / EPSILON
        expected = {(3, 0): -6 * QUARTIC * f_a, (1, 2): -2 * QUARTIC * f_a,
                    (4, 0): -6 * QUARTIC, (2, 2): -2 * QUARTIC, (0, 4): -6 * QUARTIC}
        table = vertices.radial_potential_table()
        self.assertEqual(set(table), set(expected))
        for key, value in expected.items():
            self.assertEqual(sp.simplify(table[key] - value), 0, key)

    def test_no_momentum_dependence(self):
        self.assertTrue(vertices.radial_potential_has_no_momentum())
        self.assertEqual(vertices.momentum_dependent_patterns(), [(4, 0)])
        self.assertTrue(vertices.sp_kinetic_term_gives_no_vertex())

    def test_lambda_elimination(self):
        """lambda f_a from m_sigma^2 = 2 lambda f_a^2 + eps^2 B_0 (<S^C_10> - 2<S^C_20>)."""
        s10, s20 = sp.symbols('s10 s20')
        m_sigma_squared = 2 * QUARTIC * (F_PI / EPSILON) ** 2 + B0 * EPSILON ** 2 * (s10 - 2 * s20)
        rule = cartesian.radial_lambda_times_fa()
        value = sum(coefficient * {(): 1, (('SC10',),): s10, (('SC20',),): s20}[key]
                    for key, coefficient in rule.terms.items())
        self.assertEqual(sp.simplify(value.subs(M_SIGMA ** 2, m_sigma_squared) - QUARTIC * F_PI / EPSILON), 0)


class UpOnlyTest(unittest.TestCase):
    def test_up_only_limit(self):
        self.assertTrue(all(vertices.up_only_checks().values()))
        self.assertEqual(vertices.up_only_source_patterns(), {(0, 0), (1, 0), (0, 1)})

    def test_up_only_radial_cubic(self):
        value = cartesian.up_only_substitution(vertices.closed_vertex(3, 1))
        xi2, xi3 = sp.symbols('xi2 xi3', real=True)
        expected = (TraceSum.single(('Pu', 'x2', 'x3'), -B0 * EPSILON * UP_QUARK_MASS / (2 * F_PI))
                    + TraceSum.single(('Pu', 'x3', 'x2'), -B0 * EPSILON * UP_QUARK_MASS / (2 * F_PI))
                    + TraceSum.scalar(-2 * QUARTIC * F_PI / EPSILON * xi2 * xi3))
        self.assertTrue((value - expected).is_zero())


if __name__ == '__main__':
    unittest.main()
