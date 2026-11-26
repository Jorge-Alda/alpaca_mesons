"""Unit tests for derive_sigma_pipi_loop_vertices.py (exact symbolic identities only).

Run: PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_loop_vertices
"""
import unittest

import sympy as sp

import sigma_pipi_loop_algebra as algebra
import derive_sigma_pipi_loop_vertices as vertices
from sigma_pipi_loop_algebra import B0, EPSILON, F_PI


class ConstructionTest(unittest.TestCase):
    def test_direct_expansion_equals_closed_rules(self):
        for key, (direct_ok, tree_ok) in vertices.compare_constructions().items():
            self.assertTrue(direct_ok, key)
            self.assertTrue(tree_ok, key)

    def test_no_float_in_rules(self):
        for n in range(3, vertices.MAX_LEGS + 1):
            for radial in range(n + 1):
                algebra.assert_no_float(vertices.closed_vertex(n, radial))


class RecordFormulaTest(unittest.TestCase):
    """The closed V_pot, V_kin^(4) and V_kin,a of the estimate record lack
    f_pi^{2-n}, 1/f_pi^2 and 1/f_pi^r; V(sigma;x,y) and V(sigma,sigma;x,y) are right."""

    def test_record_checks(self):
        checks = vertices.record_checks()
        for key, value in checks.items():
            if key[0] == 'V_pot':
                self.assertFalse(value['record equals rule'], key)
                self.assertTrue(value['record times f_pi^(2-n) equals rule'], key)
            elif key[0] in ('V_kin4', 'V_kin_a'):
                self.assertFalse(value['record equals rule'], key)
                self.assertTrue(list(value.values())[1], key)
            else:
                self.assertTrue(all(value.values()), key)

    def test_sigma_cubic(self):
        legs = [algebra.radial_leg({'q1': 1}), algebra.radial_leg({'q2': 1}), algebra.radial_leg({'q3': 1})]
        value = algebra.vertex(legs, algebra.symbolic_dot, eliminate_lambda=False)
        expected = (algebra.TraceSum.scalar(-6 * sp.Symbol('lambda', positive=True) * F_PI / EPSILON)
                    + algebra.TraceSum.single(('S30',), 6 * B0 * EPSILON ** 3 / F_PI))
        self.assertTrue((value - expected).is_zero())

    def test_lambda_elimination(self):
        """lambda f_a from m_sigma^2 = 2 lambda f_a^2 + B_0 eps^2(<S10> - 2<S20>)."""
        lam, s10, s20 = sp.symbols('lambda s10 s20')
        m_sigma_squared = 2 * lam * (F_PI / EPSILON) ** 2 + B0 * EPSILON ** 2 * (s10 - 2 * s20)
        rule = algebra.radial_lambda_times_fa()
        value = sum(coefficient * {(): 1, (('S10',),): s10, (('S20',),): s20}[key]
                    for key, coefficient in rule.terms.items())
        value = value.subs(algebra.M_SIGMA ** 2, m_sigma_squared)
        self.assertEqual(sp.simplify(value - lam * F_PI / EPSILON), 0)


class DimensionTest(unittest.TestCase):
    def test_dimensions(self):
        self.assertTrue(all(vertices.dimension_check().values()))


if __name__ == '__main__':
    unittest.main()
