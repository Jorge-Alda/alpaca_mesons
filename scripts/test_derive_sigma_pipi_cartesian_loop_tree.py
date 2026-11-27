"""Unit tests for derive_sigma_pipi_cartesian_loop_tree.py (exact symbolic identities only; no numerical matrices).

Run: PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_cartesian_loop_tree
"""
import unittest

import sympy as sp

import sigma_pipi_cartesian_loop_algebra as cartesian
import derive_sigma_pipi_cartesian_loop_tree as tree_script
from sigma_pipi_loop_algebra import B0, EPSILON, F_PI, QUARTIC, TraceSum, assert_no_float
from sigma_pipi_cartesian_loop_algebra import UP_QUARK_MASS


class OnePointTest(unittest.TestCase):
    def test_tadpole_conditions(self):
        functions = tree_script.one_point_functions()
        xi = sp.Symbol('xi1', real=True)
        expected = TraceSum.single(('PC00', 'x1'), B0 * F_PI) + TraceSum.single(('SC01',), B0 * EPSILON * F_PI * xi)
        self.assertTrue((functions['pseudoscalar'] - expected).is_zero())

    def test_stationarity(self):
        expected = tree_script.TRACE_SYMBOLS[('SC10',)] * EPSILON ** 2 * B0 - QUARTIC * F_PI ** 2 / EPSILON ** 2
        self.assertEqual(sp.expand(tree_script.stationarity_solution() - expected), 0)


class QuadraticTest(unittest.TestCase):
    def test_quadratic_lagrangian(self):
        for name, ok in tree_script.quadratic_lagrangian_checks().items():
            self.assertTrue(ok, name)

    def test_no_radial_alp_mixing_from_potential(self):
        """V_phi has no s p term: the radial-pseudoscalar function is pure source."""
        value = tree_script.two_point_functions()['radial-pseudoscalar']
        self.assertNotIn((), value.terms)

    def test_pole_data(self):
        for name, ok in tree_script.pole_data_checks().items():
            self.assertTrue(ok, name)

    def test_up_only_quadratic_data(self):
        data = tree_script.up_only_quadratic_data()
        self.assertEqual(sp.expand(data['A'] - B0 * UP_QUARK_MASS), 0)
        self.assertEqual(sp.expand(data['kappa_sigma'] - B0 * UP_QUARK_MASS), 0)


class CPTest(unittest.TestCase):
    def test_cp_conditions(self):
        self.assertTrue(all(tree_script.cp_checks().values()))

    def test_cp_check_is_not_empty(self):
        """With the opposite sign for P^C the check fails (the condition is not automatic)."""
        def wrong(name):
            kind, _, q = cartesian.source_indices(name)
            return (-1) ** q
        self.assertFalse(all(tree_script.cp_checks(max_legs=2, sign=wrong).values()))

    def test_odd_traces_vanish(self):
        """<S^C_{r,q}> = eta <S^C_{r,q}^T> = -<S^C_{r,q}> for odd q."""
        for q in (1, 3):
            self.assertEqual(tree_script.cp_sign(f'SC1{q}'), -1)
        self.assertEqual(tree_script.cp_sign('SC02'), 1)
        self.assertEqual(tree_script.cp_sign('PC01'), 1)


class ChargeRuleTest(unittest.TestCase):
    def test_mass_operator(self):
        difference, trace_part = tree_script.mass_operator_commutes_with_charge()
        self.assertEqual(difference, {})
        self.assertEqual(trace_part, {})


class TreeTest(unittest.TestCase):
    def test_tree(self):
        tree = tree_script.tree_amplitude()
        self.assertTrue((tree - tree_script.tree_amplitude_direct()).is_zero())
        self.assertTrue((tree_script.tree_amplitude_compact().to_trace_sum() - tree).is_zero())
        self.assertTrue(tree_script.tree_has_no_momentum())
        assert_no_float(tree)

    def test_up_only_tree(self):
        substituted, direct = tree_script.up_only_tree()
        expected = -B0 * EPSILON * UP_QUARK_MASS / F_PI
        self.assertEqual(sp.expand(substituted - expected), 0)
        self.assertEqual(sp.expand(direct - expected), 0)

    def test_assumption_table(self):
        self.assertEqual([number for number, _, _ in tree_script.ASSUMPTIONS], list(range(1, 11)))


class GeneratorTest(unittest.TestCase):
    """Completeness of the nine U(3) generators, checked once with explicit generators and
    symbolic 3x3 matrices (never printed).  It justifies the completeness rule and the
    statement that the meson matrices span all Hermitian matrices (CP conditions)."""

    @staticmethod
    def generators():
        lam = [sp.sqrt(sp.Rational(2, 3)) * sp.eye(3)]
        lam += [sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]]), sp.Matrix([[0, -sp.I, 0], [sp.I, 0, 0], [0, 0, 0]]),
                sp.Matrix([[1, 0, 0], [0, -1, 0], [0, 0, 0]]), sp.Matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]]),
                sp.Matrix([[0, 0, -sp.I], [0, 0, 0], [sp.I, 0, 0]]), sp.Matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]]),
                sp.Matrix([[0, 0, 0], [0, 0, -sp.I], [0, sp.I, 0]]),
                sp.diag(1, 1, -2) / sp.sqrt(3)]
        return lam

    def test_normalization_and_completeness(self):
        lam = self.generators()
        for i in range(9):
            for j in range(9):
                self.assertEqual(sp.simplify((lam[i] * lam[j]).trace() - 2 * (i == j)), 0)
        a = sp.Matrix(3, 3, sp.symbols('a0:9'))
        b = sp.Matrix(3, 3, sp.symbols('b0:9'))
        one = sum(((a * g).trace() * (g * b).trace() for g in lam), 0)
        two = sum(((a * g * b * g).trace() for g in lam), 0)
        self.assertEqual(sp.expand(one - 2 * (a * b).trace()), 0)
        self.assertEqual(sp.expand(two - 2 * a.trace() * b.trace()), 0)


if __name__ == '__main__':
    unittest.main()
