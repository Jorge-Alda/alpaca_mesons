"""Unit tests for derive_sigma_pipi_cartesian_loop_large_epsilon_vertices.py (exact symbolic identities only).

Run: PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_cartesian_loop_large_epsilon_vertices
"""
import unittest

import sympy as sp

import derive_sigma_pipi_cartesian_loop_large_epsilon_vertices as large
from sigma_pipi_loop_algebra import B0, EPSILON, F_PI, QUARTIC
from sigma_pipi_cartesian_loop_algebra import (ALP_MASS_COEFFICIENT as A, KAPPA_SIGMA, SINGLET_MASS_SQUARED,
                                               UP_QUARK_MASS)

P01, P01_RIGHT, C0_TRACE = large.P01, large.P01_RIGHT, large.TRACE_C_SQUARED


class RegimeTest(unittest.TestCase):
    def test_resolvent_and_self_adjointness(self):
        self.assertTrue(large.resolvent_series_check())
        self.assertTrue(large.loop_operator_self_adjoint())

    def test_normalization_is_derivative(self):
        quantities = large.pole_quantities()
        self.assertEqual(sp.expand(quantities['N'] - quantities["N from D_a'"]), 0)

    def test_root_coefficients(self):
        u = large.INVERSE_EPSILON
        c0, c1 = large.TRACE_C_POWER[:2]
        root = large.alp_pole_root()
        self.assertEqual(sp.expand(root.coeff(u, -2) - A), 0)
        self.assertEqual(sp.expand(root.coeff(u, 0) - 2 * c0 / A), 0)
        self.assertEqual(sp.expand(root.coeff(u, 2) - (2 * c1 / A ** 2 - 4 * c0 ** 2 / A ** 3)), 0)
        self.assertFalse(any(symbol.name.startswith('b_') for symbol in root.free_symbols))


class PoleTableTest(unittest.TestCase):
    def test_generic_table(self):
        table = large.pole_data_table()
        expected = {
            'm_a^2': (2, A), 'm_sigma^2': (2, KAPPA_SIGMA), 'N': (0, 1), '1/N': (0, 1),
            'M^{-1}(C)': (-2, -B0 * P01 / A),
            'N - 1': (-2, 2 * B0 ** 2 * C0_TRACE / A ** 2),
            'xi_a = 1/sqrt(N)': (0, 1),
            'meson matrix of v_a': (-1, -2 * B0 * P01 / A),
            '(i) 1_9': (0, 1),
            '(ii) 1 - 1/N': (-2, 2 * B0 ** 2 * C0_TRACE / A ** 2),
            '(iii) -2 eps M^{-1}(C)/N': (-1, 2 * B0 * P01 / A),
            '(iv) -4 eps^2 M^{-1}(C) x M^{-1}(C)/N': (-2, -4 * B0 ** 2 * P01 * P01_RIGHT / A ** 2),
            'm_a^2 - m_sigma^2 (A != kappa_sigma)': (2, A - KAPPA_SIGMA),
            'm_a^2 - m_sigma^2 (A = kappa_sigma)': (0, 2 * B0 ** 2 * C0_TRACE / A),
        }
        for name, (power, coefficient) in expected.items():
            self.assertEqual(table[name][0], power, name)
            self.assertEqual(sp.simplify(table[name][1] - coefficient), 0, name)


class UpOnlyTest(unittest.TestCase):
    def test_values(self):
        values = large.up_only_values()
        self.assertEqual(sp.expand(values['A'] - B0 * UP_QUARK_MASS), 0)
        self.assertEqual(sp.expand(values['kappa_sigma'] - B0 * UP_QUARK_MASS), 0)
        self.assertEqual(sp.expand(values['c_0'] - B0 ** 2 * UP_QUARK_MASS ** 2), 0)
        self.assertEqual(sp.expand(values['c_1'] - B0 ** 2 * UP_QUARK_MASS ** 2
                                   * (2 * B0 * UP_QUARK_MASS + SINGLET_MASS_SQUARED / 3)), 0)

    def test_root_against_exact_pole_function(self):
        self.assertTrue(large.up_only_root_check())

    def test_exact_check_is_not_empty(self):
        u = large.INVERSE_EPSILON
        root = B0 * UP_QUARK_MASS / u ** 2      # leading term only: fails at relative order u^2
        value = (u ** 2 * large.up_only_exact_alp_pole_function()).subs(large.SQUARED_MOMENTUM, root)
        self.assertNotEqual(sp.expand(sp.series(value, u, 0, 8).removeO()), 0)

    def test_mass_difference_sign(self):
        difference, (power, coefficient) = large.up_only_mass_difference()
        self.assertEqual(power, 0)
        self.assertEqual(sp.expand(coefficient - 2 * B0 * UP_QUARK_MASS), 0)
        self.assertTrue(coefficient.is_positive)
        u = large.INVERSE_EPSILON
        self.assertEqual(sp.expand(difference.coeff(u, 2) - (2 * SINGLET_MASS_SQUARED / 3 - 2 * QUARTIC * F_PI ** 2)), 0)


class VertexPowerTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.table = large.vertex_power_table()

    def test_pattern_count(self):
        self.assertEqual(len(self.table), 46)

    def test_exceptions(self):
        exceptions = {large.pattern_text(pattern) for pattern in large.exceptions_to_rule(self.table)}
        self.assertEqual(exceptions, {'sss', 'spp', 'ssss', 'sspp', 'pppp'})
        rows = {large.pattern_text(pattern): row for pattern, row in self.table.items()}
        self.assertEqual(rows['sss']['generic'], {-1, 3})
        self.assertEqual(rows['spp']['generic'], {-1, 3})
        for name in ('ssss', 'sspp', 'pppp'):
            self.assertEqual(rows[name]['generic'], {0, 4})
            self.assertEqual(rows[name]['up-only'], {0})
        self.assertEqual(rows['MMMM']['generic, loop momenta'], {0, 2})
        for name in ('MMM', 'sMMM', 'MMMMM'):
            self.assertEqual(rows[name]['generic'], set())

    def test_up_only_nonvanishing(self):
        rows = {large.pattern_text(pattern): row for pattern, row in self.table.items()}
        nonvanishing = {name for name, row in rows.items() if row['up-only']}
        self.assertEqual(nonvanishing, {'sss', 'spp', 'ssss', 'sspp', 'pppp', 'sMM', 'pMMM', 'MMMM', 'sMMMM'})

    def test_only_meson_kinetic_depends_on_momenta(self):
        for pattern, row in self.table.items():
            if large.pattern_text(pattern) != 'MMMM':
                self.assertEqual(row['generic'], row['generic, loop momenta'], pattern)

    def test_lambda_fixed(self):
        self.assertTrue(large.lambda_check())


class LegAndTreeTest(unittest.TestCase):
    def test_leg_factors(self):
        table = large.leg_factor_table()
        self.assertEqual(table['(ii) 1 - 1/N'][0], -2)
        self.assertEqual(sorted(table['(ii) 1 - 1/N'][1].values()), [-2, -1, 0])
        self.assertEqual(sorted(table['(iii) mixed'][1].values()), [-1, 0])
        self.assertEqual(sorted(table['xi_a = 1/sqrt(N)'][1].values()), [0, 1])

    def test_tree(self):
        powers, up_only = large.tree_leading()
        self.assertEqual(powers, {1})
        self.assertEqual(sp.expand(up_only + B0 * UP_QUARK_MASS * EPSILON / F_PI), 0)


if __name__ == '__main__':
    unittest.main()
