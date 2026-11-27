"""Unit tests for derive_sigma_pipi_cartesian_loop_limits.py and for the Cartesian limits record (exact symbolic identities).

Run: PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_cartesian_loop_limits
"""
import pathlib
import unittest

import sympy as sp

import derive_sigma_pipi_cartesian_loop_limits as limits
from sigma_pipi_loop_algebra import B0, F_PI, M_A
from sigma_pipi_cartesian_loop_algebra import ALP_MASS_COEFFICIENT, UP_QUARK_MASS
from test_derive_sigma_pipi_cartesian_loop_traces import FORBIDDEN_COMPONENTS, FORBIDDEN_LOOP_VALUES
from test_derive_sigma_pipi_cartesian_loop_large_epsilon_graphs import SHARE_POWERS
from test_derive_sigma_pipi_cartesian_loop_up_only_estimate import LEADING as UP_ONLY_LEADING

RECORD = (pathlib.Path(__file__).resolve().parent.parent / 'memory' / 'alp' / 'radial-mode' / 'finite-mixing'
          / 'two-pion-loop' / 'cartesian-estimate' / 'cartesian-sigma-pipi-loop-limits.md')
M_P = sp.Symbol('m_P', positive=True)


class VacuumConditionTest(unittest.TestCase):
    def test_modulus(self):
        f_a, s, p = limits.DECAY_CONSTANT, limits.RADIAL, limits.ALP_FIELD
        self.assertEqual(sp.expand(limits.modulus_expansion() - (f_a + s + p ** 2 / (2 * f_a))), 0)

    def test_shift(self):
        self.assertEqual(sp.expand(limits.vacuum_shift() - limits.ALP_FIELD_SQUARED / (2 * limits.DECAY_CONSTANT)), 0)

    def test_one_loop_relative_change(self):
        expected = -M_A ** 2 / (32 * sp.pi ** 2 * limits.DECAY_CONSTANT ** 2)
        self.assertEqual(sp.simplify(limits.relative_tree_change() - expected), 0)

    def test_all_orders_factor_is_not_the_one_loop_quantity(self):
        """Negative control: the factor with all orders in delta f_a saturates at -1 for large eps."""
        expression = limits.exact_tree_factor().subs(limits.LOOP_COUNTING, 1)
        expression = expression.subs({M_A: sp.sqrt(ALP_MASS_COEFFICIENT) / limits.U,
                                      limits.DECAY_CONSTANT: F_PI * limits.U})
        self.assertEqual(sp.limit(expression, limits.U, 0), -1)

    def test_large_epsilon(self):
        self.assertEqual(limits.leading_relative_change(),
                         (4, -ALP_MASS_COEFFICIENT / (32 * sp.pi ** 2 * F_PI ** 2)))
        self.assertEqual(limits.up_only_leading_relative_change(),
                         (4, -B0 * UP_QUARK_MASS / (32 * sp.pi ** 2 * F_PI ** 2)))

    def test_comparison_with_ratio(self):
        generic_power = max(SHARE_POWERS.values())
        self.assertEqual(limits.comparison_with_ratio(generic_power)[0], 0)
        power, coefficient = limits.up_only_leading_relative_change()
        self.assertEqual(power - 2, 2)
        self.assertEqual(sp.simplify(sp.Abs(coefficient) / sp.Abs(UP_ONLY_LEADING)
                                     - M_P ** 2 / (4 * B0 * UP_QUARK_MASS)), 0)

    def test_ratio_scale(self):
        self.assertTrue(limits.ratio_scale_identity())

    def test_up_only_d4_powers(self):
        self.assertEqual(limits.up_only_d4_powers(), {'C3': 0, 'C15': 2, 'Cpi3': 0})


class RecordTest(unittest.TestCase):
    def test_no_loop_value(self):
        text = RECORD.read_text()
        for token in FORBIDDEN_COMPONENTS + FORBIDDEN_LOOP_VALUES:
            self.assertNotIn(token, text, token)
        for token in ('A_0', 'B_0(', 'C_0', 'D_0(', 'Li_2', r'\mathrm{Li}', r'\zeta', 'Passarino', 'Veltman'):
            self.assertNotIn(token, text, token)

    def test_content(self):
        text = RECORD.read_text()
        self.assertIn('No loop integral', text)
        self.assertIn('-21/(64\\pi^2)', text)
        self.assertIn('authoritative', text)
        self.assertNotIn('action', text.lower().replace('interaction', ''))
        for statement in (r'\frac{\langle p^2\rangle}{2f_a}', r'-\frac{m_a^2}{32\pi^2f_a^2}',
                          r'-\frac{A}{32\pi^2f_\pi^2}\,\epsilon^4', r'\frac{\epsilon^2m_P^2}{4B_0m_u}'):
            self.assertIn(statement, text)


if __name__ == '__main__':
    unittest.main()
