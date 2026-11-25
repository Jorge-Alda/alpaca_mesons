"""Unit tests for derive_murayama_radial_cartesian_weights.py (task 6.1).

Run: PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_murayama_radial_cartesian_weights
"""
import itertools
import unittest

import sympy as sp

import derive_murayama_radial_cartesian_weights as weights
from derive_murayama_radial_crossing import B_ZERO, EPSILON, QUARTIC_COUPLING
from derive_murayama_radial_up_only_specialization import UP_QUARK_MASS

RADIAL_MASS_SQUARED = weights.RADIAL_MASS_SQUARED_UP_ONLY
ANGULAR_SCALE = B_ZERO * UP_QUARK_MASS * EPSILON**2   # B_0 m_u epsilon^2


class FieldMapTest(unittest.TestCase):
    def test_lagrangians_agree_under_the_field_map(self):
        potential_residual, kinetic_residual = weights.field_map_residuals()
        self.assertEqual(potential_residual, 0)
        self.assertEqual(kinetic_residual, 0)

    def test_same_quadratic_part_hence_same_projectors(self):
        hessian_difference, tadpoles = weights.quadratic_part_residual()
        self.assertTrue(hessian_difference.is_zero_matrix)
        self.assertEqual(tadpoles, [0, 0])

    def test_radial_self_couplings_and_odd_vertices(self):
        self.assertTrue(all(value == 0 for value in weights.radial_self_coupling_differences().values()))
        self.assertTrue(all(value == 0 for value in weights.odd_vertices().values()))


class VertexTest(unittest.TestCase):
    def test_radial_vertex_differs_off_shell_only(self):
        momentum_squared = sp.Symbol('q^2', real=True)
        polar, cartesian = weights.radial_two_pole_vertices(momentum_squared)
        factor = weights.common_factor((0, 1), 1)
        self.assertEqual(sp.simplify((polar - cartesian) / factor
                                     - (RADIAL_MASS_SQUARED - momentum_squared)), 0)
        self.assertEqual(sp.simplify((polar - cartesian).subs(momentum_squared, RADIAL_MASS_SQUARED)), 0)

    def test_radial_vertex_matches_closed_form_of_record(self):
        # -g_xy(q^2) of the only-up record (the record writes L = -(1/2) sigma g P P).
        momentum_squared = sp.Symbol('q^2', real=True)
        polar, _ = weights.radial_two_pole_vertices(momentum_squared)
        m_1, m_2 = weights.POLE_MASSES[:2]
        expected = -(momentum_squared - m_1 - m_2 + m_1 * m_2 / ANGULAR_SCALE) * weights.common_factor((0, 1), 1)
        self.assertEqual(sp.simplify(polar - expected), 0)

    def test_two_radial_vertex(self):
        momentum_product = sp.Symbol('p_1.p_2', real=True)
        polar, cartesian = weights.two_radial_two_pole_vertices(momentum_product)
        factor = weights.common_factor((0, 1), 2)
        self.assertEqual(sp.simplify(polar / factor + 2 * momentum_product), 0)
        self.assertEqual(sp.simplify(cartesian / factor + (RADIAL_MASS_SQUARED - ANGULAR_SCALE)), 0)

    def test_four_pole_contact_difference(self):
        masses = weights.POLE_MASSES
        pair_sum = sum(first * second for first, second in itertools.combinations(masses, 2))
        expected = (pair_sum - 2 * ANGULAR_SCALE * sum(masses) + 3 * ANGULAR_SCALE * RADIAL_MASS_SQUARED) / ANGULAR_SCALE
        difference = weights.four_pole_vertex_difference() / weights.common_factor((0, 1, 2, 3), 2)
        self.assertEqual(sp.simplify(difference - expected), 0)
        self.assertTrue(difference.has(QUARTIC_COUPLING))


class WeightTest(unittest.TestCase):
    def test_common_factor_carries_every_normalization(self):
        for name, (factor, polar, cartesian, difference) in weights.weight_comparison().items():
            for bracket in (polar, cartesian, difference):
                if bracket is not None:
                    self.assertFalse(bracket.free_symbols & set(weights.NORMALIZATIONS), name)
            self.assertNotEqual(difference, 0, name)

    def test_no_loop_quantities(self):
        rows = weights.weight_comparison()
        names = {str(symbol) for row in rows.values() for part in row if part is not None
                 for symbol in part.free_symbols}
        self.assertFalse(any('mu' == name or 'loop' in name for name in names))


class OnShellControlTest(unittest.TestCase):
    def test_four_pole_amplitude_is_field_independent(self):
        self.assertEqual(weights.four_pole_amplitude_difference(), 0)


if __name__ == '__main__':
    unittest.main()
