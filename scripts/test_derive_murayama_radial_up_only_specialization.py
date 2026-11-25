"""Unit tests for derive_murayama_radial_up_only_specialization.py (tasks 5.1, 5.2).

Run: PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_murayama_radial_up_only_specialization
"""
import unittest

import sympy as sp

import derive_murayama_radial_up_only_specialization as up_only
from derive_murayama_radial_crossing import (
    EPSILON, F_PI, B_ZERO, MAXIMUM_RATIO, QUARTIC_COUPLING, RADIAL_MASS_SQUARED,
    loop_size_estimate, reduced_contraction, LOOP_CLASSES)


class ModelInputTest(unittest.TestCase):
    def test_epsilon_from_kappa(self):
        self.assertEqual(sp.simplify(up_only.alp_decay_constant() * up_only.KAPPA_U / sp.sqrt(2)
                                     - up_only.UP_QUARK_MASS), 0)
        self.assertEqual(sp.simplify(up_only.up_only_epsilon()
                                     - F_PI * up_only.KAPPA_U / (sp.sqrt(2) * up_only.UP_QUARK_MASS)), 0)

    def test_source_coefficients_of_the_radial_record(self):
        m_u = up_only.UP_QUARK_MASS
        coefficients = up_only.up_only_source_coefficients()
        expected = {(1, 0): (m_u, 0), (0, 1): (0, -m_u), (1, 1): (0, -m_u),
                    (0, 2): (-m_u / 2, 0), (1, 2): (-m_u / 2, 0),
                    (2, 0): (0, 0), (2, 1): (0, 0), (2, 2): (0, 0)}
        for key, value in expected.items():
            self.assertEqual(coefficients[key], value, key)


class PoleDataTest(unittest.TestCase):
    def test_normalization_is_derivative_of_secular_function(self):
        self.assertEqual(up_only.normalization_is_secular_derivative(), 0)

    def test_up_entry_fixed_by_pole_equation(self):
        self.assertEqual(up_only.up_entry_pole_identity(), 0)


class TreeSpecializationTest(unittest.TestCase):
    def test_radial_couplings_reduce_to_closed_form(self):
        for pair, residual in up_only.radial_coupling_residuals().items():
            self.assertEqual(residual, 0, pair)

    def test_general_tree_reduces_to_up_only_tree(self):
        self.assertEqual(up_only.tree_substitution_residual(), 0)

    def test_on_shell_closed_form_coupling_record_value(self):
        # q^2 = m_sigma^2 with x = y: g_xx = [m_sigma^2 - 2 m^2 + m^4/(B_0 m_u eps^2)]/(f_a N)
        coupling = up_only.closed_form_radial_coupling(1, 1, RADIAL_MASS_SQUARED)
        mass = up_only.general.leg_pole_mass_squared(1)
        expected = (sp.Symbol('t_sigma', nonnegative=True)
                    * (RADIAL_MASS_SQUARED - 2 * mass + mass**2 / (B_ZERO * up_only.UP_QUARK_MASS * EPSILON**2))
                    * EPSILON / F_PI / up_only.leg_normalization_up_only(1))
        self.assertEqual(sp.simplify(coupling - expected), 0)


class CrossingSpecializationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.equation, cls.ratio = up_only.up_only_crossing()

    def test_depends_on_kappa_and_lambda_symbolically(self):
        symbols = self.equation.free_symbols
        self.assertIn(up_only.KAPPA_U, symbols)
        self.assertIn(QUARTIC_COUPLING, symbols)
        self.assertIn(MAXIMUM_RATIO, symbols)
        self.assertNotIn(EPSILON, symbols)
        self.assertNotIn(RADIAL_MASS_SQUARED, symbols)

    def test_contractions_remain_undefined(self):
        for loop_class in LOOP_CLASSES:
            if loop_class == 'R5b':
                continue
            contraction = reduced_contraction(loop_class).subs(EPSILON, up_only.up_only_epsilon())
            self.assertTrue(self.equation.has(contraction), loop_class)

    def test_radial_self_energy_of_alp_leg(self):
        external = reduced_contraction('R5a').subs(EPSILON, up_only.up_only_epsilon())
        internal = reduced_contraction('R5b').subs(EPSILON, up_only.up_only_epsilon())
        tadpole = reduced_contraction('R2').subs(EPSILON, up_only.up_only_epsilon())
        self.assertFalse(self.equation.has(internal))
        magnitudes = [term for term in self.equation.atoms(sp.Abs) if term.has(external)]
        self.assertEqual(len(magnitudes), 1)
        self.assertTrue(magnitudes[0].has(tadpole))

    def test_radial_self_energy_parts(self):
        parts = up_only.up_only_radial_alp_self_energy_parts()
        normalization = up_only.SECULAR_DERIVATIVE(up_only.general.ALP_POLE_MASS_SQUARED)
        self.assertEqual(parts['residue part'], 1 / normalization)
        self.assertEqual(parts['mass shift (absorbed into m_a)'],
                         up_only.general.ALP_POLE_MASS_SQUARED / normalization)
        self.assertEqual(parts['remainder beyond the residue'], 0)

    def test_estimate_has_secular_derivative_normalization(self):
        estimate = up_only.specialize(loop_size_estimate())
        normalization = up_only.SECULAR_DERIVATIVE(up_only.general.ALP_POLE_MASS_SQUARED)
        self.assertTrue(estimate.has(normalization))
        self.assertFalse(estimate.has(up_only.PARENT_OVERLAP_TRACE))

    def test_compact_form_agrees(self):
        compact = up_only.up_only_crossing_compact().subs(EPSILON, up_only.up_only_epsilon())
        point = {symbol: sp.Rational(3 + index, 5 + 2 * index)
                 for index, symbol in enumerate(sorted(compact.free_symbols | self.equation.free_symbols,
                                                       key=lambda symbol: symbol.name))}
        functions = {function: sp.Rational(7 + index, 3 + index) for index, function in enumerate(
            sorted(compact.atoms(sp.Function) | self.equation.atoms(sp.Function), key=str))
            if not isinstance(function, sp.Abs)}
        difference = (compact - self.equation).xreplace(functions).xreplace(point)
        self.assertLess(abs(sp.N(difference, 40)), sp.Float('1e-30'))


class ContractionEntryTest(unittest.TestCase):
    def test_invariant_entries_in_the_up_only_model(self):
        _, identities = up_only.up_only_contraction_entries()
        for name, residual in identities.items():
            self.assertEqual(residual, 0, name)


class RadialMassTest(unittest.TestCase):
    def test_closed_form(self):
        m_u, kappa = up_only.UP_QUARK_MASS, up_only.KAPPA_U
        expected = B_ZERO * F_PI**2 * kappa**2 / (2 * m_u) + 4 * QUARTIC_COUPLING * m_u**2 / kappa**2
        self.assertEqual(sp.simplify(up_only.radial_mass_squared() - expected), 0)

    def test_ratio_of_lambda_term(self):
        first, lambda_term, ratio = up_only.radial_mass_terms()
        self.assertEqual(sp.simplify(ratio - up_only.expected_radial_mass_ratio()), 0)
        self.assertFalse(first.has(QUARTIC_COUPLING))

    def test_first_term_is_bare_angular_curvature(self):
        # B_0 m_u epsilon^2 at epsilon(kappa_u): the same scale as in the pole equation.
        first, _, _ = up_only.radial_mass_terms()
        self.assertEqual(sp.simplify(first - B_ZERO * up_only.UP_QUARK_MASS
                                     * up_only.up_only_epsilon()**2), 0)


class AlpPoleTest(unittest.TestCase):
    def test_secular_polynomial_degree_and_lambda_independence(self):
        polynomial = up_only.secular_polynomial()
        self.assertEqual(polynomial.degree(), 4)
        self.assertFalse(polynomial.as_expr().has(QUARTIC_COUPLING))
        self.assertTrue(polynomial.as_expr().has(up_only.KAPPA_U))

    def test_secular_polynomial_irreducible(self):
        self.assertTrue(up_only.secular_polynomial_is_irreducible())

    def test_roots_are_zeros_of_secular_function(self):
        # Rational check point: roots of the polynomial are zeros of F(s).
        point = {B_ZERO: sp.Rational(28, 10), F_PI: sp.Rational(92, 1000), up_only.KAPPA_U: sp.Rational(3, 10),
                 up_only.UP_QUARK_MASS: sp.Rational(22, 10000), up_only.DOWN_QUARK_MASS: sp.Rational(47, 10000),
                 up_only.STRANGE_QUARK_MASS: sp.Rational(93, 1000), up_only.SINGLET_MASS_SQUARED: sp.Rational(73, 100)}
        polynomial = sp.Poly(up_only.secular_polynomial().as_expr().subs(point), up_only.SQUARED_MOMENTUM)
        secular = up_only.secular_function(epsilon=up_only.up_only_epsilon()).subs(point)
        for root in polynomial.nroots(n=40):
            self.assertLess(abs(sp.N(secular.subs(up_only.SQUARED_MOMENTUM, root), 30)), 1e-20)

    def test_small_kappa_limit(self):
        limit = up_only.alp_like_root_small_kappa_limit()
        m_u = up_only.UP_QUARK_MASS
        inverse_up = up_only.pole_inverse({'u': 1, 'd': 0, 's': 0}, 0)['u']
        expected = EPSILON**2 * B_ZERO * m_u * (1 - 2 * B_ZERO * m_u * inverse_up)
        self.assertEqual(sp.simplify(limit - expected), 0)


if __name__ == '__main__':
    unittest.main()
