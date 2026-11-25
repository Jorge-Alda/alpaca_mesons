"""Unit tests for derive_murayama_radial_crossing.py.

Run: PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_murayama_radial_crossing
"""
import unittest

import sympy as sp

from derive_murayama_radial_crossing import (
    ALP_ONLY_CLASSES, ALP_POLE_MASS_SQUARED, DAUGHTERS, EPSILON, LOOP_CLASSES, MAXIMUM_RATIO,
    PAIR_INVARIANT_12, PAIR_INVARIANT_23, PARENT_OVERLAP_TRACE, QUARTIC_COUPLING,
    RADIAL_ALP_SELF_ENERGY_CLASSES, RADIAL_COUPLING_SWITCH, RADIAL_MASS_SQUARED,
    REFERENCE_MASS_SQUARED, alp_only_check, alp_only_limit, alp_pole_normalization,
    crossing_equation, daughter_pole_mass_squared, exact_tree, frame_invariance_check,
    leading_numerator_truncation, leading_relative_power, loop_factor, loop_size_estimate,
    pair_invariants, quartic_coupling_at_fixed_radial_mass, radial_alp_self_energy_parts,
    radial_exchange, reduced_contraction, single_line_ratio_comparison)


class CrossingEquationTest(unittest.TestCase):
    def test_equation_is_ratio_equal_to_threshold(self):
        equation, ratio = crossing_equation()
        denominator = sp.Abs(exact_tree())
        numerator = loop_size_estimate()
        self.assertEqual(ratio * denominator - numerator, 0)
        self.assertEqual(equation - (numerator - MAXIMUM_RATIO * denominator), 0)
        self.assertIn(MAXIMUM_RATIO, equation.free_symbols)
        self.assertFalse(any(symbol.name in ('mu', r'\mu') for symbol in equation.free_symbols))

    def test_fixed_dalitz_invariants(self):
        invariants = pair_invariants()
        total = sum(invariants.values())
        expected = ALP_POLE_MASS_SQUARED + sum(daughter_pole_mass_squared(n) for n in DAUGHTERS)
        self.assertEqual(sp.simplify(total - expected), 0)
        self.assertEqual(invariants[(1, 2)], PAIR_INVARIANT_12)
        self.assertEqual(invariants[(2, 3)], PAIR_INVARIANT_23)


class NumeratorScopeTest(unittest.TestCase):
    def test_meson_only_loops_are_excluded(self):
        estimate = loop_size_estimate()
        self.assertNotIn('1', LOOP_CLASSES)
        self.assertFalse(estimate.has(sp.Function('reduced_contraction_class_1', real=True)(EPSILON)))
        self.assertFalse(any(function.func.__name__ == 'reduced_contraction_class_1'
                             for function in estimate.atoms(sp.Function)))

    def test_every_class_with_an_alp_like_or_radial_line_enters(self):
        estimate = loop_size_estimate()
        for loop_class in LOOP_CLASSES:
            self.assertTrue(estimate.has(reduced_contraction(loop_class)), loop_class)

    def test_numerator_is_not_truncated(self):
        estimate = alp_only_limit(loop_size_estimate())
        self.assertTrue(estimate.has(reduced_contraction('4')))
        self.assertTrue(estimate.has(PARENT_OVERLAP_TRACE))
        self.assertNotEqual(sp.simplify(estimate - leading_numerator_truncation(estimate)), 0)


class AlpOnlyLimitTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        (cls.truncation_residual, cls.full_difference, cls.split_residual,
         cls.terms) = alp_only_check()

    def test_record_equation_after_its_truncation(self):
        self.assertEqual(self.truncation_residual, 0)

    def test_difference_is_exactly_the_listed_terms(self):
        self.assertNotEqual(self.full_difference, 0)
        self.assertFalse(self.full_difference.has(MAXIMUM_RATIO))
        self.assertEqual(self.split_residual, 0)

    def test_each_listed_term_is_removed_by_the_truncation(self):
        for loop, origin, term in self.terms:
            self.assertEqual(sp.simplify(leading_numerator_truncation(term)), 0, (loop, origin))
            self.assertEqual(leading_relative_power(term), 4, (loop, origin))

    def test_listed_terms(self):
        loops = {loop for loop, _, _ in self.terms}
        self.assertIn('loop with two internal ALP-like propagators', loops)
        two_line_terms = [origin for loop, origin, _ in self.terms
                          if loop == 'loop with two internal ALP-like propagators']
        self.assertEqual(two_line_terms, ['whole contribution, absent from the record'])
        # Two terms (normalization, epsilon dependence of the contraction) for each
        # of the three moments present in the record, one for the two-line loop.
        self.assertEqual(len(self.terms), 7)

    def test_single_internal_line_matches_criterion_script(self):
        self.assertEqual(single_line_ratio_comparison(), 0)

    def test_truncation_agrees_with_limit_for_explicit_contractions(self):
        explicit = {
            '2': 1 + EPSILON, '6': -3 + EPSILON**2, '3': 2 * sp.cos(EPSILON), '4': 5 + EPSILON,
            '5': -1 + EPSILON, 'R2': sp.Rational(1, 2) - EPSILON, 'R3': 7 * sp.exp(EPSILON),
            'R4a': 3 + EPSILON, 'R4b': -2 + EPSILON, 'R5a': 4 - EPSILON**2,
            'R5b': sp.Rational(-5, 2) + EPSILON,
        }
        replacements = {reduced_contraction(name): value for name, value in explicit.items()}
        estimate = loop_size_estimate().subs(replacements).subs(RADIAL_COUPLING_SWITCH, 1)
        truncated = leading_numerator_truncation(estimate)
        limit = sp.limit(estimate / EPSILON**3, EPSILON, 0, '+')
        self.assertEqual(sp.simplify(truncated / EPSILON**3 - limit), 0)


class RadialSwitchTest(unittest.TestCase):
    def test_switch_removes_radial_loops_and_exchange(self):
        estimate = loop_size_estimate()
        self.assertTrue(estimate.has(RADIAL_COUPLING_SWITCH))
        limit_estimate = alp_only_limit(estimate)
        for name in LOOP_CLASSES:
            self.assertEqual(limit_estimate.has(reduced_contraction(name)), name in ALP_ONLY_CLASSES)
        self.assertEqual(alp_only_limit(radial_exchange()), 0)
        self.assertTrue(exact_tree().has(RADIAL_MASS_SQUARED))
        self.assertFalse(alp_only_limit(exact_tree()).has(RADIAL_MASS_SQUARED))

    def test_switch_multiplies_each_radial_coupling(self):
        exchange = radial_exchange()
        self.assertFalse(sp.expand(exchange / RADIAL_COUPLING_SWITCH**2).has(RADIAL_COUPLING_SWITCH))

    def test_large_radial_mass_is_not_the_alp_only_limit(self):
        self.assertEqual(sp.limit(radial_exchange(), RADIAL_MASS_SQUARED, sp.oo), 0)
        self.assertEqual(sp.limit(loop_factor('radial_tadpole'), RADIAL_MASS_SQUARED, sp.oo), sp.oo)

    def test_quartic_coupling_at_fixed_radial_mass(self):
        relation, quartic = quartic_coupling_at_fixed_radial_mass()
        self.assertEqual(sp.simplify(relation.rhs.subs(QUARTIC_COUPLING, quartic) - relation.lhs), 0)


class RadialAlpSelfEnergyTest(unittest.TestCase):
    """Closed sigma line on the derivative sigma-sigma-a-a vertex."""

    def test_power_scale_and_normalization(self):
        for name in RADIAL_ALP_SELF_ENERGY_CLASSES:
            moment, power, alp_lines, radial, _ = LOOP_CLASSES[name]
            self.assertEqual(power, 2, name)
            self.assertEqual(alp_lines, 0, name)
            self.assertTrue(radial, name)
            self.assertEqual(REFERENCE_MASS_SQUARED[moment], RADIAL_MASS_SQUARED, name)

    def test_combined_with_radial_tadpole_in_one_magnitude(self):
        estimate = loop_size_estimate()
        magnitudes = [term for term in estimate.atoms(sp.Abs) if term.has(reduced_contraction('R2'))]
        self.assertEqual(len(magnitudes), 1)
        for name in RADIAL_ALP_SELF_ENERGY_CLASSES:
            self.assertTrue(magnitudes[0].has(reduced_contraction(name)), name)
        # Weight of the shared magnitude: epsilon^3 m_sigma^2/(16 pi^2 f_pi^2) / sqrt(N), times t_sigma.
        placeholder = sp.Symbol('shared_magnitude')
        coefficient = sp.diff(estimate.xreplace({magnitudes[0]: placeholder}), placeholder)
        expected = (RADIAL_COUPLING_SWITCH * EPSILON**3 * loop_factor('radial_tadpole')
                    / sp.sqrt(alp_pole_normalization()))
        self.assertEqual(sp.simplify(coefficient - expected), 0)

    def test_removed_in_the_alp_only_limit(self):
        limit_estimate = alp_only_limit(loop_size_estimate())
        for name in RADIAL_ALP_SELF_ENERGY_CLASSES:
            self.assertFalse(limit_estimate.has(reduced_contraction(name)), name)

    def test_mass_shift_absorbed_and_only_residue_enters(self):
        momentum_squared = sp.Symbol('p^2', positive=True)
        component, pole = sp.symbols('v_a^2 m^2', positive=True)
        parts = radial_alp_self_energy_parts(momentum_squared, component, pole)
        self.assertEqual(parts['self-energy'], component * momentum_squared)
        self.assertEqual(parts['mass shift (absorbed into m_a)'], component * pole)
        self.assertEqual(parts['residue part'], component)
        self.assertEqual(parts['remainder beyond the residue'], 0)
        # Internal line: (Sigma(q^2) - Sigma(m^2))/(q^2 - m^2) is the residue part.
        ratio = (parts['self-energy'] - parts['mass shift (absorbed into m_a)']) / (momentum_squared - pole)
        self.assertEqual(sp.simplify(ratio - parts['residue part']), 0)


class FrameInvarianceTest(unittest.TestCase):
    def test_ratio_unchanged(self):
        dressed_difference, dressed_zero, theta = frame_invariance_check()
        self.assertEqual(theta['original'], 0)
        self.assertEqual(theta['dressed_after'], 0)
        self.assertEqual(dressed_difference, 0)
        self.assertTrue(dressed_zero)


if __name__ == '__main__':
    unittest.main()
