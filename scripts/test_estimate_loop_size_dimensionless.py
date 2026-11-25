"""Unit tests for estimate_loop_size_dimensionless.py.

No loop is computed: the tests check units, frame independence, the tree and the conversion of
the vertex products from the f_pi = 1 evaluator.

Run: PYTHONPATH=scripts .venv/bin/python -m unittest test_estimate_loop_size_dimensionless
"""
import unittest

import numpy as np

import estimate_loop_size_dimensionless as dimensionless
import evaluate_explicit_estimate as explicit
from estimate_loop_size_dimensionless import axially_redefined_frame

KAPPA = 0.3


def ratios_by_term(estimate):
    terms, tree = estimate.diagram_terms()
    return {key: value / tree for key, value in terms.items()}, sum(terms.values()) / tree


class UnitInvarianceTest(unittest.TestCase):
    """Every dimensionful input times c (c = 2, 0.5), kappa_u and lambda fixed: R unchanged."""

    @classmethod
    def setUpClass(cls):
        cls.reference_terms, cls.reference_total = ratios_by_term(dimensionless.DimensionfulEstimate(KAPPA))

    def test_ratio_and_every_term_unchanged(self):
        for factor in (2., 0.5):
            inputs = dimensionless.DEFAULT_INPUTS.rescaled(factor)
            estimate = dimensionless.DimensionfulEstimate(KAPPA, inputs)
            reference = dimensionless.DimensionfulEstimate(KAPPA)
            with self.subTest(factor=factor, quantity='epsilon'):
                self.assertAlmostEqual(estimate.epsilon / reference.epsilon, 1., places=14)
            with self.subTest(factor=factor, quantity='m_a / f_pi'):
                self.assertAlmostEqual(estimate.alp_mass_squared / reference.alp_mass_squared, factor ** 2, places=11)
            terms, total = ratios_by_term(estimate)
            with self.subTest(factor=factor, quantity='R'):
                self.assertLess(abs(abs(total) / abs(self.reference_total) - 1), 1e-9)
            for key, value in terms.items():
                with self.subTest(factor=factor, term=key):
                    reference_value = self.reference_terms[key]
                    self.assertLess(abs(value - reference_value), 1e-9 * max(abs(reference_value), 1e-12))


class GenericUnitInvarianceTest(unittest.TestCase):
    """Generic factors c = 3, 0.37, 7.1 (not powers of two, so rounding is not exact).

    Every dimensionful input (f_pi, m_u, m_d, m_s, B0 times c; m_1^2 times c^2) is rescaled
    through PhysicalInputs.rescaled; kappa_u and the dimensionless lambda are fixed.  Derived
    scales must follow: m_a^2, m_pi0^2, m_pi+^2, m_sigma^2, s_uv times c^2.  R and every
    (class, diagram labels) term divided by M0 must be unchanged to 1e-9 relative.
    """

    FACTORS = (3., 0.37, 7.1)
    KAPPAS = (0.3, 1.0)

    @classmethod
    def setUpClass(cls):
        cls.references = {}
        for kappa_u in cls.KAPPAS:
            estimate = dimensionless.DimensionfulEstimate(kappa_u)
            cls.references[kappa_u] = (estimate, *ratios_by_term(estimate))

    def test_ratio_and_every_term_unchanged(self):
        for kappa_u in self.KAPPAS:
            reference, reference_terms, reference_total = self.references[kappa_u]
            for factor in self.FACTORS:
                estimate = dimensionless.DimensionfulEstimate(kappa_u, dimensionless.DEFAULT_INPUTS.rescaled(factor))
                with self.subTest(kappa_u=kappa_u, factor=factor, quantity='epsilon'):
                    self.assertLess(abs(estimate.epsilon / reference.epsilon - 1), 1e-12)
                scales = {
                    'm_a^2': (estimate.alp_mass_squared, reference.alp_mass_squared),
                    'm_pi0^2': (estimate.neutral_pion_mass_squared, reference.neutral_pion_mass_squared),
                    'm_pi+^2': (estimate.charged_pion_mass_squared, reference.charged_pion_mass_squared),
                    'm_sigma^2': (estimate.radial_mass_squared, reference.radial_mass_squared),
                    's_uv': (estimate.pair_invariant, reference.pair_invariant),
                }
                for name, (value, reference_value) in scales.items():
                    with self.subTest(kappa_u=kappa_u, factor=factor, quantity=name):
                        self.assertLess(abs(value / (factor ** 2 * reference_value) - 1), 1e-10)
                terms, total = ratios_by_term(estimate)
                with self.subTest(kappa_u=kappa_u, factor=factor, quantity='R'):
                    self.assertLess(abs(abs(total) / abs(reference_total) - 1), 1e-9)
                self.assertEqual(set(terms), set(reference_terms))
                for key, value in terms.items():
                    with self.subTest(kappa_u=kappa_u, factor=factor, term=key):
                        reference_value = reference_terms[key]
                        self.assertLess(abs(value - reference_value), 1e-9 * max(abs(reference_value), 1e-12))


class FrameInvarianceTest(unittest.TestCase):
    """Constant axial redefinition U = xi U' xi, kappa_ax = diag(0.2, -0.2, 0)."""

    def test_estimate_frame_independent(self):
        inputs = dimensionless.DEFAULT_INPUTS
        base = dimensionless.polar_frame(inputs)
        rotated = axially_redefined_frame(base, np.diag([0.2, -0.2, 0.]))
        self.assertGreater(max(np.max(np.abs(base.sources[k] - rotated.sources[k])) for k in base.sources), 1e-4)
        first = dimensionless.DimensionfulEstimate(KAPPA, inputs, frame=base)
        second = dimensionless.DimensionfulEstimate(KAPPA, inputs, frame=rotated)
        self.assertAlmostEqual(second.alp_mass_squared / first.alp_mass_squared, 1., places=12)
        before, total_before = ratios_by_term(first)
        after, total_after = ratios_by_term(second)
        self.assertLess(abs(total_after - total_before), 1e-9 * abs(total_before))
        for key, value in before.items():
            with self.subTest(term=key):
                self.assertLess(abs(after[key] - value), 1e-9 * max(abs(value), 1e-12))


class AgreementWithExplicitEvaluatorTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.new = dimensionless.DimensionfulEstimate(KAPPA)
        cls.old = explicit.ExplicitEstimate(KAPPA)

    def test_tree(self):
        for kappa_u in (KAPPA, 1.0):
            with self.subTest(kappa_u=kappa_u):
                new = dimensionless.DimensionfulEstimate(kappa_u).tree()[0]
                old = explicit.ExplicitEstimate(kappa_u).tree()[0]
                self.assertLess(abs(new - old) / abs(old), 1e-10)

    def test_cycles_are_unit_conversions(self):
        """A_d (GeV) = f_pi^(2r - 4) A_d (f_pi = 1): the physical vertex product has dimension 2r - 4."""
        unit = explicit.PION_DECAY_CONSTANT
        for label, _, blocks, lines in explicit.DIAGRAMS:
            new = self.new.cycle(blocks, lines, self.new.loop_scale_squared(lines))
            old = self.old.cycle(blocks, lines, self.old.loop_scale_squared(lines))
            with self.subTest(label=label, blocks=blocks, lines=lines):
                self.assertLess(abs(new - unit ** (2 * len(lines) - 4) * old), 1e-9 * max(abs(new), 1e-30))

    def test_two_propagator_terms_from_units_of_f_pi(self):
        """Class 2 has r = 2 only, so (m_a^2)^0: C_2/M0 = sum_d S_d A_d/(16 pi^2 M0), with A_d and M0
        taken from the f_pi = 1 evaluator (both dimensionless there and here)."""
        new_classes, new_tree = self.new.contributions()
        old_tree = self.old.tree()[0]
        expected = sum(weight * dimensionless.LOOP_FACTOR
                       * self.old.cycle(blocks, lines, self.old.loop_scale_squared(lines))
                       for label, weight, blocks, lines in explicit.DIAGRAMS
                       if explicit.diagram_class(label) == '2') / old_tree
        self.assertLess(abs(new_classes['2'] / new_tree - expected), 1e-9 * abs(expected))

    def test_one_propagator_terms_from_units_of_f_pi(self):
        """Classes 1, 7 and 8 (r = 1): C/M0 = (m^2/f_pi^2)/(16 pi^2) x (cycles of class 1, or residues
        of class 7 or 8)/M0, with the cycles, residues, scale m^2/f_pi^2 and M0 taken from the
        f_pi = 1 evaluator."""
        new_classes, new_tree = self.new.contributions()
        old_tree = self.old.tree()[0]
        cycles = sum(weight * self.old.cycle(blocks, lines, self.old.loop_scale_squared(lines))
                     for label, weight, blocks, lines in explicit.DIAGRAMS if explicit.diagram_class(label) == '1')
        one_propagator = {'1': (('A',), cycles)}
        one_propagator.update(self.old.external_residues(old_tree))
        for name, (lines, value) in one_propagator.items():
            expected = dimensionless.LOOP_FACTOR * self.old.loop_scale_squared(lines) * value / old_tree
            with self.subTest(class_name=name):
                self.assertLess(abs(new_classes[name] / new_tree - expected), 1e-9 * abs(expected))

    def test_radial_slope_from_vertex(self):
        """eps^2 xi^2/f_pi^2 equals the p^2 slope of (1/2) V(l_e(p), l_e(-p), l_sigma, l_sigma)."""
        estimate = self.new
        sigma = explicit.Leg(np.zeros((3, 3), complex), 0j, 1. + 0j, explicit.LOOP)
        for name in ('A', '0'):
            leg = estimate.legs[name]
            legs = [explicit.Leg(leg.meson, leg.alp, 0j, 'p'), explicit.Leg(leg.meson, leg.alp, 0j, 'q'), sigma, sigma]
            values = []
            for momentum_squared in (0., 1.):
                extra = {('p', 'q'): -momentum_squared, ('q', 'p'): -momentum_squared}
                values.append(0.5 * estimate.vertex(legs, estimate.dot_rule(estimate.radial_mass_squared, extra)))
            with self.subTest(name=name):
                closed = estimate.neutral_leg_radial_slope(name)
                self.assertLess(abs(values[1] - values[0] - closed), 1e-9 * abs(closed))


class DecisionsTest(unittest.TestCase):
    """Class-7 residues on the a, pi+, pi-, pi0 legs; sigma^3 diagrams in physical units."""

    @classmethod
    def setUpClass(cls):
        cls.estimate = dimensionless.DimensionfulEstimate(KAPPA)
        cls.terms, cls.tree = cls.estimate.diagram_terms()

    def test_charged_pion_slope_closed_form(self):
        """pi+ slope d/dp^2 [(1/2) V(l_+(p), l_-(-p), l_A, l_A)] in GeV^-2 from the physical vertex
        equals -eps^2 [(M^-1)_11 - (M^-1)_22]^2 / (3 N(m_a^2) f_pi^2), M^-1 = M^{-1}_{m_a^2}(C)."""
        for kappa_u in (0.3, 1.0):
            e = dimensionless.DimensionfulEstimate(kappa_u)
            inverse = e.pole_operator_inverse(e.alp_mass_squared)
            closed = (-(e.epsilon ** 2) * (inverse[0, 0] - inverse[1, 1]) ** 2
                      / (3 * e.normalization(e.alp_mass_squared) * e.pion_decay_constant ** 2))
            for first, second in (('+', '-'), ('-', '+')):
                with self.subTest(kappa_u=kappa_u, legs=(first, second)):
                    self.assertLess(abs(e.closed_alp_line_slope(first, second) - closed), 1e-9 * abs(closed))

    def test_class_seven_residue_includes_charged_pions(self):
        """Residues 7.1-7.4 = m_a^2/(16 pi^2) x (1/2)(slope_a + slope_+ + slope_- + slope_0) M0, and the
        a and pi0 slopes vanish, so they sum to m_a^2/(16 pi^2) x slope_+ x M0.  Class 1 (diagrams
        1.1-1.3) holds the cycles with a closed A line."""
        e = self.estimate
        groups, tree = e.group_terms()
        cycles = sum(w * dimensionless.LOOP_FACTOR * e.alp_mass_squared * e.cycle(b, l, e.alp_mass_squared)
                     for label, w, b, l in explicit.DIAGRAMS if explicit.diagram_class(label) == '1')
        self.assertLess(abs(groups[('1', '1.1-1.3', 1, 'cycle')] - cycles), 1e-12 * abs(cycles))
        expected = dimensionless.LOOP_FACTOR * e.alp_mass_squared * e.closed_alp_line_slope('+', '-') * tree
        self.assertLess(abs(groups[('7', '7.1-7.4', 1, 'residue')] - expected), 1e-12 * abs(expected))
        for name in ('A', '0'):
            with self.subTest(leg=name):
                self.assertLess(abs(e.closed_alp_line_slope(name, name)),
                                1e-10 * abs(e.closed_alp_line_slope('+', '-')))
        self.assertLess(abs(self.terms[('1', '1.1-1.3')] - cycles), 1e-12 * abs(cycles))
        self.assertLess(abs(self.terms[('7', '7.1-7.4')] - expected), 1e-12 * abs(expected))

    def test_table_labels(self):
        """Printed table: every class 1-8 and its diagram labels, 44 cycles and 6 residues."""
        groups, _ = self.estimate.group_terms()
        self.assertEqual([labels for _, labels, _, _ in groups], [
            '1.1-1.3', '2.1-2.5', '3.1-3.3', '4.1-4.6', '4.7-4.14', '5.1-5.4', '5.5-5.12',
            '6.1', '6.2-6.4', '6.5-6.6', '6.7', '7.1-7.4', '8.1-8.2'])
        self.assertEqual(len(explicit.DIAGRAMS), 44)

    def test_table_entries_sum_to_terms(self):
        groups, _ = self.estimate.group_terms()
        self.assertAlmostEqual(abs(sum(groups.values()) - sum(self.terms.values())), 0., places=14)

    def test_sigma_cubed_vertex_in_physical_units(self):
        e = self.estimate
        sigma = explicit.Leg(np.zeros((3, 3), complex), 0j, 1. + 0j, explicit.LOOP)
        value = e.vertex([sigma, sigma, sigma], e.dot_rule(e.radial_mass_squared))
        self.assertLess(abs(value + 6 * e.quartic_coupling * e.alp_decay_constant), 1e-12 * abs(value))

    def test_radial_self_energy_closed_form(self):
        """Diagram 6.7 = (1/2)/(16 pi^2) rho_beta rho_gamma V_sigma3^2 (r = 2, no power of m_sigma)."""
        e = self.estimate
        cubic = -6 * e.quartic_coupling * e.alp_decay_constant
        expected = 0.5 * dimensionless.LOOP_FACTOR * e.legs['beta'].radial * e.legs['gamma'].radial * cubic ** 2
        self.assertLess(abs(self.terms[('6', '6.7')] - expected), 1e-9 * abs(expected))

    def test_radial_two_line_bubble_closed_form(self):
        """Diagram 6.1 = (1/2)/(16 pi^2) rho_beta V_sigma3 V(A, pi0, sigma, sigma),
        V(A, pi0, sigma, sigma) = -2 eps^2 xi_A xi_0 k_A.k_0 / f_pi^2."""
        e = self.estimate
        cubic = -6 * e.quartic_coupling * e.alp_decay_constant
        quartic = (-2 * e.epsilon ** 2 * e.legs['A'].alp * e.legs['0'].alp * e.external_dots[('A', '0')]
                   / e.pion_decay_constant ** 2)
        expected = 0.5 * dimensionless.LOOP_FACTOR * e.legs['beta'].radial * cubic * quartic
        self.assertLess(abs(self.terms[('6', '6.1')] - expected), 1e-9 * abs(expected))


if __name__ == '__main__':
    unittest.main()
