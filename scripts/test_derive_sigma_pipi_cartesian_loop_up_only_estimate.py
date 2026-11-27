"""Unit tests for derive_sigma_pipi_cartesian_loop_up_only_estimate.py (exact symbolic identities; no numerical matrices).

Run: PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_cartesian_loop_up_only_estimate
(about six minutes).
"""
import pathlib
import unittest

import sympy as sp

import derive_sigma_pipi_cartesian_loop_up_only_estimate as up
import derive_sigma_pipi_cartesian_loop_large_epsilon_vertices as vertices_script
from sigma_pipi_loop_algebra import B0, EPSILON, F_PI, M_A, M_MESON, M_SIGMA, NORMALIZATION, TraceSum
from sigma_pipi_cartesian_loop_algebra import UP_QUARK_MASS
from test_derive_sigma_pipi_cartesian_loop_traces import FORBIDDEN_COMPONENTS, FORBIDDEN_LOOP_VALUES

RECORD = (pathlib.Path(__file__).resolve().parent.parent / 'memory' / 'alp' / 'radial-mode' / 'finite-mixing'
          / 'two-pion-loop' / 'cartesian-estimate' / 'cartesian-sigma-pipi-loop-up-only-estimate.md')
LEADING = -B0 ** 2 * UP_QUARK_MASS ** 2 / (8 * sp.pi ** 2 * F_PI ** 2 * M_MESON ** 2)
ZERO_GRAPHS = ['C1', 'C10', 'C11', 'C14', 'C19', 'C20', 'C21', 'C39', 'C41', 'C42', 'C43', 'C44', 'C45', 'C46',
               'C47', 'C48', 'C7', 'C8', 'C9', 'Cpi7', 'Cpi9', 'Csigma5', 'Csigma8']


class PoleDataTest(unittest.TestCase):
    def test_inverse_against_trace_series(self):
        """<C M_s^{-1}(C)> = -sum c_k/s^(k+1) and <(M_s^{-1}(C))^2> = sum (n+1) c_n/s^(n+2), with the c_k of the
        generic trace sums in the up-only limit (independent of the solution on span(P_u, 1))."""
        values = vertices_script.up_only_values()
        c = [values[f'c_{k}'] for k in range(3)]
        w = sp.Symbol('w', positive=True)              # w = 1/s
        first = sp.series(up.trace_c_inverse().subs(up.SQUARED, 1 / w), w, 0, 4).removeO()
        self.assertEqual(sp.expand(first + sum(c[k] * w ** (k + 1) for k in range(3))), 0)
        second = sp.series(up.trace_inverse_squared().subs(up.SQUARED, 1 / w), w, 0, 5).removeO()
        self.assertEqual(sp.expand(second - sum((n + 1) * c[n] * w ** (n + 2) for n in range(3))), 0)

    def test_pole_function(self):
        reference = up.U ** 2 * vertices_script.up_only_exact_alp_pole_function()
        self.assertEqual(sp.simplify(reference - up.alp_pole_function()), 0)

    def test_root(self):
        self.assertGreaterEqual(up.root_residual(), 2 * up.ROOT_ORDER + 4)
        values = vertices_script.up_only_values()
        generic = vertices_script.alp_pole_root().subs(
            {vertices_script.ALP_MASS_COEFFICIENT: values['A'], vertices_script.TRACE_C_POWER[0]: values['c_0'],
             vertices_script.TRACE_C_POWER[1]: values['c_1'], vertices_script.TRACE_C_POWER[2]: values['c_2']})
        difference = sp.expand(generic - up.alp_mass_squared())
        self.assertTrue(all(term.as_coeff_exponent(up.U)[1] > 4 for term in sp.Add.make_args(difference))
                        or difference == 0)

    def test_mass_difference(self):
        self.assertEqual(up.mass_difference(), (0, 2 * B0 * UP_QUARK_MASS))


class TraceTest(unittest.TestCase):
    def test_traces(self):
        self.assertEqual(up.up_only_trace(()), 3)
        self.assertEqual(up.up_only_trace(('l+', 'l-')), 2)
        self.assertEqual(up.up_only_trace(('Pu', 'l+', 'l-')), 2)
        self.assertEqual(up.up_only_trace(('Pu', 'l-', 'l+')), 0)
        self.assertEqual(up.up_only_trace(('l+', 'Pu', 'l-')), 0)
        with self.assertRaises(up.UnreducedTrace):
            up.up_only_trace(('l+',))

    def test_tree(self):
        tree = TraceSum.single(('SC10', 'l+', 'l-')) + TraceSum.single(('SC10', 'l-', 'l+'))
        value = up.up_only_scalar(tree.scale(-B0 * EPSILON / (2 * F_PI)))
        self.assertEqual(sp.expand(value - up.TREE_UP_ONLY), 0)

    def test_inverse_letter(self):
        value = up.up_only_scalar(TraceSum.single(('MinvC', 'l+', 'l-')))
        self.assertEqual(sp.expand(value - 2 * up.PROJECTOR_COEFFICIENT - 2 * up.IDENTITY_COEFFICIENT), 0)


class LeadingTest(unittest.TestCase):
    def test_total(self):
        for variant in up.VARIANTS:
            power, coefficient, labels = up.total_leading(variant)
            self.assertEqual(power, 2)
            self.assertEqual(sp.simplify(coefficient - LEADING), 0)
            self.assertEqual(labels, ['Csigma9'])
            self.assertTrue(up.exactness_check(variant))
            self.assertTrue(up.sum_of_graphs_check(variant))

    def test_zero_graphs(self):
        for variant in up.VARIANTS:
            self.assertEqual(up.zero_graphs(variant), ZERO_GRAPHS)

    def test_series_method_against_direct_series(self):
        """The truncated-series arithmetic gives the same leading term as a direct sympy series of the fully
        substituted exact expression, written as numerator over denominator (C2 and Csigma9)."""
        projector_part, identity_part = up.inverse_components()
        root = up.alp_mass_squared()
        u = up.U
        for label in ('C2', 'Csigma9'):
            value = up.replace_maximum(up.graph_scalars('r = 1')[label], 'r = 1')
            value = value.subs(NORMALIZATION, 1 + 2 * EPSILON ** 2 * up.trace_inverse_squared())
            value = value.subs({up.PROJECTOR_COEFFICIENT: projector_part, up.IDENTITY_COEFFICIENT: identity_part}).subs(up.SQUARED, M_A ** 2)
            value = value.subs({M_A: sp.sqrt(root), M_SIGMA: sp.sqrt(up.radial_mass_squared()), EPSILON: 1 / u})
            numerator, denominator = (sp.expand(part) for part in sp.fraction(sp.together(value)))
            low_numerator = min(term.as_coeff_exponent(u)[1] for term in sp.Add.make_args(numerator))
            low_denominator = min(term.as_coeff_exponent(u)[1] for term in sp.Add.make_args(denominator))
            lead_numerator = sum(term for term in sp.Add.make_args(numerator)
                                 if term.as_coeff_exponent(u)[1] == low_numerator) / u ** low_numerator
            lead_denominator = sum(term for term in sp.Add.make_args(denominator)
                                   if term.as_coeff_exponent(u)[1] == low_denominator) / u ** low_denominator
            power, coefficient = up.graph_leading('r = 1')[label]
            self.assertEqual(-(low_numerator - low_denominator), power, label)
            self.assertEqual(sp.simplify(lead_numerator / lead_denominator - coefficient), 0, label)

    def test_generic_comparison(self):
        for variant in up.VARIANTS:
            comparison = up.generic_comparison(variant)
            self.assertEqual(comparison['total'], (4, 0))
            self.assertTrue(all(comparison['graphs'].values()))
            self.assertIn('Csigma9', comparison['graphs'])

    def test_generic_products_contain_a_vanishing_letter(self):
        vanishing = {'SC30', 'SC12', 'SC21', 'SC02', 'SC11', 'SC20'}
        for variant in up.VARIANTS:
            _, numerator, _ = up.large.total_leading(variant)
            for key in numerator.terms:
                self.assertTrue(any(name in vanishing for word in key for name in word), key)


class RecordTest(unittest.TestCase):
    def test_record_equals_script_output(self):
        self.assertIn(up.markdown(), RECORD.read_text())

    def test_record_content(self):
        text = RECORD.read_text()
        for token in FORBIDDEN_COMPONENTS + FORBIDDEN_LOOP_VALUES:
            self.assertNotIn(token, text, token)
        self.assertIn('No loop integral', text)
        self.assertIn('-21/(64\\pi^2)', text)
        self.assertIn('authoritative', text)
        self.assertNotIn('action', text.lower().replace('interaction', ''))
        for shorthand in ('\\alpha', '\\beta'):          # no notation beyond the project records
            self.assertNotIn(shorthand, text, shorthand)


if __name__ == '__main__':
    unittest.main()
