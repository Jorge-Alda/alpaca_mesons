"""Unit tests for verify_alp_loop_observable_selection.py (exact arithmetic only; no floats).

Run with:
    .venv/bin/python -m unittest discover -s scripts -p test_verify_alp_loop_observable_selection.py
"""

import sys
import unittest
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))

import verify_alp_loop_observable_candidates as candidates  # noqa: E402
import verify_alp_loop_observable_selection as selection  # noqa: E402

PAIRS = selection.ALP_THRESHOLD_PAIRS
CONDITIONS = ('P1_isovector_zero', 'P1_P3_isovector_zero', 'isospin_limit')


def isospin_symbol(stem, kind):
    return sp.Symbol(f'{stem}_{kind}', real=True)


P1p, P1m = isospin_symbol('P_hat_1', 'plus'), isospin_symbol('P_hat_1', 'minus')
S2p, S2m = isospin_symbol('S_hat_2', 'plus'), isospin_symbol('S_hat_2', 'minus')
P3m = isospin_symbol('P_hat_3', 'minus')
P1u, P1d, P1s = candidates.pseudoscalar_invariant(1)
S2u, S2d, _ = candidates.scalar_invariant(2)
P3u, P3d, _ = candidates.pseudoscalar_invariant(3)
Theta1 = candidates.topological_invariant(1)

# isoscalar combinations that appear as factors
ETA_ISOSCALAR = 15704*P1p - 1170*P1s + 71*Theta1
PION_QUADRATIC = (3538294003790*P1m**2 + 1415250938790*P1p**2 - 7726132440000*P1p*P1s
                  + 468850772000*P1p*Theta1 + 287811225000*P1s**2 - 34930935000*P1s*Theta1
                  - 387966787117*S2p + 1059870250*Theta1**2)
ETA_ISOVECTOR = (2196220*P1m*P1p - 3510000*P1m*P1s + 213000*P1m*Theta1 - 172753*S2m)
EXPECTED_LEADING = {
    ('pi0', 'alp'): sp.Rational(-60937500000, 242192415118049501)*P1m*PION_QUADRATIC,
    ('eta', 'alp'): sp.Rational(-822656250000, 37725572390564191)*ETA_ISOSCALAR*ETA_ISOVECTOR,
    ('eta_prime', 'alp'): sp.Rational(-507812500000, 37725572390564191)*ETA_ISOSCALAR*ETA_ISOVECTOR,
    ('alp', 'alp'): sp.Rational(-1015625000, 7227893319804338706056544509)*PION_QUADRATIC*(
        52774554632800*P1m**3 + 50975930288400*P1m*P1p**2 - 324152189010000*P1m*P1p*P1s
        + 19670773863000*P1m*P1p*Theta1 + 259030102500000*P1m*P1s**2 - 31437841500000*P1m*P1s*Theta1
        - 11639003613510*P1m*S2p + 953883225000*P1m*Theta1**2 - 4702884287010*P1p*S2m
        + 25497565411500*P1s*S2m - 865464371261*P3m - 1547288157450*S2m*Theta1),
}


class NeutralKineticTermTests(unittest.TestCase):
    """S1."""

    def test_closed_form_kinetic_term_is_free(self):
        self.assertEqual(selection.neutral_kinetic_term_check(), sp.zeros(2, 2))

    def test_builder_has_no_neutral_derivative_vertices(self):
        for cp_conserving, entry in selection.builder_neutral_derivative_vertices().items():
            self.assertEqual(entry['derivative_quartic'], 0, cp_conserving)
            self.assertTrue(entry['kinetic_constant'], cp_conserving)
            self.assertTrue(all(value == [] for value in entry['flags'].values()), cp_conserving)


class ContactStructureTests(unittest.TestCase):
    """S2."""

    @classmethod
    def setUpClass(cls):
        cls.result = selection.contact_structure_check()

    def test_degree_one_is_constant_on_shell(self):
        self.assertEqual(self.result['linear_derivatives'], (0, 0))

    def test_degree_two_in_span_of_one_and_sum_of_squares(self):
        self.assertEqual(self.result['quadratic_remainder_derivatives'], (0, 0))
        self.assertNotEqual(self.result['quadratic_beta'], 0)

    def test_power_sums_on_plane(self):
        e1, e2, e3 = self.result['elementary_symbols']
        self.assertEqual(self.result['power_sum_2'], e1**2 - 2*e2)
        self.assertEqual(self.result['power_sum_3'], e1**3 - 3*e1*e2 + 3*e3)
        s12, s13, _ = selection.mandelstam
        expected = sp.factor((s12 - s13)*(selection.plane_sum - 2*s12 - s13)*(selection.plane_sum - s12 - 2*s13))
        self.assertEqual(sp.expand(self.result['jacobian_e2_e3'] - expected), 0)


class BubbleFunctionTests(unittest.TestCase):
    """S3."""

    def test_denominator_decomposition(self):
        self.assertEqual(selection.bubble_denominator_identity(), 0)

    def test_taylor_coefficients_at_zero_positive(self):
        m = sp.Symbol('m', positive=True)
        result = selection.bubble_taylor_coefficients_at_zero()
        self.assertEqual(result[2]['equal_mass'], 1/(60*m**4))
        self.assertEqual(result[3]['equal_mass'], 1/(420*m**6))
        for n in (2, 3):
            self.assertEqual(result[n]['integrand_minus_expected'], 0)
            self.assertIs(result[n]['t_integrand_positive'], True)

    def test_e3_coefficient_is_half_third_derivative(self):
        coefficient, taylor = selection.e3_coefficient_from_expansion(6)
        self.assertEqual(coefficient, taylor[3]/2)

    def test_third_derivative_positive_at_symmetric_point(self):
        result = selection.bubble_third_derivative_check()
        self.assertEqual(result['third_minus_expected'], 0)
        for label in ('alp_meson', 'alp_alp'):
            self.assertIs(result[label]['threshold_distance_positive'], True, label)
            self.assertIs(result[label]['t_integrand_positive'], True, label)
            self.assertEqual(result[label]['t_form_minus_manifest'], 0, label)
            self.assertEqual(result[label]['positive_form_minus_third'], 0, label)

    def test_subtraction_scale_dependence_is_constant(self):
        result = selection.subtraction_scale_dependence()
        self.assertEqual(result['integrand'], 1)
        self.assertEqual(result['integral'], 1)


class AlpThresholdBubbleTests(unittest.TestCase):
    """S4."""

    @classmethod
    def setUpClass(cls):
        cls.result = selection.alp_threshold_bubbles()

    def test_leading_powers(self):
        powers = {k: v['leading_power'] for k, v in self.result['bubbles'].items()}
        self.assertEqual(powers, {('pi0', 'alp'): 3, ('eta', 'alp'): 3, ('eta_prime', 'alp'): 3, ('alp', 'alp'): 5})

    def test_leading_coefficients(self):
        for key, expected in EXPECTED_LEADING.items():
            coefficient = selection.isospin_form(self.result['bubbles'][key]['leading_coefficient'])
            self.assertEqual(sp.expand(coefficient - expected), 0, key)

    def test_symbols(self):
        common = sorted([P1d, P1s, P1u, S2d, S2u, Theta1], key=str)
        for key in PAIRS[:3]:
            self.assertEqual(self.result['bubbles'][key]['symbols'], common)
        self.assertEqual(self.result['bubbles'][('alp', 'alp')]['symbols'], sorted(common + [P3d, P3u], key=str))

    def test_tree_under_conditions(self):
        tree = self.result['tree']
        self.assertEqual(candidates.leading_power(tree['full']), 1)
        for name in CONDITIONS:
            self.assertFalse(tree[name], name)          # m_u = m_d test inputs

    def test_survival_without_isovector_p1(self):
        bubbles = self.result['bubbles']
        for name in ('P1_isovector_zero', 'P1_P3_isovector_zero'):
            self.assertFalse(bubbles[('pi0', 'alp')][name]['series'], name)
            for key in (('eta', 'alp'), ('eta_prime', 'alp'), ('alp', 'alp')):
                self.assertEqual(bubbles[key][name]['leading_power'], bubbles[key]['leading_power'], (key, name))
        # (eta, alp): isovector S_hat_2 times the isoscalar combination with P_hat_1u = P_hat_1d
        reduced = bubbles[('eta', 'alp')]['P1_isovector_zero']['coefficient_at_original_power'].as_expr()
        expected = sp.Rational(-822656250000, 37725572390564191)*(31408*P1d - 1170*P1s + 71*Theta1)*(-172753)*(S2u - S2d)
        self.assertEqual(sp.expand(reduced - expected), 0)
        # (alp, alp) with P_hat_1 and P_hat_3 isovector zero: proportional to S_hat_2 isovector
        reduced = bubbles[('alp', 'alp')]['P1_P3_isovector_zero']['coefficient_at_original_power']
        self.assertTrue(reduced)
        self.assertFalse(reduced.compose(candidates.generator(S2u), candidates.generator(S2d)))

    def test_isospin_limit(self):
        for key in PAIRS:
            self.assertFalse(self.result['bubbles'][key]['isospin_limit']['series'], key)


class FrameTests(unittest.TestCase):
    """S5."""

    def test_kappa_frame_equality_through_order_five(self):
        result = selection.alp_threshold_bubbles_frame_check('symbolic', 5)
        for key in PAIRS:
            self.assertFalse(result['difference'][key], key)
            self.assertFalse(result['dressed_minus_closed_form'][key], key)
            self.assertTrue(result['coefficients']['kappa'][key], key)


if __name__ == '__main__':
    unittest.main()
