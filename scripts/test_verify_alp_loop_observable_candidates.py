"""Unit tests for verify_alp_loop_observable_candidates.py (exact arithmetic only; no floats).

Run with:
    .venv/bin/python -m unittest discover -s scripts -p test_verify_alp_loop_observable_candidates.py
"""

import sys
import unittest
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))

import verify_alp_loop_observable_candidates as candidates  # noqa: E402

lp = candidates.leading_power
coefficient = candidates.epsilon_coefficient
NEUTRAL_MESONS = ('pi0', 'eta', 'eta_prime')
CHARGED_EXTERNAL = {'ud': 4, 'us': 6, 'ds': 8}


def summed(coefficients, label):
    return coefficients['tadpole'][label] + coefficients['one_point'][label]


class InfrastructureTests(unittest.TestCase):
    def test_unmixed_poles_rational(self):
        poles = candidates.unmixed_neutral_poles(candidates.meson_sector_inputs())
        self.assertEqual(poles['pi0'][:3], ([1, -1, 0, 0], sp.Rational(1, 50), 1))
        self.assertEqual(poles['eta'][1], sp.Rational(109, 450))
        self.assertEqual(poles['eta_prime'][1], sp.Rational(97, 100))
        self.assertEqual(poles['alp'][1:], (sp.Rational(3, 5), 1))

    def test_closed_form_potential_equals_builder(self):
        """Closed-form neutral dressed potential = builder (symbolic invariants and epsilon), degree <= 4."""
        self.assertFalse(candidates.closed_form_against_builder())

    def test_perturbative_pole_solves_pencil(self):
        """(M^2 - s K) u = 0 and u^T K u = n0 through the truncation order, kappa frame (K depends on epsilon)."""
        order = candidates.TWO_POINT_ORDER
        vertices = candidates.frame_vertices('kappa', 'ud')
        size = candidates.NEUTRAL_SIZE
        for name, (pole, vector, norm) in candidates.frame_poles('kappa').items():
            for i in range(size):
                residual = sum((candidates.truncate((vertices.mass[i][j] - pole*vertices.kinetic[i][j])*vector[j], order)
                                for j in range(size)), candidates.ZERO)
                self.assertFalse(candidates.truncate(residual, order), (name, i))
            norm_series = sum((candidates.truncate(vertices.kinetic[i][j]*vector[i]*vector[j], order)
                               for i in range(size) for j in range(size)), candidates.ZERO)
            self.assertEqual(candidates.truncate(norm_series, order), candidates.rational(norm))

    def test_neutral_kappa_map_is_linear(self):
        """Kappa frame (symbolic kappa) = dressed frame after N_q -> N_q - 2 eps kappa_q a, exactly."""
        self.assertEqual(candidates.kappa_map_linearity_check(), [])

    def test_no_charged_or_derivative_cubic_vertices(self):
        for frame in ('dressed', 'kappa'):
            for pair in candidates.CHARGED_PAIRS:
                vertices = candidates.frame_vertices(frame, pair)
                self.assertEqual(vertices.flags['derivative_cubic'], [])
                self.assertTrue(all(max(key) < candidates.NEUTRAL_SIZE for key in vertices.cubic))


class MesonPoleMassTests(unittest.TestCase):
    """C1."""

    def test_tree_shift_power_two(self):
        tree = candidates.tree_pole_check()
        for name in NEUTRAL_MESONS:
            self.assertEqual(lp(tree[name]), 2)
            self.assertFalse(coefficient(tree[name], 3))

    def test_meson_lines(self):
        """eps^0 part nonzero (summed over labels), eps^1 zero, eps-dependent part from eps^2.
        Per label: kaon lines give no eps^0 part to pi0, pi0 lines none to kaons, and the charged-pion
        pair none at all to the charged pion (eps-independent and zero on shell)."""
        expected = {
            'pi0': {'pi0': (0, 2), 'eta': (0, 2), 'eta_prime': (0, 2), 'ud': (0, 2), 'us': (2, 2), 'ds': (2, 2)},
            'eta': {label: (0, 2) for label in ('pi0', 'eta', 'eta_prime', 'ud', 'us', 'ds')},
            'eta_prime': {label: (0, 2) for label in ('pi0', 'eta', 'eta_prime', 'ud', 'us', 'ds')},
            4: {'pi0': (0, 2), 'eta': (0, 2), 'eta_prime': (0, 2), 'ud': (None, None)},
            6: {'pi0': (2, 2), 'eta': (0, 2), 'eta_prime': (0, 2), 'us': (None, None)},
            8: {'pi0': (2, 2), 'eta': (0, 2), 'eta_prime': (0, 2), 'ds': (None, None)},
        }
        for external, labels in expected.items():
            self.assertEqual(candidates.pole_mass_summary(external)['meson_lines'], labels, external)
            coefficients = candidates.pole_mass_coefficients('dressed', external)
            for label in labels:
                self.assertFalse(coefficient(summed(coefficients, label), 1))

    def test_alp_line_power_two(self):
        for external in NEUTRAL_MESONS + tuple(CHARGED_EXTERNAL.values()):
            summary = candidates.pole_mass_summary(external)
            self.assertEqual(summary['alp_line'], 2, external)
            self.assertEqual(summary['alp_line_tadpole_cp_even'], 2, external)

    def test_bubbles_and_one_point(self):
        for external in NEUTRAL_MESONS:
            summary = candidates.pole_mass_summary(external)
            self.assertEqual(summary['bubble_alp_alp'], 4)
            self.assertEqual(summary['bubble_alp_meson'], 6)
            self.assertEqual(summary['bubble_meson_meson'], 8)
            self.assertIsNone(summary['one_point'])              # zero through eps^4
            self.assertTrue(summary['cubic_vertices_cp_even_zero'])

    def test_alp_pair_bubble_without_p2_theta2_is_not_zero(self):
        """With P_hat_2 = Theta_hat_2 = 0 the vertex V(pi0, A, A) starts at eps^4 with
        (84500/667) tr(S_hat_3) (P_hat_1u - P_hat_1d): b_{pi0,aa} ~ eps^8, not identically zero."""
        self.assertEqual(candidates.pole_mass_summary('pi0')['bubble_alp_alp_without_P2_Theta2'], 8)
        vertex = candidates.pole_mass_coefficients('dressed', 'pi0')['bubble_vertex'][('alp', 'alp')]
        removed = list(candidates.pseudoscalar_invariant(2)) + [candidates.topological_invariant(2)]
        g = candidates.generator
        trace = sum((g(v) for v in candidates.scalar_invariant(3)), candidates.ZERO)
        difference = g(candidates.pseudoscalar_invariant(1)[0]) - g(candidates.pseudoscalar_invariant(1)[1])
        self.assertEqual(coefficient(candidates.set_to_zero(vertex, removed), 4),
                         candidates.rational(sp.Rational(84500, 667))*trace*difference)

    def test_quark_mass_shift_form_neutral(self):
        result = candidates.neutral_quark_mass_shift_check()
        self.assertFalse(result['eta_residual'])
        self.assertFalse(result['sum_ud_minus_closed_form'])
        self.assertFalse(result['strange_minus_closed_form'])

    def test_charged_remainder(self):
        for pair in candidates.CHARGED_PAIRS:
            result = candidates.charged_remainder_check(pair)
            self.assertFalse(result['difference'], pair)
            if pair != 'ud':
                self.assertTrue(result['remainder'])

    def test_axial_on_shell_equality(self):
        for external in NEUTRAL_MESONS + tuple(CHARGED_EXTERNAL.values()):
            differences = candidates.frame_comparison(external)
            self.assertTrue(all(not value for value in differences.values()), external)


class MixingEntryTests(unittest.TestCase):
    """C2."""

    def test_tree_entry_and_frame_difference(self):
        result = candidates.mixing_entry_tree()
        self.assertEqual(lp(result['dressed']), 1)
        self.assertFalse(result['difference_minus_closed_form'])
        self.assertFalse(result['at_unmixed_mass'])
        self.assertEqual(lp(result['at_pole']), 3)
        g = candidates.generator
        kappa_difference = g(candidates.kappa_symbols[0]) - g(candidates.kappa_symbols[1])
        self.assertEqual(coefficient(result['at_pole'], 3), -kappa_difference*coefficient(result['pole_shift'], 2))

    def test_one_loop_entry(self):
        result = candidates.mixing_entry_one_loop()
        for frame in ('dressed', 'kappa'):
            powers = {label: lp(value) for label, value in result['entries'][frame].items()}
            self.assertEqual(powers, {'pi0': 1, 'eta': 1, 'eta_prime': 1, 'ud': 1, 'us': 1, 'ds': 1, 'alp': 3})
        differences = {label: lp(value) for label, value in result['on_shell_difference'].items()}
        self.assertEqual(differences, {'pi0': 1, 'eta': 1, 'eta_prime': 1, 'ud': 1, 'us': 1, 'ds': 1, 'alp': 3})
        powers = [candidates.product_power(*value) for value in result['bubble_powers'].values()]
        self.assertEqual(candidates.lowest(powers), 5)


class AlpPoleMassTests(unittest.TestCase):
    """C3."""

    def test_tree_and_loop_powers(self):
        self.assertEqual(lp(candidates.tree_pole_check()['alp']), 2)
        summary = candidates.pole_mass_summary('alp')
        self.assertEqual({label: value[0] for label, value in summary['meson_lines'].items()},
                         {label: 2 for label in ('pi0', 'eta', 'eta_prime', 'ud', 'us', 'ds')})
        self.assertEqual(summary['alp_line'], 4)
        self.assertEqual(summary['alp_line_tadpole_cp_even'], 4)
        self.assertEqual(summary['one_point'], 4)
        self.assertEqual(summary['bubble_alp_meson'], 4)
        self.assertEqual(summary['bubble_alp_alp'], 6)
        self.assertEqual(summary['bubble_meson_meson'], 6)
        self.assertTrue(summary['cubic_vertices_cp_even_zero'])
        one_point = candidates.pole_mass_coefficients('dressed', 'alp')['one_point']['alp']
        self.assertFalse(candidates.set_to_zero(one_point, candidates.CP_ODD_INVARIANTS))

    def test_axial_on_shell_equality(self):
        self.assertTrue(all(not value for value in candidates.frame_comparison('alp').values()))

    def test_absorption_by_bare_mass(self):
        result = candidates.alp_mass_absorption_check()
        self.assertFalse(result['difference'])
        self.assertEqual(coefficient(result['derivative'], 0), candidates.ONE)


class AmplitudeTests(unittest.TestCase):
    """C4."""

    def test_tree_powers_and_conjecture(self):
        result = candidates.tree_amplitude_conjecture_check()
        amplitude = result['amplitude']
        self.assertEqual(lp(amplitude), 1)
        for power in (0, 2, 4):
            self.assertFalse(coefficient(amplitude, power))
        self.assertTrue(coefficient(amplitude, 3))
        self.assertFalse(result['first_order_minus_conjecture'])
        self.assertFalse(result['first_order_minus_reduced_form'])
        self.assertEqual(result['alp_vector_minus_two_y'], [candidates.ZERO]*3)

    def test_tree_frame_equality(self):
        result = candidates.tree_amplitude_frame_check()
        self.assertFalse(result['difference'])
        self.assertEqual(result['amplitudes']['dressed'], candidates.tree_amplitude())

    def test_one_loop_powers(self):
        result = candidates.amplitude_one_loop()
        self.assertEqual({r: lp(v) for r, v in result['tadpole'].items()},
                         {'pi0': 1, 'eta': 1, 'eta_prime': 1, 'alp': 3})
        bubbles = {k: lp(v) for k, v in result['bubble'].items()}
        self.assertEqual(bubbles[('alp', 'alp')], 5)
        self.assertEqual([bubbles[(r, 'alp')] for r in NEUTRAL_MESONS], [3, 3, 3])
        self.assertEqual([bubbles[(r, r)] for r in NEUTRAL_MESONS], [1, 1, 1])
        self.assertEqual({r: lp(v) for r, v in result['external'].items()},
                         {'pi0': 1, 'eta': 1, 'eta_prime': 1, 'alp': 3})

    def test_source_shift_identity(self):
        result = candidates.source_shift_identity()
        self.assertFalse(result['residual'])
        self.assertTrue(result['alp_three_meson_terms'])


if __name__ == '__main__':
    unittest.main()
