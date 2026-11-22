"""Unit tests for derive_alp_invariant_loop_crossing.py.

Run with:
    .venv/bin/python -m unittest discover -s scripts -p test_derive_alp_invariant_loop_crossing.py
"""

import math
import sys
import unittest
from functools import lru_cache
from pathlib import Path

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))

import derive_alp_invariant_loop_crossing as crossing  # noqa: E402


@lru_cache(maxsize=None)
def vertices(epsilon, generator=(0, 0, 0)):
    parameters = crossing.test_parameters(epsilon)
    terms = crossing.build_lagrangian_exact(parameters, first_frame_generator=generator)
    return crossing.extract_vertices(terms)


@lru_cache(maxsize=1)
def symbolic_frame_result():
    return crossing.symbolic_frame_difference()


def neutral_pole_indices(pole_vectors):
    alp_index = crossing.alp_like_pole(pole_vectors)
    return alp_index, [k for k in range(crossing.NEUTRAL_SECTOR_SIZE) if k != alp_index]


FINITE_MIXING = sp.Rational(1, 3)
SMALL = (sp.Rational(1, 20), sp.Rational(1, 40))


class LagrangianStructureTests(unittest.TestCase):
    def test_exact_builder_matches_symbolic_builder(self):
        """Sparse exact builder (N_q = sqrt2 P_q) equals the sympy builder on a field subset."""
        parameters = crossing.test_parameters(sp.Rational(1, 4))
        generator = crossing.EXAMPLE_FRAME_GENERATOR
        active = ('P_u', 'P_d', 'P_s', 'a', 'P_1')
        symbolic = crossing.build_lagrangian(parameters, first_frame_generator=generator, active_fields=active)
        exact_terms = crossing.build_lagrangian_exact(parameters, first_frame_generator=generator)
        size = len(crossing.FIELD_NAMES)
        active_indices = [crossing.FIELD_NAMES.index(name) for name in active]
        rebuilt = 0
        for exponents, coefficient in exact_terms:
            used = [i % size for i, power in enumerate(exponents) if power]
            if any(index not in active_indices for index in used):
                continue
            term = sp.Rational(coefficient.numerator, coefficient.denominator)
            for index, power in enumerate(exponents):
                if power:
                    variable = (crossing.meson_fields if index < size else crossing.derivative_fields)[index % size]
                    scale = 1/sp.sqrt(2) if index % size < 3 else 1      # N_q = sqrt2 P_q
                    term *= (variable/scale)**power
            rebuilt += term
        all_fields = list(crossing.meson_fields) + list(crossing.derivative_fields)
        constant = symbolic.subs({variable: 0 for variable in all_fields})   # dropped by the exact builder
        self.assertEqual(sp.expand(rebuilt - (symbolic - constant)), 0)

    def test_no_linear_terms_and_no_derivative_cubic_vertices(self):
        for generator in ((0, 0, 0), crossing.EXAMPLE_FRAME_GENERATOR):
            flags = vertices(FINITE_MIXING, generator)['flags']
            self.assertEqual(flags['linear_terms'], [])
            self.assertEqual(flags['field_derivative_quadratic'], [])
            self.assertEqual(flags['derivative_cubic'], [])
            self.assertEqual(flags['odd_derivative_quartic'], [])

    def test_charged_masses_and_no_charged_mixing(self):
        parameters = crossing.test_parameters(FINITE_MIXING)
        masses, _ = crossing.tree_spectrum(vertices(FINITE_MIXING))
        b0 = float(parameters['B0'])
        quark = [float(m) for m in parameters['quark_masses']]
        expected = [b0*(quark[0] + quark[1])]*2 + [b0*(quark[0] + quark[2])]*2 + [b0*(quark[1] + quark[2])]*2
        np.testing.assert_allclose(masses[4:], expected, rtol=1e-13)


class FrameInvarianceTests(unittest.TestCase):
    def test_pole_masses_equal_in_both_frames(self):
        dressed, _ = crossing.tree_spectrum(vertices(FINITE_MIXING))
        framed, _ = crossing.tree_spectrum(vertices(FINITE_MIXING, crossing.EXAMPLE_FRAME_GENERATOR))
        np.testing.assert_allclose(dressed, framed, rtol=1e-12)

    def test_on_shell_coefficients_equal_per_loop_function(self):
        """Finite mixing (epsilon = 1/3): every loop-function coefficient of every pole shift agrees."""
        dressed = vertices(FINITE_MIXING)
        framed = vertices(FINITE_MIXING, crossing.EXAMPLE_FRAME_GENERATOR)
        for pole in range(len(crossing.FIELD_NAMES)):
            first = crossing.pole_shift_coefficients(dressed, pole)
            second = crossing.pole_shift_coefficients(framed, pole)
            scale = np.max(np.abs(first['tadpole']))
            np.testing.assert_allclose(first['tadpole'] + first['one_point'],
                                       second['tadpole'] + second['one_point'], atol=1e-12*scale)
            for key, value in first['bubble'].items():
                self.assertAlmostEqual(value, second['bubble'][key], delta=1e-12*scale)

    def test_frame_difference_is_inverse_propagator_times_matrix(self):
        """Exact: T^-T delta'(s) T^-1 - delta(s) = D(s) X + X^T D(s), X linear in charged tadpoles."""
        result = symbolic_frame_result()
        self.assertEqual(len(result['solution']), 1)
        self.assertTrue(result['neutral_difference_free_of_G_H'])
        self.assertTrue(any(entry != 0 for entry in result['difference']))

    def test_off_shell_entry_differs(self):
        """The off-shell self-energy is frame dependent: some entry differs at s = 0."""
        difference = symbolic_frame_result()['difference']
        values = [entry.subs(sp.Symbol('s'), 0) for entry in difference]
        self.assertTrue(any(value != 0 for value in values))

    def test_cubic_vertices_congruent(self):
        result = symbolic_frame_result()
        self.assertEqual(crossing.cubic_congruence_residual(result['dressed'], result['framed'],
                                                            result['transform']), 0.0)

    def test_second_order_generator_changes_vertex_by_inverse_propagator(self):
        result = crossing.second_order_frame_vertex_check()
        self.assertEqual(len(result['solution']), 1)
        solution = result['solution'][0]
        generator = (sp.Rational(1, 2), sp.Rational(-3, 4), sp.Rational(1, 5))
        ratios = [sp.simplify(solution[sp.Symbol(name)]/value)
                  for name, value in zip(('w_u', 'w_d', 'w_s'), generator)]
        self.assertTrue(all(sp.simplify(r - ratios[0]) == 0 for r in ratios))
        self.assertTrue(any(entry.has(sp.Symbol('s')) for entry in result['vertex_difference']))

    def test_linear_spurion_tadpole_combination_is_not_invariant(self):
        result = crossing.raw_tadpole_coefficient_shift()
        self.assertTrue(result['sum_equals_invariant'])
        self.assertNotEqual(sp.simplify(result['linear_spurion_minus_invariant']), 0)


class ClosedFormTests(unittest.TestCase):
    def residuals(self, closed_form):
        values = []
        for epsilon in SMALL:
            values.append(closed_form(epsilon))
        return values

    def assert_relative_scaling(self, pairs, exponent):
        """Relative residual |exact - closed|/|exact| decreases by 2^exponent when epsilon halves."""
        relative = [abs(exact - closed)/abs(exact) for exact, closed in pairs]
        ratio = relative[0]/relative[1]
        self.assertGreater(ratio, 2**exponent*0.8)
        self.assertLess(ratio, 2**exponent*1.25)

    def alp_line_pairs(self, meson):
        pairs = []
        for epsilon in SMALL:
            parameters = crossing.test_parameters(epsilon)
            data = vertices(epsilon)
            masses, vectors = crossing.tree_spectrum(data)
            alp_index, neutral = neutral_pole_indices(vectors)
            if meson == 'pion-like':
                _, flavour = crossing.neutral_meson_basis(parameters)
                external, label = neutral[0], flavour[0]
            else:
                external, label = crossing.CHARGED_PAIRS[meson][0], meson
            coefficients = crossing.pole_shift_coefficients(data, external)
            exact = coefficients['tadpole'][alp_index] + coefficients['one_point'][alp_index]
            pairs.append((exact, crossing.alp_line_coefficient(parameters, label, masses[alp_index])))
        return pairs

    def test_alp_line_neutral(self):
        self.assert_relative_scaling(self.alp_line_pairs('pion-like'), 2)

    def test_alp_line_charged(self):
        for pair in crossing.CHARGED_PAIRS:
            self.assert_relative_scaling(self.alp_line_pairs(pair), 2)

    def test_charged_symbolic_reduction_to_pole_aligned_form(self):
        self.assertEqual(crossing.derive_charged_alp_line_coefficient()['residual'], 0)

    def test_alp_pair_bubble(self):
        pairs = []
        for epsilon in SMALL:
            parameters = crossing.test_parameters(epsilon)
            data = vertices(epsilon)
            _, vectors = crossing.tree_spectrum(data)
            alp_index, neutral = neutral_pole_indices(vectors)
            _, flavour = crossing.neutral_meson_basis(parameters)
            for n in range(3):
                exact = crossing.pole_shift_coefficients(data, neutral[n])['bubble'][(alp_index, alp_index)]
                pairs.append((exact, crossing.alp_pair_bubble_coefficient(parameters, flavour[n])))
        for n in range(3):
            self.assert_relative_scaling([pairs[n], pairs[3 + n]], 2)

    def test_meson_lines_through_second_order(self):
        maxima = []
        for epsilon in SMALL:
            parameters = crossing.test_parameters(epsilon)
            data = vertices(epsilon)
            masses, vectors = crossing.tree_spectrum(data)
            alp_index, neutral = neutral_pole_indices(vectors)
            residual = 0.0
            for n in range(3):
                coefficients = crossing.pole_shift_coefficients(data, neutral[n])
                total = coefficients['tadpole'] + coefficients['one_point']
                closed_neutral, closed_charged = crossing.meson_line_coefficients(
                    parameters, n, masses[alp_index], masses[neutral[n]])
                for p in range(3):
                    residual = max(residual, abs(total[neutral[p]] - closed_neutral[p]))
                for pair, (first, second) in crossing.CHARGED_PAIRS.items():
                    residual = max(residual, abs(total[first] + total[second] - closed_charged[pair]))
            maxima.append(residual)
        ratio = maxima[0]/maxima[1]
        self.assertGreater(ratio, 16*0.75)
        self.assertLess(ratio, 16*1.3)

    def test_one_point_contributions_start_beyond_second_order(self):
        """Loop-induced one-point functions do not enter the meson poles at O(epsilon^2)."""
        sizes = []
        for epsilon in SMALL:
            data = vertices(epsilon)
            _, vectors = crossing.tree_spectrum(data)
            _, neutral = neutral_pole_indices(vectors)
            sizes.append(max(np.max(np.abs(crossing.pole_shift_coefficients(data, pole)['one_point']))
                             for pole in neutral))
        self.assertGreater(sizes[0]/sizes[1], 2**4*0.8)


class LoopFunctionAndTableTests(unittest.TestCase):
    def test_bubble_integral_against_closed_form(self):
        for s_value, mass in ((0.0182, 0.1), (0.0182, 1.0), (0.3, 0.5)):
            numerical = crossing.loop_function_bubble(s_value, mass, mass, 0.5929)
            self.assertAlmostEqual(numerical.real, crossing.equal_mass_bubble_below_threshold(s_value, mass, 0.5929),
                                   places=9)
            self.assertEqual(numerical.imag, 0.0)

    def test_bubble_imaginary_part_above_threshold(self):
        s_value, mass = 1.0, 0.3
        value = crossing.loop_function_bubble(s_value, mass, mass, 0.5929)
        self.assertAlmostEqual(value.imag, math.pi*math.sqrt(1 - 4*mass**2/s_value), places=6)

    def test_tadpole_zero(self):
        self.assertAlmostEqual(crossing.loop_function_tadpole(math.e*0.5929, 0.5929), 0.0, places=15)

    def test_table_crossings_solve_ratio_equal_one(self):
        b0, rows = crossing.table_rows()
        self.assertAlmostEqual(b0, crossing.NEUTRAL_PION_MASS**2/(crossing.QUARK_MASSES[0] + crossing.QUARK_MASSES[1]))
        f = crossing.PION_DECAY_CONSTANT
        pion_squared = crossing.NEUTRAL_PION_MASS**2
        for row in rows:
            epsilon = row['alp_line_crossing_with_logarithm']
            self.assertAlmostEqual(epsilon**2*abs(row['A0'])/(2*f**2), 1.0, places=12)
            epsilon = row['alp_pair_crossing_with_logarithm']
            self.assertAlmostEqual(epsilon**4*pion_squared*abs(row['Bbar'])/(32*math.pi**2*f**2), 1.0, places=6)


if __name__ == '__main__':
    unittest.main()
