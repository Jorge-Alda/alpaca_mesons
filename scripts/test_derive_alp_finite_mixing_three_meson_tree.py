"""Unit tests for derive_alp_finite_mixing_three_meson_tree.py.

Run: .venv/bin/python -m unittest discover -s scripts -p test_derive_alp_finite_mixing_three_meson_tree.py
All checks are exact (SymPy) unless a decimal comparison is stated.
"""
import random
import unittest

import sympy as sp

from derive_alp_finite_mixing_three_meson_tree import (
    ALP_POLE_SQUARED_MASS, DOWN_STRANGE_BASIS, KINEMATIC_POINT, NEUTRAL_BASIS, OFFSETS, alp_inverse_propagator,
    alp_like_leg, alp_mass_entry, amplitude_series, amplitudes_at_both_signs, compare_amplitude, cubic_vertex,
    dressed_source_arrays, inverse_pole_operator, is_exactly_zero, leading_coefficient, matrix_C, mixed_meson_leg,
    normalization, overlap_argument_difference_series, pole_operator, quartic_contact_classes, solver_theory,
    tuned_point, unmixed_neutral_matrix, unmixed_neutral_squared_masses, verification_inputs)
from derive_alp_lo_joint_vertices import vertex

M_A2 = ALP_POLE_SQUARED_MASS
S12, S23 = KINEMATIC_POINT


def random_hermitian(generator, down_strange=True):
    """Hermitian matrix commuting with the quark charge (diagonal plus a d-s entry)."""
    matrix = sp.diag(*[sp.Rational(generator.randint(-4, 4), generator.randint(1, 3)) for _ in range(3)])
    if down_strange:
        entry = sp.Rational(generator.randint(-3, 3), 2) + sp.I * sp.Rational(generator.randint(-3, 3), 4)
        matrix[1, 2], matrix[2, 1] = entry, sp.conjugate(entry)
    return matrix


def momentum_polynomial_to_products(expression, symbols, products):
    total = 0
    for powers, coefficient in sp.Poly(expression, *symbols).terms():
        legs = [i for i, power in enumerate(powers) for _ in range(power)]
        total += coefficient * (products[legs[0]][legs[1]] if legs else 1)
    return sp.expand(total)


class VerticesAndQuadraticBlocks(unittest.TestCase):
    """The closed vertex classes, M^2_aa and the mixing entries equal the CCWZ Lagrangian."""

    @classmethod
    def setUpClass(cls):
        cls.invariants = verification_inputs(down_strange=True)
        cls.eps = sp.Rational(2, 7)
        cls.sources = dressed_source_arrays(cls.invariants)
        generator = random.Random(11)
        cls.legs = [(random_hermitian(generator), sp.Rational(generator.randint(-3, 3), generator.randint(1, 3)))
                    for _ in range(4)]

    def ccwz(self, legs):
        symbols = sp.symbols(f'q0:{len(legs)}')
        expression = vertex([leg[0] for leg in legs], [leg[1] for leg in legs], symbols, *self.sources,
                            self.invariants['f_pi'], self.invariants['f_pi'] / self.eps,
                            self.invariants['anomaly_mass_squared'], 0)
        return expression, symbols

    def test_quartic_contact_classes(self):
        products = [[sp.Symbol(f'pp{min(i, j)}{max(i, j)}') for j in range(4)] for i in range(4)]
        expression, symbols = self.ccwz(self.legs)
        expected = momentum_polynomial_to_products(expression, symbols, products)
        classes = quartic_contact_classes(self.invariants, self.eps, self.legs, products)
        self.assertEqual(sp.expand(expected - sum(classes.values())), 0)
        for name, value in classes.items():
            self.assertNotEqual(sp.expand(value), 0, name)

    def test_cubic_vertex_and_no_three_meson_vertex(self):
        expression, _ = self.ccwz(self.legs[:3])
        self.assertEqual(sp.expand(expression - cubic_vertex(self.invariants, self.eps, self.legs[:3])), 0)
        mesons_only = [(leg[0], 0) for leg in self.legs[:3]]
        self.assertEqual(self.ccwz(mesons_only)[0], 0)

    def test_alp_mass_entry_and_mixing(self):
        eps, bare = sp.symbols('epsilon bare', positive=True)
        kinetic, mass_squared = solver_theory(self.invariants, eps, bare, DOWN_STRANGE_BASIS)
        self.assertEqual(sp.expand(mass_squared[-1, -1] - alp_mass_entry(self.invariants, eps, bare)), 0)
        source = matrix_C(self.invariants)
        for k, basis in enumerate(DOWN_STRANGE_BASIS):
            self.assertEqual(sp.simplify(mass_squared[k, -1] + eps * (source * basis).trace()), 0)
        self.assertEqual(kinetic, sp.diag(sp.Rational(1, 2), sp.Rational(1, 2), sp.Rational(1, 2), 1, 1, 1))

    def test_inverse_operator_and_normalization_derivative(self):
        generator = random.Random(5)
        target = random_hermitian(generator)
        s = sp.Rational(1, 3)
        self.assertEqual(sp.simplify(pole_operator(inverse_pole_operator(target, self.invariants, s),
                                                   self.invariants, s) - target), sp.zeros(3))
        variable = sp.Symbol('s')
        derivative = sp.diff(alp_inverse_propagator(self.invariants, self.eps, 1, variable), variable).subs(variable, s)
        self.assertEqual(sp.simplify(derivative - normalization(self.invariants, self.eps, s)), 0)


class ExactComparisonWithSolver(unittest.TestCase):
    """Closed formula equals the solver-based amplitude exactly at tuned rational points."""

    @classmethod
    def setUpClass(cls):
        cls.invariants = verification_inputs()
        cls.unmixed_mass = unmixed_neutral_squared_masses(cls.invariants)[0]
        cls.neutral = [compare_amplitude(cls.invariants, M_A2, cls.unmixed_mass + offset, S12, S23)
                       for offset in OFFSETS]
        cls.charged = compare_amplitude(cls.invariants, M_A2, cls.unmixed_mass + OFFSETS[1], S12, S23, 'charged')
        cls.down_strange = compare_amplitude(verification_inputs(down_strange=True), M_A2,
                                             cls.unmixed_mass + OFFSETS[1], S12, S23, 'down_strange')

    def test_unmixed_masses_rational(self):
        self.assertEqual(unmixed_neutral_squared_masses(self.invariants),
                         [sp.Rational(3, 50), sp.Rational(7, 25), sp.Rational(41, 50)])

    def check(self, result):
        self.assertTrue(result.eps_squared.is_Rational and result.eps_squared > 0)
        self.assertTrue(result.bare_alp_mass_squared.is_Rational)
        self.assertTrue(is_exactly_zero(result.closed - result.solver))
        self.assertTrue(result.rational_form.is_Rational, result.rational_form)
        self.assertTrue(is_exactly_zero(result.closed * sp.sqrt(result.normalization_product)
                                        - result.rational_form))

    def test_three_neutral_daughters_at_three_eps(self):
        self.assertEqual(len({result.eps_squared for result in self.neutral}), 3)
        for result in self.neutral:
            self.check(result)
            self.assertTrue(all(pole > 0 for pole in result.spectrum))
            for name, value in result.pieces.items():
                if name != 'kinetic':  # commuting diagonal matrices: no kinetic four-meson term
                    self.assertNotEqual(sp.N(value), 0, name)

    def test_charged_channel(self):
        self.check(self.charged)
        self.assertNotEqual(sp.N(self.charged.pieces['kinetic']), 0)
        self.assertNotEqual(sp.N(self.charged.pieces['exchange A n1 | n2 n3']), 0)

    def test_down_strange_channel(self):
        self.check(self.down_strange)
        self.assertNotEqual(sp.N(self.down_strange.pieces['exchange A n1 | n2 n3']), 0)

    def test_leg_normalization(self):
        result = self.neutral[0]
        eps = sp.sqrt(result.eps_squared)
        unmixed = unmixed_neutral_matrix(self.invariants, self.unmixed_mass)
        for leg in (alp_like_leg(self.invariants, eps, M_A2),
                    mixed_meson_leg(self.invariants, eps, result.daughter_squared_mass, unmixed)):
            self.assertEqual(sp.simplify((leg[0] * leg[0]).trace() / 2 + leg[1]**2), 1)

    def test_odd_in_eps(self):
        plus, minus = amplitudes_at_both_signs(self.invariants, M_A2, self.unmixed_mass + OFFSETS[1], S12, S23)
        self.assertTrue(is_exactly_zero(plus + minus))
        self.assertNotEqual(sp.N(plus), 0)


class ExpansionInEps(unittest.TestCase):
    """Series at fixed invariants, m_A^2, s_12 and s_23."""

    @classmethod
    def setUpClass(cls):
        cls.invariants = verification_inputs()
        cls.coefficients, _ = amplitude_series(cls.invariants, M_A2, 0, S12, S23)
        unmixed_mass = unmixed_neutral_squared_masses(cls.invariants)[0]
        unmixed = unmixed_neutral_matrix(cls.invariants, unmixed_mass)
        mandelstam = {(1, 2): S12, (2, 3): S23, (1, 3): M_A2 + 3 * unmixed_mass - S12 - S23}
        cls.leading = leading_coefficient(cls.invariants, M_A2, [unmixed] * 3, [unmixed_mass] * 3, mandelstam)
        cls.unmixed, cls.unmixed_mass = unmixed, unmixed_mass

    def test_only_odd_powers(self):
        for k in (0, 2, 4):
            self.assertEqual(sp.simplify(self.coefficients[k]), 0)
        self.assertNotEqual(sp.N(self.coefficients[3]), 0)

    def test_first_order_coefficient(self):
        self.assertTrue(is_exactly_zero(self.coefficients[1] - self.leading))
        # diagonal invariants: (B_0/f_pi^2) <(2 Y m_q - hat P_1) Phi_0^3>
        Y = inverse_pole_operator(matrix_C(self.invariants), self.invariants, M_A2)
        reduced = self.invariants['B_0'] / self.invariants['f_pi']**2 * (
            (2 * Y * self.invariants['S'][0] - self.invariants['P'][1]) * self.unmixed**3).trace()
        self.assertTrue(is_exactly_zero(reduced - self.leading))

    def test_remainder_approaches_third_order_coefficient(self):
        differences, eps_squared_values = [], []
        for offset in OFFSETS:
            result = compare_amplitude(self.invariants, M_A2, self.unmixed_mass + offset, S12, S23)
            eps = sp.sqrt(result.eps_squared)
            remainder = sp.N((result.closed - eps * self.leading) / eps**3 - self.coefficients[3], 30)
            differences.append(remainder)
            eps_squared_values.append(sp.N(result.eps_squared, 30))
        # the remainder minus the eps^3 coefficient is O(eps^2): ratio to eps^2 stays bounded and similar
        ratios = [d / e for d, e in zip(differences, eps_squared_values)]
        self.assertTrue(abs(differences[0]) < abs(differences[1]) < abs(differences[2]))
        self.assertLess(abs(ratios[0] - ratios[1]) / abs(ratios[0]), sp.Rational(1, 10))

    def test_overlap_argument_shift_is_third_order(self):
        series = overlap_argument_difference_series(self.invariants, M_A2, 0, S12, S23)
        for k in (0, 1, 2):
            self.assertEqual(sp.simplify(series[k]), 0)
        self.assertNotEqual(sp.N(series[3]), 0)

    def test_tuned_point_is_a_double_root(self):
        target = self.unmixed_mass + OFFSETS[0]
        eps_squared, bare = tuned_point(self.invariants, M_A2, target)
        eps = sp.sqrt(eps_squared)
        for pole in (M_A2, target):
            self.assertEqual(sp.simplify(alp_inverse_propagator(self.invariants, eps, bare, pole)), 0)


if __name__ == '__main__':
    unittest.main()
