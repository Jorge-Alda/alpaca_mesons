"""Tests for the generic finite-mixing a_phys -> P_i P_j P_k amplitude.

Run: .venv/bin/python -m unittest discover -s scripts -p test_derive_alp_finite_mixing_three_pseudoscalars.py
"""
import itertools
import unittest

import sympy as sp

from derive_alp_finite_mixing_three_pseudoscalars import (
    B_0, anticommutator, component_quadratic_blocks, derive, f_a, f_pi,
    field_basis_directions, generic_amplitude, hermitian_unit_basis,
    linear_source_ccwz_derivatives,
    aligned_weak_residuals, numerical_route_comparison, permuted_inputs, record_formula_checks,
    identical_daughter_residuals,
    three_pion_substitution_residuals, transformation_checks, weak_mixing_residuals)
from symbolic_utils import u3_generators

NUMBERS = {f_pi: 1, f_a: 3, B_0: 2}


def exact_setup(generators=None):
    """Exact rational sources, quadratic blocks and field basis (no physical meaning)."""
    quark_masses = sp.diag(1, 2, 5)
    axial = sp.Matrix([[1, 1, 0], [1, 0, 2], [0, 2, -1]])
    pseudoscalar = sp.Matrix([[0, 1, 1], [1, 2, 0], [1, 0, 1]])
    vector = sp.Matrix([[1, 0, 1], [0, -1, 1], [1, 1, 0]])
    scalar = sp.Matrix([[0, 2, 0], [2, 1, 1], [0, 1, 2]])
    plus, minus = linear_source_ccwz_derivatives(axial, pseudoscalar, vector, scalar, quark_masses)
    plus = [sp.Matrix(matrix).subs(NUMBERS) for matrix in plus]
    minus = [sp.Matrix(matrix).subs(NUMBERS) for matrix in minus]
    coupling = (-B_0*anticommutator(quark_masses, axial - pseudoscalar) - 2*sp.eye(3)).subs(NUMBERS)
    generators = hermitian_unit_basis() if generators is None else generators
    kinetic, mass_squared = component_quadratic_blocks(generators, quark_masses, sp.Rational(9, 2),
                                                       coupling, sp.Rational(17, 2))
    return plus, minus, field_basis_directions(generators), kinetic, mass_squared.subs(NUMBERS)


def directions():
    """Three different daughter directions with nonzero ALP components and a parent."""
    parent = (sp.Matrix([[1, 2, 0], [2, -1, 1], [0, 1, 3]])/7, 1)
    first = (sp.Matrix([[1, 0, 0], [0, -1, 0], [0, 0, 0]]), sp.Rational(1, 5))
    second = (sp.Matrix([[0, 1, 1], [1, 0, 0], [1, 0, 0]]), sp.Rational(-1, 3))
    third = (sp.Matrix([[0, 0, 0], [0, 1, 2], [0, 2, -1]]), sp.Rational(2, 7))
    return parent, [first, second, third]


class GenericAmplitudeStructure(unittest.TestCase):
    def setUp(self):
        self.plus, self.minus, self.basis, self.kinetic, self.mass_squared = exact_setup()
        self.parent, self.daughters = directions()
        self.parent_mass_sq = sp.Rational(40)
        self.masses = {1: sp.Rational(3), 2: sp.Rational(5), 3: sp.Rational(7)}
        s12, s13 = sp.Rational(19, 2), sp.Rational(27, 2)
        self.invariants = {frozenset({1, 2}): s12, frozenset({1, 3}): s13,
                           frozenset({2, 3}): self.parent_mass_sq + sum(self.masses.values()) - s12 - s13}

    def amplitude(self, daughters, masses, invariants):
        total, contact, exchanges = generic_amplitude(
            self.parent, daughters, self.parent_mass_sq, masses, invariants,
            self.plus, self.minus, self.basis, self.kinetic, self.mass_squared)
        return sp.nsimplify(sp.expand(total.subs(NUMBERS))), contact, exchanges

    def test_daughter_permutation_symmetry(self):
        reference, contact, exchanges = self.amplitude(self.daughters, self.masses, self.invariants)
        self.assertNotEqual(reference, 0)
        self.assertTrue(all(sp.expand(value.subs(NUMBERS)) != 0 for value in exchanges.values()))
        for order in itertools.permutations((1, 2, 3)):
            permuted = permuted_inputs(self.daughters, self.masses, self.invariants, order)
            value, _, _ = self.amplitude(*permuted)
            self.assertEqual(sp.simplify(value - reference), 0, order)

    def test_internal_basis_independence(self):
        reference = self.amplitude(self.daughters, self.masses, self.invariants)[0]
        _, _, basis, kinetic, mass_squared = exact_setup(u3_generators())
        self.assertEqual(kinetic, sp.eye(10))
        total = generic_amplitude(self.parent, self.daughters, self.parent_mass_sq, self.masses,
                                  self.invariants, self.plus, self.minus, basis, kinetic,
                                  mass_squared)[0]
        self.assertEqual(sp.simplify(sp.expand(total.subs(NUMBERS)) - reference), 0)

    def test_daughter_alp_component_enters_contact_and_exchange(self):
        component = sp.Symbol('x_1')
        daughters = [(self.daughters[0][0], component)] + self.daughters[1:]
        _, contact, exchanges = self.amplitude(daughters, self.masses, self.invariants)
        self.assertNotEqual(sp.diff(contact.subs(NUMBERS), component), 0)
        self.assertNotEqual(sp.diff(sum(exchanges.values()).subs(NUMBERS), component), 0)

    def test_identical_daughters_use_same_column_without_factorial(self):
        same = [self.daughters[0]]*3
        masses = {r: self.masses[1] for r in (1, 2, 3)}
        s12, s13 = sp.Rational(23, 2), sp.Rational(25, 2)
        invariants = {frozenset({1, 2}): s12, frozenset({1, 3}): s13,
                      frozenset({2, 3}): self.parent_mass_sq + 3*self.masses[1] - s12 - s13}
        value, _, _ = self.amplitude(same, masses, invariants)
        for order in itertools.permutations((1, 2, 3)):
            permuted = self.amplitude(*permuted_inputs(same, masses, invariants, order))[0]
            self.assertEqual(sp.simplify(permuted - value), 0)


class DerivationPieces(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.derived = derive()

    def test_kinetic_contact_reduction(self):
        self.assertEqual(sp.expand(self.derived['kinetic_contact_summed']
                                   - self.derived['kinetic_contact_reduced']), 0)

    def test_power_counting(self):
        index_level = self.derived['index_level']
        self.assertEqual(index_level['series'][0], 0)
        self.assertEqual(index_level['series'][2], 0)
        for order in range(3):
            self.assertEqual(index_level['exchange_series'][order], 0)
        self.assertEqual(sp.expand(index_level['series'][1] - index_level['expected_leading']), 0)
        self.assertNotEqual(index_level['series'][1], 0)

    def test_daughter_alp_components_in_exact_index_amplitude(self):
        index_level = self.derived['index_level']
        xi = sp.Symbol('xi_1')
        self.assertNotEqual(sp.diff(index_level['contact'], xi), 0)
        self.assertNotEqual(sp.diff(index_level['exchange'], xi), 0)

    def test_symmetry_factors(self):
        self.assertEqual(self.derived['symmetry_factors'],
                         {'iii': (6, 6), 'iij': (2, 2), 'ijk': (1, 1)})


class Reproduction(unittest.TestCase):
    def test_three_pion_record_section_4(self):
        residuals = three_pion_substitution_residuals()
        self.assertEqual(residuals['000'][0], 0)
        self.assertEqual(residuals['0+-'][0], 0)
        # Every exchange is nonzero except (a,pi0)|(pi+,pi-), which vanishes for
        # linear sources because diag([c_q^VS, m_q]) = 0.
        self.assertTrue(all(residuals['000'][1].values()))
        self.assertEqual(residuals['0+-'][1], {(1, 2, 3): False, (2, 1, 3): True, (3, 1, 2): True})

    def test_weak_mixing_against_small_mixing_forms(self):
        residuals = weak_mixing_residuals()
        for name, residual in residuals.items():
            if 'nonzero' in name:
                self.assertIs(residual, True, name)
            else:
                self.assertEqual(residual, 0, name)


class PoleAlignedRoute(unittest.TestCase):
    """Route 2: transformation of the notes, transformed vertices, comparison with route 1."""

    @classmethod
    def setUpClass(cls):
        cls.checks = transformation_checks()

    def test_mixing_elimination(self):
        self.assertTrue(self.checks['aligned blocks = T^T (K, M^2) T (congruence of the single-pole record)'])
        self.assertTrue(self.checks["meson-ALP entry of m_a^2 K' - M'^2 vanishes at the parent pole"])
        self.assertTrue(self.checks["K'_aa = N_a = 1 + 2 (f_pi/f_a)^2 <Y^2>"])
        self.assertTrue(self.checks['derivative-free blocks = three-pion record section 2'])

    def test_aligned_sources(self):
        self.assertTrue(self.checks['first aligned derivatives: d chi_+ unchanged, d chi_- = 4iB_0{c_P_hat,m_q}/f_a'])
        self.assertTrue(self.checks["aligned d_a^n chi_+- = 4B_0 S_n'/f_a^n, 4iB_0 P_n'/f_a^n (Taylor law), orders 0-4"])
        self.assertTrue(self.checks['dressed first order = first-order record: S_1 = -i[c_VS,m_q], P_1 = -{c_AP,m_q}'])
        self.assertTrue(self.checks['dressed second order = second-order record (linear sources)'])
        self.assertTrue(self.checks["S_1' = widehat S_1, P_1' = widehat P_1 + {c_A_hat,m_q}"])
        self.assertTrue(self.checks["S_2', P_2' explicit forms; second-order closed form maps them back to widehat S_2, P_2"])

    def test_chi_plus_contact_piece(self):
        self.assertTrue(self.checks['chi_+ contact (momentum-independent part) = V4_chi with aligned derivatives'])

    def test_full_contact_and_cubic(self):
        self.assertTrue(self.checks['full aligned contact = joint CCWZ expansion'])
        self.assertTrue(self.checks['aligned cubic vertex = derivative-free cubic on shifted directions'])

    def test_weak_limit_without_on_shell_condition(self):
        self.assertEqual(aligned_weak_residuals(), 0)

    def test_routes_agree_on_shell_only(self):
        route_one, route_two, contact_difference, columns = numerical_route_comparison()
        self.assertLess(abs(route_one - route_two), 1e-10*abs(route_one))
        self.assertLess(abs(contact_difference), 1e-10*abs(route_one))
        self.assertTrue(all(columns.values()), columns)
        self.assertGreater(abs(numerical_route_comparison(off_shell_shift=0.37)[2]), 1e-6)


class CompactRecordFormulas(unittest.TestCase):
    """The record's compact quartic, cubic and amplitude formulas equal the vertex functions."""

    def test_compact_quartic_cubic_amplitude(self):
        results = record_formula_checks()
        self.assertEqual(len(results), 8)
        self.assertTrue(all(results.values()), [name for name, ok in results.items() if not ok])

    def test_identical_daughter_formulas(self):
        results = identical_daughter_residuals()
        self.assertEqual(len(results), 4)
        self.assertTrue(all(results.values()), [name for name, ok in results.items() if not ok])


if __name__ == '__main__':
    unittest.main()
