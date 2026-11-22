"""Unit tests for derive_alp_dressed_vertex_power_counting.py.

Run: .venv/bin/python -m unittest discover -s scripts -p test_derive_alp_dressed_vertex_power_counting.py
All checks are exact (SymPy rationals, Gaussian rationals and square roots); the inputs are
test values only, not physical values.
"""
import itertools
import unittest

import sympy as sp

from derive_alp_dressed_vertex_power_counting import (
    ALP_LABEL, ALP_POLE_SQUARED_MASS, EPSILON, MESON_FIELD_NAMES, MOMENTUM_K, MOMENTUM_P, QUARK_CHARGE,
    SIXTH_ORDER_FIELDS, TEST_EPSILON, alp_pole_end_scaling, alp_pole_overlap, alp_two_meson_cubic_expected,
    alp_two_meson_cubic_table, closed_alp_line_full, closed_alp_line_potential,
    closed_alp_line_potential_matrix_form, closed_line_scalar_source_piece,
    closed_line_scalar_source_piece_down_strange, closed_line_scalar_source_piece_down_strange_pair,
    closed_line_single_meson, closed_unmixed_meson_line, commutator, commuting_class_a_form,
    down_strange_overlap_entry, first_order_scalar_source, joint_vertex_cubic_residuals,
    nonzero_alp_two_meson_pairs, without_down_strange_pseudoscalar,
    compare_with_builder, compare_with_direct_expansion, compare_with_joint_vertex, cubic_vertices,
    flavour_generators, homogeneity_with_joint_vertex, kaon_like_unmixed, kinetic_and_topological_structure,
    kinetic_closed_line_form, l1_mass_vertex, meson_pole_alp_component_series, normalized_contraction_residual,
    on_shell_closed_line_form, one_end_bracket_residual, potential_vertex_pieces, representative_odd_even_entries,
    rotated_scalar_source, rotated_source_prediction, single_meson_check_with_builder, single_meson_combination,
    tadpole_kinematics, test_invariants, topological_monomials_with_two_or_more_mesons, unit_matrix)
from derive_alp_finite_mixing_three_meson_tree import GENERATORS as GELL_MANN_GENERATORS
from derive_alp_finite_mixing_three_meson_tree import (
    matrix_C, pole_operator, unmixed_neutral_matrix, unmixed_neutral_squared_masses)

GENERATORS = flavour_generators()
REPRESENTATIVE_LEGS = (0, 2, 3, 7, 8, ALP_LABEL)      # E_uu, E_ss, lambda_1, lambda_6, lambda_7, ALP
DIAGONAL_NEUTRAL = (3, 8, 0)                          # lambda_3, lambda_8, lambda_0 (Gell-Mann indices)
DOWN_STRANGE_GENERATORS = (6, 7)                      # lambda_6, lambda_7


def off_diagonal_entries(matrix):
    return [matrix[i, j] for i in range(3) for j in range(3) if i != j]


class FirstOrderScalarSource(unittest.TestCase):
    """S_1 = -i[c_q^{VS}, m_q] with c_q^{VS} commuting with the quark charge."""

    def test_structure(self):
        for down_strange in (False, True):
            invariants = test_invariants(down_strange)
            quark_mass = invariants['S'][0]
            vector_scalar, S_1 = first_order_scalar_source(quark_mass, down_strange)
            self.assertEqual(invariants['S'][1], S_1)
            self.assertEqual(invariants['vector_scalar_invariant'], vector_scalar)
            self.assertEqual(commutator(vector_scalar, QUARK_CHARGE), sp.zeros(3))
            self.assertEqual(vector_scalar, vector_scalar.H)
            self.assertEqual(S_1, S_1.H)
            self.assertEqual([S_1[q, q] for q in range(3)], [0, 0, 0])
            self.assertEqual(commutator(S_1, QUARK_CHARGE), sp.zeros(3))
            expected = sp.zeros(3)
            expected[1, 2] = sp.expand(-sp.I * vector_scalar[1, 2] * (quark_mass[2, 2] - quark_mass[1, 1]))
            expected[2, 1] = sp.conjugate(expected[1, 2])
            self.assertEqual(S_1, expected)
            if down_strange:
                self.assertNotEqual(sp.re(S_1[1, 2]), 0)
                self.assertNotEqual(sp.im(S_1[1, 2]), 0)
            else:
                self.assertEqual(S_1, sp.zeros(3))



class InvariantVertexFormula(unittest.TestCase):
    """(a) The closed invariant formula against the existing builders and the L1 code."""

    def test_builder_every_monomial_up_to_fourth_order(self):
        compared, failures, unmatched = compare_with_builder()
        self.assertEqual(compared, 1000)
        self.assertEqual(failures, [])
        self.assertEqual(unmatched, [])

    def test_single_meson_vertex_combination_against_builder(self):
        residuals = single_meson_check_with_builder(max_alps=3)
        self.assertEqual(len(residuals), 4 * 9)
        self.assertEqual(set(residuals.values()), {0})

    def test_single_meson_vertex_is_c_at_first_order(self):
        for invariants in (test_invariants(False), test_invariants(True)):
            C = matrix_C(invariants)
            for generator in GENERATORS:
                self.assertEqual(single_meson_combination(1, generator, invariants, EPSILON),
                                 sp.expand(EPSILON * (generator * C).trace()))
                for n_a in range(4):
                    pieces = potential_vertex_pieces(n_a, [generator], invariants, EPSILON)
                    self.assertEqual(sp.expand(pieces['mass'] + pieces['topological']),
                                     single_meson_combination(n_a, generator, invariants, EPSILON))

    def test_down_strange_invariants_against_joint_vertex_code(self):
        compared, failures = compare_with_joint_vertex()
        self.assertGreater(compared, 100)
        self.assertEqual(failures, [])

    def test_mass_piece_equals_l1_u_form(self):
        invariants = test_invariants(True, True)
        extra = [sp.sqrt(2) * unit_matrix(0, 1), sp.sqrt(2) * unit_matrix(2, 1)]   # non-Hermitian legs
        pool = GENERATORS + extra
        cases = [(n_a, combination) for n_a in range(4) for n_p in range(1, 5)
                 for combination in itertools.combinations_with_replacement((1, 3, 7, 8, 9, 10), n_p)
                 if n_a + n_p <= 5][::7]
        self.assertGreater(len(cases), 20)
        for n_a, combination in cases:
            legs = [pool[i] for i in combination]
            pieces = potential_vertex_pieces(n_a, legs, invariants, EPSILON)
            self.assertEqual(sp.simplify(pieces['mass'] - l1_mass_vertex(n_a, legs, invariants, EPSILON)), 0)

    def test_direct_expansion_all_fields_down_strange_fourth_order(self):
        compared, non_homogeneous, failures = compare_with_direct_expansion(test_invariants(True), MESON_FIELD_NAMES, 4)
        self.assertEqual(compared, 1000)
        self.assertEqual(non_homogeneous, [])
        self.assertEqual(failures, [])


class PowersOfEpsilon(unittest.TestCase):
    """(b) Homogeneity of degree n_a in eps at fixed f."""

    def test_joint_vertex_code_with_symbolic_epsilon(self):
        cases = [(1, (7, 7)), (2, (1, 8)), (1, (3, 4, 7)), (3, (7,)), (4, ()), (2, (0, 7)), (0, (7, 8, 1, 1))]
        leftovers = homogeneity_with_joint_vertex(cases)
        self.assertEqual(set(leftovers.values()), {0})

    def test_formula_is_homogeneous(self):
        invariants = test_invariants(True, True)
        for n_a in range(5):
            for combination in itertools.combinations_with_replacement((0, 3, 7, 8), 2):
                value = sum(potential_vertex_pieces(n_a, [GENERATORS[i] for i in combination], invariants,
                                                    EPSILON, include_bare_alp_mass=False).values())
                self.assertEqual(sp.diff(sp.expand(value / EPSILON**n_a), EPSILON), 0)


class KineticAndTopologicalTerms(unittest.TestCase):
    """(c) No ALP vertex from the kinetic term; topological term at most quadratic in mesons."""

    def test_builder_structure(self):
        kinetic_alp, topological_extra = kinetic_and_topological_structure()
        self.assertEqual(kinetic_alp, [])
        self.assertEqual(topological_extra, [])

    def test_topological_term_to_sixth_order(self):
        extra = topological_monomials_with_two_or_more_mesons(test_invariants(True, True), SIXTH_ORDER_FIELDS, 6)
        self.assertEqual(extra, [])


class SixthOrderVertices(unittest.TestCase):
    """(d) n_a + n_P <= 6 against a direct Taylor expansion (d-s invariants, d-s pair, P_1, neutrals)."""

    def test_sixth_order(self):
        compared, non_homogeneous, failures = compare_with_direct_expansion(
            test_invariants(True, True), SIXTH_ORDER_FIELDS, 6)
        self.assertEqual(compared, 1715)
        self.assertEqual(non_homogeneous, [])
        self.assertEqual(failures, [])


class CubicVertices(unittest.TestCase):
    """(e) Cubic vertices exist for general invariants, not in the c_q^A, c_q^P representative."""

    def test_cubic_vertices_nonzero(self):
        for invariants in (test_invariants(False), test_invariants(True)):
            cubic = cubic_vertices(invariants, EPSILON)
            self.assertTrue(any(v != 0 for v in cubic['Paa'].values()))
            self.assertNotEqual(cubic['aaa'], 0)
        # V_PPa comes from S_1 only: zero for S_1 = 0, nonzero for (E_dd, lambda_6) with d-s S_1;
        # zero for E_uu^2 and lambda_6^2 (diagonal anticommutator, zero-diagonal S_1)
        self.assertEqual(set(cubic_vertices(test_invariants(False), EPSILON)['PPa'].values()), {0})
        down_strange = cubic_vertices(test_invariants(True), EPSILON)['PPa']
        self.assertNotEqual(down_strange[(1, 7)], 0)
        self.assertEqual(down_strange[(0, 0)], 0)
        self.assertEqual(down_strange[(7, 7)], 0)
        invariants = test_invariants(False)
        f, B_0 = invariants['f_pi'], invariants['B_0']
        m0 = invariants['anomaly_mass_squared'] / 3
        self.assertEqual(cubic_vertices(invariants, EPSILON)['aaa'],
                         sp.expand(EPSILON**3 / f * (B_0 * invariants['S'][3].trace()
                                                     - 6 * m0 * invariants['theta'][1] * invariants['theta'][2])))

    def test_alp_two_meson_pairs(self):
        invariants = test_invariants(True)
        expected_pairs = sorted([(diagonal, d_s) if diagonal < d_s else (d_s, diagonal)
                                 for diagonal in DIAGONAL_NEUTRAL for d_s in DOWN_STRANGE_GENERATORS]
                                + [(i, j) for i in (1, 2) for j in (4, 5)])
        self.assertEqual(len(expected_pairs), 10)
        self.assertEqual(nonzero_alp_two_meson_pairs(invariants), expected_pairs)
        self.assertEqual(nonzero_alp_two_meson_pairs(test_invariants(False)), [])
        table = alp_two_meson_cubic_table(invariants, EPSILON)
        self.assertEqual(len(table), 45)
        for (i, j), value in table.items():
            self.assertEqual(sp.simplify(value - alp_two_meson_cubic_expected(
                GELL_MANN_GENERATORS[i], GELL_MANN_GENERATORS[j], invariants, EPSILON)), 0)

    def test_alp_two_meson_pairs_against_joint_vertex_code(self):
        residuals = joint_vertex_cubic_residuals(test_invariants(True))
        self.assertEqual(len(residuals), 45)
        self.assertEqual(set(residuals.values()), {0})

    def test_representative_has_no_odd_vertices(self):
        for n, (scalar, pseudoscalar) in representative_odd_even_entries().items():
            if n % 2:
                self.assertEqual(scalar, sp.zeros(3))
                self.assertNotEqual(pseudoscalar, sp.zeros(3))
            else:
                self.assertEqual(pseudoscalar, sp.zeros(3))
                self.assertNotEqual(scalar, sp.zeros(3))


class AlpLikePoleEnds(unittest.TestCase):
    """(f) Each end on the ALP-like pole costs eps / sqrt(N)."""

    def test_one_end(self):
        invariants = test_invariants(True)
        sets = [others for size in range(2, 4)                     # interaction vertices: >= 3 legs
                for others in itertools.combinations_with_replacement(REPRESENTATIVE_LEGS, size)]
        result = alp_pole_end_scaling(invariants, sets, 1)
        self.assertEqual({derivative for _, derivative in result.values()}, {0})
        self.assertTrue(any(value != 0 for value, _ in result.values()))

    def test_two_ends(self):
        invariants = test_invariants(True)
        sets = [others for size in range(1, 3)                     # interaction vertices: >= 3 legs
                for others in itertools.combinations_with_replacement(REPRESENTATIVE_LEGS, size)]
        result = alp_pole_end_scaling(invariants, sets, 2)
        self.assertEqual({derivative for _, derivative in result.values()}, {0})

    def test_bracket_and_normalization(self):
        invariants = test_invariants(True)
        for others in ((3, 7), (7, ALP_LABEL), (0, 8, 8)):
            self.assertEqual(one_end_bracket_residual(invariants, others), 0)
        for others, ends in (((3, 7), 1), ((7,), 2), ((ALP_LABEL, 0), 2)):
            self.assertEqual(normalized_contraction_residual(invariants, others, ends), 0)


class ClosedAlpLine(unittest.TestCase):
    """(f) Closed ALP-like line at the quartic vertex with two meson legs."""

    def test_commuting_class_a_combination(self):
        invariants = test_invariants(False)
        Y = alp_pole_overlap(invariants)
        lam_3 = sp.diag(1, -1, 0)
        for generator in GENERATORS[:3] + [lam_3]:
            value = closed_alp_line_potential(generator, generator, invariants)
            self.assertEqual(sp.simplify(value - commuting_class_a_form(generator, invariants, Y)), 0)
        for a, b in itertools.combinations(GENERATORS[:3], 2):
            value = closed_alp_line_potential(a, b, invariants)
            self.assertEqual(sp.simplify(value - rotated_source_prediction(a, b, invariants, Y)), 0)

    def test_rotated_source_fails_and_matrix_form_holds(self):
        invariants = test_invariants(True)
        Y = alp_pole_overlap(invariants)
        rotated_failures = 0
        for i, j in itertools.combinations_with_replacement(range(9), 2):
            a, b = GENERATORS[i], GENERATORS[j]
            value = closed_alp_line_potential(a, b, invariants)
            self.assertEqual(sp.simplify(value - closed_alp_line_potential_matrix_form(a, b, invariants, Y)), 0)
            if sp.simplify(value - rotated_source_prediction(a, b, invariants, Y)) != 0:
                rotated_failures += 1
        self.assertGreater(rotated_failures, 0)

    def test_matrix_form_reduces_to_rotated_source_for_commuting(self):
        invariants = test_invariants(False)
        Y = alp_pole_overlap(invariants)
        S2_prime = rotated_scalar_source(invariants, Y)
        S, P = invariants['S'], invariants['P']
        self.assertEqual(S2_prime - (S[2] + 4 * Y * P[1] - 4 * Y * Y * S[0]), sp.zeros(3))

    def test_kinetic_part_and_on_shell_form(self):
        invariants = test_invariants(True)
        Y = alp_pole_overlap(invariants)
        for i, j in ((1, 7), (3, 8), (7, 7)):
            a, b = GENERATORS[i], GENERATORS[j]
            full = closed_alp_line_full(a, b, invariants)
            potential = full.subs({MOMENTUM_P: 0, MOMENTUM_K: 0})
            self.assertEqual(sp.simplify(potential - closed_alp_line_potential(a, b, invariants)), 0)
            self.assertEqual(sp.simplify(full - potential - kinetic_closed_line_form(a, b, invariants, Y)), 0)
            self.assertEqual(sp.simplify(tadpole_kinematics(full) - on_shell_closed_line_form(a, b, invariants, Y)), 0)

    def test_on_shell_charged_meson_gives_class_a(self):
        for invariants in (test_invariants(False), test_invariants(True)):
            Y = alp_pole_overlap(invariants)
            f, B_0, m = invariants['f_pi'], invariants['B_0'], invariants['S'][0]
            generator = GENERATORS[3]                                   # lambda_1, mass B_0 (m_u + m_d)
            on_shell = on_shell_closed_line_form(generator, generator, invariants, Y).subs(
                MOMENTUM_P**2, B_0 * (m[0, 0] + m[1, 1]))
            class_a = (-B_0 / f**2 * (generator**2 * rotated_scalar_source(invariants, Y)).trace()
                       + ALP_POLE_SQUARED_MASS / f**2 * (commutator(Y, generator)**2).trace())
            self.assertEqual(sp.simplify(on_shell - class_a), 0)
        # Y solves M_{m_A^2}(Y) = C exactly (input of the on-shell reduction)
        invariants = test_invariants(True)
        self.assertEqual((pole_operator(alp_pole_overlap(invariants), invariants, ALP_POLE_SQUARED_MASS)
                          - matrix_C(invariants)).applyfunc(sp.simplify), sp.zeros(3))

    def test_single_meson_closed_line_rotated_pseudoscalar(self):
        invariants = test_invariants(True)
        for generator in GENERATORS + list(GELL_MANN_GENERATORS):
            value, rotated = closed_line_single_meson(generator, invariants)
            self.assertEqual(sp.simplify(value - rotated), 0)

    def test_scalar_source_piece_diagonal_neutral(self):
        both = test_invariants(True)
        diagonal_pseudoscalar = without_down_strange_pseudoscalar(both)
        self.assertEqual(off_diagonal_entries(diagonal_pseudoscalar['P'][1]), [0] * 6)
        self.assertNotEqual(off_diagonal_entries(diagonal_pseudoscalar['S'][1]), [0] * 6)
        for invariants, nonzero in ((both, True), (diagonal_pseudoscalar, False), (test_invariants(False), False)):
            Y = alp_pole_overlap(invariants)
            self.assertEqual(sp.simplify(Y[1, 2] - down_strange_overlap_entry(invariants)), 0)
            for index in DIAGONAL_NEUTRAL:
                generator = GELL_MANN_GENERATORS[index]
                piece = closed_line_scalar_source_piece(generator, invariants, Y)
                self.assertEqual(sp.simplify(piece - closed_line_scalar_source_piece_down_strange(
                    generator, invariants, Y)), 0)
                if nonzero:
                    self.assertNotEqual(sp.simplify(piece), 0)
                else:
                    self.assertEqual(sp.simplify(piece), 0)
        self.assertEqual(alp_pole_overlap(diagonal_pseudoscalar)[1, 2], 0)
        self.assertNotEqual(alp_pole_overlap(both)[1, 2], 0)

    def test_scalar_source_piece_down_strange_generators(self):
        for invariants in (test_invariants(True), without_down_strange_pseudoscalar(test_invariants(True))):
            Y = alp_pole_overlap(invariants)
            for index in DOWN_STRANGE_GENERATORS:
                generator = GELL_MANN_GENERATORS[index]
                piece = closed_line_scalar_source_piece(generator, invariants, Y)
                self.assertNotEqual(sp.simplify(piece), 0)
                self.assertEqual(sp.simplify(piece - closed_line_scalar_source_piece_down_strange_pair(
                    generator, invariants, Y)), 0)


class ClosedUnmixedMesonLine(unittest.TestCase):
    """Closed unmixed meson line on the ALP-two-meson vertex: zero with zero-diagonal S_1."""

    def test_flavour_matrix_is_diagonal_and_trace_vanishes(self):
        invariants = test_invariants(True)
        flavour, trace, eigen_traces = closed_unmixed_meson_line(invariants)
        self.assertEqual(off_diagonal_entries(flavour), [0] * 6)
        self.assertNotEqual([flavour[q, q] for q in range(3)], [0, 0, 0])
        self.assertEqual(trace, 0)
        self.assertEqual(len(eigen_traces), 9)
        self.assertEqual(set(eigen_traces), {0})


class MesonLikePoles(unittest.TestCase):
    """(g) ALP component of the meson-like poles is eps <Phi_{n,0} C>/(m_A^2 - m_{n,0}^2) + O(eps^3)."""

    def test_neutral_poles(self):
        invariants = test_invariants(False)
        for mass in unmixed_neutral_squared_masses(invariants):
            matrix = unmixed_neutral_matrix(invariants, mass)
            coefficients, expected = meson_pole_alp_component_series(invariants, mass, matrix)
            self.assertEqual(coefficients[0], 0)
            self.assertEqual(coefficients[2], 0)
            self.assertNotEqual(expected, 0)
            self.assertEqual(sp.radsimp(coefficients[1] - expected), 0)

    def test_kaon_like_pole_down_strange(self):
        invariants = test_invariants(True)
        mass, matrix = kaon_like_unmixed(invariants)
        coefficients, expected = meson_pole_alp_component_series(invariants, mass, matrix)
        self.assertEqual(coefficients[0], 0)
        self.assertNotEqual(expected, 0)
        self.assertEqual(sp.simplify(coefficients[1] - expected), 0)


if __name__ == '__main__':
    unittest.main()
