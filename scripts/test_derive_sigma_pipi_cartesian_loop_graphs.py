"""Unit tests for derive_sigma_pipi_cartesian_loop_graphs.py (exact symbolic identities only; no numerical matrices).

Run: PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_cartesian_loop_graphs
"""
import unittest

import sympy as sp

import derive_sigma_pipi_cartesian_loop_graphs as graphs
from sigma_pipi_loop_algebra import EPSILON, M_A, M_PI, MOMENTUM_SQUARED, TraceSum, assert_no_float


def labels(rows):
    return [row['label'] for row in rows]


class EnumerationTest(unittest.TestCase):
    def test_counts(self):
        rows = graphs.table('one_particle_irreducible')
        self.assertEqual(len(rows), 48)
        counts = graphs.class_counts(rows)
        self.assertEqual({name: entry['graphs'] for name, entry in counts.items()},
                         {'tadpole': 3, 'bubble': 18, 'triangle': 27})
        self.assertEqual({name: entry['kept'] for name, entry in counts.items()},
                         {'tadpole': 3, 'bubble': 12, 'triangle': 11})
        self.assertEqual({name: entry['vanishing identically'] for name, entry in counts.items()},
                         {'tadpole': 0, 'bubble': 6, 'triangle': 16})
        for group in ('sigma', 'pi'):
            self.assertEqual(len(graphs.table(group)), 9)
        self.assertEqual(labels(graphs.kept(graphs.table('sigma'))), [f'Csigma{n}' for n in range(4, 10)])
        self.assertEqual(labels(graphs.kept(graphs.table('pi'))), ['Cpi2', 'Cpi3', 'Cpi7', 'Cpi8', 'Cpi9'])

    def test_no_box_and_no_two_point_vertex(self):
        result = graphs.no_box_and_no_two_point_vertex()
        self.assertEqual(result['box graphs'], 0)
        self.assertEqual(result['distributions with four vertices'], 0)
        self.assertEqual(result['(2,5) graphs'], 0)
        self.assertGreater(result['configurations with a two-leg vertex (not vertices)'], 0)
        for group in graphs.GROUPS:
            for row in graphs.table(group):
                self.assertTrue(all(count >= 3 for count in graphs.vertex_leg_counts(row)))

    def test_symmetry_factors(self):
        """S = 1/2 for a tadpole and for a bubble with two lines of the same type, else 1."""
        for group in graphs.GROUPS:
            for row in graphs.table(group):
                r = len(row['blocks'])
                expected = sp.Rational(1, 2) if r == 1 or (r == 2 and row['lines'][0] == row['lines'][1]) else 1
                self.assertEqual(row['symmetry'], expected, row['label'])

    def test_vanishing_reasons(self):
        """Every identically vanishing graph has a vertex with one charged pion and no generator end,
        which vanishes by the charge rule (its trace has charge +-1)."""
        for group in graphs.GROUPS:
            for row in graphs.table(group):
                if row['status'] != 'vanishes identically':
                    continue
                self.assertTrue(row['reasons'])
                for reason in row['reasons']:
                    self.assertTrue(reason.startswith('vertex '), reason)
                    self.assertNotIn('generator', reason)

    def test_external_kinematics(self):
        """p_sigma + p_1 + p_2 = 0 with the exact products."""
        dot = graphs.loop_dot(M_A)
        total = {'ps': 1, 'p1': 1, 'p2': 1}
        for label in ('ps', 'p1', 'p2'):
            self.assertEqual(sp.expand(dot(total, {label: 1})), 0)


class LinePieceTest(unittest.TestCase):
    def test_four_piece_split_equals_residue(self):
        """(i)+(ii)+(iii)+(iv) = (1_9 + 1) - v_a v_a^T on every meson line, every graph, every routing."""
        for group in graphs.GROUPS:
            for row in graphs.table(group):
                if 'M' not in row['lines']:
                    continue
                three = [graphs.THREE_PIECE_RESIDUE if line == 'M' else graphs.pieces_of(line) for line in row['lines']]
                for line, orientation in graphs.routings(row):
                    four_value = graphs.routing_value(row, group, line, orientation)
                    three_value = graphs.routing_value(row, group, line, orientation, piece_lists=three)
                    self.assertTrue((four_value - three_value).is_zero(), (row['label'], line, orientation))

    def test_alp_like_line_plus_remainder_is_bare_p_line(self):
        """With the common scale m_a and a common routing: ALP-like line + pieces (ii)+(iii)+(iv)
        of a meson line = bare p line (per-term scales of the other lines kept)."""
        for group in graphs.GROUPS:
            for row in graphs.table(group):
                for position, line_type in enumerate(row['lines']):
                    if line_type != 'A':
                        continue
                    unit = {**row, 'symmetry': 1}

                    def lists(replacement):
                        return [replacement if j == position else graphs.pieces_of(t) for j, t in enumerate(row['lines'])]
                    for line, orientation in graphs.routings(row):
                        alp_like = graphs.routing_value(unit, group, line, orientation, piece_lists=lists((graphs.ALP_LIKE_LINE,)))
                        remainder = graphs.routing_value(unit, group, line, orientation,
                                                         piece_lists=lists(graphs.ALP_LIKE_REMAINDER))
                        bare = graphs.routing_value(unit, group, line, orientation, piece_lists=lists((graphs.BARE_P_LINE,)))
                        self.assertTrue((alp_like + remainder - bare).is_zero(), (row['label'], position))

    def test_closed_alp_like_and_meson_lines_on_pion_leg(self):
        """dZ(Cpi2) + dZ(Cpi3) = 0 exactly: with the common scale m_a they form the bare p line,
        whose vertex (pi+, pi-, p, p) has no p^2 dependence."""
        total = graphs.row_by_label('Cpi2')['routing values'][0] + graphs.row_by_label('Cpi3')['routing values'][0]
        self.assertTrue(total.is_zero())


class RoutingTest(unittest.TestCase):
    def test_graphs_without_kinetic_vertex_are_routing_independent(self):
        found = []
        for group in graphs.GROUPS:
            for row in graphs.kept(graphs.table(group)):
                values = row['routing values']
                independent = all((value - values[0]).is_zero() for value in values)
                if not graphs.can_contain_kinetic_vertex(row, group):
                    self.assertTrue(independent, row['label'])
                if not independent:
                    found.append(row['label'])
        self.assertEqual(found, ['C12', 'C13', 'C15'])

    def test_number_of_routings(self):
        for group in graphs.GROUPS:
            for row in graphs.kept(graphs.table(group)):
                self.assertEqual(len(row['routing values']), 2 * len(row['blocks']))


class D4Test(unittest.TestCase):
    def test_excluded_terms(self):
        excluded = [row['label'] for group in graphs.GROUPS for row in graphs.table(group)
                    if row['D4 status'] == 'excluded (D4)']
        self.assertEqual(excluded, ['C3', 'C15', 'Cpi3'])
        zero = [row['label'] for group in graphs.GROUPS for row in graphs.table(group)
                if row['D4 status'] == 'all-generator-sum term is zero']
        self.assertEqual(zero, ['C9', 'C21', 'C48', 'Csigma3', 'Cpi9'])

    def test_criterion_and_exact_removal(self):
        for label in ('C3', 'C15', 'Cpi3'):
            row = graphs.row_by_label(label)
            group = row['group']
            full = graphs.routing_values(row, group, 'all')
            lambda_kept = graphs.routing_values(row, group, 'D4', eliminate_lambda=False)
            for d4, d4_lambda, kept, total in zip(row['D4 term'], lambda_kept, row['routing values'], full):
                self.assertFalse(d4.is_zero())
                self.assertTrue(graphs.satisfies_d4_criterion(d4_lambda, group))
                self.assertTrue((d4 - d4_lambda).is_zero())
                self.assertTrue((kept + d4 - total).is_zero())

    def test_criterion_is_not_empty(self):
        """The all-generator-sum term of the sigma self-energy Csigma9 has two s ends (eps^2): kept."""
        row = graphs.row_by_label('Csigma9')
        self.assertEqual(row['D4 status'], 'all-generator-sum term kept (criterion fails)')
        d4 = graphs.routing_values(row, 'sigma', 'D4')
        self.assertFalse(d4[0].is_zero())
        for value in d4[0].terms.values():
            self.assertEqual(sp.Poly(sp.expand(value), EPSILON).monoms(), [(2,)])
        self.assertFalse(graphs.satisfies_d4_criterion(d4[0], 'sigma'))


class QuarticPotentialTest(unittest.TestCase):
    def test_quartic_vertices_only_in_sigma_tadpoles(self):
        found = graphs.quartic_radial_potential_vertices()
        self.assertEqual(sorted({label for _, label, _ in found}), ['Csigma1', 'Csigma2', 'Csigma3'])

    def test_cubic_vertices_in_kept_graphs(self):
        self.assertEqual(graphs.kept_graphs_with_lambda(),
                         ['C10', 'C12', 'C13', 'C15', 'C40', 'C44', 'C45', 'C47', 'C48',
                          'Csigma4', 'Csigma6', 'Csigma7', 'Csigma9'])

    def test_sigma_tadpoles_have_no_p2_dependence(self):
        for label in ('Csigma1', 'Csigma2', 'Csigma3'):
            row = graphs.row_by_label(label)
            parts = graphs.lambda_part(row, 'sigma')
            self.assertTrue(any(part.terms for part in parts))
            for part in parts:
                for value in part.terms.values():
                    self.assertFalse(sp.sympify(value).has(MOMENTUM_SQUARED))
            self.assertEqual(row['status'], 'mass shift only (dZ = 0)')


class ExternalLegTest(unittest.TestCase):
    def test_pi_plus_equals_pi_minus(self):
        for row in graphs.kept(graphs.table('pi')):
            for index, (line, orientation) in enumerate(graphs.routings(row)):
                total = TraceSum()
                for choice, scale, value in graphs.term_values(row, graphs.setup('pi-'), line, orientation):
                    if row['selection'] == 'without D4' and graphs.is_all_generator_sum(choice):
                        continue
                    total = total + graphs.residue_correction(row, value, scale, M_PI)
                self.assertTrue((total - row['routing values'][index]).is_zero(), row['label'])

    def test_exact_coefficients(self):
        for group in graphs.GROUPS:
            for row in graphs.kept(graphs.table(group)):
                for value in row['routing values']:
                    assert_no_float(value)
                    for coefficient in value.terms.values():
                        self.assertFalse(sp.sympify(coefficient).has(MOMENTUM_SQUARED))


class CompletenessTest(unittest.TestCase):
    """U(3) completeness, checked once with explicit generators and symbolic 3x3 matrices (never printed)."""

    def test_completeness(self):
        lam = [sp.sqrt(sp.Rational(2, 3)) * sp.eye(3),
               sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]]), sp.Matrix([[0, -sp.I, 0], [sp.I, 0, 0], [0, 0, 0]]),
               sp.Matrix([[1, 0, 0], [0, -1, 0], [0, 0, 0]]), sp.Matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]]),
               sp.Matrix([[0, 0, -sp.I], [0, 0, 0], [sp.I, 0, 0]]), sp.Matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]]),
               sp.Matrix([[0, 0, 0], [0, 0, -sp.I], [0, sp.I, 0]]), sp.diag(1, 1, -2) / sp.sqrt(3)]
        for i in range(9):
            for j in range(9):
                self.assertEqual(sp.simplify((lam[i] * lam[j]).trace() - 2 * (i == j)), 0)
        a = sp.Matrix(3, 3, sp.symbols('a0:9'))
        b = sp.Matrix(3, 3, sp.symbols('b0:9'))
        self.assertEqual(sp.expand(sum(((a * g).trace() * (g * b).trace() for g in lam), 0) - 2 * (a * b).trace()), 0)
        self.assertEqual(sp.expand(sum(((a * g * b * g).trace() for g in lam), 0) - 2 * a.trace() * b.trace()), 0)

    def test_completeness_rule_on_words(self):
        """The word rule used by the scripts: <A G B G> -> 2<A><B>, <A G><G B> -> 2<A B>."""
        import sigma_pipi_loop_algebra as generic
        single = generic.completeness(TraceSum.single(('A', 'G', 'B', 'G')), 'G')
        self.assertTrue((single - TraceSum({(('A',), ('B',)): 2})).is_zero())
        double = generic.completeness(TraceSum({(('A', 'G'), ('G', 'B')): 1}), 'G')
        self.assertTrue((double - TraceSum.single(('A', 'B'), 2)).is_zero())


if __name__ == '__main__':
    unittest.main()
