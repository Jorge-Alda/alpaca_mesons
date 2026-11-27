"""Unit tests for derive_sigma_pipi_cartesian_loop_large_epsilon_graphs.py (exact symbolic identities; no numerical matrices).

Run: PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_cartesian_loop_large_epsilon_graphs
(about four minutes).
"""
import pathlib
import unittest

import sympy as sp

import derive_sigma_pipi_cartesian_loop_large_epsilon_graphs as large
from sigma_pipi_loop_algebra import F_PI, TraceSum
from test_derive_sigma_pipi_cartesian_loop_traces import FORBIDDEN_COMPONENTS, FORBIDDEN_LOOP_VALUES

RECORD = (pathlib.Path(__file__).resolve().parent.parent / 'memory' / 'alp' / 'radial-mode' / 'finite-mixing'
          / 'two-pion-loop' / 'cartesian-estimate' / 'cartesian-sigma-pipi-loop-large-epsilon-graphs.md')
REGIONS = ('r < 1', 'r > 1', 'r = 1')
SHARE_POWERS = {'C1': 4, 'C2': 4, 'C3': 2, 'C7': 2, 'C8': 2, 'C9': 0, 'C10': 4, 'C11': 4, 'C12': 4, 'C13': 2,
                'C14': 2, 'C15': 2, 'C19': 2, 'C20': 2, 'C21': 0, 'C30': 0, 'C39': 0, 'C40': 2, 'C41': 2, 'C42': 0,
                'C43': 2, 'C44': 2, 'C45': 0, 'C46': 0, 'C47': 0, 'C48': -2, 'Csigma4': 4, 'Csigma5': 4,
                'Csigma6': 4, 'Csigma7': 2, 'Csigma8': 2, 'Csigma9': 2, 'Cpi2': 0, 'Cpi3': 0, 'Cpi7': 0, 'Cpi8': 0,
                'Cpi9': -2}
LEADING_GRAPHS = ['C1', 'C10', 'C11', 'C12', 'C2', 'Csigma4', 'Csigma5', 'Csigma6']


def anticommutator_trace(letter):
    """<X {lambda_pi+, lambda_pi-}> as a TraceSum."""
    return TraceSum.single((letter, 'l+', 'l-')) + TraceSum.single((letter, 'l-', 'l+'))


class LaurentTest(unittest.TestCase):
    def test_inverse_of_infinite_series_is_cut_and_exact_below_the_cut(self):
        u = large.U
        series = large.Laurent({0: 1, 2: 3})          # 1 + 3 u^2, exact
        inverse = series.inverse()
        self.assertEqual(inverse.order, large.INVERSE_DEPTH)
        expected = sp.series(1 / (1 + 3 * u ** 2), u, 0, large.INVERSE_DEPTH).removeO()
        self.assertEqual(sp.expand(sum(v * u ** p for p, v in inverse.terms.items()) - expected), 0)
        product = series * inverse
        self.assertEqual(product.terms, {0: 1})

    def test_inverse_with_negative_lowest_power(self):
        series = large.Laurent({-2: large.A, 2: 5})    # A/u^2 + 5 u^2
        inverse = series.inverse()
        self.assertEqual(inverse.lowest(), 2)
        self.assertEqual((series * inverse).terms, {0: 1})


class RegionTest(unittest.TestCase):
    def test_regions(self):
        merged, boundaries = large.regions()
        self.assertEqual([(interval['lower'], interval['upper']) for interval in merged],
                         [(0, 1), (1, sp.oo)])
        self.assertTrue(all(interval['constant'] for interval in merged))
        # three points per interval between consecutive comparison roots
        self.assertEqual(len(merged[0]['points']), 3 * (1 + len(merged[0]['interior roots'])))
        self.assertEqual(len(merged[1]['points']), 3 * (1 + len(merged[1]['interior roots'])))
        self.assertEqual(set(boundaries), set(large.comparison_roots()) | {sp.Integer(1)})

    def test_routing_choice(self):
        merged, boundaries = large.regions()
        below = {index for index, _ in merged[0]['choices'].values()}
        above = {index for index, _ in merged[1]['choices'].values()}
        self.assertEqual(below, {2})
        self.assertEqual(above, {3})
        self.assertEqual({index for index, _ in boundaries[sp.Integer(1)].values()}, {2})
        self.assertFalse(any(marked for interval in merged for _, marked in interval['choices'].values()))
        self.assertEqual(set(large.routing_dependent_products()), {'C12', 'C13', 'C15'})

    def test_stage_b_branches_agree(self):
        merged, boundaries = large.regions()
        for point in [point for interval in merged for point in interval['points']] + list(boundaries):
            self.assertTrue(large.stage_b_branch_agrees(point), point)

    def test_mass_conditions(self):
        (power, coefficient), _ = large.large_epsilon_conditions(large.REGION_BELOW)
        self.assertEqual((power, sp.expand(coefficient - (large.A - large.KAPPA))), (2, 0))
        (power, coefficient), _ = large.large_epsilon_conditions(large.REGION_BOUNDARY)
        self.assertEqual(power, 0)
        self.assertEqual(sp.expand(coefficient - 2 * large.TRACE_C[0] / large.A), 0)


class GraphTest(unittest.TestCase):
    def test_share_powers(self):
        for name in REGIONS:
            powers = {label: data['contribution'][0] for label, data in large.graph_results(name).items()}
            self.assertEqual(powers, SHARE_POWERS, name)

    def test_hand_derived_tadpoles(self):
        """C1 (s tadpole on V(pi+, pi-, s, s, s) = -3 B_0 eps^3/f_pi^3 <S^C_30 {.}>, size m_sigma^2/(32 pi^2)) and
        C2 (ALP-like tadpole, leading part -B_0 eps^3/f_pi^3 <S^C_12 {.}>/N, size m_a^2/(32 pi^2)):
        X/M_0 = X (-2 f_pi/(B_0 eps))/T_0 with m_sigma^2 -> eps^2 kappa_sigma, m_a^2 -> eps^2 A, N -> 1."""
        expected = {'C1': anticommutator_trace('SC30').scale(3 * large.KAPPA / (16 * sp.pi ** 2 * F_PI ** 2)),
                    'C2': anticommutator_trace('SC12').scale(large.A / (16 * sp.pi ** 2 * F_PI ** 2))}
        for label, value in expected.items():
            power, numerator = large.graph_results('r < 1')[label]['contribution']
            self.assertEqual(power, 4)
            self.assertTrue((numerator - value).is_zero(), label)

    def test_total(self):
        for name in REGIONS:
            power, numerator, labels = large.total_leading(name)
            self.assertEqual(power, 4)
            self.assertEqual(labels, LEADING_GRAPHS)
            self.assertTrue(large.graphs_sum_check(name))
            self.assertEqual(large.common_roots_in_ratio(numerator), [])

    def test_leading_numerator_letters(self):
        """The leading numerators contain no P^C_01 and no S^C_00: the meson components of the ALP-like lines
        (which carry M^{-1}_{m_a^2}(C), written with P^C_01 and S^C_00) do not reach the leading power, so the
        four-meson kinetic vertex (the only momentum-dependent vertex) does not enter it."""
        allowed = {'SC02', 'SC10', 'SC11', 'SC12', 'SC20', 'SC21', 'SC30', 'l+', 'l-'}
        for name in REGIONS:
            _, numerator, _ = large.total_leading(name)
            self.assertEqual(numerator.letters() - allowed, set(), name)
            for label in LEADING_GRAPHS:
                self.assertEqual(large.graph_results(name)[label]['contribution'][1].letters() - allowed, set())

    def test_regions_differ_only_where_expected(self):
        """r > 1 differs from r < 1 only by max(m_a, m_sigma) and the routing; the S^C_21 S^C_21 term of the
        leading numerator has 1/A for r < 1 and 1/kappa_sigma for r > 1."""
        key = (('SC10', 'l+', 'l-'), ('SC21',), ('SC21',))
        below = large.total_leading('r < 1')[1].terms[key]
        above = large.total_leading('r > 1')[1].terms[key]
        self.assertEqual(sp.simplify(below * large.A - above * large.KAPPA), 0)

    def test_d4_powers(self):
        self.assertEqual(large.d4_powers(), {'C3': 0, 'C15': 2, 'Cpi3': 0})

    def test_groups(self):
        for name in REGIONS:
            groups = large.same_power_groups(name)
            self.assertEqual(sorted(groups), [-2, 0, 2, 4])
            self.assertEqual(groups[4][0], LEADING_GRAPHS)
            self.assertTrue(all(status in ('single graph', 'same sign', 'could cancel', 'depends on sources')
                                for _, status, _ in groups.values()))
            self.assertIn(('Cpi2', 'Cpi3'), large.exactly_opposite_pairs(name))


class RecordTest(unittest.TestCase):
    def test_record_equals_script_output(self):
        text = RECORD.read_text()
        blocks = [large.markdown_graph_table(), large.markdown_graph_coefficients(), large.markdown_d4()]
        for name in large.region_names() + ['r = 1']:
            blocks += [large.markdown_leading(name), large.markdown_contributions(name), large.markdown_groups(name)]
        for block in blocks:
            self.assertIn(block, text)

    def test_record_content(self):
        text = RECORD.read_text()
        for token in FORBIDDEN_COMPONENTS + FORBIDDEN_LOOP_VALUES:
            self.assertNotIn(token, text, token)
        self.assertIn('No loop integral', text)
        self.assertIn('authoritative', text)
        self.assertNotIn('action', text.lower().replace('interaction', ''))
        for shorthand in ('T_0', '\\mathcal N'):          # no notation beyond the project records
            self.assertNotIn(shorthand, text, shorthand)


if __name__ == '__main__':
    unittest.main()
