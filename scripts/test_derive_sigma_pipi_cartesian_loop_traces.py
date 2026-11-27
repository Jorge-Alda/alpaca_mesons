"""Unit tests for derive_sigma_pipi_cartesian_loop_traces.py (exact symbolic identities only; no numerical matrices).

Run: PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_cartesian_loop_traces
"""
import pathlib
import unittest

import sympy as sp

import sigma_pipi_cartesian_loop_algebra as cartesian
import derive_sigma_pipi_cartesian_loop_traces as traces
from sigma_pipi_loop_algebra import B0, F_PI, M_A, M_MESON, M_PI, M_SIGMA, assert_no_float

RECORD_FOLDER = (pathlib.Path(__file__).resolve().parent.parent / 'memory' / 'alp' / 'radial-mode' / 'finite-mixing'
                 / 'two-pion-loop' / 'cartesian-estimate')
RECORDS = ('cartesian-sigma-pipi-loop-graph-selection.md', 'cartesian-sigma-pipi-loop-graph-list.md',
           'cartesian-sigma-pipi-loop-traces-1pi.md', 'cartesian-sigma-pipi-loop-traces-external-legs.md')
W = sp.Symbol('w', positive=True)


class ConstructionTest(unittest.TestCase):
    def test_words_equal_compact_trees(self):
        for entry in traces.entries():
            self.assertTrue(entry['constructions agree'], entry['label'])

    def test_exact(self):
        for entry in traces.entries():
            for value in entry['words']:
                assert_no_float(value)


class RoutingMaximumTest(unittest.TestCase):
    def test_routing_dependent_graphs(self):
        dependent = {entry['label']: entry['maximum']['counts']['mass-dependent']
                     for entry in traces.entries() if entry['common routing'] is None}
        self.assertEqual(dependent, {'C12': 18, 'C13': 34, 'C15': 34})
        for entry in traces.entries():
            self.assertEqual(entry['maximum']['counts']['numerical'], 0)
            self.assertEqual(entry['maximum']['marked'], [])

    def test_branches(self):
        """Routing 3 for m_sigma^2 - m_pi^2 <= s^2 and routing 4 above, s = m_a or max(m_P, m_a)."""
        for label in ('C12', 'C13', 'C15'):
            entry = traces.entry_by_label(label)
            self.assertEqual(set(entry['segments']), {M_A} if label == 'C12' else {M_A, sp.Max(M_MESON, M_A)})
            for scale, (conditions, indices) in entry['segments'].items():
                combination = M_SIGMA ** 2 - M_PI ** 2
                self.assertEqual(conditions, (sp.Le(combination, scale ** 2), sp.Gt(combination, scale ** 2)))
                self.assertEqual(indices, (2, 3))
            self.assertEqual(len(entry['maximum']['marked boundary']), entry['maximum']['counts']['mass-dependent'])

    def test_routing_pattern(self):
        """Every routing-dependent product has C_t^rho = c_t (a_rho s^2 - (m_sigma^2 - m_pi^2)) with
        (a_1, ..., a_4) = (2, 0, 3, -1); the sum over the routings vanishes at m_sigma^2 - m_pi^2 = s^2."""
        for entry in traces.entries():
            for key, (kind, info) in entry['maximum']['choices'].items():
                if kind != 'mass-dependent':
                    continue
                self.assertEqual(info['gamma'], -1)
                factor = -info['slopes'][0]
                self.assertEqual(list(info['offsets']), [2 * factor, 0, 3 * factor, -factor])
                self.assertEqual(list(info['slopes']), [-factor] * 4)

    def test_maximum_is_largest(self):
        """|C_t^{rho*}| >= |C_t^rho| for every rho and every product.  For mass-dependent products the
        coefficients are first checked to be g (a_rho + b_rho w) with w = (m_sigma^2 - m_pi^2)/s^2,
        then the inequality is proved on each branch with an exact univariate inequality."""
        for entry in traces.entries():
            for key, (kind, info) in entry['maximum']['choices'].items():
                coefficients = [sp.expand(value.terms.get(key, 0)) for value in entry['words']]
                if kind == 'independent':
                    continue
                self.assertEqual(kind, 'mass-dependent')
                scale, gamma = info['scale'], info['gamma']
                substitution = {M_SIGMA: sp.sqrt(W * scale ** 2 - gamma * M_PI ** 2)}
                reduced = [sp.expand(c.subs(substitution)) for c in coefficients]
                reference = sp.cancel(reduced[0] / (info['offsets'][0] + info['slopes'][0] * W))
                for value, offset, slope in zip(reduced, info['offsets'], info['slopes']):
                    self.assertEqual(sp.cancel(value - reference * (offset + slope * W)), 0)
                functions = [offset + slope * W for offset, slope in zip(info['offsets'], info['slopes'])]
                for segment in info['segments']:
                    chosen = functions[segment['index']]
                    interval = sp.Interval(segment['lower'], segment['upper'], not segment['lower closed'],
                                           not segment['upper closed'])
                    for function in functions:
                        allowed = sp.solve_univariate_inequality(chosen ** 2 - function ** 2 >= 0, W, relational=False)
                        self.assertTrue(interval.is_subset(allowed), (entry['label'], key))

    def test_displayed_form(self):
        for entry in traces.entries():
            self.assertTrue(traces.displayed_equals_maximum(entry), entry['label'])

    def test_maximum_of_routing_independent_graph(self):
        for entry in traces.entries():
            if entry['common routing'] is not None:
                self.assertTrue((entry['maximum']['value'] - entry['words'][0]).is_zero())


def rescale_expression(expression, scale):
    """B_0, f_pi and the masses -> c x; Max(c x, c y) -> c Max(x, y)."""
    value = sp.sympify(expression).subs({B0: scale * B0, F_PI: scale * F_PI, M_SIGMA: scale * M_SIGMA,
                                         M_A: scale * M_A, M_MESON: scale * M_MESON, M_PI: scale * M_PI},
                                        simultaneous=True)
    return value.replace(lambda node: isinstance(node, sp.Max),
                         lambda node: scale * sp.Max(*[sp.cancel(argument / scale) for argument in node.args]))


class DimensionTest(unittest.TestCase):
    def test_rescaling(self):
        """X -> c^{d_X} X with d = 1 for B_0, f_pi, the masses and each source letter, d = 0 for eps, N,
        lambda_pi, M^{-1}_{m_a^2}(C): X^max scales as c^1, dZ^max as c^0; Piecewise conditions are
        homogeneous of degree 2."""
        scale = sp.Symbol('c', positive=True)
        for entry in traces.entries():
            target = 1 if entry['group'] == 'one_particle_irreducible' else 0
            for key, coefficient in entry['maximum']['value'].terms.items():
                power = target - cartesian.source_count(key)
                coefficient = sp.piecewise_fold(sp.sympify(coefficient))
                branches = coefficient.args if isinstance(coefficient, sp.Piecewise) else ((coefficient, sp.true),)
                for piece, condition in branches:
                    self.assertEqual(sp.simplify(rescale_expression(piece, scale) - scale ** power * piece), 0)
                    for relation in condition.atoms(sp.core.relational.Relational):
                        difference = relation.lhs - relation.rhs
                        self.assertEqual(sp.simplify(rescale_expression(difference, scale) - scale ** 2 * difference), 0)


class ExternalLegTest(unittest.TestCase):
    def test_no_mass_difference_factor_in_pion_residue(self):
        """dZ_pi contains no factor m_a^2 - m_sigma^2 (in any graph or routing, nor in the maximum)."""
        for entry in traces.entries():
            if entry['group'] != 'pi':
                continue
            for value in list(entry['words']) + [entry['maximum']['value']]:
                for coefficient in value.terms.values():
                    numerator, denominator = sp.fraction(sp.factor(sp.together(coefficient)))
                    for part in (numerator, denominator):
                        for factor, _ in sp.factor_list(part)[1]:
                            self.assertFalse(sp.expand(factor - (M_A ** 2 - M_SIGMA ** 2)) == 0)
                            self.assertFalse(sp.expand(factor - (M_A - M_SIGMA)) == 0)
                            self.assertFalse(sp.expand(factor + (M_A - M_SIGMA)) == 0)

    def test_pi_plus_equals_pi_minus(self):
        minus = traces.pi_minus_residue_corrections()
        for entry in traces.entries():
            if entry['group'] == 'pi':
                self.assertTrue((minus[entry['label']] - entry['words'][0]).is_zero())

    def test_estimate_inputs(self):
        inputs = traces.estimate_inputs()
        self.assertEqual(len(inputs['one particle irreducible']), 26)
        total = sum((entry['maximum']['value'] for entry in traces.entries() if entry['group'] == 'pi'),
                    start=type(inputs['delta Z pi'])())
        self.assertTrue((total - inputs['delta Z pi']).is_zero())


FORBIDDEN_COMPONENTS = ('E_{', 'E_1', 'E_2', 'E_3', 'E_i', '_{ij}', '_{ab}', r'\lambda_i', r'\lambda_j', r'\lambda_{i',
                        r'\lambda_{j', r'\lambda_a', 'G1', 'G2', "G'")
FORBIDDEN_LOOP_VALUES = (r'\log', r'\ln', 'A_0(', 'C_0(', 'B_0(p', 'B_0(m', 'B_{0}(', 'C_{0}(', r'\overline{', 'Passarino')


class OutputTest(unittest.TestCase):
    def test_markdown_has_no_component_notation(self):
        for entry in traces.entries():
            text = traces.entry_markdown(entry)
            self.assertTrue(traces.no_component_notation(text))
            for token in FORBIDDEN_COMPONENTS + FORBIDDEN_LOOP_VALUES:
                self.assertNotIn(token, text)

    def test_records_equal_script_output(self):
        """The result blocks of the trace records and the tables of the graph-list record are the
        current --markdown output of the scripts."""
        import derive_sigma_pipi_cartesian_loop_graphs as graphs
        one = (RECORD_FOLDER / RECORDS[2]).read_text()
        legs = (RECORD_FOLDER / RECORDS[3]).read_text()
        for entry in traces.entries():
            block = traces.entry_heading(entry) + '\n' + traces.entry_markdown(entry)
            self.assertIn(block, one if entry['group'] == 'one_particle_irreducible' else legs, entry['label'])
        graph_list = (RECORD_FOLDER / RECORDS[1]).read_text()
        for group in graphs.GROUPS:
            self.assertIn(graphs.markdown_table(group), graph_list)

    def test_records(self):
        """The trace and graph records contain no component notation, no E_ij and no loop value."""
        for name in RECORDS:
            path = RECORD_FOLDER / name
            if not path.exists():
                self.fail(f'missing record {name}')
            text = path.read_text()
            for token in FORBIDDEN_COMPONENTS + FORBIDDEN_LOOP_VALUES:
                self.assertNotIn(token, text, (name, token))
            self.assertIn('No loop integral', text)
            self.assertNotIn('action', text.lower().replace('interaction', ''))


if __name__ == '__main__':
    unittest.main()
