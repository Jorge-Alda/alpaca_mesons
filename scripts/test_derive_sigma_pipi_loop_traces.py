"""Unit tests for derive_sigma_pipi_loop_traces.py and the trace algebra of
sigma_pipi_loop_algebra.py (exact symbolic identities only).

Run: PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_loop_traces

The only explicit matrices are the U(3) generators and two symbolic 3x3
matrices in CompletenessTest, used to confirm the completeness identities
once; they are never printed.
"""
import itertools
import unittest

import sympy as sp

import sigma_pipi_loop_algebra as algebra
import derive_sigma_pipi_loop_traces as traces
from sigma_pipi_loop_algebra import CompactSum, TraceSum, letter


def _generators():
    """Normalized U(3) generators lambda_0, ..., lambda_8 (internal to this test)."""
    gell_mann = [
        sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]]), sp.Matrix([[0, -sp.I, 0], [sp.I, 0, 0], [0, 0, 0]]),
        sp.Matrix([[1, 0, 0], [0, -1, 0], [0, 0, 0]]), sp.Matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]]),
        sp.Matrix([[0, 0, -sp.I], [0, 0, 0], [sp.I, 0, 0]]), sp.Matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]]),
        sp.Matrix([[0, 0, 0], [0, 0, -sp.I], [0, sp.I, 0]]),
        sp.Matrix([[1, 0, 0], [0, 1, 0], [0, 0, -2]]) / sp.sqrt(3)]
    return [sp.sqrt(sp.Rational(2, 3)) * sp.eye(3)] + gell_mann


class CompletenessTest(unittest.TestCase):
    def test_u3_completeness(self):
        generators = _generators()
        first = sp.Matrix(3, 3, sp.symbols('a0:9'))
        second = sp.Matrix(3, 3, sp.symbols('b0:9'))
        for i, j in itertools.product(range(9), repeat=2):
            self.assertEqual(sp.simplify((generators[i] * generators[j]).trace()), 2 if i == j else 0)
        self.assertEqual(generators[0], sp.sqrt(sp.Rational(2, 3)) * sp.eye(3))
        same_trace = sp.expand(sum((first * g * second * g).trace() for g in generators))
        self.assertEqual(sp.expand(same_trace - 2 * first.trace() * second.trace()), 0)
        two_traces = sp.expand(sum((first * g).trace() * (g * second).trace() for g in generators))
        self.assertEqual(sp.expand(two_traces - 2 * (first * second).trace()), 0)
        octet = sp.expand(sum((first * g * second * g).trace() for g in generators[1:]))
        self.assertNotEqual(sp.expand(octet - 2 * first.trace() * second.trace()), 0)


class TreeCompletenessTest(unittest.TestCase):
    """Completeness on trees equals completeness on words, for generic letters."""

    def check(self, traces_list):
        compact = CompactSum()
        compact.add_term(tuple(traces_list), 1)
        words = algebra.completeness(compact.to_trace_sum(), 'G1')
        trees = algebra.compact_completeness(compact, 'G1')
        self.assertTrue((trees.to_trace_sum() - words).is_zero())

    def test_two_traces(self):
        self.check([('P', letter('X'), ('S', letter('G1'), letter('Y'), letter('Z'))),
                    ('P', letter('W'), ('S', letter('G1'), letter('V')))])

    def test_commutators(self):
        self.check([('P', ('C', letter('G1'), letter('Y')), ('C', letter('Z'), letter('W'))),
                    ('P', letter('V'), ('A', letter('G1'), letter('U')))])

    def test_same_trace(self):
        self.check([('P', letter('X'), ('S', letter('G1'), letter('G1'), letter('Y'), letter('Z')))])
        self.check([('P', ('C', letter('G1'), letter('Y')), ('C', letter('G1'), letter('W')))])

    def test_cut_and_split(self):
        tree = ('P', letter('X'), ('S', letter('G1'), letter('Y'), letter('Z')))
        rest = algebra.cut(tree, algebra.UNIT, 'G1')
        self.assertTrue((algebra.trace_of_tree(tree)
                         - algebra.trace_of_tree(('P', letter('G1'), rest))).is_zero())
        pieces = algebra.split(tree, 'G1')
        total = TraceSum()
        for k, alpha, beta in pieces:
            total = total + algebra.trace_of_tree(('P', alpha, letter('G1'), beta)).scale(k)
        self.assertTrue((total - algebra.trace_of_tree(tree)).is_zero())


ROUTING_DEPENDENT = {'bubble {pi+,pi-}-A-{sigma}-A-', 'bubble {pi+,pi-}-A-{sigma}-M-', 'bubble {pi+,pi-}-M-{sigma}-M-',
                     'triangle {pi+}-M-{pi-}-A-{sigma}-A-', 'triangle {pi+}-M-{pi-}-M-{sigma}-A-',
                     'triangle {pi+}-M-{pi-}-A-{sigma}-M-', 'triangle {pi+}-M-{pi-}-M-{sigma}-M-',
                     'sigma self-energy bubble {in}-A-{out}-A-', 'sigma self-energy bubble {in}-A-{out}-M-',
                     'sigma self-energy bubble {in}-M-{out}-M-'}

D4_REMOVED = {'tadpole {pi+,pi-,sigma}-M-', 'bubble {pi+,pi-}-M-{sigma}-M-',
              'pi+ self-energy tadpole {in,out}-M-', 'pi- self-energy tadpole {in,out}-M-'}


class KeptTracesTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.entries = traces.kept_graph_traces()
        cls.by_label = {entry['label']: entry for entry in cls.entries}

    def test_kept_set(self):
        counts = {group: sum(1 for entry in self.entries if entry['group'] == group) for group in traces.GROUPS}
        self.assertEqual(counts, {'one_particle_irreducible': 26, 'sigma': 6, 'pi+': 5, 'pi-': 5})
        removed = {entry['label'] for entry in self.entries if entry['excluded'] is not None}
        self.assertEqual(removed, D4_REMOVED)

    def test_constructions_agree(self):
        """Words and compact trees agree in every routing (exact)."""
        for entry in self.entries:
            self.assertTrue(entry['constructions_agree'], entry['label'])
            algebra.assert_no_float(entry['words'])

    def test_d4_removal_is_exact(self):
        """For the graphs with a D4 term: kept value = full value - excluded term, routing by routing."""
        for group, (legs, basis_dot) in traces.setups().items():
            for entry in self.entries:
                if entry['group'] != group or entry['excluded'] is None:
                    continue
                full = traces.graphs.routing_values(entry, legs, basis_dot, group, pieces='all')
                for total, bare, rest in zip(full, entry['excluded_routings'], entry['routing_values']):
                    self.assertTrue((total - bare - rest).is_zero(), entry['label'])
                    self.assertFalse(bare.is_zero(), entry['label'])

    def test_no_generators_and_compact_output(self):
        for entry in self.entries:
            self.assertFalse(any(algebra.is_generator(name) for name in entry['words'].letters()))
            for part in traces.displayed(entry):
                self.assertTrue(traces.no_component_notation(algebra.compact_latex(part)))
        for blocks in traces.markdown_blocks(self.entries).values():
            for text in blocks:
                self.assertTrue(traces.no_component_notation(text))

    def test_routing_dependence(self):
        dependent = {entry['label'] for entry in self.entries
                     if entry['spread']['numerical'] or entry['spread']['mass-dependent']}
        self.assertEqual(dependent, ROUTING_DEPENDENT)
        without_common = {entry['label'] for entry in self.entries if entry['common routing'] is None}
        self.assertEqual(without_common, {label for label in ROUTING_DEPENDENT if 'self-energy' not in label})
        for label in ROUTING_DEPENDENT:
            if 'self-energy' in label:
                self.assertEqual(self.by_label[label]['common routing'], 1, label)
        spread = self.by_label['bubble {pi+,pi-}-A-{sigma}-A-']['spread']
        self.assertEqual((spread['numerical'], spread['mass-dependent']), (20, 12))
        self.assertEqual(spread['minimum ratio'], sp.Rational(-1, 3))
        self.assertEqual(self.by_label['sigma self-energy bubble {in}-A-{out}-A-']['spread']['minimum ratio'], 0)

    def test_displayed_form_reproduces_the_maximum(self):
        """Base routing (trees) plus correction (words) = termwise maximum."""
        for entry in self.entries:
            base, correction = traces.displayed(entry)
            difference = (base.to_trace_sum() + correction.to_trace_sum() - entry['words']).simplified()
            self.assertTrue(all(traces.graphs.expression_is_zero(value) for value in difference.terms.values()),
                            entry['label'])

    def test_known_closed_forms(self):
        sixteen_pi_squared = 16 * sp.pi ** 2
        tadpole = self.by_label['tadpole {pi+,pi-,sigma}-S-']['words']
        expected = (TraceSum.single(('S30', 'l+', 'l-')) + TraceSum.single(('S30', 'l-', 'l+'))).scale(
            -3 * algebra.B0 * algebra.EPSILON ** 3 / algebra.F_PI ** 3 * algebra.M_SIGMA ** 2 / (2 * sixteen_pi_squared))
        self.assertTrue((tadpole - expected).is_zero())
        # pi+ residue from the closed ALP-like line: d V/d p^2 = <[l_a, l+][l_a, l-]>/(12 f^2), l_a = 2 eps Minv/sqrt(N)
        residue = self.by_label['pi+ self-energy tadpole {in,out}-A-']['residue']
        commutators = algebra.trace_of_tree(('P', ('C', letter('Minv'), letter('l+')),
                                             ('C', letter('Minv'), letter('l-'))))
        expected = commutators.scale(-sp.Rational(1, 2) * algebra.M_A ** 2 / sixteen_pi_squared
                                     * 4 * algebra.EPSILON ** 2 / (12 * algebra.NORMALIZATION * algebra.F_PI ** 2))
        self.assertTrue((residue - expected).is_zero())
        # closed meson line on the five-leg vertex with the D4 term removed (scale m_a): ALP-field entry and
        # ALP-like pole only; together with the closed ALP-like line it is the bare ALP-field line
        meson_tadpole = self.by_label['tadpole {pi+,pi-,sigma}-M-']['words']
        alp_tadpole = self.by_label['tadpole {pi+,pi-,sigma}-A-']['words']
        field_entry = (TraceSum.single(('S12', 'l+', 'l-')) + TraceSum.single(('S12', 'l-', 'l+'))).scale(
            -algebra.B0 * algebra.EPSILON ** 3 / algebra.F_PI ** 3 * algebra.M_A ** 2 / (2 * sixteen_pi_squared))
        self.assertTrue((meson_tadpole + alp_tadpole - field_entry).is_zero())

    def test_sources(self):
        self.assertEqual(traces.sources_in(self.entries),
                         ['P01', 'P11', 'P20', 'S00', 'S01', 'S02', 'S10', 'S11', 'S12', 'S20', 'S21', 'S30'])


if __name__ == '__main__':
    unittest.main()
