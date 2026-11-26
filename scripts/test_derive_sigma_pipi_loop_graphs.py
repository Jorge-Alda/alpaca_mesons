"""Unit tests for derive_sigma_pipi_loop_graphs.py (exact symbolic identities only).

Run: PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_loop_graphs
"""
import unittest

import sympy as sp

import sigma_pipi_loop_algebra as algebra
import derive_sigma_pipi_loop_graphs as graphs

VANISHING_ONE_PARTICLE_IRREDUCIBLE = {
    'bubble {pi+,sigma}-S-{pi-}-S-', 'bubble {pi+,sigma}-A-{pi-}-S-', 'bubble {pi+,sigma}-A-{pi-}-A-',
    'bubble {pi+}-S-{pi-,sigma}-S-', 'bubble {pi+}-A-{pi-,sigma}-S-', 'bubble {pi+}-A-{pi-,sigma}-A-',
    'triangle {pi+}-S-{pi-}-S-{sigma}-S-', 'triangle {pi+}-S-{pi-}-A-{sigma}-S-',
    'triangle {pi+}-S-{pi-}-M-{sigma}-S-', 'triangle {pi+}-A-{pi-}-S-{sigma}-S-',
    'triangle {pi+}-A-{pi-}-A-{sigma}-S-', 'triangle {pi+}-A-{pi-}-M-{sigma}-S-',
    'triangle {pi+}-S-{pi-}-S-{sigma}-A-', 'triangle {pi+}-S-{pi-}-A-{sigma}-A-',
    'triangle {pi+}-S-{pi-}-M-{sigma}-A-', 'triangle {pi+}-A-{pi-}-S-{sigma}-A-',
    'triangle {pi+}-A-{pi-}-A-{sigma}-A-', 'triangle {pi+}-A-{pi-}-M-{sigma}-A-',
    'triangle {pi+}-S-{pi-}-S-{sigma}-M-', 'triangle {pi+}-S-{pi-}-A-{sigma}-M-',
    'triangle {pi+}-A-{pi-}-S-{sigma}-M-', 'triangle {pi+}-A-{pi-}-A-{sigma}-M-'}

EXCLUDED_PIECES = {'tadpole {pi+,pi-,sigma}-M-', 'bubble {pi+,pi-}-M-{sigma}-M-',
                   'pi+ self-energy tadpole {in,out}-M-', 'pi- self-energy tadpole {in,out}-M-'}


class GraphSelectionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = graphs.summary()
        cls.all_rows = [row for group in ('one_particle_irreducible', 'sigma', 'pi+', 'pi-')
                        for row in cls.data[group]]

    def test_counts(self):
        rows = self.data['one_particle_irreducible']
        self.assertEqual(len(rows), 48)
        counts = {kind: sum(1 for row in rows if graphs.topology(row) == kind)
                  for kind in ('tadpole', 'bubble', 'triangle')}
        self.assertEqual(counts, {'tadpole': 3, 'bubble': 18, 'triangle': 27})
        self.assertEqual(len(graphs.kept(rows)), 26)
        self.assertEqual({kind: len(graphs.kept(self.data[kind])) for kind in ('sigma', 'pi+', 'pi-')},
                         {'sigma': 6, 'pi+': 5, 'pi-': 5})

    def test_no_box_and_no_two_point_vertex(self):
        self.assertEqual(self.data['surjections_three'][4], 0)
        self.assertEqual(self.data['surjections_three'][5], 0)
        self.assertEqual(self.data['surjections_two'][3], 0)
        for row in self.data['one_particle_irreducible']:
            self.assertLessEqual(len(row['blocks']), 3)
            self.assertTrue(all(len(block) >= 1 for block in row['blocks']), row['label'])

    def test_symmetry_factors(self):
        for row in self.data['one_particle_irreducible'] + self.data['sigma']:
            lines = row['lines']
            if len(lines) == 1 or (len(lines) == 2 and lines[0] == lines[1]):
                self.assertEqual(row['symmetry'], sp.Rational(1, 2), row['label'])
            else:
                self.assertEqual(row['symmetry'], 1, row['label'])

    def test_no_series_drop(self):
        self.assertTrue(graphs.no_series_drop(self.data))
        for row in self.all_rows:
            if not row['status'].startswith('kept'):
                self.assertIn(row['status'], ('vanishes identically', 'no residue correction (exact: dZ = 0)'))

    def test_vanishing_set_is_exact(self):
        legs, basis_dot = graphs.three_point_setup()
        vanishing = {row['label'] for row in self.data['one_particle_irreducible']
                     if row['status'] == 'vanishes identically'}
        self.assertEqual(vanishing, VANISHING_ONE_PARTICLE_IRREDUCIBLE)
        for row in self.data['one_particle_irreducible']:
            if row['label'] in vanishing:
                for start, sign in graphs.all_routings(row):
                    self.assertTrue(graphs.graph_value(row, legs, basis_dot, start=start,
                                                       orientation=sign).is_zero())

    def test_graphs_with_radial_and_alp_lines_kept(self):
        both = {row['label'] for row in self.all_rows
                if {'S', 'A'} <= set(row['lines']) and row in graphs.kept([row])}
        self.assertEqual(both, {'bubble {pi+,pi-}-A-{sigma}-S-', 'triangle {pi+}-M-{pi-}-A-{sigma}-S-',
                                'triangle {pi+}-M-{pi-}-S-{sigma}-A-', 'sigma self-energy bubble {in}-A-{out}-S-'})

    def test_term_scales(self):
        """D3: m_term = largest line scale; m_P only for the bare generator sum."""
        m_p, m_a, m_sigma = graphs.M_P, algebra.M_A, algebra.M_SIGMA
        choice = [(1, None, 'G1', graphs.GENERATOR_SUM), (1, None, None, graphs.RADIAL_LINE)]
        self.assertEqual(graphs.term_scale(choice), sp.Max(m_p, m_sigma))
        choice = [(1, None, None, graphs.ALP_FIELD), (1, None, None, graphs.ALP_LIKE_LINE)]
        self.assertEqual(graphs.term_scale(choice), m_a)
        scales = {row['label']: set(graphs.possible_scales(row)) for row in self.all_rows if row in graphs.kept([row])}
        expected = {
            'tadpole {pi+,pi-,sigma}-S-': {m_sigma}, 'tadpole {pi+,pi-,sigma}-A-': {m_a},
            'tadpole {pi+,pi-,sigma}-M-': {m_a},                                   # D4 term (m_P) removed
            'bubble {pi+,sigma}-M-{pi-}-S-': {sp.Max(m_p, m_sigma)},               # only the generator sum
            'bubble {pi+,pi-}-M-{sigma}-S-': {sp.Max(m_a, m_sigma)},               # no generator sum
            'bubble {pi+,sigma}-M-{pi-}-M-': {sp.Max(m_p, m_a)},                   # all-generator term zero
            'bubble {pi+,pi-}-A-{sigma}-A-': {m_a},
            'bubble {pi+,pi-}-A-{sigma}-M-': {m_a, sp.Max(m_p, m_a)},
            'bubble {pi+,pi-}-M-{sigma}-M-': {m_a, sp.Max(m_p, m_a)},             # D4 term (m_P) removed
            'triangle {pi+}-M-{pi-}-A-{sigma}-S-': {sp.Max(m_p, m_a, m_sigma)},
            'sigma self-energy bubble {in}-M-{out}-S-': {sp.Max(m_p, m_sigma), sp.Max(m_a, m_sigma)},
            'sigma self-energy bubble {in}-M-{out}-M-': {m_p, m_a, sp.Max(m_p, m_a)},   # kept: two sigma ends
            'pi+ self-energy tadpole {in,out}-M-': {m_a}}
        for label, value in expected.items():
            self.assertEqual(scales[label], value, label)
        for row in self.all_rows:
            if row['excluded'] is not None:
                self.assertNotIn(m_p, graphs.possible_scales(row), row['label'])

    def test_excluded_pieces(self):
        excluded = {row['label'] for row in self.all_rows if row['excluded'] is not None}
        self.assertEqual(excluded, EXCLUDED_PIECES)
        for row in self.all_rows:
            if row['excluded'] is not None:
                self.assertTrue(graphs.is_meson_only(row))
                self.assertFalse(row['excluded'].is_zero())
                for value in row['excluded_routings']:
                    self.assertTrue(graphs.satisfies_restated_d4(value, row['group']), row['label'])
                    self.assertFalse(any(graphs.is_alp_sector_letter(name) for name in value.letters()))
                    for coefficient in value.terms.values():
                        self.assertEqual(sp.sympify(coefficient).free_symbols & {algebra.M_A, algebra.NORMALIZATION},
                                         set())
                self.assertEqual(row['pieces'], 'rest')
        sigma_meson = next(row for row in self.data['sigma'] if row['label'].endswith('{in}-M-{out}-M-'))
        self.assertEqual(sigma_meson['excluded_status'], 'kept (more than one sigma end)')
        self.assertEqual(sigma_meson['pieces'], 'all')

    def test_alp_sector_letters(self):
        for name in ('Minv', 'P01', 'S01', 'S12', 'P11', 'S21'):
            self.assertTrue(graphs.is_alp_sector_letter(name), name)
        for name in ('S00', 'S10', 'S20', 'S30', 'P20', 'l+', 'l-'):
            self.assertFalse(graphs.is_alp_sector_letter(name), name)

    def test_split_is_exact(self):
        """all = bare + rest for every meson-only graph and routing (size-rule values)."""
        legs, basis_dot = graphs.three_point_setup()
        for row in self.data['one_particle_irreducible']:
            if not graphs.is_meson_only(row):
                continue
            pieces = {kind: graphs.routing_values(row, legs, basis_dot, 'one_particle_irreducible', pieces=kind)
                      for kind in ('all', 'bare', 'rest')}
            for total, bare, rest in zip(pieces['all'], pieces['bare'], pieces['rest']):
                self.assertTrue((total - bare - rest).is_zero(), row['label'])


class RoutingMaximumTest(unittest.TestCase):
    """D12: termwise routing maximum (exact)."""

    @classmethod
    def setUpClass(cls):
        cls.data = graphs.summary()
        cls.kept_rows = [row for group in ('one_particle_irreducible', 'sigma', 'pi+', 'pi-')
                         for row in graphs.kept(cls.data[group])]

    def test_routing_independent_graphs_unchanged(self):
        """(a) A graph whose products do not depend on the routing enters with its
        single-routing value."""
        independent = 0
        for row in self.kept_rows:
            counts = row['maximum']['counts']
            if counts['numerical'] == 0 and counts['mass-dependent'] == 0:
                independent += 1
                for value in row['routing_values']:
                    self.assertTrue((row['value'] - value).is_zero(), row['label'])
        self.assertEqual(independent, 42 - 10)

    def test_largest_magnitude(self):
        """(b) |C_t^{rho*}| >= |C_t^rho| for every routing and product: exactly for
        numerical ratios, per mass branch for the mass-dependent ones (checked on
        the original coefficients with m_sigma^2 = w s^2 - gamma m_pi^2)."""
        w = graphs.BRANCH_VARIABLE
        mass_dependent = 0
        for row in self.kept_rows:
            for key, (kind, info) in row['maximum']['choices'].items():
                coefficients = [sp.expand(value.terms.get(key, 0)) for value in row['routing_values']]
                if kind == 'independent':
                    continue
                if kind == 'numerical':
                    chosen = coefficients[info]
                    for coefficient in coefficients:
                        ratio = sp.cancel(coefficient / chosen)
                        self.assertFalse(ratio.free_symbols)
                        self.assertLessEqual(abs(ratio), 1, (row['label'], key))
                    continue
                mass_dependent += 1
                scale, gamma = info['scale'], info['gamma']
                common = sp.gcd_list([c for c in coefficients if c != 0])
                pion = tree_mass = graphs.tree_script.M_PI
                placeholder = sp.Dummy('s', positive=True)
                functions = []
                for coefficient in coefficients:
                    reduced = sp.cancel((coefficient / common).subs(scale, placeholder)
                                        .subs(algebra.M_SIGMA, sp.sqrt(w * placeholder ** 2 - gamma * pion ** 2)))
                    reduced = sp.cancel(reduced.subs(placeholder, 1))
                    self.assertFalse(reduced.has(tree_mass), (row['label'], key))
                    functions.append(sp.expand(reduced))
                for segment in info['segments']:
                    chosen = functions[segment['index']]
                    lower, upper = segment['lower'], segment['upper']
                    probe = lower + 1 if upper == sp.oo else (lower + upper) / 2
                    sign = 1 if chosen.subs(w, probe) > 0 else -1
                    roots = [root for root in sp.solve(chosen, w) if lower < root < upper] if chosen.has(w) else []
                    self.assertEqual(roots, [], (row['label'], key))
                    for function in functions:
                        for other in (function, -function):
                            difference = sp.expand(sign * chosen - other)   # affine in w, must be >= 0
                            self.assertGreaterEqual(difference.subs(w, lower), 0)
                            if upper == sp.oo:
                                self.assertGreaterEqual(difference.coeff(w, 1), 0)
                            else:
                                self.assertGreaterEqual(difference.subs(w, upper), 0)
        self.assertEqual(mass_dependent, 80)

    def test_branch_patterns(self):
        """The three branch patterns of the mass-dependent products."""
        patterns = set()
        for row in self.kept_rows:
            for kind, info in row['maximum']['choices'].values():
                if kind == 'mass-dependent':
                    patterns.add((info['gamma'], tuple((segment['lower'], segment['upper'], segment['index'])
                                                      for segment in info['segments'])))
        self.assertEqual(patterns, {(-1, ((0, 1, 2), (1, sp.oo, 3))),
                                    (-1, ((0, 2, 2), (2, 5, 3), (5, sp.oo, 2))),
                                    (-2, ((0, 2, 0), (2, sp.oo, 1)))})

    def test_ties(self):
        self.assertEqual(sum(len(row['maximum']['marked']) for row in self.kept_rows), 0)
        self.assertEqual(sum(len(row['maximum']['marked boundary']) for row in self.kept_rows), 59)
        self.assertEqual(graphs.tie_choice([1, -1], [0, 1]), (0, True))
        self.assertEqual(graphs.tie_choice([1, -1, 1], [0, 1]), (0, False))
        self.assertEqual(graphs.tie_choice([1, -1, -1], [0, 1]), (1, False))
        self.assertEqual(graphs.largest_numerical([sp.Rational(1, 2), -2, 1]), (1, False))

    def test_half_line(self):
        """|1 - w| against |w/2| on w > 0: break points 2/3 and 2."""
        segments = graphs.largest_on_half_line([1, 0], [-1, sp.Rational(1, 2)])
        self.assertEqual([(s['lower'], s['upper'], s['index']) for s in segments],
                         [(0, sp.Rational(2, 3), 0), (sp.Rational(2, 3), 2, 1), (2, sp.oo, 0)])

    def test_closed_lines_cancel_in_charged_residue(self):
        """(c) dZ of the closed ALP-like line plus dZ of the closed meson line (D4
        term removed) vanishes exactly, so dZ_{pi+-} is the sum of the three bubbles."""
        for kind in ('pi+', 'pi-'):
            rows = {row['label'].split('self-energy ')[1]: row for row in graphs.kept(self.data[kind])}
            closed = rows['tadpole {in,out}-A-']['value'] + rows['tadpole {in,out}-M-']['value']
            self.assertTrue(closed.is_zero(), kind)

    def test_field_basis_identity(self):
        """(d) With a common routing and the per-term scales, an ALP-like line plus
        pieces (ii)+(iii)+(iv) of a meson line equals the bare ALP-field line
        (all three have the line scale m_a)."""
        checked = 0
        for group, externals in (('one_particle_irreducible', ('sigma', 'pi+', 'pi-')), ('sigma', ('in', 'out')),
                                 ('pi+', ('in', 'out'))):
            legs, basis_dot = graphs.three_point_setup() if group == 'one_particle_irreducible' \
                else graphs.two_point_setup(group)
            graph_list, _ = graphs.enumerate_cycle_graphs(externals)
            for graph in graph_list:
                for position, line in enumerate(graph['lines']):
                    if line != 'A':
                        continue
                    def replaced(new):
                        lines = list(graph['lines'])
                        lines[position] = new
                        return {**graph, 'lines': tuple(lines)}
                    not_bare = (lambda choice, i=position: choice[i][3] != graphs.GENERATOR_SUM)
                    for start, orientation in graphs.all_routings(graph):
                        options = dict(start=start, orientation=orientation)
                        alp = graphs.graph_value(graph, legs, basis_dot, **options)
                        meson = graphs.graph_value(replaced('M'), legs, basis_dot, choice_filter=not_bare,
                                                   **options)
                        field = graphs.graph_value(replaced('F'), legs, basis_dot, **options)
                        self.assertTrue((alp + meson - field).is_zero(), (graph['lines'], position))
                        checked += 1
        self.assertGreater(checked, 50)


if __name__ == '__main__':
    unittest.main()
