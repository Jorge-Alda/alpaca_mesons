"""Unit tests for derive_sigma_pipi_loop_large_epsilon.py (exact symbolic identities only).

Run: PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_loop_large_epsilon
(about two minutes: the large-eps results of all kept graphs are computed once).
"""
import unittest

import sympy as sp

import sigma_pipi_loop_algebra as algebra
import derive_sigma_pipi_loop_large_epsilon as large
from sigma_pipi_loop_algebra import B0, EPSILON, TraceSum

Q = large.TRACE_MINV_SQUARED
T = large.INVERSE_EPSILON_SQUARED
A = large.ALP_MASS_COEFFICIENT
C0, C1 = large.TRACE_C_POWER[0], large.TRACE_C_POWER[1]


class QuadraticLagrangianTest(unittest.TestCase):
    def test_operator_identities(self):
        for name, value in large.quadratic_form_checks().items():
            self.assertTrue(value, name)

    def test_stationarity(self):
        self.assertTrue(large.stationarity_check())

    def test_pole_function(self):
        _, pole_function = large.alp_inverse_propagator()
        expected = (large.SQUARED_MOMENTUM + 2 * B0 * EPSILON ** 2 * large.TRACE_S02
                    + 2 * EPSILON ** 2 * large.TRACE_C_MINV)
        self.assertEqual(sp.expand(pole_function - expected), 0)

    def test_normalization_is_norm_of_pole_vector(self):
        self.assertEqual(sp.expand(large.pole_vector_norm() - (1 + 2 * EPSILON ** 2 * Q)), 0)


class ResidueIdentityTest(unittest.TestCase):
    def test_residue_identities(self):
        for name, value in large.residue_identities().items():
            self.assertTrue(value, name)

    def test_exact_zero_with_normalization(self):
        n = algebra.NORMALIZATION
        zero = {EPSILON / n - EPSILON / n, n - 1 - 2 * EPSILON ** 2 * Q, 1 / n - 1 / (1 + 2 * EPSILON ** 2 * Q),
                2 * EPSILON / sp.sqrt(n) * (n - 1 - 2 * EPSILON ** 2 * Q)}
        nonzero = {(1 - 1 / n) * EPSILON ** 2, 1 + EPSILON / sp.sqrt(n), EPSILON ** 3 / n ** 2}
        for value in zero:
            self.assertTrue(algebra.coefficient_is_exactly_zero(value), value)
        for value in nonzero:
            self.assertFalse(algebra.coefficient_is_exactly_zero(value), value)
        self.assertTrue(algebra.is_exactly_zero(TraceSum.single(('Minv', 'Minv'), 2 * EPSILON ** 2)
                                                + TraceSum.scalar(1 - n)))


class AlpPoleSeriesTest(unittest.TestCase):
    def test_rewriting(self):
        self.assertTrue(large.pole_function_with_alp_mass_entry())

    def test_resolvent_series(self):
        for order in range(5):
            self.assertTrue(large.resolvent_series_check(order))

    def test_root_and_normalization(self):
        result = large.alp_pole_series(order=2)
        self.assertEqual(sp.expand(result['root'] - (A / T + 2 * C0 / A + 2 * C1 * T / A ** 2
                                                     - 4 * C0 ** 2 * T / A ** 3)), 0)
        self.assertEqual(sp.expand(result['N from <(M^{-1})^2>'] - result["N from D_a'"]), 0)
        self.assertEqual(sp.expand(result['N from <(M^{-1})^2>'] - (1 + 2 * C0 * T / A ** 2)), 0)
        self.assertEqual(sp.expand(result['coefficient of C in M^{-1}(C) (leading)'] + T / A), 0)

    def test_up_only(self):
        result = large.up_only_check()
        m_u, m1 = large.UP_QUARK_MASS, large.SINGLET_MASS_SQUARED
        self.assertEqual(result['A'], B0 * m_u)
        self.assertEqual(sp.expand(result['root'] - (B0 * m_u / T + 2 * B0 * m_u + 2 * m1 * T / 3)), 0)
        self.assertEqual(sp.expand(result['N'] - (1 + 2 * T)), 0)


# ---------------------------------------------------------------------------
# Every quantity of the estimate at large eps
# ---------------------------------------------------------------------------
_RESULTS = {}


def large_epsilon_results():
    if not _RESULTS:
        _RESULTS.update(large.all_results())
    return _RESULTS


KAPPA = large.KAPPA_SIGMA
P01P01 = (('P01', 'P01'),)


class TruncatedSeriesTest(unittest.TestCase):
    def test_series_identities(self):
        for order in (2, 3):
            checks = large.series_identities(order)
            self.assertTrue(checks['pole condition'])
            self.assertTrue(checks['defining relation of M^{-1}(C)'])
            self.assertTrue(checks['root agrees with alp_pole_series(order=2)'])

    def test_precision_of_products(self):
        exact = large.TruncatedSeries.monomial(1, -2)
        truncated = large.TruncatedSeries.monomial(1, 2, precision=4)
        self.assertEqual((exact * truncated).precision, 2)
        self.assertEqual(large.one_plus_power(truncated, -1).precision, 4)

    def test_cancellation_raises(self):
        def build(ctx):
            return ctx.alp_mass_squared - ctx.alp_mass_squared
        with self.assertRaises(large.InsufficientOrder):
            large.leading_of(build)

    def test_mass_trace_rewriting(self):
        rewritten = large.substitute_mass_traces(TraceSum.single(('S10',)) + TraceSum.single(('S02',)))
        expected = (TraceSum.scalar(KAPPA / B0 - A / (2 * B0)) + TraceSum.single(('S20',), 2))
        self.assertTrue((rewritten - expected).is_zero())


class PoleDataTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rows = large_epsilon_results()['pole']

    def test_powers(self):
        expected = {'m_a^2': 2, 'm_sigma^2': 2, 'm_a^2 - m_sigma^2': 2, 'M^{-1}_{m_a^2}(C)': -2, 'N(m_a^2)': 0,
                    'N - 1': -2, '1/N': 0, 'xi_a = 1/sqrt(N)': 0, 'v_a, meson matrix 2 eps M^{-1}/sqrt(N)': -1,
                    'v_a, ALP-field component 1/sqrt(N)': 0, 'residue (i): bare generator sum 1_9': 0,
                    'residue (ii): ALP-field entry 1 - 1/N': -2, 'residue (iii): mixed entry -2 eps M^{-1}/N': -1,
                    'residue (iv): meson-block pole -4 eps^2 M^{-1} (x) M^{-1}/N': -2,
                    'm_a^2 - m_sigma^2 (A = kappa_sigma)': 0, 'm_a^2 - m_sigma^2 (up-only)': 0}
        self.assertEqual({name: row['power'] for name, row in self.rows.items()}, expected)

    def test_coefficients(self):
        def coefficient(name):
            return self.rows[name]['coefficient']
        self.assertTrue((coefficient('m_a^2 - m_sigma^2') - TraceSum.scalar(A - KAPPA)).is_zero())
        self.assertTrue((coefficient('M^{-1}_{m_a^2}(C)') - TraceSum.single(('P01', 'Y'), -B0 / A)).is_zero())
        self.assertTrue((coefficient('residue (ii): ALP-field entry 1 - 1/N')
                         - TraceSum({P01P01: 2 * B0 ** 2 / A ** 2})).is_zero())
        self.assertTrue((coefficient('residue (iii): mixed entry -2 eps M^{-1}/N')
                         - TraceSum.single(('P01', 'Y'), 2 * B0 / A)).is_zero())
        self.assertTrue((coefficient('residue (iv): meson-block pole -4 eps^2 M^{-1} (x) M^{-1}/N')
                         - TraceSum({(('P01', 'Y1'), ('P01', 'Y2')): -4 * B0 ** 2 / A ** 2})).is_zero())
        self.assertTrue((coefficient('m_a^2 - m_sigma^2 (A = kappa_sigma)')
                         - TraceSum({P01P01: 2 * B0 ** 2 / A})).is_zero())

    def test_up_only_mass_difference(self):
        row = self.rows['m_a^2 - m_sigma^2 (up-only)']
        u = large.INVERSE_EPSILON
        expected = (2 * B0 * large.UP_QUARK_MASS
                    + (2 * large.SINGLET_MASS_SQUARED / 3 - 2 * large.QUARTIC * large.F_PI ** 2) * u ** 2)
        self.assertEqual(sp.expand(row['up-only series'] - expected), 0)


class VertexPowerTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.table = large_epsilon_results()['vertices']

    def test_count(self):
        self.assertEqual(len(self.table), 46)

    def test_lambda_pitfall(self):
        radial3, radial4 = ('radial',) * 3, ('radial',) * 4
        self.assertEqual(self.table[radial3][('generic', 'fixed momenta')], [-1, 3])
        self.assertEqual(self.table[radial3][('r <= 1', 'fixed momenta')], [-1])
        self.assertEqual(self.table[radial4][('generic', 'fixed momenta')], [0, 4])
        self.assertEqual(self.table[radial4][('r <= 1', 'fixed momenta')], [0])
        with_lambda, eliminated = large.radial_cubic_vertex_forms()
        self.assertTrue((with_lambda - eliminated).is_zero())
        self.assertTrue((with_lambda - TraceSum.single(('S30',), 6 * B0 * EPSILON ** 3 / large.F_PI)
                         - TraceSum.scalar(-6 * large.QUARTIC * large.F_PI / EPSILON)).is_zero())

    def test_kinetic_vertices_with_loop_momenta(self):
        self.assertEqual(self.table[('radial', 'ALP field', 'ALP field')][('generic', 'fixed momenta')], [1, 3])
        self.assertEqual(self.table[('radial', 'ALP field', 'ALP field')][('generic', 'loop momenta')], [3])
        self.assertEqual(self.table[('radial', 'radial', 'ALP field', 'ALP field')][('generic', 'loop momenta')], [4])
        self.assertEqual(self.table[('meson',) * 4][('generic', 'loop momenta')], [0, 2])

    def test_one_power_per_radial_or_alp_field_end(self):
        """Generic sources, fixed momenta: every nonvanishing pattern carries one power
        of eps per radial end and per ALP-field end, except the four patterns with
        lambda or with the kinetic sigma-a-a terms."""
        exceptions = {}
        vanishing = set()
        for kinds, row in self.table.items():
            exponents = row[('generic', 'fixed momenta')]
            if not exponents:
                vanishing.add(kinds)
            elif exponents != [sum(1 for kind in kinds if kind != 'meson')]:
                exceptions[kinds] = exponents
        self.assertEqual(exceptions, {('radial',) * 3: [-1, 3], ('radial',) * 4: [0, 4],
                                      ('radial', 'ALP field', 'ALP field'): [1, 3],
                                      ('radial', 'radial', 'ALP field', 'ALP field'): [2, 4]})
        self.assertEqual(vanishing, {('meson',) * 3, ('radial', 'meson', 'meson', 'meson'), ('meson',) * 5})


class LegFactorTest(unittest.TestCase):
    def test_net_powers(self):
        rows = {row['name']: row for row in large_epsilon_results()['legs']}
        self.assertEqual(rows['ALP-field component of an ALP-like end, xi_a']['net'],
                         {'1 at source, 0 at kinetic vertices': 1, '0 at source, 1 at kinetic vertices': 0})
        self.assertEqual(rows['meson-line residue (ii): ALP-field entry 1 - 1/N']['net'],
                         {'2 at source, 0 at kinetic vertices': 0, '1 at source, 1 at kinetic vertices': -1,
                          '0 at source, 2 at kinetic vertices': -2})
        self.assertEqual(rows['meson-line residue (iii): mixed entry']['net'],
                         {'1 at source, 0 at kinetic vertices': 0, '0 at source, 1 at kinetic vertices': -1})
        self.assertEqual(rows['meson-line residue (iv): meson-block pole']['power'], -2)
        self.assertEqual(rows['meson component of v_a']['power'], -1)


GRAPH_POWERS = {
    # name: (T, T/M0, contribution/M0) generic, then r <= 1; None = exact zero
    '1': ((3, 2, 4), None), '2': ((3, 2, 4), (3, 2, 4)), '3': ((1, 0, 2), (1, 0, 2)),
    '4': ((3, 2, 2), None), '5': ((3, 2, 2), (3, 2, 2)), '6': ((1, 0, 0), (1, 0, 0)),
    '7': ((3, 2, 2), None), '8': ((3, 2, 2), (3, 2, 2)), '9': ((1, 0, 0), (1, 0, 0)),
    '10': ((5, 4, 4), None), '11': ((5, 4, 4), None), '12': ((3, 2, 2), None),
    '13': ((5, 4, 4), (5, 4, 4)), '14': ((3, 2, 2), (3, 2, 2)), '15': ((3, 2, 2), (3, 2, 2)),
    '16': ((5, 4, 2), (1, 0, -2)), '17': ((5, 4, 2), None), '18': ((3, 2, 0), None),
    '19': ((5, 4, 2), None), '20': ((5, 4, 2), (5, 4, 2)), '21': ((3, 2, 0), (3, 2, 0)),
    '22': ((3, 2, 0), (3, 2, 0)), '23': ((3, 2, 0), (3, 2, 0)), '24': ((3, 2, 0), None),
    '25': ((3, 2, 0), (3, 2, 0)), '26': ((1, 0, -2), (1, 0, -2)),
}
SELF_ENERGY_POWERS = {
    # name: (T, dZ) generic, then r <= 1
    'sigma1': ((6, 4), (-2, -4)), 'sigma2': ((6, 4), None), 'sigma3': ((4, 2), None),
    'sigma4': ((6, 4), (6, 4)), 'sigma5': ((4, 2), (4, 2)), 'sigma6': ((2, 2), (2, 2)),
    'pi1': ((2, 0), (2, 0)), 'pi2': ((0, 0), (0, 0)), 'pi3': ((2, 0), (2, 0)),
    'pi4': ((2, 0), (2, 0)), 'pi5': ((0, -2), (0, -2)),
}


class GraphPowerTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.results = large_epsilon_results()

    def rows(self, region_index=0):
        result = {}
        for row in self.results['graphs'][region_index]:
            power = (lambda value: None if value is None else value['power'])
            if row['group'] == 'one_particle_irreducible':
                values = (power(row['T']), row['T/M0'], row['contribution/M0'])
            else:
                values = (power(row['T']), power(row['dZ']))
            result[(row['name'], row['sources'])] = None if values[0] is None else values
        return result

    def test_tree(self):
        self.assertEqual(self.results['tree']['power'], 1)
        expected = (TraceSum.single(('S10', 'l+', 'l-'), -B0 / (2 * large.F_PI))
                    + TraceSum.single(('S10', 'l-', 'l+'), -B0 / (2 * large.F_PI)))
        self.assertTrue((self.results['tree']['coefficient'] - expected).is_zero())

    def test_regions(self):
        """Regions of kappa_sigma/A in which every routing choice is constant."""
        bounds = [(region['lower'], region['upper']) for region in self.results['regions']]
        self.assertEqual(bounds, [(0, 1), (1, 2), (2, 5), (5, sp.oo)])
        self.assertEqual([region['branch'] for region in self.results['regions']],
                         [large.BRANCHES[0]] + [large.BRANCHES[1]] * 3)
        for region in self.results['regions']:
            self.assertEqual({name for name, (choices, _) in region['choices'].items() if choices},
                             {'13', '14', '15', '20', '21', '25', '26', 'sigma4', 'sigma5', 'sigma6'})
            self.assertFalse(any(marks for _, marks in region['choices'].values()))

    def test_graph_powers(self):
        for index in range(len(self.results['regions'])):
            rows = self.rows(index)
            for name, (generic, restricted) in {**GRAPH_POWERS, **SELF_ENERGY_POWERS}.items():
                self.assertEqual(rows[(name, 'generic')], generic, name)
                self.assertEqual(rows[(name, 'r <= 1')], restricted, name)

    def test_counts(self):
        names = {row['name'] for row in self.results['graphs'][0]}
        self.assertEqual(len([n for n in names if n.isdigit()]), 26)
        self.assertEqual(len([n for n in names if n.startswith('sigma')]), 6)
        self.assertEqual(len([n for n in names if n.startswith('pi') and not n.startswith('pi-')]), 5)

    def test_leading_coefficients_free_of_pion_mass(self):
        for rows in self.results['graphs']:
            for row in rows:
                for quantity in ('T', 'contribution', 'dZ'):
                    if row.get(quantity):
                        for value in row[quantity]['coefficient'].terms.values():
                            self.assertFalse(sp.sympify(value).has(algebra.M_PI), row['name'])

    def test_lambda_vertices_cancel_exactly(self):
        # graphs with V(sigma^3): no single trace <S_20> left from the elimination of lambda
        for row in self.results['graphs'][0]:
            if row['name'] in ('10', '16', 'sigma1'):
                for quantity in ('T', 'contribution', 'dZ'):
                    if row.get(quantity):
                        self.assertFalse(any(('S20',) in key for key in row[quantity]['coefficient'].terms))

    def test_sigma_tadpole_from_lambda_in_restricted_sources(self):
        row = next(r for r in self.results['graphs'][0] if r['name'] == 'sigma1' and r['sources'] == 'r <= 1')
        expected = TraceSum.scalar(36 * large.QUARTIC ** 2 * large.F_PI ** 2)
        self.assertTrue((row['T']['coefficient'] - expected).is_zero())

    def test_routing_maximum_coefficients(self):
        """Graph 13 and sigma4 at the leading power with the routing maximum."""
        rows = {row['name']: row for row in self.results['graphs'][0] if row['sources'] == 'generic'}
        pi2 = sp.pi ** 2
        expected = (TraceSum({(('S12',), algebra.canonical_word(('S02', 'l+', 'l-'))): -B0 ** 2 / (16 * pi2 * large.F_PI ** 3),
                              (('S12',), algebra.canonical_word(('S02', 'l-', 'l+'))): -B0 ** 2 / (16 * pi2 * large.F_PI ** 3)})
                    + TraceSum.single(('S02', 'l+', 'l-'), -3 * A * B0 / (16 * pi2 * large.F_PI ** 3))
                    + TraceSum.single(('S02', 'l-', 'l+'), -3 * A * B0 / (16 * pi2 * large.F_PI ** 3)))
        self.assertEqual(rows['13']['contribution']['power'], 5)
        self.assertTrue((rows['13']['contribution']['coefficient'] - expected).is_zero())
        expected = (TraceSum({(('S12',), ('S12',)): -B0 ** 2 / (8 * pi2 * A * large.F_PI ** 2)})
                    + TraceSum.single(('S12',), -B0 / (2 * pi2 * large.F_PI ** 2))
                    + TraceSum.scalar(-A / (2 * pi2 * large.F_PI ** 2)))
        self.assertTrue((rows['sigma4']['dZ']['coefficient'] - expected).is_zero())


class TotalEstimateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.results = large_epsilon_results()
        cls.total = cls.results['total']

    def test_ratio_power_and_graphs(self):
        for region_total in self.total:
            generic = region_total['generic']
            restricted = region_total['r <= 1']
            self.assertEqual(generic['R power'], 4)
            self.assertEqual(restricted['R power'], 4)
            self.assertEqual(set(generic['attaining']),
                             {'1', '2', '10', '11', '13', 'sigma1', 'sigma2', 'sigma4'})
            self.assertEqual(set(restricted['attaining']), {'2', '13', 'sigma4'})

    def test_generic_ratio_coefficient(self):
        """M_1^est/M_0 -> eps^4 [...]/(16 pi^2 f_pi^2), region kappa_sigma < A."""
        def trace(word):
            return sp.Symbol(r'\langle ' + large.word_latex(word) + r'\rangle')
        anticommutator = ('{+-}',)
        s02, s12, s11, s20, s30, s10 = (trace(('S02',) + anticommutator), trace(('S12',) + anticommutator),
                                        trace(('S11',) + anticommutator), trace(('S20',) + anticommutator),
                                        trace(('S30',) + anticommutator), trace(('S10',) + anticommutator))
        lone_s12, lone_s21, lone_s30 = trace(('S12',)), trace(('S21',)), trace(('S30',))
        kappa = KAPPA
        bracket = ((6 * A * s02 + A * s12 + 2 * B0 * s02 * lone_s12 + 2 * B0 * s11 * lone_s21
                    + 6 * B0 * s20 * lone_s30 + 3 * kappa * s30) / s10
                   - 4 * A - 4 * B0 * lone_s12 - 9 * B0 ** 2 * lone_s30 ** 2 / kappa
                   - B0 ** 2 * lone_s12 ** 2 / A - 2 * B0 ** 2 * lone_s21 ** 2 / A)
        expected = bracket / (16 * sp.pi ** 2 * large.F_PI ** 2)
        self.assertEqual(sp.simplify(self.total[0]['generic']['R coefficient'] - expected), 0)
        other = expected + (2 * B0 ** 2 * lone_s21 ** 2 / A - 2 * B0 ** 2 * lone_s21 ** 2 / kappa) / (
            16 * sp.pi ** 2 * large.F_PI ** 2)
        for region_total in self.total[1:]:
            self.assertEqual(sp.simplify(region_total['generic']['R coefficient'] - other), 0)

    def test_up_only_ratio(self):
        for region_total in self.total:
            value = large.up_only_ratio(region_total['r <= 1'])
            expected = -21 * B0 * large.UP_QUARK_MASS / (64 * sp.pi ** 2 * large.F_PI ** 2)
            self.assertEqual(sp.simplify(value - expected), 0)
            self.assertEqual(set(large.up_only_contributors(region_total['r <= 1'])), {'2', '13', 'sigma4'})

    def test_regions_differ_only_by_sigma2(self):
        one = self.total[0]['generic']['by graph']
        for region_total in self.total[1:]:
            two = region_total['generic']['by graph']
            different = {name for name in one if not (one[name] - two[name]).is_zero()}
            self.assertEqual(different, {'sigma2'})


SURVIVORS = {
    ('3', 0): ('ii', 'iii', 'iv'), ('4', 0): ('i',), ('5', 1): ('i',), ('6', 0): ('i', 'ii'),
    ('6', 1): ('i', 'ii'), ('7', 0): ('i',), ('8', 1): ('i',), ('9', 0): ('i', 'ii'), ('9', 1): ('i', 'ii'),
    ('12', 0): ('ii', 'iii'), ('14', 1): ('i', 'ii', 'iii', 'iv'), ('15', 0): ('i', 'ii', 'iii', 'iv'),
    ('15', 1): ('i', 'ii', 'iii', 'iv'), ('16', 0): ('i',), ('17', 0): ('i',), ('18', 0): ('i',),
    ('18', 1): ('ii', 'iii'), ('19', 0): ('i',), ('20', 0): ('i',), ('21', 0): ('i',), ('21', 1): ('ii', 'iii'),
    ('22', 1): ('i',), ('22', 2): ('i',), ('23', 1): ('i',), ('23', 2): ('i',), ('24', 0): ('i',),
    ('24', 2): ('ii', 'iii'), ('25', 0): ('i',), ('25', 2): ('ii', 'iii'), ('26', 0): ('i', 'ii'),
    ('26', 1): ('i', 'ii', 'iii'), ('26', 2): ('i', 'ii', 'iii'), ('sigma3', 0): ('i', 'ii', 'iii', 'iv'),
    ('sigma5', 1): ('i', 'ii', 'iii', 'iv'), ('sigma6', 0): ('i', 'ii', 'iii', 'iv'),
    ('sigma6', 1): ('i', 'ii', 'iii', 'iv'), ('pi2', 0): ('iv',), ('pi3', 0): ('i',), ('pi4', 1): ('i',),
    ('pi5', 0): ('i', 'ii'), ('pi5', 1): ('i', 'ii'),
}


class SurvivalTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rows = large_epsilon_results()['survival']

    def test_four_piece_decomposition_is_exact(self):
        self.assertTrue(all(row['decomposition exact'] for row in self.rows))

    def test_survivors(self):
        found = {(row['name'], line['line']): tuple(name.split()[0].strip('()') for name in line['survivors'])
                 for row in self.rows for line in row['lines']}
        self.assertEqual(found, SURVIVORS)

    def test_no_line_with_only_the_alp_field_entry(self):
        self.assertFalse(any(line['only ALP-field entry'] for row in self.rows for line in row['lines']))

    def test_all_generator_sum_term_vanishes_at_three_meson_vertices(self):
        combinations = {row['name']: row['combinations'] for row in self.rows}
        for name in ('6', '9', '26', 'pi5'):
            bare = tuple(large.RESIDUE_PIECES[0] for _ in next(iter(combinations[name])))
            self.assertNotIn(bare, combinations[name], name)

    def test_alp_field_entry_and_pole_cancel_by_two_powers(self):
        # entry 1 alone and -v_a v_a^T alone are each eps^2 larger than the combined ALP-field entry 1 - 1/N
        for row in self.rows:
            for line in row['lines']:
                piece = line['pieces'][large.RESIDUE_PIECES[1]]
                if not piece['nonzero'] or line['pieces'][large.RESIDUE_PIECES[2]]['nonzero']:
                    continue
                separate = line['separate pieces']
                for power in separate.values():
                    self.assertEqual(power, piece['power'] + 2, (row['name'], line['line']))

    def test_excluded_terms(self):
        excluded = {name: (value['measure']['power'], value['relative'])
                    for name, value in large_epsilon_results()['D4 terms'].items()}
        self.assertEqual(excluded, {'3': (1, 0), '15': (3, 2), 'pi2': (0, 0)})
        for name, value in large_epsilon_results()['D4 terms'].items():
            self.assertTrue(value['criterion'], name)
            self.assertFalse(any(letter in ('Minv', 'P01') or (letter[0] in 'SP' and letter[1:].isdigit()
                                                                 and int(letter[2:]) >= 1)
                                 for letter in value['letters']), name)
            self.assertNotIn('m_a', value['scales'], name)


if __name__ == '__main__':
    unittest.main()
