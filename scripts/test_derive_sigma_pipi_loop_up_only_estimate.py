"""Unit tests for derive_sigma_pipi_loop_up_only_estimate.py (exact symbolic identities only).

Run: PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_loop_up_only_estimate
(about one minute: the kept graphs and their routing choices are built once).
"""
import unittest

import sympy as sp

import derive_sigma_pipi_loop_up_only_estimate as up_only
from sigma_pipi_loop_algebra import B0, EPSILON, F_PI, TraceSum

M_U = up_only.UP_QUARK_MASS
A = up_only.ALP_MASS_COEFFICIENT
KAPPA = up_only.KAPPA_SIGMA
PI2 = sp.pi ** 2
_RESULTS = {}


def results():
    if not _RESULTS:
        _RESULTS.update(up_only.all_results())
    return _RESULTS


def symbol(text):
    return sp.Symbol('<' + text + '>')


class SourceTest(unittest.TestCase):
    """(a) Up-only sources and the reduction of traces."""

    def test_letters(self):
        expected = {'P01': (-M_U, 'Pu'), 'P11': (-M_U, 'Pu'), 'S02': (-M_U / 2, 'Pu'), 'S10': (M_U, 'Pu'),
                    'S12': (-M_U / 2, 'Pu'), 'P03': (M_U / 6, 'Pu'), 'S00': (1, 'S00'),
                    'S01': None, 'S11': None, 'P02': None, 'P00': None, 'P10': None, 'S20': None, 'S21': None,
                    'P20': None, 'S30': None}
        for name, image in expected.items():
            self.assertEqual(up_only.up_only_letter(name), image, name)

    def test_trace_reduction(self):
        self.assertEqual(up_only.up_only_trace(('Pu', 'l+', 'l-')), 2)
        self.assertEqual(up_only.up_only_trace(('Pu', 'l-', 'l+')), 0)
        self.assertEqual(up_only.up_only_trace(('S00', 'Pu')), M_U)
        self.assertEqual(up_only.up_only_trace(('Pu', 'l+', 'Pu', 'l-')), 0)
        self.assertEqual(up_only.up_only_trace(('l+', 'l-', 'l+', 'l-')), 4)
        self.assertEqual(up_only.up_only_trace(()), 3)
        with self.assertRaises(up_only.UnreducedTrace):
            up_only.up_only_trace(('S00', 'Pd'))

    def test_letters_of_the_kept_graphs(self):
        self.assertEqual(results()['letters'],
                         ['P01', 'P11', 'P20', 'S00', 'S01', 'S02', 'S10', 'S11', 'S12', 'S20', 'S21', 'S30'])


class DerivedInputTest(unittest.TestCase):
    """(b) A, kappa_sigma, M_0, trace ratios and single traces."""

    @classmethod
    def setUpClass(cls):
        cls.inputs = up_only.derived_inputs(results())

    def test_masses_and_tree(self):
        self.assertEqual(self.inputs['A = -2 B_0 <S02>'], B0 * M_U)
        self.assertEqual(self.inputs['kappa_sigma = B_0 (<S10> - 2 <S20>)'], B0 * M_U)
        self.assertEqual(sp.simplify(self.inputs['M_0'] + B0 * M_U * EPSILON / F_PI), 0)
        self.assertEqual(self.inputs['<S10 {l+,l-}>'], 2 * M_U)

    def test_ratios_and_traces(self):
        self.assertEqual(self.inputs['ratios'], {('S02',): sp.Rational(-1, 2), ('S11',): 0,
                                                 ('S12',): sp.Rational(-1, 2), ('S20',): 0, ('S30',): 0})
        self.assertEqual(self.inputs['single traces'], {('S12',): -M_U / 2, ('S21',): 0, ('S30',): 0})


class ContributionTest(unittest.TestCase):
    """(c) Generic eps^4 terms graph by graph, their up-only values and the other graphs."""

    def test_generic_coefficients(self):
        s10, s02, s12 = symbol('S10 {l+,l-}'), symbol('S02 {l+,l-}'), symbol('S12 {l+,l-}')
        lone_s12 = symbol('S12')
        expected = {
            '2': A * s12 / (16 * PI2 * F_PI ** 2 * s10),
            '13': (B0 * lone_s12 + 3 * A) * s02 / (8 * PI2 * F_PI ** 2 * s10),
            'sigma4': -(B0 ** 2 * lone_s12 ** 2 / A + 4 * B0 * lone_s12 + 4 * A) / (16 * PI2 * F_PI ** 2),
        }
        for label in up_only.RATIO_POINTS:
            contributions = results()['contributions'][label]
            self.assertEqual(set(contributions), {'1', '2', '10', '11', '13', 'sigma1', 'sigma2', 'sigma4'})
            for name, value in expected.items():
                self.assertEqual(sp.simplify(contributions[name]['generic'] - value), 0, (label, name))

    def test_up_only_values(self):
        expected = {'2': sp.Rational(-1, 32), '13': sp.Rational(-5, 32), 'sigma4': sp.Rational(-9, 64),
                    '1': 0, '10': 0, '11': 0, 'sigma1': 0, 'sigma2': 0}
        for label in up_only.RATIO_POINTS:
            contributions = results()['contributions'][label]
            for name, value in expected.items():
                self.assertEqual(sp.simplify(contributions[name]['up-only'] - value * up_only.UNIT_OF_RATIO), 0,
                                 (label, name))

    def test_other_graphs_below_fourth_power(self):
        nonzero = {'3': 2, '14': 2, '15': 2, '16': -2, '22': 0, 'sigma1': -4, 'sigma5': 2, 'sigma6': 2,
                   'pi1': 0, 'pi2': 0, 'pi3': 0, 'pi-1': 0, 'pi-2': 0, 'pi-3': 0, '2': 4, '13': 4, 'sigma4': 4}
        for label in up_only.RATIO_POINTS:
            powers = results()['powers'][label]
            self.assertEqual(len(powers), 26 + 6 + 5 + 5)
            found = {name: power['power'] for name, power in powers.items() if power['power'] is not None}
            self.assertEqual(found, nonzero, label)
            for name, power in powers.items():
                if name not in nonzero:
                    self.assertTrue(power['note'].startswith('zero (exact'), (label, name))


class TotalTest(unittest.TestCase):
    """(d) and (e): the sum, R, the branch of max(m_sigma, m_a) and the full total."""

    def test_sum_and_ratio(self):
        expected = -21 * B0 * M_U / (64 * PI2 * F_PI ** 2)
        for label in up_only.RATIO_POINTS:
            value = up_only.total_value(results()['contributions'][label])
            self.assertEqual(sp.simplify(value - expected), 0, label)
            self.assertEqual(sp.simplify(sp.Abs(value) + expected), 0, label)

    def test_branch_of_the_larger_mass(self):
        mass = results()['mass difference']
        self.assertEqual(mass['power'], 0)
        self.assertTrue((mass['coefficient'] - TraceSum.scalar(2 * B0 * M_U)).is_zero())
        generic = results()['mass difference generic']
        self.assertEqual(generic['power'], 2)
        self.assertTrue((generic['coefficient'] - TraceSum.scalar(A - KAPPA)).is_zero())
        self.assertEqual([region['branch'] for region in results()['regions'].values()],
                         list(up_only.large.BRANCHES))

    def test_full_total(self):
        expected = -21 * B0 * M_U / (64 * PI2 * F_PI ** 2)
        for label in up_only.RATIO_POINTS:
            total = results()['totals'][label]['r <= 1']
            self.assertEqual(total['R power'], 4)
            self.assertEqual(sp.simplify(up_only.large.up_only_ratio(total) - expected), 0, label)
            self.assertEqual(set(up_only.large.up_only_contributors(total)), {'2', '13', 'sigma4'})


if __name__ == '__main__':
    unittest.main()
