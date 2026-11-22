"""Unit tests for derive_alp_loop_crossing_table.py (units GeV).

Run: .venv/bin/python -m unittest discover -s scripts -p test_derive_alp_loop_crossing_table.py
"""
import unittest

import numpy as np
import sympy as sp
from scipy.integrate import quad

from derive_alp_loop_crossing_table import (
    ALP_MASSES, PHYSICAL_INPUTS, QUARK_MASS_SETS, bubble_function_at_threshold,
    crossing_alp_line_in_alp_entry, crossing_alp_line_in_meson_entry, crossing_alp_meson_bubble,
    crossing_alp_pair_bubble, crossing_table, crossing_untuned_alp_line_in_meson_entry,
    low_energy_constant_B0, meson_line_squared_masses, meson_loop_ratios,
    ratio_alp_line_in_alp_entry, ratio_alp_line_in_meson_entry, ratio_alp_meson_bubble,
    ratio_alp_pair_bubble, ratio_untuned_alp_line_in_meson_entry, source_coefficients,
    untuned_alp_squared_mass)
from derive_alp_loop_power_counting import (
    bubble_to_tree_ratios, loop_function_A0, regime_map, unmixed_closed_forms)

F_PI = PHYSICAL_INPUTS['pion_decay_constant']
SCALE = PHYSICAL_INPUTS['renormalization_scale']
PRIMARY = QUARK_MASS_SETS['primary']


def primary_squared_masses():
    return meson_line_squared_masses(PRIMARY, low_energy_constant_B0(PRIMARY))


class BubbleFunctionTests(unittest.TestCase):
    def test_threshold_closed_form_against_integral(self):
        for first, second in ((0.1, 0.1078), (1.0, 0.7083), (2.0, 2.0), (0.3, 0.3)):
            s = (first+second)**2
            integrand = lambda x: -np.log(abs(x*first**2+(1-x)*second**2-x*(1-x)*s)/SCALE**2)
            numerical, _ = quad(integrand, 0, 1, points=[second/(first+second)], limit=200)
            self.assertAlmostEqual(bubble_function_at_threshold(first, second, SCALE),
                                   numerical, places=8)

    def test_equal_mass_threshold(self):
        for mass in (0.1, 0.5, 2.0):
            self.assertAlmostEqual(bubble_function_at_threshold(mass, mass, SCALE),
                                   2-np.log(mass**2/SCALE**2), places=12)


class ClosedFormTests(unittest.TestCase):
    """Crossings against the regime_map and unmixed closed forms of the power-counting script."""

    def setUp(self):
        _, self.crossings, self.symbols = regime_map()

    def test_class_i_against_regime_map(self):
        sym = self.symbols
        for coupling in (1, 3):
            for alp_mass in ALP_MASSES:
                expected = float(self.crossings['alp_line_meson_entry_fixed_mass'].subs(
                    {sym['f']: F_PI, sym['kappa']: coupling**2, sym['m_a']: alp_mass}))
                self.assertAlmostEqual(crossing_alp_line_in_meson_entry(
                    coupling, 0.0, alp_mass, F_PI, SCALE, False), expected, places=10)

    def test_class_i_against_unmixed_ratio_and_A0(self):
        forms = unmixed_closed_forms()
        sym = forms['symbols']
        ratio = forms['alp_line_in_meson_entry_ratio'][(0, 0)]
        for alp_mass in (0.3, 2.0):
            loop = loop_function_A0(alp_mass**2, SCALE**2)
            value = float(ratio.subs({sym['axial'][0]: 3, sym['pseudoscalar'][0]: 0,
                                      sym['f_a']: F_PI/0.7}).replace(
                sym['A0'], sp.Lambda(sp.Symbol('y'), loop)))
            self.assertAlmostEqual(ratio_alp_line_in_meson_entry(
                0.7, 3, 0.0, alp_mass, F_PI, SCALE, True), abs(value), places=10)

    def test_bubbles_against_regime_map(self):
        sym = self.symbols
        squared = primary_squared_masses()
        for coupling in (1, 3):
            for q in (0, 2):
                common = {sym['f']: F_PI, sym['kappa']: coupling**2, sym['c_AP']: coupling,
                          sym['m_qq']: np.sqrt(squared[(q, q)]), sym['B0_loop']: 1}
                self.assertAlmostEqual(crossing_alp_meson_bubble(
                    coupling, 0.0, squared[(q, q)], 0.5, F_PI, SCALE, False),
                    float(self.crossings['bubble_aP'].subs(common)), places=10)
                self.assertAlmostEqual(crossing_alp_pair_bubble(
                    coupling, 0.0, squared[(q, q)], 0.5, F_PI, SCALE, False),
                    float(self.crossings['bubble_aa'].subs(common)), places=10)

    def test_bubble_ratios_against_power_counting(self):
        ratios = bubble_to_tree_ratios()
        sym = ratios['symbols']
        squared = primary_squared_masses()
        B0 = low_energy_constant_B0(PRIMARY)
        bubble = bubble_function_at_threshold(0.3, np.sqrt(squared[(2, 2)]), SCALE)
        values = {sym['f']: F_PI, sym['f_a']: F_PI/0.4, sym['B_0']: B0, sym['m_q']: PRIMARY[2],
                  sym['c_AP']: 3, sym['kappa']: 9, sym['B0_loop']: bubble}
        self.assertAlmostEqual(ratio_alp_meson_bubble(0.4, 3, 0.0, squared[(2, 2)], 0.3, F_PI,
                                                      SCALE, True),
                               float(ratios['aP'].subs(values)), places=10)
        values[sym['B0_loop']] = bubble_function_at_threshold(0.3, 0.3, SCALE)
        self.assertAlmostEqual(ratio_alp_pair_bubble(0.4, 3, 0.0, squared[(2, 2)], 0.3, F_PI,
                                                     SCALE, True),
                               float(ratios['aa'].subs(values)), places=10)

    def test_class_iv_against_unmixed_alp_entry(self):
        forms = unmixed_closed_forms()
        sym = forms['symbols']
        B0 = low_energy_constant_B0(PRIMARY)
        squared = primary_squared_masses()
        for convention in ('universal', 'up_only'):
            axial, pseudoscalar = source_coefficients(convention, 3)
            substitution = {sym['f']: F_PI, sym['f_a']: F_PI/0.5, sym['B_0']: B0}
            substitution.update(dict(zip(sym['masses'], PRIMARY)))
            substitution.update(dict(zip(sym['axial'], axial)))
            substitution.update(dict(zip(sym['pseudoscalar'], pseudoscalar)))
            alp_mass = 0.3
            loop = loop_function_A0(alp_mass**2, SCALE**2)
            correction = float(forms['alp_line_in_alp_entry'].subs(substitution).replace(
                sym['A0'], sp.Lambda(sp.Symbol('y'), loop)))
            self.assertAlmostEqual(ratio_alp_line_in_alp_entry(
                0.5, axial, pseudoscalar, squared, alp_mass, F_PI, SCALE, True),
                abs(correction)/alp_mass**2, places=10)

    def test_untuned_against_regime_map(self):
        sym = self.symbols
        squared = primary_squared_masses()
        for coupling in (1, 3):
            axial, pseudoscalar = source_coefficients('universal', coupling)
            total = sum(squared[(q, q)] for q in range(3))
            expected = float(self.crossings['alp_line_meson_entry_quark_mass'].subs(
                {sym['f']: F_PI, sym['kappa']: coupling**2,
                 sym['sum_kappa_m2']: coupling**2*total}))
            self.assertAlmostEqual(crossing_untuned_alp_line_in_meson_entry(
                0, axial, pseudoscalar, squared, F_PI, SCALE, False), expected, places=10)
            # ALP line in M^2_aa equals the untuned case (i) for universal couplings
            self.assertAlmostEqual(crossing_alp_line_in_alp_entry(
                axial, pseudoscalar, squared, 0.7, F_PI, SCALE, False), expected, places=10)


class CrossingTests(unittest.TestCase):
    def test_ratio_is_one_at_crossing(self):
        squared = primary_squared_masses()
        for coupling in (1, 3):
            for keep in (False, True):
                axial, pseudoscalar = source_coefficients('universal', coupling)
                for alp_mass in ALP_MASSES:
                    r = crossing_alp_line_in_meson_entry(coupling, 0.0, alp_mass, F_PI, SCALE, keep)
                    self.assertAlmostEqual(ratio_alp_line_in_meson_entry(
                        r, coupling, 0.0, alp_mass, F_PI, SCALE, keep), 1.0, places=10)
                    for q in (0, 2):
                        r = crossing_alp_meson_bubble(coupling, 0.0, squared[(q, q)], alp_mass,
                                                      F_PI, SCALE, keep)
                        self.assertAlmostEqual(ratio_alp_meson_bubble(
                            r, coupling, 0.0, squared[(q, q)], alp_mass, F_PI, SCALE, keep),
                            1.0, places=10)
                        r = crossing_alp_pair_bubble(coupling, 0.0, squared[(q, q)], alp_mass,
                                                     F_PI, SCALE, keep)
                        self.assertAlmostEqual(ratio_alp_pair_bubble(
                            r, coupling, 0.0, squared[(q, q)], alp_mass, F_PI, SCALE, keep),
                            1.0, places=10)
                    r = crossing_alp_line_in_alp_entry(axial, pseudoscalar, squared, alp_mass,
                                                       F_PI, SCALE, keep)
                    self.assertAlmostEqual(ratio_alp_line_in_alp_entry(
                        r, axial, pseudoscalar, squared, alp_mass, F_PI, SCALE, keep),
                        1.0, places=10)
                r = crossing_untuned_alp_line_in_meson_entry(0, axial, pseudoscalar, squared,
                                                             F_PI, SCALE, keep)
                self.assertAlmostEqual(ratio_untuned_alp_line_in_meson_entry(
                    r, 0, axial, pseudoscalar, squared, F_PI, SCALE, keep), 1.0, places=8)
                if keep:  # smallest root: ratio below one just before it
                    self.assertLess(ratio_untuned_alp_line_in_meson_entry(
                        0.99*r, 0, axial, pseudoscalar, squared, F_PI, SCALE, keep), 1.0)

    def test_class_i_independent_of_B0(self):
        reference = crossing_table(PRIMARY, 'universal', 3, True)
        for masses in (QUARK_MASS_SETS['alternative'], (0.001, 0.009, 0.2)):
            rows = crossing_table(masses, 'universal', 3, True)
            for label in ('(i) ALP line, uu entry', '(i) ALP line, ss entry'):
                np.testing.assert_allclose(rows[label], reference[label], rtol=1e-14)

    def test_up_only_switches_off_strange_classes(self):
        rows = crossing_table(PRIMARY, 'up_only', 1, False)
        for label in ('(i) ALP line, ss entry', '(ii) aP bubble, s', '(iii) aa bubble, s'):
            self.assertTrue(all(np.isinf(v) for v in rows[label]))
        np.testing.assert_allclose(rows['(iv) ALP line in M2_aa'], rows['(iii) aa bubble, u'])

    def test_physical_inputs(self):
        B0 = low_energy_constant_B0(PRIMARY)
        self.assertAlmostEqual(B0, 0.13498**2/0.0069, places=12)
        squared = primary_squared_masses()
        self.assertAlmostEqual(squared[(0, 0)]+squared[(1, 1)], 2*0.13498**2, places=12)
        self.assertAlmostEqual(squared[(0, 1)], 0.13498**2, places=12)

    def test_untuned_alp_mass(self):
        squared = primary_squared_masses()
        axial, pseudoscalar = source_coefficients('universal', 1)
        self.assertAlmostEqual(untuned_alp_squared_mass(1.0, axial, pseudoscalar, squared),
                               2*sum(squared[(q, q)] for q in range(3)), places=12)

    def test_meson_loops_independent_of_decay_constant_ratio(self):
        ratios = meson_loop_ratios(0.13498**2, F_PI, SCALE)
        self.assertAlmostEqual(ratios['size_estimate'],
                               0.13498**2/(16*np.pi**2*F_PI**2), places=14)
        self.assertLess(ratios['bubble_at_threshold'], 1.0)


if __name__ == '__main__':
    unittest.main()
