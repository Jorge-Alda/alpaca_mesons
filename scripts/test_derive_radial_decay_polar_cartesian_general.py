"""Unit tests: polar and Cartesian sigma -> P_i P_j vertices for general sources (GeV units).

Run: .venv/bin/python -m unittest scripts.test_derive_radial_decay_polar_cartesian_general
"""
import unittest

import numpy as np
import sympy as sp

from scripts.derive_radial_decay_polar_cartesian_general import (
    B_0, f_a, f_pi, quartic, angular_row, compare_vertices, derive_coordinate_map,
    derive_general_cubic_lagrangians, derive_vertex_difference, direct_on_shell_vertex,
    first_and_second_derivatives, only_up_parameters, pole_equation, pole_normalization,
    pole_shift_matrix, pole_vectors, random_general_parameters, trace)


class SymbolicDerivationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = derive_general_cubic_lagrangians()
        cls.trace_scalar = cls.data['trace_scalar']

    def test_coordinate_map_and_flat_metric(self):
        sigma_series, angular_series, metric_difference = derive_coordinate_map()
        s, p = sp.symbols('s p', real=True)
        self.assertEqual(sp.expand(sigma_series - (s + p ** 2 / (2 * f_a) - s * p ** 2 / (2 * f_a ** 2))), 0)
        self.assertEqual(sp.expand(angular_series - (p - s * p / f_a + s ** 2 * p / f_a ** 2 - p ** 3 / (3 * f_a ** 2))), 0)
        self.assertEqual(metric_difference, sp.zeros(2, 2))

    def test_stationarity_and_radial_mass(self):
        t = self.trace_scalar
        self.assertEqual(sp.simplify(self.data['mass_phi_squared'] - (B_0 * f_pi ** 2 * t[(1, 0)] / f_a ** 2 - quartic * f_a ** 2)), 0)
        expected = 2 * quartic * f_a ** 2 + B_0 * f_pi ** 2 / f_a ** 2 * (t[(1, 0)] - 2 * t[(2, 0)])
        self.assertEqual(sp.simplify(self.data['radial_mass_squared'] - expected), 0)
        self.assertEqual(sp.simplify(self.data['cartesian_radial_mass_squared'] - expected), 0)

    def test_radial_mixing_requires_s11_and_p10(self):
        # <S_01> also mixes s and p in Cartesian coordinates; it vanishes by angular stationarity.
        expected = B_0 * f_pi ** 2 * (self.trace_scalar[(1, 1)] - self.trace_scalar[(0, 1)]) / f_a ** 2
        self.assertEqual(sp.simplify(self.data['cartesian_s_p_mixing'] - expected), 0)
        nonzero = {k: v for k, v in self.data['cartesian_s_meson_mixing'].items() if v != 0}
        self.assertEqual(list(nonzero), [(1, 0)])

    def test_general_cartesian_cubic_coefficients(self):
        t = self.trace_scalar
        expected_s_p_p = B_0 * f_pi ** 2 / f_a ** 3 * (-t[(1, 0)] / 2 + t[(2, 0)] - 2 * t[(0, 2)] + t[(1, 2)]) - quartic * f_a
        self.assertEqual(sp.simplify(self.data['cartesian_s_p_p'] - expected_s_p_p), 0)
        s_p_meson = {k: v for k, v in self.data['cartesian_s_p_meson'].items() if v != 0}
        self.assertEqual(s_p_meson, {(0, 1): -B_0 * f_pi / f_a ** 2, (1, 1): B_0 * f_pi / f_a ** 2})
        s_meson_pair = {k: v for k, v in self.data['cartesian_s_meson_pair'].items() if v != 0}
        self.assertEqual(s_meson_pair, {(1, 0): -B_0 / (2 * f_a)})

    def test_polar_cubic_coefficients(self):
        self.assertEqual(sp.simplify(self.data['polar_sigma_a_a'] - B_0 * f_pi ** 2 * self.trace_scalar[(1, 2)] / f_a ** 3), 0)
        polar_a_meson = {k: v for k, v in self.data['polar_sigma_a_meson'].items() if v != 0}
        self.assertEqual(polar_a_meson, {(1, 1): B_0 * f_pi / f_a ** 2})

    def test_difference_vanishes_by_pole_and_radial_mass_relations(self):
        *_, reduced = derive_vertex_difference()
        self.assertEqual(reduced, 0)


class NumericalGeneralSourcesTest(unittest.TestCase):
    def test_angular_row_equals_pole_equation_off_shell(self):
        parameters = random_general_parameters()
        for momentum_squared in (0.05, 0.4, 1.3):
            self.assertAlmostEqual(angular_row(parameters, momentum_squared), pole_equation(parameters, momentum_squared), places=12)

    def test_random_sources(self):
        for seed in (2026, 7, 11):
            with self.subTest(seed=seed):
                results = compare_vertices(random_general_parameters(seed))
                self.assertTrue(np.all(results['masses_squared'] > 0))
                self.assertGreater(results['radial_mass_squared'], 0)
                self.assertGreaterEqual(len(results['mixed_poles']), 5)
                self.assertLess(results['tadpole'], 1e-14)
                self.assertLess(results['hessian_difference'], 1e-14)
                self.assertLess(results['radial_mixing'], 1e-14)
                self.assertLess(results['radial_mass_difference'], 1e-14)
                self.assertLess(results['pole_equation'], 1e-9)
                self.assertLess(results['pole_direction'], 1e-9)
                self.assertLess(results['normalization'], 1e-9)
                self.assertLess(results['polar_direct_vs_cartesian_direct'], 1e-12)
                self.assertLess(results['polar_direct_vs_trace'], 1e-10)
                self.assertLess(results['cartesian_direct_vs_trace'], 1e-10)

    def test_only_up_model(self):
        parameters = only_up_parameters()
        sources = parameters['sources']
        self.assertTrue(np.allclose(sources['P'][(1, 1)] - sources['P'][(0, 1)], 0))
        self.assertAlmostEqual(trace(sources['S'][(1, 0)]) - 2 * trace(sources['S'].get((2, 0), np.zeros((3, 3))))
                               + 4 * trace(sources['S'][(0, 2)]) - 2 * trace(sources['S'][(1, 2)]), 0, places=15)
        results = compare_vertices(parameters)
        self.assertLess(results['polar_direct_vs_cartesian_direct'], 1e-12)
        self.assertLess(results['polar_direct_vs_trace'], 1e-9)
        self.assertLess(results['cartesian_direct_vs_trace'], 1e-9)

    def test_general_terms_are_needed(self):
        """Dropping the P_11 - P_01 and trace terms (only-up form) fails for general sources."""
        parameters = random_general_parameters()
        results = compare_vertices(parameters)
        B0, fpi, fa = parameters['B_0'], parameters['f_pi'], parameters['f_a']
        masses_squared = results['masses_squared']
        i, j = results['mixed_poles'][0], results['mixed_poles'][-1]
        omega_i = pole_shift_matrix(parameters, masses_squared[i])
        omega_j = pole_shift_matrix(parameters, masses_squared[j])
        only_up_form = ((2 * parameters['lambda'] * fa
                         + B0 / (2 * fa) * trace(parameters['sources']['S'][(1, 0)] @ (omega_i @ omega_j + omega_j @ omega_i)))
                        / np.sqrt(pole_normalization(parameters, masses_squared[i]) * pole_normalization(parameters, masses_squared[j])))
        direct = results['tables']['cartesian_direct'][i, j]
        self.assertGreater(abs(only_up_form - direct) / results['vertex_scale'], 1e-3)

    def test_stationarity_is_needed(self):
        parameters = random_general_parameters()
        parameters['m_phi_squared'] += 0.05
        vertices = {}
        for coordinates in ('polar', 'cartesian'):
            _, hessian = first_and_second_derivatives(parameters, coordinates)
            masses_squared, vectors = pole_vectors(hessian)
            vertices[coordinates] = direct_on_shell_vertex(parameters, coordinates, -hessian[0, 0], masses_squared, vectors, 0, 0)
        self.assertGreater(abs(vertices['polar'] - vertices['cartesian']), 1e-3 * abs(vertices['polar']))


if __name__ == '__main__':
    unittest.main()
