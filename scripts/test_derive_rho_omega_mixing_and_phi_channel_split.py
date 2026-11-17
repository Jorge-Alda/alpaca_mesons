"""Tests of the rho0-omega mixing and phi -> K Kbar charged/neutral split derivation.

Run: .venv/bin/python -m unittest scripts.test_derive_rho_omega_mixing_and_phi_channel_split
"""
import math
import sys
import unittest
from pathlib import Path

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import derive_rho_omega_mixing_and_phi_channel_split as derivation
from implementation import vector_sector

d = derivation


class SymbolicDerivation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = derivation.derive()

    def test_isospin_limit_reproduces_mass_records(self):
        self.assertTrue(all(value == 0 for value in self.result['isospin_limit_differences']))

    def test_off_diagonal_entries(self):
        self.assertEqual(sp.simplify(self.result['rho_omega_ideal']
                                     - (d.x_u - d.x_d)*(-2*d.mass_splitting_coefficient
                                                        + d.double_trace_mass_coefficient/2)), 0)
        self.assertEqual(sp.simplify(self.result['rho_phi_ideal']
                                     - d.double_trace_mass_coefficient*(d.x_u - d.x_d)/(2*sp.sqrt(2))), 0)
        # Neither entry depends on Delta_V or on x_s.
        for entry in (self.result['rho_omega_ideal'], self.result['rho_phi_ideal']):
            self.assertNotIn(d.singlet_shift, entry.free_symbols)
            self.assertNotIn(d.x_s, entry.free_symbols)

    def test_kstar_identity_and_rho_degeneracy(self):
        self.assertEqual(self.result['kstar_identity'], 0)
        self.assertEqual(self.result['rho_degeneracy'], 0)
        self.assertTrue(all(value == 0 for value in self.result['diagonal_shifts']))

    def test_pole_component(self):
        component, (a, b, delta) = derivation.first_order_pole_component()
        self.assertEqual(sp.simplify(component - delta/(b - a)), 0)

    def test_residue_isospin_breaking(self):
        breaking = self.result['residue_isospin_breaking']
        difference = d.x_u - d.x_d
        self.assertEqual(sp.simplify(breaking['rho'] - d.z_V*difference/(3*sp.sqrt(2))), 0)
        self.assertEqual(sp.simplify(breaking['omega_I'] - (d.z_V + d.y_V)*difference/sp.sqrt(2)), 0)
        self.assertEqual(sp.simplify(breaking['phi_I'] - d.y_V*difference/2), 0)
        # F_V gives no isospin breaking (<Q> structure is spurion independent).
        for value in breaking.values():
            self.assertNotIn(d.F_V, value.free_symbols)

    def test_phi_traces(self):
        self.assertEqual(self.result['phi_sum_first_order'], 0)
        c = self.result['component_symbol']
        charged, neutral = self.result['phi_traces']
        self.assertEqual(sp.simplify(charged - neutral - sp.sqrt(2)*c), 0)

    def test_residues_match_implementation(self):
        values = {d.x_u: 0.026272, d.x_d: 0.04696, d.x_s: 0.943684,
                  d.F_V: 0.169, d.z_V: -0.171, d.y_V: 0.07}
        generators = [np.array(sp.N(g), dtype=float) for g in
                      (d.GENERATOR_RHO, d.GENERATOR_OMEGA_IDEAL, d.GENERATOR_PHI_IDEAL)]
        implemented = vector_sector.source_residues(generators, np.diag([0.026272, 0.04696, 0.943684]),
                                                    0.169, -0.171, 0.07)
        symbolic = [float(self.result['residues'][name].subs(values)) for name in ('rho', 'omega_I', 'phi_I')]
        np.testing.assert_allclose(symbolic, implemented, rtol=0, atol=1e-14)


class EndpointNumbers(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = derivation.derive()
        cls.inputs = derivation.endpoint_inputs()
        cls.numbers = derivation.evaluate(cls.result, cls.inputs)

    def test_first_order_matches_exact_complex_eigenvectors(self):
        masses = dict(zip(('rho', 'omega', 'phi'), self.inputs['sector'].masses))
        totals = self.inputs['totals']
        diagonal = [masses['rho']**2 - 1j*masses['rho']*totals['rho0'],
                    masses['omega']**2 - 1j*masses['omega']*totals['omega'],
                    masses['phi']**2 - 1j*masses['phi']*totals['phi']]
        matrix = np.diag(diagonal).astype(complex)
        matrix[0, 1] = matrix[1, 0] = self.numbers['rho_omega_physical_gev2']
        matrix[0, 2] = matrix[2, 0] = self.numbers['rho_phi_physical_gev2']
        eigenvalues, eigenvectors = np.linalg.eig(matrix)
        for column, key in ((1, 'rho_in_omega'), (2, 'rho_in_phi')):
            index = int(np.argmin(abs(eigenvalues - diagonal[column])))
            exact = eigenvectors[0, index]/eigenvectors[column, index]
            first_order = complex(*self.numbers[key])
            self.assertLess(abs(exact - first_order), 0.05*abs(first_order))

    def test_real_diagonalization_is_not_perturbative_for_rho_omega(self):
        self.assertGreater(abs(self.numbers['rho_in_omega_zero_width']), 0.3)
        self.assertLess(self.numbers['rho_in_omega_modulus'], 0.06)

    def test_kstar_identity_numerically(self):
        kappa = self.inputs['kappa']
        xu, xd, _ = self.inputs['spurion']
        self.assertAlmostEqual(self.numbers['rho_omega_ideal_gev2'] - kappa*(xu - xd)/2,
                               self.numbers['kstar_plus_minus_zero_mass_squared_gev2'], places=14)

    def test_phi_split_and_rotated_residuals(self):
        self.assertAlmostEqual(self.numbers['phi_ratio_change'], self.numbers['phi_ratio_change_linear'],
                               delta=0.05*abs(self.numbers['phi_ratio_change_linear']))
        rotated = np.array(self.numbers['phi_shifted_symmetric_residuals'])
        cholesky = np.array(self.numbers['phi_shifted_cholesky_residuals'])
        self.assertAlmostEqual(rotated @ rotated, cholesky @ cholesky, places=10)
        covariance = np.array([[2.5e-5, -1.44e-5], [-1.44e-5, 1.6e-5]])
        difference = np.array(self.numbers['phi_shifted_predictions']) - np.array([0.491, 0.339])
        self.assertAlmostEqual(difference @ np.linalg.solve(covariance, difference), rotated @ rotated, places=8)
        self.assertAlmostEqual(self.numbers['phi_endpoint_chi2'], 17.878514018128385, places=8)

    def test_projection_remainders(self):
        decomposition = self.numbers['phi_residual_decomposition']
        weight = np.linalg.inv(np.array([[2.5e-5, -1.44e-5], [-1.44e-5, 1.6e-5]]))
        predicted = np.array([0.5084953, 0.3355123])
        difference = predicted - np.array([0.491, 0.339])
        # Brute-force one-parameter minimum along the common-normalization direction.
        scales = np.linspace(-0.05, 0.05, 200001)
        values = [(difference + s*predicted) @ weight @ (difference + s*predicted) for s in scales[::100]]
        self.assertAlmostEqual(min(values), decomposition['remaining_after_common_normalization'], delta=1e-3)
        self.assertGreater(decomposition['removed_by_common_normalization'], decomposition['removed_by_split'])
        two_directions = np.column_stack((predicted, [predicted[0], -predicted[1]]))
        coefficients = np.linalg.solve(two_directions.T @ weight @ two_directions, -two_directions.T @ weight @ difference)
        remainder = difference + two_directions @ coefficients
        self.assertLess(remainder @ weight @ remainder, 1e-20)

    def test_omega_pipi_order_of_magnitude(self):
        for value in self.numbers['omega_pipi_branching'].values():
            self.assertGreater(value, 0.0153)
            self.assertLess(value, 0.05)


if __name__ == '__main__':
    unittest.main()
