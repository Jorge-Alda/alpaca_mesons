"""Tests of the double-trace source derivation for vector masses and the KK pi current.

Run: .venv/bin/python -m unittest scripts.test_derive_double_trace_vector_mass_and_kkpi
"""
import sys
import unittest
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.linalg import expm

sys.path.insert(0, str(Path(__file__).resolve().parent))
import derive_double_trace_vector_mass_and_kkpi as derivation
from implementation import vector_sector

CHARGE = np.diag([2/3, -1/3, -1/3])


class SymbolicDerivation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = derivation.derive()

    def test_goldstone_expansion_has_only_even_orders(self):
        coefficients = self.result['expansion_coefficients']
        self.assertEqual(coefficients[1], 0)
        self.assertEqual(coefficients[3], 0)
        self.assertNotEqual(coefficients[2], 0)
        self.assertEqual(self.result['vacuum_check'], 0)

    def test_no_quadratic_vector_term(self):
        self.assertEqual(self.result['vector_degree'], 1)
        self.assertTrue(self.result['hessian_is_zero'])
        self.assertTrue(self.result['lagrangian_vanishes_without_photon'])

    def test_only_vector_trace_couples(self):
        for coefficient in self.result['off_diagonal_vector_coefficients']:
            self.assertEqual(sp.expand(coefficient), 0)
        for value in self.result['traceless_derivatives'].values():
            self.assertEqual(value, 0)

    def test_residue_shifts(self):
        shift = self.result['residue_shift']
        y, xl, xs = derivation.y_V, derivation.x_l, derivation.x_s
        self.assertEqual(shift['rho'], 0)
        self.assertEqual(sp.simplify(shift['omega_I'] - sp.sqrt(2)*y*(xl - xs)/3), 0)
        self.assertEqual(sp.simplify(shift['phi_I'] - y*(xl - xs)/3), 0)

    def test_kkpi_substitution_reproduces_insertions(self):
        for value in self.result['kkpi']['differences'].values():
            self.assertEqual(value, 0)

    def test_leading_electromagnetic_shift(self):
        shift = self.result['electromagnetic_mass_shift']
        self.assertEqual(shift['leading_shift'], [2*derivation.e**2*shift['residue']**2])

    def test_vector_pseudoscalar_bilinear_and_weak_current(self):
        self.assertEqual(self.result['vector_to_pseudoscalar'], 0)
        self.assertEqual(self.result['weak_current_trace'], 0)


class NumericalCrossChecks(unittest.TestCase):
    def test_residues_agree_with_vector_sector_module(self):
        rng = np.random.default_rng(7)
        FV, z, y, xl, xs = rng.normal(size=5)
        generators = [np.diag([1, -1, 0])/np.sqrt(2), np.diag([1, 1, 0])/np.sqrt(2), np.diag([0, 0, 1.])]
        module = vector_sector.source_residues(generators, np.diag([xl, xl, xs]), FV, z, y)
        symbolic = derivation.photon_residues()
        substitutions = {derivation.F_V: FV, derivation.z_V: z, derivation.y_V: y,
                         derivation.x_l: xl, derivation.x_s: xs}
        expected = [float(symbolic[name].subs(substitutions)) for name in ('rho', 'omega_I', 'phi_I')]
        np.testing.assert_allclose(module, expected, rtol=1e-12, atol=1e-14)

    def test_source_trace_is_even_in_goldstones_to_all_orders(self):
        rng = np.random.default_rng(11)
        matrix = rng.normal(size=(3, 3)) + 1j*rng.normal(size=(3, 3))
        phi = (matrix + matrix.conj().T)/2
        chi = np.diag(rng.uniform(.01, 1, size=3))

        def source_trace(field):
            u = expm(1j*field/(2*0.0924))
            ud = u.conj().T
            chi_plus = ud@chi@ud + u@chi@u
            f_plus = -(u@CHARGE@ud + ud@CHARGE@u)
            return np.trace(chi_plus@f_plus)
        for scale in (.01, .05, .2):
            self.assertAlmostEqual(abs(source_trace(scale*phi) - source_trace(-scale*phi)), 0, places=12)


if __name__ == '__main__':
    unittest.main()
