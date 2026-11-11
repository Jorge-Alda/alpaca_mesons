"""Unit tests of the axial-vector nonet flavor conventions.

Run: .venv/bin/python -m unittest scripts.implementation.test_axial_nonet_conventions
"""
import unittest
import sympy as sp
from scripts.derive_axial_nonet_conventions import (
    alda_exponent_comparison, axial_field_symbols, derive, gell_mann_matrices,
    ideal_generators, mixing_rotation, nonet_generators, trace_gram_matrix)


class AxialNonetConventionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = derive()
        cls.fields = axial_field_symbols()

    def test_gell_mann_normalization(self):
        lambdas = gell_mann_matrices()
        self.assertEqual(trace_gram_matrix(lambdas), 2 * sp.eye(8))
        for lam in lambdas:
            self.assertEqual(sp.trace(lam), 0)
            self.assertEqual(lam.H, lam)

    def test_generator_normalization_and_orthogonality(self):
        self.assertEqual(self.result['generator_gram'], sp.eye(9))
        self.assertEqual(sp.trace(nonet_generators()[0] ** 2), 1)

    def test_completeness_relation(self):
        # sum_i (T_i)_{ab} (T_i)_{cd} = delta_ad delta_bc for a U(3) nonet.
        generators = nonet_generators()
        for a in range(3):
            for b in range(3):
                for c in range(3):
                    for d in range(3):
                        total = sp.simplify(sum(T[a, b] * T[c, d] for T in generators))
                        self.assertEqual(total, int(a == d and b == c))

    def test_explicit_matrix_entries(self):
        f = self.fields
        A = self.result['explicit_matrix']
        self.assertEqual(sp.simplify(A[0, 0] - (f['a1_0'] + f['f1_N']) / sp.sqrt(2)), 0)
        self.assertEqual(sp.simplify(A[1, 1] - (f['f1_N'] - f['a1_0']) / sp.sqrt(2)), 0)
        self.assertEqual(A[2, 2], f['f1_S'])
        self.assertEqual((A[0, 1], A[1, 0]), (f['a1_p'], f['a1_m']))
        self.assertEqual((A[0, 2], A[2, 0]), (f['K1_p'], f['K1_m']))
        self.assertEqual((A[1, 2], A[2, 1]), (f['K1_0'], f['K1bar_0']))

    def test_matrix_from_generators(self):
        self.assertEqual(self.result['matrix_from_generators_minus_explicit'], sp.zeros(3))

    def test_ideal_basis(self):
        T_N, T_S = ideal_generators()
        self.assertEqual(self.result['ideal_gram'], sp.eye(2))
        self.assertEqual(self.result['ideal_projection_minus_explicit'], sp.zeros(3))
        # Positive strange phase: both ideal states overlap positively with the singlet.
        self.assertEqual(self.result['ideal_singlet_overlaps'],
                         sp.Matrix([sp.sqrt(2) / sp.sqrt(3), 1 / sp.sqrt(3)]))
        # Ideal states are orthogonal to the isovector and kaon generators.
        for T in nonet_generators()[1:8]:
            self.assertEqual(sp.trace(T * T_N), 0)
            self.assertEqual(sp.trace(T * T_S), 0)

    def test_rotation(self):
        self.assertEqual(self.result['rotation_orthogonality'], sp.zeros(2))
        self.assertEqual(self.result['rotation_determinant'], 1)
        self.assertEqual(mixing_rotation(0), sp.eye(2))
        self.assertEqual(self.result['physical_gram'], sp.eye(2))
        self.assertEqual(self.result['isoscalar_rotation_invariance'], sp.zeros(3))

    def test_physical_generators_explicit(self):
        theta = sp.Symbol('theta', real=True)
        T_light, T_heavy = self.result['physical_generators']
        expected_light = sp.diag(sp.cos(theta) / sp.sqrt(2), sp.cos(theta) / sp.sqrt(2), -sp.sin(theta))
        expected_heavy = sp.diag(sp.sin(theta) / sp.sqrt(2), sp.sin(theta) / sp.sqrt(2), sp.cos(theta))
        self.assertEqual((T_light - expected_light).applyfunc(sp.simplify), sp.zeros(3))
        self.assertEqual((T_heavy - expected_heavy).applyfunc(sp.simplify), sp.zeros(3))

    def test_alda_dictionary(self):
        comparison = alda_exponent_comparison()
        self.assertEqual(comparison['lambda_gram'], 2 * sp.eye(9))
        self.assertEqual(comparison['u_exponent_difference'], sp.zeros(3))
        self.assertEqual(comparison['big_u_exponent_difference'], sp.zeros(3))


if __name__ == '__main__':
    unittest.main()
