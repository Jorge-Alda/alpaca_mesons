#!/usr/bin/env python3
"""Weak ALP invariant for diagonal quark couplings, in natural units.

Assume an inert remaining weak current, three light flavours, w_0(2,3)=1,
and real diagonal entries of the Hermitian ALP quark coupling matrices.
The first-order coefficient w_1 remains a general complex 3 by 3 matrix.
Run: .venv/bin/python scripts/weak_partonic_diagonal_couplings.py
Test: .venv/bin/python scripts/weak_partonic_diagonal_couplings.py --test
"""
import sys
import unittest
import sympy as sp

weak_alp_matrix = sp.MatrixSymbol('w_1', 3, 3)
vector_up, vector_down, vector_strange = sp.symbols('c_u^V c_d^V c_s^V', real=True)
axial_up, axial_down, axial_strange = sp.symbols('c_u^A c_d^A c_s^A', real=True)


def expressions():
    vector_coupling = sp.diag(vector_up, vector_down, vector_strange)
    axial_coupling = sp.diag(axial_up, axial_down, axial_strange)
    weak_matrix = sp.zeros(3)
    weak_matrix[1, 2] = 1
    commutator = (vector_coupling - axial_coupling) * weak_matrix - weak_matrix * (vector_coupling - axial_coupling)
    invariant = weak_alp_matrix.as_explicit() - sp.I * commutator
    conjugate_invariant = invariant.adjoint()
    return vector_coupling, axial_coupling, weak_matrix, sp.simplify(commutator), sp.simplify(invariant), sp.simplify(conjugate_invariant)


class DiagonalCouplingTests(unittest.TestCase):
    def test_commutator_and_entries(self):
        _, _, _, commutator, invariant, _ = expressions()
        expected = sp.zeros(3)
        expected[1, 2] = vector_down - axial_down - vector_strange + axial_strange
        self.assertEqual(commutator, expected)
        self.assertEqual(sp.simplify(invariant - weak_alp_matrix.as_explicit() + sp.I * expected), sp.zeros(3))

    def test_conjugate_sign(self):
        _, _, _, _, _, conjugate = expressions()
        self.assertEqual(sp.simplify(conjugate[2, 1] - sp.conjugate(weak_alp_matrix[1, 2]) - sp.I * (vector_down - axial_down - vector_strange + axial_strange)), 0)

    def test_universal_left_couplings(self):
        _, _, _, commutator, _, _ = expressions()
        universal = {vector_strange: vector_down, axial_strange: axial_down}
        self.assertEqual(commutator.subs(universal), sp.zeros(3))


def main():
    print('Natural units; first order in a/f_a; diagonal Hermitian c_q^V and c_q^A; w_1 unrestricted.')
    for label, expression in zip(('Vector coupling', 'Axial coupling', 'Weak matrix w_0', 'Commutator', 'Weak invariant', 'Adjoint invariant'), expressions()):
        print(f'{label}:\n{sp.latex(expression)}\n')


if __name__ == '__main__':
    if '--test' in sys.argv:
        sys.argv.remove('--test')
        unittest.main()
    else:
        main()
