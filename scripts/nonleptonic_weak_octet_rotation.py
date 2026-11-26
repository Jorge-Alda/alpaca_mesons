#!/usr/bin/env python3
"""First-order rotation and covariant flavor contraction of a weak operator.

Run: .venv/bin/python scripts/nonleptonic_weak_octet_rotation.py
Test: .venv/bin/python scripts/nonleptonic_weak_octet_rotation.py --test
Natural units; x=a/f_a; current matrices in the (u,d,s) basis.
"""
import sys
import unittest
import sympy as sp


def calculate():
    x = sp.Symbol('a/f_a', real=True)
    rotation = sp.Matrix(3, 3, lambda row, col: sp.Symbol(f't_{row+1}{col+1}'))
    rotation_matrix = sp.eye(3) + sp.I*x*rotation
    inverse_rotation_matrix = sp.eye(3) - sp.I*x*rotation
    first_current = sp.zeros(3)
    first_current[0, 2] = 1
    second_current = sp.zeros(3)
    second_current[1, 0] = 1
    first_rotated = rotation_matrix * first_current * inverse_rotation_matrix
    second_rotated = rotation_matrix * second_current * inverse_rotation_matrix
    # SymPy's Kronecker product represents both current flavor indices without
    # introducing a separate symbolic tensor algebra.
    four_quark_tensor = sp.kronecker_product(first_rotated, second_rotated)
    four_quark_linear = four_quark_tensor.diff(x).subs(x, 0).applyfunc(sp.expand)
    flavor_contraction = second_rotated * first_rotated
    octet_at_zero = flavor_contraction.subs(x, 0)
    octet_linear = flavor_contraction.diff(x).subs(x, 0).applyfunc(sp.expand)
    diagonal_rotation = {rotation[i, j]: (sp.Symbol(f't_{i+1}', real=True) if i == j else 0)
                         for i in range(3) for j in range(3)}
    return (rotation, first_current, second_current, four_quark_linear, octet_at_zero,
            octet_linear, diagonal_rotation)


def print_calculation():
    (rotation, first_current, second_current, four_quark_linear, octet_at_zero,
     octet_linear, diagonal_rotation) = calculate()
    print('Natural units; x=a/f_a; linear order in x; Wilson coefficients held fixed.')
    for label, expression in [('First current matrix', first_current),
                              ('Second current matrix', second_current),
                              ('Coefficient of x in the four-quark tensor', four_quark_linear),
                              ('Contracted flavor matrix at x=0', octet_at_zero),
                              ('Coefficient of x in contracted flavor matrix', octet_linear),
                              ('Diagonal rotation coefficient', octet_linear.subs(diagonal_rotation))]:
        print(f'{label}:\n{sp.latex(expression)}\n')
    print('The contraction is an equivariant octet map, not the complete weak Hamiltonian or a unique matching normalization.')


class RotationTests(unittest.TestCase):
    def test_generic_contraction(self):
        rotation, _, _, _, octet, linear, _ = calculate()
        self.assertEqual(linear, (sp.I*(rotation*octet-octet*rotation)).applyfunc(sp.expand))

    def test_diagonal_contraction(self):
        _, _, _, _, octet, linear, diagonal = calculate()
        self.assertEqual(linear.subs(diagonal),
                         (sp.I*(sp.Symbol('t_2', real=True)-sp.Symbol('t_3', real=True))*octet).applyfunc(sp.expand))

    def test_tensor_both_currents(self):
        rotation, first, second, linear, _, _, _ = calculate()
        expected = sp.I*(sp.kronecker_product(rotation*first-first*rotation, second)
                         + sp.kronecker_product(first, rotation*second-second*rotation))
        self.assertEqual(linear.applyfunc(sp.expand), expected.applyfunc(sp.expand))


if __name__ == '__main__':
    if '--test' in sys.argv:
        sys.argv.remove('--test')
        unittest.main()
    else:
        print_calculation()
