#!/usr/bin/env python3
"""Derive a left-handed weak coefficient invariant at first order in a/f_a.

Run: .venv/bin/python scripts/weak_partonic_invariant.py
Test: .venv/bin/python scripts/weak_partonic_invariant.py --test
Natural units. The weak current outside the displayed bilinear is inert.
The weak matrices need not be Hermitian; the rotation generator is Hermitian.
"""
import sys
import unittest
import sympy as sp


a_over_f = sp.Symbol('a/f_a', real=True)
weak_matrix = sp.Symbol('w', commutative=False)
weak_alp_matrix = sp.Symbol('w^(1)', commutative=False)
left_rotation = sp.Symbol('delta-kappa', commutative=False)
vector_coupling = sp.Symbol('c_q^V', commutative=False)
axial_coupling = sp.Symbol('c_q^A', commutative=False)
scalar_coupling = sp.Symbol('c_q^S', commutative=False)
meson_mixing = sp.Symbol('M_a^-1(C)', commutative=False)


def expressions():
    weak_source = weak_matrix + a_over_f * weak_alp_matrix
    rotated_source = (1 + sp.I * a_over_f * left_rotation) * weak_source * (1 - sp.I * a_over_f * left_rotation)
    rotated_linear = sp.diff(sp.expand(rotated_source), a_over_f).subs(a_over_f, 0).expand()
    left_coupling = vector_coupling - axial_coupling
    rotated_left_coupling = left_coupling + left_rotation
    invariant = weak_alp_matrix - sp.I * (left_coupling * weak_matrix - weak_matrix * left_coupling)
    rotated_invariant = sp.expand(rotated_linear - sp.I * (rotated_left_coupling * weak_matrix - weak_matrix * rotated_left_coupling))
    no_mixing_linear = rotated_linear.subs(left_rotation, -scalar_coupling + axial_coupling + meson_mixing).expand()
    return weak_source, rotated_linear, invariant.expand(), rotated_invariant, no_mixing_linear


class InvarianceTests(unittest.TestCase):
    def test_invariant_under_rotation(self):
        _, _, invariant, rotated_invariant, _ = expressions()
        self.assertEqual(sp.expand(rotated_invariant - invariant), 0)

    def test_initially_alp_independent(self):
        _, _, _, _, no_mixing_linear = expressions()
        generator = axial_coupling - scalar_coupling + meson_mixing
        expected = sp.I * (generator * weak_matrix - weak_matrix * generator)
        self.assertEqual(sp.expand(no_mixing_linear.subs(weak_alp_matrix, 0) - expected), 0)


def main():
    print('Natural units; order a/f_a; weak current outside bilinear is inert.')
    for label, value in zip(('Initial source', 'Rotated linear coefficient', 'Invariant coefficient', 'Rotated invariant', 'Linear coefficient in no-mixing variables'), expressions()):
        print(f'{label}:\n{sp.latex(value)}\n')


if __name__ == '__main__':
    if '--test' in sys.argv:
        sys.argv.remove('--test')
        unittest.main()
    else:
        main()
