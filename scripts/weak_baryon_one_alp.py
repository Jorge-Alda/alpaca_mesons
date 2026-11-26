#!/usr/bin/env python3
"""Expand the weak meson-baryon Lagrangian to first order in a/f_a.

Run: .venv/bin/python scripts/weak_baryon_one_alp.py
Test: .venv/bin/python scripts/weak_baryon_one_alp.py --test
All matrices are abstract three-flavour matrices, and the weak spurion is
assumed ALP-independent in the starting parametrization. The left-handed
rotation, its Hermitian conjugate, and u are expanded to first order.
Natural units; a/f_a, u and flavour matrices are dimensionless.
"""
import sys
import unittest

import sympy as sp


def show(label, expression):
    print(f"{label}:\n{sp.latex(expression)}\n")


def calculate():
    n = 3
    identity = sp.Identity(n)
    a_over_fa = sp.Symbol('a/f_a', real=True)
    d_w, f_w = sp.symbols('d_w f_w', real=True)
    w = sp.MatrixSymbol('w', n, n)
    w_dagger = sp.MatrixSymbol('w^dagger', n, n)
    lambda_b1_dagger = sp.MatrixSymbol('lambda_B1^dagger', n, n)
    lambda_b2 = sp.MatrixSymbol('lambda_B2', n, n)
    delta = sp.MatrixSymbol('delta', n, n)
    kappa = sp.MatrixSymbol('kappa', n, n)
    cq_a = sp.MatrixSymbol('c_q^A', n, n)
    cq_s = sp.MatrixSymbol('c_q^S', n, n)
    mass_inverse = sp.MatrixSymbol('M_a^-1(C)', n, n)

    left = identity + sp.I * a_over_fa * (delta - kappa)
    left_dagger = identity - sp.I * a_over_fa * (delta - kappa)
    # u = 1 once the no-mixing variables are chosen and meson fields vanish.
    dressed = left * (w + w_dagger) * left_dagger
    dressed_linear = sp.diff(dressed, a_over_fa).subs(a_over_fa, 0).doit().expand()
    dressed_no_mixing = dressed_linear.subs({delta: -cq_s, kappa: -cq_a - mass_inverse}).doit().expand()
    lagrangian_coefficient = (
        d_w * sp.Trace(lambda_b1_dagger * (dressed_no_mixing * lambda_b2 + lambda_b2 * dressed_no_mixing))
        + f_w * sp.Trace(lambda_b1_dagger * (dressed_no_mixing * lambda_b2 - lambda_b2 * dressed_no_mixing))
    ).expand().doit()
    return locals()


class TestWeakBaryonOneAlp(unittest.TestCase):
    def test_left_rotation_and_no_mixing(self):
        e = calculate()
        matrix = e['cq_a'] - e['cq_s'] + e['mass_inverse']
        weak_matrix = e['w'] + e['w_dagger']
        expected = sp.I * (matrix * weak_matrix - weak_matrix * matrix)
        self.assertEqual((e['dressed_no_mixing'] - expected).doit().expand(), sp.ZeroMatrix(3, 3))

    def test_vanishing_for_identity_rotation(self):
        e = calculate()
        actual = e['dressed_no_mixing'].subs({e['cq_s']: e['cq_a'] + e['mass_inverse']}).doit().expand()
        self.assertEqual(actual, sp.ZeroMatrix(3, 3))

    def test_generic_baryon_projection(self):
        e = calculate()
        self.assertTrue(e['lagrangian_coefficient'].has(e['lambda_b1_dagger'], e['lambda_b2'], e['d_w'], e['f_w']))


def main():
    e = calculate()
    print('Natural units; first order in a/f_a; no meson fields; weak spurion initially ALP-independent.')
    show('Left-handed rotation, first order', e['left'])
    show('Weak matrix after the left-handed rotation', e['dressed'])
    show('Coefficient of a/f_a before choosing field variables', e['dressed_linear'])
    show('Coefficient of a/f_a in the no-mixing parametrization', e['dressed_no_mixing'])
    show('Coefficient of (a/f_a) bar(B_1) B_2', e['lagrangian_coefficient'])


if __name__ == '__main__':
    if '--test' in sys.argv:
        sys.argv.remove('--test')
        unittest.main()
    else:
        main()
