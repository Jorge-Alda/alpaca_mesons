#!/usr/bin/env python3
"""Verify mixed finite tadpole contractions; run .venv/bin/python scripts/verify_mixed_finite_contractions.py.

Natural units, Q=sK-M^2+Sigma; one loop and first order in 1/f_a.
The symmetric coincident propagator is arbitrary, including nonet mixing.
"""
import unittest
import sympy as sp
from derive_nlo_loop_vertices import (leading_lagrangian, coefficient_of_orders,
    F_PI, F_A, ALP_FIELD, ALP_DERIVATIVE)
from derive_nlo_quadratic_chipt import trace

external, first, second, chi, axial, pseudoscalar = sp.symbols(
    'external first second chi_0 c_A c_P', commutative=False)
background, loop_first, loop_second = sp.symbols('background loop_first loop_second')


def derive_comparison():
    field = background*external + loop_first*first + loop_second*second
    fields = {'Phi':field, 'partial_mu_Phi':background*external, 'chi_0':chi,
              'c_q^A':axial, 'c_q^P':pseudoscalar,
              'c_q^V':sp.Integer(0), 'c_q^S':sp.Integer(0)}
    vertex = coefficient_of_orders(leading_lagrangian(fields, trace), 3, 1)
    coefficient = sp.expand(vertex).coeff(background).coeff(loop_first).coeff(loop_second)
    axial_words = first*second*external + external*first*second - 2*second*external*first
    axial_swapped = second*first*external + external*second*first - 2*first*external*second
    words = external*first*second + first*external*second + first*second*external
    swapped = external*second*first + second*external*first + second*first*external
    expected_axial = trace(axial*(axial_words+axial_swapped))/(6*F_PI*F_A)
    expected_pseudo = -trace(pseudoscalar*(chi*(words+swapped)+(words+swapped)*chi))/(12*F_PI*F_A)
    return (sp.expand(coefficient).coeff(ALP_DERIVATIVE), expected_axial,
            sp.expand(coefficient).coeff(ALP_FIELD), expected_pseudo)

class FiniteContractionTests(unittest.TestCase):
    def test_all_symmetric_contractions(self):
        derived_axial, expected_axial, derived_pseudo, expected_pseudo = derive_comparison()
        self.assertEqual(sp.expand(derived_axial-expected_axial), 0)
        self.assertEqual(sp.expand(derived_pseudo-expected_pseudo), 0)

if __name__ == '__main__':
    values = derive_comparison()
    for label, value in zip(('Axial from vertex', 'Axial from contraction',
                             'Pseudoscalar from vertex', 'Pseudoscalar from contraction'), values):
        print(label + ': ' + sp.latex(value))
    unittest.main(argv=[__file__], exit=True)
