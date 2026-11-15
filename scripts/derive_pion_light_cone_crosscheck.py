#!/usr/bin/env python3
"""Pion-source check of the ALP light-cone derivation, natural units; run --test."""
import sys
import unittest

import sympy as sp
import derive_alp_light_cone_tff as alp

F_pi = sp.symbols('F_pi', positive=True)
quark_charges = sp.diag(sp.Rational(2, 3), -sp.Rational(1, 3), -sp.Rational(1, 3))
pion_axial_currents = sp.diag(F_pi, -F_pi, 0)
charge_weighted_current = sp.trace(quark_charges**2 * pion_axial_currents)
prefactor = 2 * charge_weighted_current
equal_result = sp.simplify(prefactor * alp.equal_integral)
single_result = sp.simplify(prefactor * alp.singly_integral)
asymptotic_single_result = single_result.subs(alp.a_2, 0)


class TestPionCrosscheck(unittest.TestCase):
    def test_flavor_trace(self):
        self.assertEqual(charge_weighted_current, F_pi / 3)

    def test_symmetric(self):
        self.assertEqual(equal_result, 2 * F_pi / (3 * alp.Q_1_squared))

    def test_single_asymptotic(self):
        self.assertEqual(asymptotic_single_result, 2 * F_pi / alp.Q_1_squared)

    def test_single_shape_dependence(self):
        self.assertEqual(single_result, 2 * F_pi * (1 + alp.a_2) / alp.Q_1_squared)


if __name__ == '__main__':
    if '--test' in sys.argv:
        unittest.main(argv=[sys.argv[0]])
    else:
        for label, expression in [('pion axial currents', pion_axial_currents),
                                  ('charge-weighted axial current', charge_weighted_current),
                                  ('symmetric limit', equal_result),
                                  ('singly virtual limit', single_result),
                                  ('asymptotic shape', asymptotic_single_result)]:
            print(f'{label}:\n{sp.latex(expression)}\n')
