#!/usr/bin/env python3
"""Leading-twist ALP light-current contribution; natural units. Run --test."""
import sys
import unittest
import sympy as sp

u = sp.symbols('u', real=True)
a_2 = sp.symbols('a_2', real=True)
Q_1_squared, Q_2_squared = sp.symbols('Q_1_squared Q_2_squared', positive=True)
f_pi, f_a = sp.symbols('f_pi f_a', positive=True)
current_matrix_element = sp.symbols('current_matrix_element', real=True)
phi = 6*u*(1-u)*(1+a_2*sp.gegenbauer(2, sp.Rational(3, 2), 2*u-1))
denominator = u*Q_1_squared+(1-u)*Q_2_squared
normalization = sp.integrate(phi, (u, 0, 1))
equal_integral = sp.integrate(phi / Q_1_squared, (u, 0, 1))
singly_integral = sp.integrate(sp.cancel(phi/(u*Q_1_squared)), (u, 0, 1))
coefficient = 4*f_pi**2*current_matrix_element/f_a

class TestLightCone(unittest.TestCase):
    def test_normalization(self): self.assertEqual(normalization, 1)
    def test_bose(self): self.assertEqual(sp.simplify(phi-phi.subs(u, 1-u)), 0)
    def test_equal(self): self.assertEqual(sp.simplify(coefficient*equal_integral-coefficient/Q_1_squared), 0)
    def test_single(self): self.assertEqual(sp.simplify(singly_integral-3*(1+a_2)/Q_1_squared), 0)
    def test_asymptotic(self): self.assertEqual(sp.simplify(phi.subs(a_2, 0)-6*u*(1-u)), 0)
    def test_ward(self):
        epsilon = sp.LeviCivita
        for nu in range(4):
            self.assertEqual(sum(epsilon(mu,nu,rho,sigma)*sp.Symbol('q1_'+str(mu))*sp.Symbol('q1_'+str(rho))*sp.Symbol('q2_'+str(sigma)) for mu in range(4) for rho in range(4) for sigma in range(4)), 0)

if __name__ == '__main__':
    if '--test' in sys.argv:
        unittest.main(argv=[sys.argv[0]])
    else:
        for label, value in [('distribution amplitude', phi), ('normalization', normalization), ('hard denominator', denominator), ('equal virtuality integral', equal_integral), ('singly virtual integral', singly_integral), ('source coefficient', coefficient)]:
            print(f'{label}:\n{sp.latex(value)}\n')
