#!/usr/bin/env python3
"""Tree-level propagator mass term; .venv/bin/python scripts/derive_alp_finite_mass_hard_tff.py --test. Natural units."""
import sys
import unittest
import sympy as sp

u, Q1, Q2, m_a, a2, t = sp.symbols('u Q_1_squared Q_2_squared m_a a_2 t', positive=True)
phi = 6*u*(1-u)*(1+a2*sp.gegenbauer(2, sp.Rational(3,2), 2*u-1))
denominator = u*Q1+(1-u)*Q2+u*(1-u)*m_a**2
expansion = sp.series(1/denominator, m_a, 0, 4).removeO().expand()
equal = sp.integrate(phi*expansion.subs(Q2,Q1), (u,0,1)).simplify()
single = sp.integrate(sp.cancel(phi*expansion.subs(Q2,0)), (u,0,1)).simplify()

class Checks(unittest.TestCase):
    def test_massless(self):
        self.assertEqual(sp.simplify(expansion.subs(m_a,0)-1/(u*Q1+(1-u)*Q2)),0)
    def test_equal(self):
        self.assertEqual(sp.simplify(equal-(1/Q1-m_a**2*(7-3*a2)/(35*Q1**2))),0)
    def test_single(self):
        self.assertEqual(sp.simplify(single-sp.integrate(sp.cancel(phi/(u*Q1)-phi*u*(1-u)*m_a**2/(u*Q1)**2),(u,0,1))),0)
    def test_scaling(self):
        self.assertEqual(sp.limit(t*equal.subs(Q1,t),t,sp.oo),1)
    def test_exchange(self):
        self.assertEqual(sp.simplify(expansion-expansion.xreplace({Q1:Q2,Q2:Q1,u:1-u})),0)

if __name__=='__main__':
    if '--test' in sys.argv:
        sys.argv.remove('--test'); unittest.main()
    else:
        for label, value in [('distribution amplitude',phi),('propagator denominator',denominator),('mass expansion',expansion),('equal spacelike virtualities',equal),('one real photon',single)]:
            print(label+':\n'+sp.latex(value)+'\n')
