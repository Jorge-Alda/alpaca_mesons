"""Checks of the LO radial vertices and on-shell field-coordinate equivalence."""
import unittest
import sympy as sp
from scripts.derive_murayama_sigma_pp import calculate, polar_check


class SigmaVerticesTest(unittest.TestCase):
    def test_vertices(self):
        (B_0,m_u,f_pi,f_phi,lam), _, neutral, charged, kaon = calculate()
        direction = sp.Matrix([1, 1/sp.sqrt(3), sp.sqrt(sp.Rational(2,3)), 0])
        angular = sp.Matrix([0,0,0,1])
        expected = B_0*m_u/f_phi*direction*direction.T+2*lam*f_phi*angular*angular.T
        self.assertEqual(sp.simplify(neutral-expected), sp.zeros(4))
        self.assertEqual(charged, B_0*m_u/f_phi)
        self.assertEqual(kaon, charged)

    def test_polar_cartesian(self):
        self.assertEqual(polar_check()[2], 0)


if __name__ == '__main__':
    unittest.main()
