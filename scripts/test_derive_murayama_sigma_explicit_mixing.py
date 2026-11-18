"""Checks of the explicit pole columns and sigma vertex."""
import unittest
import sympy as sp
from scripts.derive_murayama_sigma_explicit_mixing import derivation


class ExplicitMixingTest(unittest.TestCase):
    def test_pole_column(self):
        d = derivation()
        residual = (d['mass_matrix'] - d['s']*sp.eye(4))*d['unnormalized']
        self.assertEqual([sp.simplify(x) for x in residual[:3]], [0, 0, 0])
        self.assertEqual(sp.simplify(residual[3] - d['pole_polynomial']), 0)

    def test_explicit_physical_poles(self):
        d = derivation()
        substitutions = {d['B0']: sp.Rational(5,2), d['mu']: sp.Rational(2,1000),
                         d['md']: sp.Rational(5,1000), d['ms']: sp.Rational(95,1000),
                         d['m1']: sp.Rational(9,10), d['fpi']: sp.Rational(92,1000),
                         d['fphi']: sp.Rational(1,5)}
        mass_matrix = d['mass_matrix'].subs(substitutions)
        for root in sp.nroots(d['pole_polynomial'].subs(substitutions)):
            column = d['unnormalized'].subs(substitutions).subs(d['s'], root).evalf(15)
            residual = (mass_matrix - root*sp.eye(4))*column
            self.assertLess(max(abs(complex(x)) for x in residual), 1e-12)


if __name__ == '__main__':
    unittest.main()
