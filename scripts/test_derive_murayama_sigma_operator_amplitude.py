"""Tests for the pole operator and radial on-shell vertex."""
import unittest
import sympy as sp
from scripts.derive_murayama_sigma_operator_amplitude import derive, derive_dimensionful_coefficients


class OperatorAmplitudeTests(unittest.TestCase):
    def test_dimensionful_source_coefficients(self):
        _, coefficients = derive_dimensionful_coefficients()
        symbols = {str(symbol): symbol for matrix in coefficients.values() for symbol in matrix.free_symbols}
        mu, md, ms, fa = [symbols[name] for name in ('m_u', 'm_d', 'm_s', 'f_a')]
        self.assertEqual(coefficients['S00'], sp.diag(mu, md, ms))
        self.assertEqual(coefficients['S10'], sp.diag(mu/fa, 0, 0))
        self.assertEqual(coefficients['P01'], sp.diag(-mu/fa, 0, 0))
        self.assertEqual(coefficients['P11'], sp.diag(-mu/fa**2, 0, 0))
        self.assertEqual(coefficients['S02'], sp.diag(-mu/(2*fa**2), 0, 0))
        self.assertEqual(coefficients['S12'], sp.diag(-mu/(2*fa**3), 0, 0))

    def test_dimensionful_contraction(self):
        _, _, projected, norm_first, norm_second = derive()
        symbols = {str(symbol): symbol for symbol in projected.free_symbols | norm_first.free_symbols | norm_second.free_symbols}
        B0, mu, fpi, fa = [symbols[name] for name in ('B_0', 'm_u', 'f_pi', 'f_a')]
        mi2, mj2, msigma2 = [symbols[name] for name in ('m_i_squared', 'm_j_squared', 'm_sigma_squared')]
        first = sp.diag(*(symbols['A_i_'+q]/fa for q in ('u', 'd', 's')))
        second = sp.diag(*(symbols['A_j_'+q]/fa for q in ('u', 'd', 's')))
        scalar10 = sp.diag(mu/fa, 0, 0)
        pseudoscalar11 = sp.diag(-mu/fa**2, 0, 0)
        scalar12 = sp.diag(-mu/(2*fa**3), 0, 0)
        ordinary = (2*B0**3*fpi**2*sp.trace(scalar10*(first*second+second*first))
                    - 2*B0**2*fpi**2*sp.trace(pseudoscalar11*(first+second))
                    - 2*B0*fpi**2*sp.trace(scalar12)
                    + (msigma2-mi2-mj2)/fa)
        self.assertEqual(sp.simplify(projected-ordinary), 0)

    def test_normalization_and_vertex(self):
        first, second, projected, norm_first, norm_second = derive()
        symbols = {str(symbol): symbol for symbol in projected.free_symbols | norm_first.free_symbols | norm_second.free_symbols}
        B0, mu, fpi, fa = [symbols[name] for name in ('B_0', 'm_u', 'f_pi', 'f_a')]
        first_up, second_up = [symbols[name] for name in ('A_i_u', 'A_j_u')]
        mi2, mj2, msigma2 = [symbols[name] for name in ('m_i_squared', 'm_j_squared', 'm_sigma_squared')]
        expected = B0*mu*fpi**2/fa**3*(1+2*B0*first_up)*(1+2*B0*second_up)+(msigma2-mi2-mj2)/fa
        self.assertEqual(sp.simplify(projected-expected), 0)
        self.assertEqual(sp.simplify(norm_first-(1+2*B0**2*fpi**2/fa**2*sum(symbols['A_i_'+q]**2 for q in ('u','d','s')))), 0)
        self.assertEqual(sp.simplify(norm_second-(1+2*B0**2*fpi**2/fa**2*sum(symbols['A_j_'+q]**2 for q in ('u','d','s')))), 0)

    def test_operator_and_pole_equations(self):
        B0, fpi, fa, s, m1 = sp.symbols('B_0 f_pi f_a s m_1', real=True, nonzero=True)
        masses = sp.diag(1, 2, 3)
        pseudoscalar_source = sp.diag(-1, 0, 0)
        # Exact rational point keeps the verification short and deterministic.
        substitutions = {B0: 1, fpi: 2, fa: 3, m1: 4, s: 5}
        operator = sp.diag(*(2*B0*masses[i,i]-s for i in range(3)))+sp.ones(3)*m1**2/3
        inverse_source = operator.inv()*sp.Matrix(list(pseudoscalar_source.diagonal()))
        inverse_source = inverse_source.subs(substitutions)
        operator = operator.subs(substitutions)
        self.assertEqual(operator*inverse_source, sp.Matrix(list(pseudoscalar_source.diagonal())))
        meson_block = operator
        mixed = -2*sp.sqrt(2)/3*sp.Matrix([1, 0, 0])
        meson_component = -sp.sqrt(2)*sp.Rational(2,3)*inverse_source
        self.assertEqual(sp.simplify(meson_block*meson_component+mixed), sp.zeros(3,1))

if __name__ == '__main__':
    unittest.main()
