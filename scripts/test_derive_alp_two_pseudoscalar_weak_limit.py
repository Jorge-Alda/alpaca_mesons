import unittest
import sympy as sp
from derive_alp_two_pseudoscalar_weak_limit import derive, derive_projector_and_pole_shift, derive_kaon_pion_comparison
from derive_alp_two_pseudoscalar_weak_limit import eom_identity_left, eom_identity_right


class WeakMixingLimit(unittest.TestCase):
    def test_general_mass_operator_identity(self):
        self.assertEqual(sp.expand(eom_identity_left-eom_identity_right), sp.zeros(3))

    def test_previous_kaon_pion_result(self):
        _, actual, earlier = derive_kaon_pion_comparison()
        m_d, m_s, B_0, f_a, real_c, imag_c = sp.symbols('m_d m_s B_0 f_a c_real c_imag', real=True)
        expected = B_0*(m_s-m_d)*(real_c+sp.I*imag_c)/(sp.sqrt(2)*f_a)
        self.assertEqual(sp.simplify(actual-earlier), 0)
        self.assertEqual(sp.simplify(actual-expected), 0)

    def test_projector_residue_and_meson_pole_shift(self):
        projector, source, residue, shift = derive_projector_and_pole_shift()
        self.assertEqual(sp.simplify(projector*projector-projector), sp.zeros(2))
        self.assertEqual(sp.simplify(projector*source-residue), sp.zeros(2, 1))
        m1, ma, c1 = sp.symbols('m_1_squared m_a_squared c_1', real=True)
        self.assertEqual(sp.simplify(shift-2*c1**2/(m1-ma)), 0)

    def test_leading_coefficient_and_remainder(self):
        _, amplitude, leading = derive()
        epsilon, f_pi, B_0 = sp.symbols('epsilon f_pi B_0', nonzero=True, real=True)
        scalar = sp.MatrixSymbol('hat_S_1', 3, 3)
        meson_i = sp.MatrixSymbol('Phi_i', 3, 3)
        meson_j = sp.MatrixSymbol('Phi_j', 3, 3)
        expected = -B_0*epsilon*(sp.Trace(scalar*meson_i*meson_j)+sp.Trace(scalar*meson_j*meson_i))/(2*f_pi)
        self.assertEqual(sp.expand(leading-expected), 0)
        self.assertEqual(sp.expand(amplitude).coeff(epsilon, 2), 0)
        self.assertNotEqual(sp.expand(amplitude).coeff(epsilon, 3), 0)


if __name__ == '__main__':
    unittest.main()
