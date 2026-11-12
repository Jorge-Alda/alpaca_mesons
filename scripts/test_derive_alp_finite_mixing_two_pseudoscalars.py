import unittest
import sympy as sp
from derive_alp_finite_mixing_two_pseudoscalars import derive, projected_vertex
from derive_alp_charged_u3_vector_redefinition import derive as charged_check


class FiniteMixingTwoPseudoscalars(unittest.TestCase):
    def test_bose_and_all_alp_components(self):
        tensor, parent, first, second, amplitude, weak = derive()
        self.assertEqual(sp.expand(amplitude-projected_vertex(tensor,parent,second,first)),0)
        self.assertNotEqual(sp.diff(amplitude, parent[2], first[2], second[2]),0)
        self.assertNotEqual(sp.diff(amplitude, parent[0], first[2], second[2]),0)
        self.assertNotEqual(sp.diff(amplitude, parent[2], first[2], second[0]),0)

    def test_leading_weak_mixing(self):
        tensor,parent,first,second,amplitude,weak = derive()
        expected = sp.Symbol("epsilon", real=True)*parent[2]*sum(tensor[2,i,j]*first[i]*second[j]
                                  for i in range(2) for j in range(2))
        self.assertEqual(sp.expand(weak-expected),0)

    def test_charged_vector_redefinition(self):
        contact, amplitude, on_shell_contact = charged_check()
        self.assertNotEqual(contact,0)
        self.assertNotEqual(amplitude,0)
        self.assertEqual(on_shell_contact,0)



class NeutralKaonWeakMixing(unittest.TestCase):
    def test_same_lo_mass_operator_and_complex_conjugate_source(self):
        from derive_alp_flavour_violating_two_pseudoscalars import (
            extracted_channel, masses, B_0, m_d, m_s, cV_ds, cS_ds, f_a)
        # Isospin-symmetric LO masses: m_K^2-m_pi^2 = B_0(m_s-m_d).
        kaon = extracted_channel('K0','pi0')
        antikaon = extracted_channel('K0bar','pi0')
        self.assertEqual(sp.simplify(antikaon + sp.conjugate(kaon)), 0)
        lo = sp.expand(kaon).subs(masses['K0']**2,
                                  masses['pi0']**2+B_0*(m_s-m_d))
        expected = B_0*(m_s-m_d)*(cV_ds-cS_ds)/(sp.sqrt(2)*f_a)
        self.assertEqual(sp.simplify(lo-expected), 0)

class SourceConsistentWeakMixing(unittest.TestCase):
    def test_ccwz_quadratic_and_vertex_share_sources(self):
        from derive_alp_lo_joint_vertices import quadratic_blocks
        from derive_alp_pole_projected_amplitudes import ccwz_vertex
        epsilon = sp.Symbol('epsilon', positive=True)
        zero = sp.zeros(3)
        generator = sp.diag(1, 0, -1)
        sources = ([sp.diag(2,4,6), sp.diag(1+sp.I,-1,0), zero, zero],
                   [zero]*4, [zero]*4, [0]*4, sp.S.One, 1/epsilon, 0, sp.S(12))
        kinetic, mass = quadratic_blocks([generator], *sources)
        self.assertEqual(kinetic, sp.eye(2))
        parent_shift, daughter_shift = sp.symbols('parent_shift daughter_shift')
        parent = sp.Matrix([epsilon*parent_shift,1])
        daughter = sp.Matrix([1,epsilon*daughter_shift])
        parent_equation = sp.expand(((12*kinetic-mass)*parent)[0]).coeff(epsilon)
        daughter_equation = sp.expand(((4*kinetic-mass)*daughter)[1]).coeff(epsilon)
        shifts = sp.solve([parent_equation,daughter_equation], [parent_shift,daughter_shift])
        self.assertTrue(shifts[parent_shift] != 0 and shifts[daughter_shift] != 0)
        parent, daughter = parent.subs(shifts), daughter.subs(shifts)
        gram = sp.Matrix([[12,-6,-6],[-6,4,2],[-6,2,4]])
        momenta = [sp.eye(3)[:,i] for i in range(3)]
        exact = ccwz_vertex([generator],sources,[parent,daughter,daughter],momenta,gram)
        bare = ccwz_vertex([generator],sources,[sp.Matrix([0,1]),sp.Matrix([1,0]),sp.Matrix([1,0])],momenta,gram)
        self.assertNotEqual(bare,0)
        self.assertEqual(sp.series(exact-bare,epsilon,0,2).removeO(),0)

    def test_same_source_finite_poles_approach_weak_limit(self):
        from derive_alp_finite_mixing_two_pseudoscalars import projected_ccwz_vertex
        from derive_alp_lo_joint_vertices import quadratic_blocks
        zero = sp.zeros(3)
        generator = sp.diag(1,0,-1)
        residuals = []
        for epsilon in (sp.Rational(1,5), sp.Rational(1,10)):
            sources = ([sp.diag(2,4,6),sp.diag(1+sp.I,-1,0),zero,zero],
                       [zero]*4,[zero]*4,[0]*4,sp.S.One,1/epsilon,0,sp.S(12))
            kinetic,mass = quadratic_blocks([generator],*sources)
            self.assertEqual(mass,sp.Matrix([[4,-epsilon/2],[-epsilon/2,12]]))
            splitting = sp.sqrt(64+epsilon**2)
            daughter_mass = (16-splitting)/2
            parent_mass = (16+splitting)/2
            overlap = -epsilon/(8+splitting)
            norm = sp.sqrt(1+overlap**2)
            parent = sp.Matrix([overlap,1])/norm
            daughter = sp.Matrix([1,-overlap])/norm
            gram = sp.Matrix([[parent_mass,-parent_mass/2,-parent_mass/2],
                              [-parent_mass/2,daughter_mass,(parent_mass-2*daughter_mass)/2],
                              [-parent_mass/2,(parent_mass-2*daughter_mass)/2,daughter_mass]])
            exact = projected_ccwz_vertex([generator],sources,parent,daughter,daughter,gram)
            bare = projected_ccwz_vertex([generator],sources,sp.Matrix([0,1]),
                                         sp.Matrix([1,0]),sp.Matrix([1,0]),gram)
            residuals.append(abs(float(sp.N(exact-bare))))
        self.assertGreater(residuals[0], 0)
        self.assertAlmostEqual(residuals[0]/residuals[1],8,delta=0.02)

class ChargedComplexVectorSource(unittest.TestCase):
    def test_ccwz_charged_amplitude_with_complex_ds_source(self):
        from derive_alp_finite_mixing_two_pseudoscalars import projected_ccwz_vertex
        zero=sp.zeros(3)
        real,imag,overlap=sp.symbols('c_ds_real c_ds_imag v_a', real=True)
        vector=sp.Matrix([[0,0,0],[0,0,real+sp.I*imag],[0,real-sp.I*imag,0]])
        matrix_unit=lambda i,j:sp.eye(3)[:,i]*sp.eye(3)[j,:]
        basis=[sp.sqrt(2)*matrix_unit(2,0),sp.sqrt(2)*matrix_unit(0,1)]
        sources=([zero]*4,[vector,zero,zero,zero],[vector,zero,zero,zero],
                 [0]*4,sp.S.One,sp.S.One,0,0)
        gram=sp.Matrix([[9,-5,-4],[-5,3,2],[-4,2,2]])
        amplitude=sp.I*projected_ccwz_vertex(basis,sources,sp.Matrix([0,0,overlap]),
                                            sp.Matrix([1,0,0]),sp.Matrix([0,1,0]),gram)
        # m_K^2-m_pi^2=3-2=1=B_0(m_s-m_d), f_a=1; scalar source is zero.
        self.assertEqual(sp.expand(amplitude+overlap*(real+sp.I*imag)),0)
