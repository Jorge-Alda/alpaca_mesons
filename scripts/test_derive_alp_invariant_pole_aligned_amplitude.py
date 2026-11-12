"""Checks for the dressed-source expression with all pole admixtures."""
import unittest
import sympy as sp
from derive_alp_explicit_neutral_amplitude import aligned_sources, aligned_cubic_terms, cubic_amplitude, overlap_flavor_matrix
from derive_alp_invariant_pole_aligned_amplitude import derive_dressed_contact

class InvariantAmplitudeTest(unittest.TestCase):
    def test_invariant_source_and_inverse(self):
        B0, fpi, fa, m0sq, squared_mass = map(sp.Rational, [2,3,7,5,11])
        mass=sp.diag(1,2,4)
        axial=sp.Matrix([[1,sp.I,0],[-sp.I,2,1],[0,1,3]])
        pseudoscalar=sp.Matrix([[3,1,0],[1,0,sp.I],[0,-sp.I,2]])
        theta=sp.Rational(2,5)
        invariant=B0*(pseudoscalar-mass*axial-axial*mass)-m0sq*(theta+sp.trace(axial))*sp.eye(3)/3
        overlap=overlap_flavor_matrix(squared_mass, list((2*B0*mass).diagonal()),
            2*B0*sp.I*pseudoscalar, axial, -axial, theta, fpi, fa, m0sq)
        dressed_overlap=overlap-2*fpi*axial/fa
        meson_operator=B0*(mass*dressed_overlap+dressed_overlap*mass)+m0sq*sp.trace(dressed_overlap)*sp.eye(3)/3-squared_mass*dressed_overlap
        self.assertEqual(sp.simplify(meson_operator-2*fpi*invariant/fa),sp.zeros(3))
        shift=sp.Matrix([[0,1,0],[1,1,0],[0,0,0]])
        shifted_pseudoscalar=pseudoscalar+mass*shift+shift*mass
        shifted_theta=theta-sp.trace(shift)
        shifted_axial=axial+shift
        shifted=B0*(shifted_pseudoscalar-mass*shifted_axial-shifted_axial*mass)-m0sq*(shifted_theta+sp.trace(shifted_axial))*sp.eye(3)/3
        self.assertEqual(shifted,invariant)

    def test_dressed_sources_and_pole_alignment(self):
        zero=sp.zeros(3)
        first=sp.Matrix([[1,1,0],[1,0,sp.I],[0,-sp.I,2]])
        second=sp.Matrix([[0,sp.I,1],[-sp.I,2,0],[1,0,1]])
        parent=sp.diag(1,2,-1)
        scalar1=sp.Matrix([[2,1,0],[1,-1,sp.I],[0,-sp.I,3]])
        scalar3=sp.diag(2,1,-1)
        pseudoscalar2=sp.Matrix([[1,sp.I,0],[-sp.I,0,1],[0,1,2]])
        B0=sp.Rational(3,2); fpi=sp.Rational(2); fa=sp.Rational(5); m0sq=sp.Rational(7)
        theta1=sp.Rational(2,3); theta2=sp.Rational(-1,4)
        sa,s1,s2=map(sp.Rational,[3,2,1])
        gram=sp.Matrix([[sa,(s2-sa-s1)/2,(s1-sa-s2)/2],[(s2-sa-s1)/2,s1,(sa-s1-s2)/2],[(s1-sa-s2)/2,(sa-s1-s2)/2,s2]])
        chi=[2*B0*sp.diag(1,2,3),2*B0*scalar1,2*B0*sp.I*pseudoscalar2,2*B0*scalar3]
        current=[zero,zero]
        theta=[0,theta1,theta2]
        matrices=[2*fpi*item/fa for item in (parent,first,second)]
        direct=cubic_amplitude(matrices,[1,1,1],gram,chi,current,current,theta,fpi,fa,m0sq)
        pairs=first*second+second*first+parent*second+second*parent+parent*first+first*parent
        total=parent+first+second
        expected=fpi**2/fa**3*(-2*B0*sp.trace(scalar1*pairs)+2*B0*sp.trace(pseudoscalar2*total)-2*m0sq*theta2*sp.trace(total)/3+B0*sp.trace(scalar3)-2*m0sq*theta1*theta2)
        self.assertEqual(sp.simplify(direct-expected),0)
        symbolic=derive_dressed_contact()[1]
        replacements={"M_a_inverse_C":parent,"M_1_inverse_C":first,"M_2_inverse_C":second,
                      "hat_S_1":scalar1,"hat_S_3":scalar3,"hat_P_2":pseudoscalar2,
                      "f_pi":fpi,"f_a":fa,"B_0":B0,"m_0":sp.sqrt(m0sq),"N_f":sp.Integer(3),
                      "hat_Theta_1":theta1,"hat_Theta_2":theta2}
        symbolic=symbolic.xreplace({symbol:replacements[symbol.name] for symbol in symbolic.free_symbols})
        self.assertEqual(sp.simplify(symbolic.doit()-expected),0)
        mass,right,left,transformed_theta=aligned_sources(matrices[0],chi,current,current,theta,fpi,fa)
        aligned=sum(aligned_cubic_terms(matrices[1]-matrices[0],matrices[2]-matrices[0],1,1,gram,mass,right,left,transformed_theta,fpi,fa,m0sq))
        self.assertEqual(sp.simplify(aligned-expected),0)

if __name__ == '__main__':
    unittest.main()
