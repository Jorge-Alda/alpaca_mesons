#!/usr/bin/env python3
"""Conditional nearly symmetric ALP TFF pieces, LO in QCD; GeV units. Run with --test for checks."""
import sys
import unittest
import sympy as sp
from derive_alp_heavy_quark_triangle import derive_axial_current_matrix_element

cG, cP, cA, lightP, lightA, kappa = sp.symbols('c_G c_hP c_hA c_lightP c_lightA kappa')
cgamma, Nc, charge, fa, fpi, qsq = sp.symbols('c_gamma N_c Q_h f_a f_pi q_squared', nonzero=True)
light_trace, heavy_trace = sp.symbols('Tr_Q2_cA_light Tr_Q2_MinvC_light')
alpha_s = sp.Symbol('alpha_s', nonnegative=True)
light_rotation = sp.Symbol('Tr_Q2_kappa_light')
b4 = sp.symbols('B_4')

# c_GP is invariant and includes the heavy pseudoscalar coupling.
full_cGP = cG + lightP + cP
effective_cG = sp.solve(sp.Eq(sp.Symbol('c_G_eff') + lightP, full_cGP), sp.Symbol('c_G_eff'))[0]
effective_cGA = sp.expand(effective_cG + lightA)
full_cGA = cG + lightA + cA
# Triangle and photon contact; light axial matrix element from the local OPE.
heavy = 2*Nc*charge**2/(4*sp.pi**2*fa)*(cP+(cA-cP)*b4)
contact = (cgamma + 2*Nc*light_trace)/(4*sp.pi**2*fa)
local_ope = -4*fpi**2*heavy_trace*(1-alpha_s/sp.pi)/(fa*qsq)
form_factor = sp.expand(contact + heavy + local_ope)
rotation = {cA:cA+kappa,cP:cP+kappa,cgamma:cgamma-2*Nc*charge**2*kappa,cG:cG-kappa}

def derive_hard_remainder():
    """Leading large-spacelike-momentum terms, for p² below the pair threshold."""
    momentum_fraction, boundary = sp.symbols("x delta", positive=True)
    virtuality, mass, alp_mass_squared, cutoff_squared = sp.symbols(
        "Q_squared m_h p_squared Lambda_squared", positive=True)
    denominator = sp.Symbol("D", positive=True)  # m_h² - x(1-x)p² below threshold
    # The boundary obeys denominator/Q_squared << delta << 1.  The first
    # integral is the soft endpoint, the second the complementary hard region.
    parameter = sp.Symbol("z", positive=True)
    soft = sp.integrate(1/(denominator + parameter*virtuality), (parameter,0,boundary))
    hard = sp.integrate(1/(parameter*virtuality), (parameter,boundary,1))
    leading_integral = (sp.log(virtuality/denominator))/virtuality
    matrix_element = (sp.log(cutoff_squared/denominator)-1)/virtuality
    return sp.simplify(leading_integral-matrix_element), sp.simplify(
        sp.limit(virtuality*(soft+hard-leading_integral), virtuality, sp.oo)), virtuality, cutoff_squared

def derive_msbar_subtraction():
    """Scalar two-propagator integral in d=4-2 epsilon, MS-bar subtraction."""
    epsilon = sp.Symbol("epsilon", positive=True)
    denominator, scale_squared, virtuality = sp.symbols("D mu_squared Q_squared", positive=True)
    integral = sp.exp(sp.EulerGamma*epsilon) * sp.gamma(epsilon) * (scale_squared/denominator)**epsilon
    finite = sp.simplify(sp.limit(integral-1/epsilon, epsilon, 0))
    hard = sp.expand_log(sp.log(virtuality/denominator)-finite, force=True)
    return finite, sp.simplify(hard), scale_squared, virtuality


def heavy_current_matching():
    """Two-loop heavy-current matching; inputs are the source's operator coefficient."""
    coupling, heavy_mass, scale = sp.symbols("alpha_s m_h mu", positive=True)
    light_currents = sp.symbols("J_u J_d J_s")
    heavy_charge = sp.Symbol("Q_h", real=True)
    coefficient = (coupling/sp.pi)**2*(sp.log(heavy_mass/scale)+sp.Rational(1, 8))
    charge_weighted = sum(weight*current for weight, current in
                          zip((sp.Rational(4,9),sp.Rational(1,9),sp.Rational(1,9)), light_currents))
    return coefficient, sp.expand(charge_weighted+heavy_charge**2*coefficient*sum(light_currents)), coupling, heavy_mass, scale


def active_flavor_evolution(number_of_flavors, singlet_evolution):
    """Evolve axial currents: non-singlets unchanged, full singlet multiplied."""
    identity = sp.eye(number_of_flavors)
    singlet_projector = sp.ones(number_of_flavors) / number_of_flavors
    return identity - singlet_projector + singlet_evolution * singlet_projector


class Checks(unittest.TestCase):
    def test_threshold(self):
        self.assertEqual(sp.simplify(effective_cGA-(full_cGA-cA+cP)),0)
        self.assertEqual(sp.simplify(effective_cG.subs(rotation, simultaneous=True)-effective_cG),0)
    def test_rotation(self):
        self.assertEqual(sp.simplify(form_factor.subs(rotation, simultaneous=True)-form_factor),0)
    def test_light_rotation(self):
        shifted = form_factor.subs({light_trace: light_trace + light_rotation,
                                   cgamma: cgamma - 2*Nc*light_rotation}, simultaneous=True)
        self.assertEqual(sp.simplify(shifted - form_factor), 0)

    def test_limits(self):
        self.assertEqual(sp.simplify(heavy.subs(b4,0)-2*Nc*charge**2*cP/(4*sp.pi**2*fa)),0)
        self.assertEqual(sp.simplify(heavy.subs(b4,1)-2*Nc*charge**2*cA/(4*sp.pi**2*fa)),0)
    def test_regions(self):
        remainder, boundary_check, virtuality, cutoff_squared = derive_hard_remainder()
        self.assertEqual(boundary_check, 0)
        self.assertEqual(sp.simplify(sp.expand_log(remainder-sp.log(virtuality/cutoff_squared)/virtuality-1/virtuality, force=True)), 0)

    def test_perturbative_overlap(self):
        result = derive_axial_current_matrix_element()
        self.assertEqual(result["remainder_check"], 0)
        self.assertNotEqual(result["remainder_expected"], 0)

    def test_msbar_subtraction(self):
        finite, hard, scale_squared, virtuality = derive_msbar_subtraction()
        self.assertEqual(sp.simplify(hard-sp.log(virtuality/scale_squared)), 0)
        self.assertEqual(sp.simplify(finite-sp.log(scale_squared/sp.Symbol("D", positive=True))), 0)

    def test_heavy_current_matching(self):
        coefficient, weighted, coupling, heavy_mass, scale = heavy_current_matching()
        self.assertEqual(sp.simplify(coefficient.subs(scale, heavy_mass)-(coupling/sp.pi)**2/8), 0)
        heavy_charge = sp.Symbol("Q_h", real=True)
        for current, charge_squared in zip(sp.symbols("J_u J_d J_s"),
                                           (sp.Rational(4,9), sp.Rational(1,9), sp.Rational(1,9))):
            self.assertEqual(sp.simplify(sp.diff(weighted, current) - charge_squared
                                         - heavy_charge**2*coefficient), 0)

    def test_active_flavor_evolution(self):
        evolution, later_evolution = sp.symbols("E E_later")
        for count in (3, 4, 5):
            matrix = active_flavor_evolution(count, evolution)
            self.assertEqual(matrix * sp.ones(count, 1), evolution * sp.ones(count, 1))
            non_singlet = sp.Matrix([1, -1] + [0] * (count-2))
            self.assertEqual(matrix * non_singlet, non_singlet)
            self.assertEqual(sp.simplify(active_flavor_evolution(count, later_evolution) * matrix
                                         - active_flavor_evolution(count, evolution * later_evolution)),
                             sp.zeros(count))
        for count in (3, 4, 5):
            weights = sp.Matrix([[sp.Rational(4, 9), sp.Rational(1, 9), sp.Rational(1, 9)] + [0]*(count-3)])
            actual = weights * active_flavor_evolution(count, evolution)
            expected = weights + sp.ones(1, count)*sp.Rational(2, 3)*(evolution-1)/count
            self.assertEqual(sp.simplify(actual-expected), sp.zeros(1, count))
        matrix = active_flavor_evolution(4, evolution)
        self.assertEqual(sp.simplify(matrix[0, 3] - (evolution-1)/4), 0)

    def test_tensor(self):
        p=sp.symbols('p0:4'); k=sp.symbols('k0:4'); metric=[1,-1,-1,-1]
        tensor=sp.Matrix(4,4,lambda mu,nu: sum(sp.LeviCivita(mu,nu,r,s)*metric[r]*k[r]*metric[s]*p[s] for r in range(4) for s in range(4)))
        q1=[k[i]+p[i]/2 for i in range(4)]; q2=[-k[i]+p[i]/2 for i in range(4)]
        self.assertTrue(all(sp.expand(x)==0 for x in (sp.Matrix([metric[i]*q1[i] for i in range(4)]).T*tensor)))
        self.assertTrue(all(sp.expand(x)==0 for x in (tensor*sp.Matrix([metric[i]*q2[i] for i in range(4)]))))
        self.assertEqual(tensor.subs({k[i]:-k[i] for i in range(4)}, simultaneous=True).T,tensor)

if __name__=='__main__':
    if '--test' in sys.argv: unittest.main(argv=[sys.argv[0]])
    else:
        remainder, boundary_check, _, _ = derive_hard_remainder()
        print('Hard remainder per Feynman parameter at finite p squared:\n' + sp.latex(remainder) + '\n')
        print('Cancellation of the boundary: ' + str(boundary_check))
        threshold, weighted, _, _, _ = heavy_current_matching()
        print('Heavy current threshold coefficient:\n' + sp.latex(threshold) + '\n')
        print('Charge-weighted matched current:\n' + sp.latex(weighted) + '\n')
        finite_msbar, hard_msbar, _, _ = derive_msbar_subtraction()
        print('MS-bar finite two-point integral:\n' + sp.latex(finite_msbar) + '\n')
        print('MS-bar hard remainder per parameter:\n' + sp.latex(hard_msbar) + '\n')
        print('Four-flavor axial-current evolution:\n' + sp.latex(active_flavor_evolution(4, sp.Symbol('E'))) + '\n')
        result = derive_axial_current_matrix_element()
        print('Direct triangle minus local axial-current term (cutoff scheme):\n' + sp.latex(result['remainder_expected']) + '\n')
        for label, value in [('full c_GP', full_cGP),('effective c_G',effective_cG),('effective c_GA',effective_cGA),('heavy triangle',heavy),('conditional LO ChPT sum',form_factor)]:
            print(label+':\n'+sp.latex(value)+'\n')
