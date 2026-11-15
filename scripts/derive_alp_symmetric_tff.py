#!/usr/bin/env python3
"""Conditional nearly symmetric ALP TFF pieces, LO in QCD; GeV units. Run with --test for checks."""
import sys
import unittest
import sympy as sp

cG, cP, cA, lightP, lightA, kappa = sp.symbols('c_G c_hP c_hA c_lightP c_lightA kappa')
cgamma, Nc, charge, fa, fpi, qsq = sp.symbols('c_gamma N_c Q_h f_a f_pi q_squared', nonzero=True)
light_trace, heavy_trace = sp.symbols('Tr_Q2_cA_light Tr_Q2_MinvC_light')
b4 = sp.symbols('B_4')

# c_GP is invariant and includes the heavy pseudoscalar coupling.
full_cGP = cG + lightP + cP
effective_cG = sp.solve(sp.Eq(sp.Symbol('c_G_eff') + lightP, full_cGP), sp.Symbol('c_G_eff'))[0]
effective_cGA = sp.expand(effective_cG + lightA)
full_cGA = cG + lightA + cA
# Triangle and photon contact; light axial matrix element from the local OPE.
heavy = 2*Nc*charge**2/(4*sp.pi**2*fa)*(cP+(cA-cP)*b4)
contact = (cgamma + 2*Nc*light_trace)/(4*sp.pi**2*fa)
local_ope = -4*fpi**2*heavy_trace/(fa*qsq)
form_factor = sp.expand(contact + heavy + local_ope)
rotation = {cA:cA+kappa,cP:cP+kappa,cgamma:cgamma-2*Nc*charge**2*kappa,cG:cG-kappa}

class Checks(unittest.TestCase):
    def test_threshold(self):
        self.assertEqual(sp.simplify(effective_cGA-(full_cGA-cA+cP)),0)
        self.assertEqual(sp.simplify(effective_cG.subs(rotation, simultaneous=True)-effective_cG),0)
    def test_rotation(self):
        self.assertEqual(sp.simplify(form_factor.subs(rotation, simultaneous=True)-form_factor),0)
    def test_limits(self):
        self.assertEqual(sp.simplify(heavy.subs(b4,0)-2*Nc*charge**2*cP/(4*sp.pi**2*fa)),0)
        self.assertEqual(sp.simplify(heavy.subs(b4,1)-2*Nc*charge**2*cA/(4*sp.pi**2*fa)),0)
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
        for label, value in [('full c_GP', full_cGP),('effective c_G',effective_cG),('effective c_GA',effective_cGA),('heavy triangle',heavy),('conditional sum; overlap not established',form_factor)]:
            print(label+':\n'+sp.latex(value)+'\n')
