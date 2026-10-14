"""Unit checks of physical tensor orientation and source conversion."""
import unittest
import sympy as sp
from scripts.derive_levi_civita_convention import derive
from scripts.symbolic_utils import levi_civita_contract


class ConventionTests(unittest.TestCase):
    def test_metric_lowering(self):
        result = derive()
        self.assertEqual(result['epsilon_lower0123'], -1)
        self.assertEqual(result['canonical_basis_contraction'], -1)
        self.assertEqual(result['legacy_basis_contraction'], 1)
        self.assertEqual(result['canonical_Fdual_over_E'], -4)
        self.assertEqual(result['canonical_epsilon_dAdA_over_E'], -2)

    def test_action_and_interference(self):
        result = derive()
        self.assertEqual(result['action_invariance'], 0)
        self.assertEqual(result['total_squared_invariance'], 0)
        A, Lc, Ld, c, d = sp.symbols('A L_c L_d c d', real=True)
        self.assertEqual(result['canonical_scalar'], -A+Lc*c+Ld*d)
        self.assertNotEqual(sp.expand((-A+Lc*c+Ld*d)**2-(A+Lc*c+Ld*d)**2),0)

    def test_ward_and_antisymmetry(self):
        q, r, e, f = (sp.Matrix(x) for x in ((2,0,0,2),(2,0,0,-2),(0,1,0,0),(0,0,1,0)))
        self.assertEqual(levi_civita_contract(q,q,r,f),0)
        self.assertEqual(levi_civita_contract(e,f,q,r),-levi_civita_contract(q,e,r,f))


if __name__ == '__main__':
    unittest.main()
