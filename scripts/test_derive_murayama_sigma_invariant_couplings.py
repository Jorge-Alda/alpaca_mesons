"""Checks for only-up radial source and invariant matching."""
import unittest
import sympy as sp
from scripts.derive_murayama_sigma_invariant_couplings import derive


class TestInvariantCouplings(unittest.TestCase):
    def test_matching_and_bracket(self):
        r = derive()
        self.assertEqual(r['invariant_pseudoscalar'], -2*r['quark_mass']*r['invariant_coupling'])
        self.assertEqual(r['radial_source'], 2*r['quark_mass']*r['invariant_coupling']/r['f_a'])
        self.assertEqual(sp.simplify(r['invariant_bracket'].subs(r['c_ap'], sp.Rational(1, 2))-r['model_bracket']), 0)

    def test_first_phase_invariant_does_not_determine_radial_source(self):
        r = derive()
        self.assertEqual(r['general_phase_source'], r['mass_source'].diff(r['x']).subs({r['x']: 0, r['sigma']: 0}))
        self.assertNotEqual(r['general_radial_source'].diff(r['radial_coupling']), sp.zeros(3))


if __name__ == '__main__':
    unittest.main()
