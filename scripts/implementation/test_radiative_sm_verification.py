"""Independent regression checks of the action-derived radiative amplitude.

Run .venv/bin/python -m unittest scripts.implementation.test_radiative_sm_verification
"""
import unittest
from dataclasses import replace

import numpy as np

from .odd_parity_rcht import OddParityCouplings, VJPCouplings, VVPCouplings, ideal_vector_generators
from .radiative_odd_parity import OmegaPiGammaInputs, G_omega_pi_gamma
from .radiative_sm import RadiativeSMInputs, Q, G_radiative_sm, radiative_sm_flavor, radiative_sm_width


class RadiativeSMVerification(unittest.TestCase):
    def setUp(self):
        gen = ideal_vector_generators()
        self.x = RadiativeSMInputs(.3028, .092, .145, .770, .783, .135,
                                  gen['omega'], gen['rho0'], np.diag([.135**2]*3))
        self.g = OddParityCouplings(VJPCouplings(.01, -.02, .03, .04, .05, -.06, .07),
                                   VVPCouplings(.08, -.09, .10, .11))

    def test_published_omega_eq27_every_coupling(self):
        old = OmegaPiGammaInputs(self.x.e, self.x.F, self.x.FV, self.x.MV,
                                self.x.Mvector, self.x.mp)
        for j in range(11):
            values = np.eye(11)[j]
            g = OddParityCouplings(VJPCouplings(*values[:7]), VVPCouplings(*values[7:]))
            with self.subTest(coupling=j):
                self.assertAlmostEqual(G_radiative_sm(g, self.x).real,
                                       G_omega_pi_gamma(g, old), places=13)
                self.assertEqual(G_radiative_sm(g, self.x).imag, 0)

    def test_complete_neutral_exchange_for_all_flavor_basis_pairs(self):
        basis = [np.eye(3)[:, i:i+1] @ np.eye(3)[j:j+1] for i in range(3) for j in range(3)]
        chi = np.diag([.02, .03, .48])
        for v in basis:
            for p in basis:
                fq, fh, _, exchange = radiative_sm_flavor(replace(self.x, TV=v, TP=p, chi=chi))
                self.assertAlmostEqual(sum(w*a for w,a,_,_ in exchange), fq)
                self.assertAlmostEqual(sum(w*b for w,_,b,_ in exchange), fh)

    def test_c4_charged_kaon_and_isospin_pion_selection(self):
        v = np.zeros((3,3)); v[2,0] = 1
        p = v.T
        x = replace(self.x, TV=v, TP=p, chi=np.diag([.02,.02,.48]))
        _,_,fk,_ = radiative_sm_flavor(x)
        self.assertAlmostEqual(fk, -.46)
        pion = replace(x, TV=np.roll(v, -1, axis=0), TP=np.roll(v, -1, axis=0).T)
        self.assertEqual(radiative_sm_flavor(pion)[2], 0)

    def test_crossed_width_spin_factor_and_closed_channel(self):
        x = replace(self.x, mp=1.2)
        amp = abs(G_radiative_sm(self.g,x))**2
        expected = amp*(x.mp**2-x.Mvector**2)**3/(32*np.pi*x.mp**3)
        self.assertAlmostEqual(radiative_sm_width(self.g,x,'P_to_V'),expected)
        with self.assertRaises(ValueError):
            radiative_sm_width(self.g,self.x,'P_to_V')
        with self.assertRaises(ValueError):
            radiative_sm_width(self.g,x,'V_to_P')


if __name__ == '__main__':
    unittest.main()
