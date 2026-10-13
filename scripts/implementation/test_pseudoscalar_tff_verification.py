"""Independent normalization, symmetry, and observation checks for the TFF."""
import math
import unittest
import numpy as np
from .odd_parity_rcht import OddParityCouplings, VJPCouplings, VVPCouplings
from .pseudoscalar_tff import TFFInputs, tff_components, pseudoscalar_tff
from .odd_parity_data import load_cleo_tff_observations, load_babar_pi0_tff_observation, load_cleo_tff, load_babar_pi0_table11


class TFFVerificationTests(unittest.TestCase):
    def setUp(self):
        self.x = TFFInputs(.092, .145, .770, .135, math.sqrt(2)/3,
                           2*.135**2*math.sqrt(2)/3)
        self.c = OddParityCouplings(VJPCouplings(.01,-.02,.03,99,.04,-.05,88),
                                   VVPCouplings(.06,-.01,.02,77))

    def test_anomaly_and_published_both_real_corrections(self):
        x, c, d = self.x, self.c.vjp, self.c.vvp
        wzw, vjp, vvp = tff_components(self.c, x)
        self.assertAlmostEqual(wzw, 1/(4*math.pi**2*x.F))
        self.assertAlmostEqual(vjp, 8*math.sqrt(2)*x.FV*x.mp**2/(3*x.F*x.MV**3)
                               *(c.c1+c.c2+8*c.c3-c.c5))
        self.assertAlmostEqual(vvp, -8*x.FV**2*x.mp**2/(3*x.F*x.MV**4)
                               *(d.d1+8*d.d2-d.d3))

    def test_bose_symmetry_and_ope_homogeneous_null_directions(self):
        np.testing.assert_allclose(tff_components(self.c,self.x,-2,-5),
                                   tff_components(self.c,self.x,-5,-2), atol=1e-14)
        # Six independent null shifts of the five stated OPE equations.
        vectors = ([1,1,-.25,0,0,0,0,0,0,0,0],
                   [0,1,0,0,1,1,0,0,0,0,0],
                   [0,0,0,1,0,0,0,0,0,0,0],
                   [0,0,0,0,0,0,1,0,0,0,0],
                   [0,0,0,0,0,0,0,-8,1,0,0],
                   [0,0,0,0,0,0,0,0,0,0,1])
        baseline = pseudoscalar_tff(self.c,self.x,[0,1,5,20])
        for vector in vectors:
            shifted = self.c.as_array()+np.array(vector)*.2
            couplings = OddParityCouplings(VJPCouplings(*shifted[:7]),VVPCouplings(*shifted[7:]))
            np.testing.assert_allclose(pseudoscalar_tff(couplings,self.x,[0,1,5,20]),
                                       baseline, rtol=1e-13, atol=1e-13)

    def test_observation_units_and_separate_normalization(self):
        cleo = load_cleo_tff()
        observations = load_cleo_tff_observations({'pi0':self.x})
        self.assertEqual(len(observations),1)
        rows = cleo.rows[cleo.rows['meson']=='pi0']
        np.testing.assert_array_equal(observations[0].values, rows['q2F_GeV'])
        np.testing.assert_allclose(observations[0].covariance,
                                   np.diag(rows['stat_GeV']**2+rows['syst_GeV']**2))
        babar = load_babar_pi0_tff_observation(self.x)
        rows = load_babar_pi0_table11().rows
        np.testing.assert_allclose(babar.values,rows['q2F_MeV']*.001)
        self.assertAlmostEqual(math.sqrt(babar.covariance[-2,-2]),.0365)
        self.assertEqual(len(babar.nuisances),1)
        self.assertEqual(babar.nuisances[0].fraction,.023)
        self.assertEqual(babar.unit,'GeV')

    def test_invalid_virtualities_and_pole(self):
        for q in (-1,float('nan'),float('inf'),1+2j):
            with self.assertRaises((ValueError,TypeError)):
                pseudoscalar_tff(self.c,self.x,q)
        with self.assertRaises(ValueError):
            tff_components(self.c,self.x,self.x.MV**2,0)
