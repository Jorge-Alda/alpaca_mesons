"""Meaningful normalization and source-formula checks for the KKpi predictor.

Run: .venv/bin/python -m unittest scripts.test_kkpi_vector_current
"""

import unittest
import numpy as np

from scripts.implementation.kkpi_vector_current import (
    GEV_MINUS_TWO_TO_NB, KKPiInputs, babar_table7_prediction,
    electromagnetic_cross_section_matrix, energy_bin_matrices,
    predict_from_matrices, resonance_scalar_basis,
)


class KKPiCurrentTests(unittest.TestCase):
    def setUp(self):
        self.inputs = KKPiInputs(.49, .14, .092, .15, .07, .77,
            lambda x: 1/(.6-np.asarray(x)-.11j),
            lambda x: 1/(.61-np.asarray(x)-.01j),
            lambda x: 1/(.8-np.asarray(x)-.045j))

    def test_massless_constant_scalar_absolute_normalization(self):
        # Integral ds dt |epsilon p1 p2 p3|² = (Q²)^5/480.
        alpha, energy = .01, 2.
        result = electromagnetic_cross_section_matrix(energy, (0., 0., 0.),
            lambda q, s, t: np.ones(t.shape+(1,)), alpha, 8, 8)
        expected = GEV_MINUS_TWO_TO_NB*alpha**2*energy**4/(23040*np.pi)
        np.testing.assert_allclose(result, [[expected]], rtol=1e-13)

    def test_source_charged_polynomials_and_relative_signs(self):
        q, s, t = 2.5, .7, .9
        c = np.arange(1, 8)*.001
        d = np.arange(1, 5)*.002
        g = np.arange(1, 6)*.003
        kaon, pion = self.inputs.kaon_mass**2, self.inputs.pion_mass**2
        def direct(pair, m1, m2, m3):
            return ((c[0]-c[1]+c[4])*q-(c[0]-c[1]-c[4]+2*c[5])*pair
                    +(c[0]+c[1]+8*c[2]-c[4])*m3+8*c[3]*(m1-m2))
        def double(pair, mass):
            return d[2]*(q+pair)+(d[0]+8*d[1]-d[2])*mass
        local = ((g[0]+2*g[1]-g[2])*(s+t)-2*g[1]*(q+kaon)
                 -(g[0]-g[2])*(3*kaon+pion)+2*g[3]*(kaon+pion)+2*g[4]*kaon)
        inputs = self.inputs
        expected = -4*inputs.vector_pion_coupling/(inputs.vector_mass_normalization
            *inputs.pion_decay_constant**3)*(direct(s,kaon,kaon,pion)*inputs.omega_propagator(s)
            +direct(t,kaon,pion,kaon)*inputs.kstar_propagator(t)
            -2*inputs.vector_source_coupling/inputs.vector_pion_coupling*local*inputs.rho_propagator(q)
            -np.sqrt(2)*inputs.vector_source_coupling*inputs.vector_mass_normalization
            *inputs.rho_propagator(q)*(double(s,pion)*inputs.omega_propagator(s)
                                    +double(t,kaon)*inputs.kstar_propagator(t)))
        actual = resonance_scalar_basis(q,s,t,inputs,"charged")@np.r_[c,d,g]
        np.testing.assert_allclose(actual, expected, rtol=1e-13)

    def test_kstar_selection_removes_all_local_vppp(self):
        for mode in ("charged", "neutral"):
            selected = resonance_scalar_basis(2.5,.7,.9,self.inputs,mode,"kstar")
            np.testing.assert_array_equal(selected[11:], np.zeros(5))
        with self.assertRaisesRegex(NotImplementedError, "Dalitz-component"):
            babar_table7_prediction()

    def test_coherent_interference_and_bin_average(self):
        def scalar(q, s, t):
            return np.broadcast_to(np.array([1., 1j, 2.]),t.shape+(3,))
        matrices = energy_bin_matrices([[1.,2.]], (0.,0.,0.),scalar,.01,4,8,8)
        # F=1+2i-1=2i for coefficient vector (1,2,-.5).
        actual = predict_from_matrices([1.,2.,-.5],matrices)[0]
        # Uniform integral E^4 dE from 1 to 2 is 31/5.
        expected = 4*GEV_MINUS_TWO_TO_NB*.01**2/(23040*np.pi)*31/5
        self.assertAlmostEqual(actual/expected,1., places=13)
        self.assertGreaterEqual(np.linalg.eigvalsh(matrices[0]).min(),-1e-14)


if __name__ == "__main__":
    unittest.main()
