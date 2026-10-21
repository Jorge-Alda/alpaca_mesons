"""Independent normalization and exact Gaussian-fit checks for G_V."""
import math
import json
from pathlib import Path
import unittest
import numpy as np
from fit_vector_gv import width_coefficient,fit_squared_coupling,flavor_factors,fixed_parameters


class VectorCouplingChecks(unittest.TestCase):
    def test_equal_mass_rho_normalization(self):
        mass, pion, decay_constant = .775,.140,.093
        expected = mass**3/(48*math.pi*decay_constant**4)*(1-4*pion**2/mass**2)**1.5
        self.assertAlmostEqual(width_coefficient(mass,mass,pion,pion,2),expected)

    def test_flavor_and_threshold(self):
        factors=flavor_factors(0.)
        self.assertAlmostEqual(factors['rho0_pi+_pi-'],2.)
        self.assertAlmostEqual(factors['kstar_charged_K0_pi+'],1.)
        self.assertAlmostEqual(factors['kstar_neutral_K+_pi-'],1.)
        self.assertAlmostEqual(factors['kstar_neutral_K0_pi0'],(-.98257+.01140)**2/2)
        self.assertAlmostEqual(factors['phi_K+_K-'],1.)
        self.assertAlmostEqual(factors['phi_K0_K0bar'],1.)
        self.assertEqual(width_coefficient(.6,.9,.5,.14,1),0.)

    def test_correlated_fit_stationarity(self):
        coefficients=np.array([2.,3.]); covariance=np.array([[.000004,.000001],[.000001,.000009]])
        result=fit_squared_coupling(coefficients,np.array([.02,.033]),covariance)
        self.assertAlmostEqual(float(coefficients@np.linalg.solve(covariance,
            result['residuals'])),0.)
        for coupling in result['G_V_delta_chi2_one_GeV']:
            residual=coefficients*coupling**2-np.array([.02,.033])
            self.assertAlmostEqual(float(residual@np.linalg.solve(covariance,residual))-result['chi2'],1.)

    def test_one_width_exact(self):
        result=fit_squared_coupling([10.],[.05],np.array([[.001**2]]))
        self.assertAlmostEqual(result['G_V_GeV'],math.sqrt(.005))
        self.assertAlmostEqual(result['chi2'],0.)
        self.assertEqual(result['nominal_dof'],0)
        selected,masses,angle=fixed_parameters()
        self.assertEqual(selected['label'],'central_single_trace_negative')
        self.assertAlmostEqual(masses['kstar_charged'],selected['even_parameters_by_name']['m_Kstar_GeV'])

    def test_four_channel_artifact(self):
        root=Path(__file__).resolve().parents[1]
        result=json.loads((root/'data/vector_sector/gv_conditional_fit.json').read_text())
        self.assertEqual(result['nominal_dof'],3)
        self.assertEqual(len(result['channels']),4)
        self.assertAlmostEqual(sum(result['chi2_by_block'].values()),result['chi2'])
        phi_covariance=np.array(result['measurements']['covariance'])[2:,2:]
        self.assertAlmostEqual(phi_covariance[0,1]/math.sqrt(phi_covariance[0,0]*phi_covariance[1,1]),-.72)
        coefficients=np.array(result['prediction_coefficients_per_GeV2'])
        residuals=np.array(result['residuals'])
        covariance=np.array(result['measurements']['covariance'])
        self.assertAlmostEqual(float(coefficients@np.linalg.solve(covariance,residuals)),0.,places=6)


if __name__=='__main__':
    unittest.main()
