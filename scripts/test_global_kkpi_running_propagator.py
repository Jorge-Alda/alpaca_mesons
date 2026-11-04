"""Regression checks for the global fit's running propagators.

Run: .venv/bin/python -m unittest scripts.test_global_kkpi_running_propagator
The fit uses GeV, GeV squared, and nb.
"""
import json
import unittest
from unittest.mock import patch

import numpy as np

from scripts import fit_babar_kkpi_global_vector_sector as fit


def original_running_propagator(name,invariant,beta,mass,shared_width,reference):
    """Evaluate the running propagator as it was before reuse of the ground term."""
    denominator=1/fit.reference_propagator(name,invariant,reference,0.)
    on_shell=1/fit.reference_propagator(name,mass**2,reference,0.)
    width=-denominator.imag/mass*shared_width/(-on_shell.imag/mass)
    excited=(fit.reference_propagator(name,invariant,reference,beta)
        -fit.reference_propagator(name,invariant,reference,0.))
    return 1/(mass**2-invariant-1j*mass*width)+excited


class GlobalRunningPropagatorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with (fit.ROOT/'data/vector_sector/babar_kkpi_global_vector_sector_fit.json').open() as source:
            cls.stored=json.load(source)
        cls.candidate=np.array(cls.stored['parameter_values'])

    def test_reference_term_is_reused_only_within_one_call(self):
        model=fit.GlobalLikelihood(1,3,4)
        point,_,mass_kstar,_,totals=model.unpack(self.candidate)
        reference=dict(model.preliminary,m_rho_GeV=point[0],m_Kstar_GeV=mass_kstar,
            Gamma_total_omega_GeV=totals['omega'],Gamma_total_phi_GeV=totals['phi'])
        invariant=np.array([.25,.8,1.4,2.1])
        for name,mass,width in (
                ('rho',point[0],totals['rho0']),
                ('kstar',mass_kstar,(totals['Kstar0']+totals['Kstarplus'])/2)):
            with self.subTest(name=name):
                expected=original_running_propagator(name,invariant,(.3,-.1),mass,width,reference)
                with patch.object(fit,'reference_propagator',wraps=fit.reference_propagator) as propagator:
                    actual=fit.running_propagator(name,invariant,(.3,-.1),mass,width,reference)
                self.assertEqual(propagator.call_count,3)
                np.testing.assert_array_equal(actual,expected)
                changed=fit.running_propagator(name,invariant+.01,(.2,.1),mass,width*1.01,reference)
                expected_changed=original_running_propagator(name,invariant+.01,(.2,.1),mass,width*1.01,reference)
                np.testing.assert_array_equal(changed,expected_changed)

    def test_residuals_predictions_and_jacobian_match_original_path(self):
        for orders,block_size in (((2,16,16),None),((1,3,4),7)):
            with self.subTest(orders=orders,block_size=block_size):
                original=fit.GlobalLikelihood(*orders,bin_block_size=block_size)
                improved=fit.GlobalLikelihood(*orders,bin_block_size=block_size)
                points=[self.candidate.copy()]
                for index in range(len(self.candidate)):
                    shifted=self.candidate.copy()
                    shifted[index]+=1e-4
                    points.append(shifted)
                points.append(self.candidate.copy())
                with patch.object(fit,'running_propagator',side_effect=original_running_propagator):
                    expected=[]
                    for values in points:
                        residual=original.residual(values)
                        _,predictions=original.babar_blocks(values,original.unpack(values))
                        expected.append((residual,{name:array.copy() for name,array in predictions.items()}))
                    expected_jacobian=original.jacobian(self.candidate)
                for values,(residual,predictions) in zip(points,expected):
                    np.testing.assert_array_equal(improved.residual(values),residual)
                    _,actual_predictions=improved.babar_blocks(values,improved.unpack(values))
                    for channel in fit.CHANNELS:
                        np.testing.assert_array_equal(actual_predictions[channel],predictions[channel])
                np.testing.assert_array_equal(improved.jacobian(self.candidate),expected_jacobian)
                if orders==(2,16,16):
                    np.testing.assert_allclose(expected[0][0],self.stored['residuals'],rtol=1e-12,atol=1e-11)
                    for channel in fit.CHANNELS:
                        np.testing.assert_allclose(expected[0][1][channel],
                            self.stored['babar_predictions'][channel],rtol=1e-12,atol=1e-11)


if __name__=='__main__':
    unittest.main()
