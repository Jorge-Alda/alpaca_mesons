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


def original_reference_propagator_pair(name,invariant,beta,reference):
    """Evaluate both reference propagators separately as in the original fit."""
    return (fit.reference_propagator(name,invariant,reference,0.),
        fit.reference_propagator(name,invariant,reference,beta))


def original_jacobian(model,values,finite_difference_steps=None):
    """Use the full residual for every finite-difference column."""
    baseline=model.residual(values)
    jacobian=np.zeros((75,32))
    state=model.unpack(values)
    _,odd_derivatives,_=model.old_blocks(values,state)
    jacobian[:12,12:24]=odd_derivatives
    jacobian[28:73,12:24]=model.babar_blocks(values,state,derivatives=True)
    for index in list(range(12))+list(range(24,32)):
        step=(np.sqrt(np.finfo(float).eps)*max(1.,abs(values[index]))
            if finite_difference_steps is None else finite_difference_steps[index])
        shifted=values.copy();shifted[index]+=step
        jacobian[:,index]=(model.residual(shifted)-baseline)/step
    return jacobian


def width_optimized_jacobian(model,values):
    """Use the previous width-pull optimization for a local timing comparison."""
    baseline=model.residual(values)
    jacobian=np.zeros((75,32))
    state=model.unpack(values)
    _,odd_derivatives,_=model.old_blocks(values,state)
    jacobian[:12,12:24]=odd_derivatives
    jacobian[28:73,12:24]=model.babar_blocks(values,state,derivatives=True)
    for index in list(range(12))+list(range(24,32)):
        step=np.sqrt(np.finfo(float).eps)*max(1.,abs(values[index]))
        shifted=values.copy();shifted[index]+=step
        if index in (8,11):
            shifted_old,_,_=model.old_blocks(shifted,model.unpack(shifted))
            jacobian[:28,index]=(shifted_old-baseline[:28])/step
            jacobian[28:,index]=(baseline[28:]-baseline[28:])/step
        else:
            jacobian[:,index]=(model.residual(shifted)-baseline)/step
    return jacobian


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
                self.assertEqual(propagator.call_count,2)
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
                with (patch.object(fit,'running_propagator',side_effect=original_running_propagator),
                      patch.object(fit,'reference_propagator_pair',side_effect=original_reference_propagator_pair)):
                    expected=[]
                    for values in points:
                        residual=original.residual(values)
                        _,predictions=original.babar_blocks(values,original.unpack(values))
                        expected.append((residual,{name:array.copy() for name,array in predictions.items()}))
                    expected_jacobian=original_jacobian(original,self.candidate)
                for values,(residual,predictions) in zip(points,expected):
                    actual_residual=improved.residual(values)
                    np.testing.assert_array_equal(actual_residual.view(np.uint64),residual.view(np.uint64))
                    _,actual_predictions=improved.babar_blocks(values,improved.unpack(values))
                    for channel in fit.CHANNELS:
                        np.testing.assert_array_equal(actual_predictions[channel].view(np.uint64),
                            predictions[channel].view(np.uint64))
                np.testing.assert_array_equal(improved.jacobian(self.candidate).view(np.uint64),
                    expected_jacobian.view(np.uint64))
                if orders==(2,16,16):
                    np.testing.assert_allclose(expected[0][0],self.stored['residuals'],rtol=1e-12,atol=1e-11)
                    for channel in fit.CHANNELS:
                        np.testing.assert_allclose(expected[0][1][channel],
                            self.stored['babar_predictions'][channel],rtol=1e-12,atol=1e-11)

    def test_combined_reference_propagators_match_separate_calls(self):
        model=fit.GlobalLikelihood(1,3,4)
        point,_,mass_kstar,_,totals=model.unpack(self.candidate)
        reference=dict(model.preliminary,m_rho_GeV=point[0],m_Kstar_GeV=mass_kstar,
            Gamma_total_omega_GeV=totals['omega'],Gamma_total_phi_GeV=totals['phi'])
        for name in ('rho','kstar','omega','phi'):
            for invariant in (.8,np.array([.2,.8,1.4,2.1])):
                for beta in (0.,(.3,-.1)):
                    with self.subTest(name=name,invariant=invariant,beta=beta):
                        expected=original_reference_propagator_pair(name,invariant,beta,reference)
                        actual=fit.reference_propagator_pair(name,invariant,beta,reference)
                        for result,reference_result in zip(actual,expected):
                            np.testing.assert_array_equal(result,reference_result)

    def test_jacobian_skips_full_residual_when_only_one_block_changes(self):
        model=fit.GlobalLikelihood(1,3,4)
        with patch.object(model,'residual',wraps=model.residual) as residual:
            model.jacobian(self.candidate)
        self.assertEqual(residual.call_count,11)

    def test_single_channel_predictions_match_full_integration(self):
        for block_size in (None,7):
            with self.subTest(block_size=block_size):
                full=fit.GlobalLikelihood(1,3,4,bin_block_size=block_size)
                selected=fit.GlobalLikelihood(1,3,4,bin_block_size=block_size)
                values=self.candidate.copy()
                full_residual,full_predictions=full.babar_blocks(values,full.unpack(values))
                start=0
                for channel in fit.CHANNELS:
                    count=len(full.babar[channel]['observed'])
                    residual,predictions=selected.babar_blocks(values,selected.unpack(values),channels=(channel,))
                    np.testing.assert_array_equal(residual.view(np.uint64),full_residual[start:start+count].view(np.uint64))
                    np.testing.assert_array_equal(predictions[channel].view(np.uint64),full_predictions[channel].view(np.uint64))
                    start+=count

    def test_shifted_candidate_jacobians_match_full_residual_path(self):
        original=fit.GlobalLikelihood(1,3,4)
        improved=fit.GlobalLikelihood(1,3,4)
        for index in range(len(self.candidate)):
            with self.subTest(index=index):
                values=self.candidate.copy();values[index]+=1e-4
                expected=original_jacobian(original,values)
                actual=improved.jacobian(values)
                np.testing.assert_array_equal(actual.view(np.uint64),expected.view(np.uint64))

    def test_default_grid_channel_coordinate_jacobians_match_full_path(self):
        original=fit.GlobalLikelihood()
        improved=fit.GlobalLikelihood()
        for index in range(24,32):
            with self.subTest(index=index):
                values=self.candidate.copy();values[index]+=1e-4
                expected=original_jacobian(original,values)
                actual=improved.jacobian(values)
                np.testing.assert_array_equal(actual.view(np.uint64),expected.view(np.uint64))

    def test_beta_shift_outside_domain_uses_full_residual(self):
        original=fit.GlobalLikelihood(1,3,4)
        improved=fit.GlobalLikelihood(1,3,4)
        values=self.candidate.copy();values[26]=5.
        steps=np.sqrt(np.finfo(float).eps)*np.maximum(1.,np.abs(values))
        steps[26]=1e-4
        expected=original_jacobian(original,values,steps)
        actual=improved.jacobian(values,steps)
        np.testing.assert_array_equal(actual.view(np.uint64),expected.view(np.uint64))

    def test_invalid_and_negative_steps_match_full_residual_path(self):
        original=fit.GlobalLikelihood(1,3,4)
        improved=fit.GlobalLikelihood(1,3,4)
        steps=np.sqrt(np.finfo(float).eps)*np.maximum(1.,np.abs(self.candidate))
        steps[8]=-1e6
        steps[11]=-1e-5
        expected=original_jacobian(original,self.candidate,steps)
        actual=improved.jacobian(self.candidate,steps)
        np.testing.assert_array_equal(actual,expected)
        np.testing.assert_array_equal(np.signbit(actual[28:,11]),np.signbit(expected[28:,11]))

    def test_penalty_baseline_uses_full_residual_path(self):
        def evaluate(model,calculator):
            actual_residual=model.residual
            actual_residual(self.candidate)
            first_call=[True]
            def residual_with_penalty(values):
                if first_call[0]:
                    first_call[0]=False
                    return np.full(75,1e8+np.linalg.norm(values))
                return actual_residual(values)
            with patch.object(model,'residual',side_effect=residual_with_penalty) as residual:
                result=calculator(model,self.candidate)
            self.assertEqual(residual.call_count,21)
            return result
        original=fit.GlobalLikelihood(1,3,4)
        improved=fit.GlobalLikelihood(1,3,4)
        expected=evaluate(original,original_jacobian)
        actual=evaluate(improved,lambda model,values:model.jacobian(values))
        np.testing.assert_array_equal(actual,expected)

    def test_nonfinite_baseline_uses_full_residual_path(self):
        model=fit.GlobalLikelihood(1,3,4)
        actual_residual=model.residual
        actual_residual(self.candidate)
        first_call=[True]
        def residual_with_nonfinite(values):
            if first_call[0]:
                first_call[0]=False
                return np.full(75,np.nan)
            return actual_residual(values)
        with patch.object(model,'residual',side_effect=residual_with_nonfinite) as residual:
            model.jacobian(self.candidate)
        self.assertEqual(residual.call_count,21)


if __name__=='__main__':
    unittest.main()
