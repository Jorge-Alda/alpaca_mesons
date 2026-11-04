"""Independent quadrature comparisons for the joint KKpi likelihood.

Run: .venv/bin/python -m unittest scripts.test_kkpi_quadrature_certification
"""
import copy
import unittest
from unittest.mock import patch

import numpy as np

from scripts import fit_babar_kkpi_global_vector_sector as fit


class QuadratureComparisonTests(unittest.TestCase):
    def setUp(self):
        self.previous=dict(predictions={channel:np.array([2.,3.]) for channel in fit.CHANNELS},
            residual=np.zeros(4),jacobian=np.ones((4,32)),chi2=10.)
        self.errors={channel:np.array([.5,2.]) for channel in fit.CHANNELS}

    def compare(self,current,scales=None):
        return fit.compare_quadrature_outputs(self.previous,current,self.errors,
            np.ones(32) if scales is None else scales,fit.QUADRATURE_TOLERANCES)

    def test_unchanged_outputs_pass_every_comparison(self):
        result=self.compare(copy.deepcopy(self.previous))
        self.assertTrue(result['passed'])
        for name in ('predictions','residuals','jacobian','chi2'):
            self.assertTrue(result[name]['passed'])

    def test_each_output_can_fail_without_other_comparisons_failing(self):
        for name in ('predictions','residuals','jacobian','chi2'):
            with self.subTest(output=name):
                current=copy.deepcopy(self.previous)
                if name=='predictions':current[name][fit.CHANNELS[0]][0]+=.006
                elif name=='residuals':current['residual'][0]=.011
                elif name=='jacobian':current[name][:,7]*=1.002
                else:current['chi2']+=.011
                result=self.compare(current)
                self.assertFalse(result['passed'])
                for comparison in ('predictions','residuals','jacobian','chi2'):
                    self.assertEqual(result[comparison]['passed'],comparison!=name)
                if name=='predictions':
                    self.assertAlmostEqual(result[name]['maximum_change_in_sigma'],.012)
                if name=='jacobian':
                    self.assertEqual([column['parameter'] for column in result[name]['columns']
                        if not column['passed']],[fit.PARAMETER_NAMES[7]])

    def test_negligible_scaled_columns_report_absolute_changes(self):
        self.previous['jacobian'][:,0]=0.
        current=copy.deepcopy(self.previous)
        current['jacobian'][:,0]=1e-15
        scales=np.ones(32);scales[0]=10.
        column=self.compare(current,scales)['jacobian']['columns'][0]
        self.assertEqual(column['comparison'],'absolute')
        self.assertIsNone(column['relative_change'])
        self.assertAlmostEqual(column['absolute_change'],2e-14)
        self.assertTrue(column['passed'])

    def test_disappearing_column_uses_relative_comparison(self):
        self.previous['jacobian'][:,0]=1e-4
        current=copy.deepcopy(self.previous);current['jacobian'][:,0]=0.
        column=self.compare(current)['jacobian']['columns'][0]
        self.assertEqual(column['comparison'],'relative')
        self.assertEqual(column['relative_change'],1.)
        self.assertFalse(column['passed'])


class IndependentRefinementTests(unittest.TestCase):
    def test_early_failure_keeps_ladder_and_point_uncertified(self):
        class Model:
            phase='original'
            babar={channel:dict(errors=np.ones(2)) for channel in fit.CHANNELS}
            def residual(self,values):
                return np.array([.02 if self.phase[0]>=4 else 0.,0.,0.,0.])
            def jacobian(self,values,finite_difference_steps):return np.ones((4,32))
            def unpack(self,values):return None
            def babar_blocks(self,values,state):
                return np.zeros(4),{channel:np.ones(2) for channel in fit.CHANNELS}
        model=Model()
        with patch.object(fit,'phase_points',side_effect=lambda data,*orders:orders):
            result=fit.audit_quadrature(model,np.zeros(32),
                dict(energy=[2,4,8],pair=[3,6,12],angle=[5,10,20]))
        for name in ('energy','combined'):
            ladder=result['refinements'][name]
            self.assertFalse(ladder['comparisons'][0]['passed'])
            self.assertTrue(ladder['comparisons'][-1]['passed'])
            self.assertFalse(ladder['passed'])
        self.assertFalse(result['passed'])
        self.assertEqual(result['status'],'uncertified')

    def test_orders_vary_one_dimension_at_a_time_then_combined(self):
        class Model:
            phase='original'
            babar={channel:dict(errors=np.ones(2)) for channel in fit.CHANNELS}
            def residual(self,values):return np.zeros(4)
            def jacobian(self,values,finite_difference_steps):
                steps.append(finite_difference_steps.copy())
                return np.ones((4,32))
            def unpack(self,values):return None
            def babar_blocks(self,values,state):
                return np.zeros(4),{channel:np.ones(2) for channel in fit.CHANNELS}
        model=Model();steps=[];evaluated=[]
        def phase(data,*orders):
            evaluated.append(orders)
            return orders
        orders=dict(energy=[2,4,8],pair=[3,6,12],angle=[5,10,20])
        values=np.arange(32,dtype=float)
        with patch.object(fit,'phase_points',side_effect=phase):
            result=fit.audit_quadrature(model,values,orders)
        self.assertEqual(set(evaluated),{(2,3,5),(4,3,5),(8,3,5),(2,6,5),
            (2,12,5),(2,3,10),(2,3,20),(4,6,10),(8,12,20)})
        self.assertEqual(len(evaluated),len(set(evaluated)))
        for dimension in ('energy','pair','angle'):
            for comparison in result['refinements'][dimension]['comparisons']:
                changed=[name for name in orders if comparison['from_orders'][name]
                    !=comparison['to_orders'][name]]
                self.assertEqual(changed,[dimension])
                for name in orders:
                    if name!=dimension:
                        self.assertEqual(comparison['to_orders'][name],orders[name][0])
        for recorded in steps:
            np.testing.assert_array_equal(recorded,np.sqrt(np.finfo(float).eps)*np.maximum(1.,abs(values)))
        self.assertEqual(model.phase,'original')
        self.assertIsNone(model.babar_cache_key)
        self.assertTrue(result['passed'])
        self.assertEqual(result['status'],'certified')

    def test_real_likelihood_energy_comparison_matches_direct_evaluations(self):
        model=fit.GlobalLikelihood(1,2,2)
        values=model.initial()
        original_phase=model.phase
        residual=model.residual(values)
        _,predictions=model.babar_blocks(values,model.unpack(values))
        model.phase=fit.phase_points(model.babar,2,2,2);model.babar_cache_key=None
        refined_residual=model.residual(values)
        _,refined_predictions=model.babar_blocks(values,model.unpack(values))
        model.phase=original_phase;model.babar_cache_key=None
        result=fit.audit_quadrature(model,values,dict(energy=[1,2],pair=[2,3],angle=[2,3]))
        comparison=result['refinements']['energy']['comparisons'][0]
        expected_prediction=max(float(np.max(abs(refined_predictions[channel]-predictions[channel])
            /model.babar[channel]['errors'])) for channel in fit.CHANNELS)
        self.assertAlmostEqual(comparison['predictions']['maximum_change_in_sigma'],expected_prediction)
        self.assertAlmostEqual(comparison['residuals']['maximum_absolute_change'],
            float(np.max(abs(refined_residual-residual))))
        self.assertAlmostEqual(comparison['chi2']['absolute_change'],
            float(abs(refined_residual@refined_residual-residual@residual)))
        self.assertEqual(len(comparison['jacobian']['columns']),32)
        self.assertTrue(all(np.isfinite(column['absolute_change']) for column in comparison['jacobian']['columns']))
        self.assertIs(model.phase,original_phase)


if __name__=='__main__':unittest.main()
