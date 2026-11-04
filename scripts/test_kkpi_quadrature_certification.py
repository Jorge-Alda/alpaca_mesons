"""Independent quadrature comparisons for the joint KKpi likelihood.

Run: .venv/bin/python -m unittest scripts.test_kkpi_quadrature_certification
"""
import copy
import json
import unittest
from unittest.mock import patch

import numpy as np
from scipy.integrate import quad

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
    def test_maximum_orders_retain_all_failing_output_differences(self):
        class Model:
            phase='original'
            babar={channel:dict(errors=np.ones(2)) for channel in fit.CHANNELS}
            def residual(self,values):return np.array([float(self.phase[0])])
            def jacobian(self,values,finite_difference_steps):
                return np.full((1,32),1.+.01*self.phase[0])
            def unpack(self,values):return None
            def babar_blocks(self,values,state):
                return None,{channel:np.full(2,float(self.phase[0])) for channel in fit.CHANNELS}
        model=Model()
        with patch.object(fit,'phase_points',side_effect=lambda data,*orders:orders):
            result=fit.audit_quadrature(model,np.zeros(32),
                dict(energy=[1,2,4],pair=[2,4],angle=[2,4]))
        result=json.loads(json.dumps(result))
        self.assertEqual(result['status'],'uncertified')
        self.assertFalse(result['passed'])
        comparison=result['refinements']['energy']['comparisons'][-1]
        self.assertEqual(comparison['from_orders'],dict(energy=2,pair=2,angle=2))
        self.assertEqual(comparison['to_orders'],dict(energy=4,pair=2,angle=2))
        for output in ('predictions','residuals','jacobian','chi2'):
            self.assertFalse(comparison[output]['passed'])
        self.assertEqual(comparison['predictions']['maximum_change_in_sigma'],2.)
        self.assertEqual(comparison['residuals']['maximum_absolute_change'],2.)
        self.assertEqual(comparison['chi2']['absolute_change'],12.)
        self.assertAlmostEqual(comparison['jacobian']['columns'][0]['relative_change'],.02/1.04)
        self.assertIs(model.phase,Model.phase)

    def test_invalid_physical_points_are_unevaluated_with_reason_and_orders(self):
        model=fit.GlobalLikelihood(1,2,2)
        original_phase=model.phase
        for parameter,value in ((6,0.),(26,6.)):
            with self.subTest(parameter=fit.PARAMETER_NAMES[parameter]):
                values=model.initial();values[parameter]=value
                result=fit.audit_quadrature(model,values,
                    dict(energy=[1,2],pair=[2,3],angle=[2,3]))
                result=json.loads(json.dumps(result))
                self.assertEqual(result['status'],'unevaluated')
                self.assertFalse(result['passed'])
                self.assertIn('domain crossed',result['reason'])
                self.assertEqual(result['error_type'],'ValueError')
                self.assertEqual(result['attempted_orders'],dict(energy=1,pair=2,angle=2))
                self.assertEqual(result['refinements']['energy']['comparisons'],[])
                self.assertIs(model.phase,original_phase)
                self.assertIsNone(model.babar_cache_key)

    def test_later_evaluation_failure_preserves_finished_comparisons(self):
        class Model:
            phase='original'
            babar={channel:dict(errors=np.ones(2)) for channel in fit.CHANNELS}
            def residual(self,values):
                if self.phase[0]==4:raise FloatingPointError('nonfinite pole evaluation')
                return np.zeros(2)
            def jacobian(self,values,finite_difference_steps):return np.ones((2,32))
            def unpack(self,values):return None
            def babar_blocks(self,values,state):
                return None,{channel:np.ones(2) for channel in fit.CHANNELS}
        model=Model()
        with patch.object(fit,'phase_points',side_effect=lambda data,*orders:orders):
            result=fit.audit_quadrature(model,np.zeros(32),
                dict(energy=[1,2,4],pair=[2,4],angle=[2,4]))
        self.assertEqual(result['status'],'unevaluated')
        self.assertEqual(result['reason'],'nonfinite pole evaluation')
        self.assertEqual(result['attempted_orders'],dict(energy=4,pair=2,angle=2))
        self.assertEqual(len(result['refinements']['energy']['comparisons']),1)
        self.assertTrue(result['refinements']['energy']['comparisons'][0]['passed'])
        self.assertFalse(result['refinements']['energy']['passed'])
        self.assertEqual(model.phase,'original')

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


class PoleMappingTests(unittest.TestCase):
    def test_phi_and_both_kstar_pole_integrals_match_adaptive_dalitz_integration(self):
        energy=1.65
        data={channel:dict(edges=np.array([[energy,energy]])) for channel in fit.CHANNELS}
        mapped=fit.phase_points(data,1,64,32)
        for channel,masses in (('kkpi0',(.493677,.493677,.1349768)),
                ('kskpi',(.497611,.493677,.13957039))):
            first_mass,second_mass,third_mass=masses
            invariant_mass_squared=energy**2
            pair_lower=(first_mass+second_mass)**2
            pair_upper=(energy-third_mass)**2
            phase=mapped[channel]
            invariant_sum=invariant_mass_squared+sum(mass**2 for mass in masses)
            self.assertTrue(np.all(phase['s']>pair_lower))
            self.assertTrue(np.all(phase['s']<pair_upper))
            self.assertTrue(np.all(phase['weight']>=0))
            outer=fit.triangle(invariant_mass_squared,phase['s'],third_mass**2)
            inner=fit.triangle(phase['s'],first_mass**2,second_mass**2)
            span=np.sqrt(outer*inner)/phase['s']
            other_center=(first_mass**2+third_mass**2+(phase['s']+first_mass**2-second_mass**2)
                *(invariant_mass_squared-phase['s']-third_mass**2)/(2*phase['s']))
            center=invariant_sum-phase['s']-other_center
            self.assertTrue(np.all(phase['t']>center-span/2))
            self.assertTrue(np.all(phase['t']<center+span/2))
            coefficient=389379.3656*fit.ALPHA**2/(48*np.pi*invariant_mass_squared**3)
            coefficient*=np.exp(-2*invariant_mass_squared**3/2.5**6)*(2 if channel=='kskpi' else 1)
            for pole in ('phi','Kstar_t','Kstar_u'):
                with self.subTest(channel=channel,pole=pole):
                    pole_position,pole_scale=((1.019461**2,1.019461*.004249)
                        if pole=='phi' else (.8917**2,.8917*.05))
                    def pair_integrand(pair_invariant):
                        outer=fit.triangle(invariant_mass_squared,pair_invariant,third_mass**2)
                        inner=fit.triangle(pair_invariant,first_mass**2,second_mass**2)
                        span=np.sqrt(max(outer*inner,0))/pair_invariant
                        other_center=(first_mass**2+third_mass**2+(pair_invariant+first_mass**2-second_mass**2)
                            *(invariant_mass_squared-pair_invariant-third_mass**2)/(2*pair_invariant))
                        center=invariant_sum-pair_invariant-other_center
                        lower,upper=center-span/2,center+span/2
                        def angular_integrand(pair_t):
                            selected=(pair_invariant if pole=='phi' else pair_t if pole=='Kstar_t'
                                else invariant_sum-pair_invariant-pair_t)
                            return max(1-(2*(pair_t-center)/span)**2,0)/((selected-pole_position)**2+pole_scale**2)
                        angular_pole=(pole_position if pole=='Kstar_t'
                            else invariant_sum-pair_invariant-pole_position)
                        points=[angular_pole] if pole!='phi' and lower<angular_pole<upper else None
                        angular_integral=quad(angular_integrand,lower,upper,points=points,
                            epsabs=1e-9,epsrel=1e-9)[0]
                        return outer*inner/(16*pair_invariant)*angular_integral
                    points=[pole_position] if pole=='phi' else None
                    reference=coefficient*quad(pair_integrand,pair_lower,pair_upper,points=points,
                        epsabs=1e-9,epsrel=1e-8)[0]
                    selected=(phase['s'] if pole=='phi' else phase['t'] if pole=='Kstar_t'
                        else invariant_sum-phase['s']-phase['t'])
                    numerical=float(np.sum(phase['weight']/((selected-pole_position)**2+pole_scale**2)))
                    self.assertAlmostEqual(numerical/reference,1.,delta=1e-5)


if __name__=='__main__':unittest.main()
