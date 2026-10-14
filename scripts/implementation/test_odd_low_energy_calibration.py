"""Unit checks for exact fast traces, sign starts and measured-image gauges."""
import sys,unittest,json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
import calibrate_odd_low_energy as C

class CalibrationTests(unittest.TestCase):
    def test_signed_linear_recovery_and_phase_covariance(self):
        A=np.array([[1.,0.],[0.,2.],[1.,1.]])
        c=np.array([.4,-.2]);offset=np.array([.1,0.,.03]);target=np.abs(A@c+offset)
        starts,score,count=C.signed_starts(A,offset,target,np.full(3,.01),8)
        self.assertEqual(count,8);self.assertLess(score,1e-20)
        self.assertTrue(any(np.linalg.norm(x-c)<1e-10 for x in starts))
        flipped=C.signed_starts(-A,-offset,target,np.full(3,.01),8)
        self.assertLess(flipped[1],1e-20)
    def test_gauge_contains_image_and_detects_nuisance_lift(self):
        A=np.array([[1.,2.,0.],[0.,0.,1.]])
        gauge,rank=C.image_gauge(A,np.eye(3));self.assertEqual(rank,2)
        self.assertLess(C.containment(A,np.eye(3),gauge),1e-12)
        lift=np.vstack([A,[0.,1.,0.]])
        self.assertGreater(C.containment(lift,np.eye(3),gauge),1e-3)
    def test_fast_channels_at_displaced_scenarios(self):
        cal=json.loads((ROOT/'data/vector_sector/stage3_results.json').read_text());alda=json.loads((ROOT/'data/vector_sector/alda_arxiv_v1.json').read_text());ass=cal['assumptions'];tree=C.J.module('alda_neutral_inputs').alda_neutral_inputs(mass_scheme='tree',F=C.J.F,fpi=C.J.F)
        fast=C.FastChannels(tree,alda['physical_masses_gev'],np.sqrt(4*np.pi*ass['alpha_em']))
        for hypothesis,p in [('universal',[.775,.895,.018,.043,.16]),('double_trace',[.774,.894,.019,.044,.15,.2,-.1]),('single_trace_negative',[.776,.896,.685,-.61,.15,-.2])]:
            sector,mk=C.J.action_sector(np.array(p),hypothesis,ass)
            exact=C.J.channel_rows(sector,mk,tree,alda['physical_masses_gev'],fast.e);new=fast(sector,mk)
            for key in exact:
                np.testing.assert_allclose(new[key]['row'],exact[key]['row'],rtol=3e-12,atol=3e-12)
                np.testing.assert_allclose(new[key]['offset'],exact[key]['offset'],rtol=3e-12,atol=3e-12)
                np.testing.assert_allclose(new[key]['phase_factor'],exact[key]['phase_factor'],rtol=3e-12,atol=3e-12)
    def test_inverse_mass_root_certificate(self):
        import verify_odd_low_energy_calibration as V
        ass=json.loads((ROOT/'data/vector_sector/stage3_results.json').read_text())['assumptions']
        # Two actual stage3 roots retain the physical inverse-mass branch.
        cal=json.loads((ROOT/'data/vector_sector/stage3_results.json').read_text());pdg=json.loads((ROOT/'data/vector_sector/pdg2024.json').read_text());scenario=cal['scenarios'][0]
        roots=sorted(scenario['extended_mass_branches'],key=lambda x:abs(x['theta_deg']))
        for branch,r in zip(('near_ideal','other'),roots):
            p=[pdg['masses'][scenario['rho_input']]['value_mev']/1000,pdg['masses'][scenario['kstar_input']]['value_mev']/1000,r['delta_gev2'],r['kappa'],.16]
            self.assertTrue(V.branch_certificate(p,'universal',branch,ass)['valid'])
            other='other' if branch=='near_ideal' else 'near_ideal'
            self.assertFalse(V.branch_certificate(p,'universal',other,ass)['valid'])

    def test_local_observable_covariance_is_basis_invariant(self):
        import verify_odd_low_energy_calibration as V
        j=json.loads((ROOT/'data/vector_sector/odd_low_energy_calibration.json').read_text());r=next(x for x in j['results'] if x['model']=='original_11' and x['control']=='narrow_nonconstrained_summary_control')
        ass=json.loads((ROOT/'data/vector_sector/stage3_results.json').read_text())['assumptions'];alda=json.loads((ROOT/'data/vector_sector/alda_arxiv_v1.json').read_text());data=json.loads((ROOT/'data/odd_parity/odd_low_energy_measurements.json').read_text());tree=C.J.module('alda_neutral_inputs').alda_neutral_inputs(mass_scheme='tree',F=C.J.F,fpi=C.J.F);fast=C.FastChannels(tree,alda['physical_masses_gev'],np.sqrt(4*np.pi*ass['alpha_em']));mapping=C.model_maps()['original_11']
        left=V.local_diagnostics(r,data,ass,fast,mapping);right=V.local_diagnostics(r,data,ass,fast,mapping@np.diag(np.arange(1,12)))
        np.testing.assert_allclose(left['local_prediction_covariance'],right['local_prediction_covariance'],rtol=2e-8,atol=1e-20)
        self.assertLess(left['derivative_step_relative_change'],1e-6)
        self.assertEqual([x['joint_nuisance_identifiable'] for x in left['held_out']],[x['joint_nuisance_identifiable'] for x in right['held_out']])

    def test_symbolic_width_derivatives(self):
        deriv=C.symbolic_derivation();self.assertIn('odd_derivative',deriv);self.assertIn('total_derivative',deriv)

if __name__=='__main__':unittest.main()
