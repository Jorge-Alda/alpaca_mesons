"""Unit validation of independent-source ALP target projection adapter."""
import sys
import json
import unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
import derive_odd_low_energy_jacobians as J
API=J.module('odd_bounded_amplitudes')
ALP=J.module('odd_alp_target_null')
import audit_odd_alp_target_nulls as AUDIT


class TargetTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ald=json.loads((ROOT/'data/vector_sector/alda_2026nrd_arxiv_v1.json').read_text())
        cls.tree=J.module('alda_2026nrd_neutral_inputs').alda_2026nrd_neutral_inputs(mass_scheme='tree')
        ass=json.loads((ROOT/'data/vector_sector/stage3_results.json').read_text())['assumptions']
        cls.sector,_=J.action_sector(np.array([.775,.895,.018,.043,.16]),'universal',ass)
        cls.e=np.sqrt(4*np.pi*ass['alpha_em'])
    def test_sm_recovery_with_matching_source_masses(self):
        for x in self.tree.values():
            U=-np.sqrt(2)*x['TP']/x['F']
            Y=np.sqrt(2)*(x['TP']@x['chi']+x['chi']@x['TP'])/x['F']
            rows=ALP.target_rows(x['chi'],U,Y,x['mp']**2,self.sector,self.e)
            args=J.module('vector_sector_amplitudes').VectorAmplitudeInputs(
                J.F,J.OPERATOR_MV,x['mp'],x['TP'],x['chi'],self.sector)
            np.testing.assert_allclose(rows['gammagamma'],API.tff_gradients(args),atol=2e-13)
            for i,name in enumerate(('rho','omega','phi')):
                np.testing.assert_allclose(rows['Vgamma:'+name],API.radiative_gradients(
                    args,self.sector.generators[i],self.sector.masses[i],self.e),atol=2e-13)
    def test_uv_linearity_eom_and_scale(self):
        chi=self.tree['pi0']['chi'];m0=self.ald['m0_squared_gev2']
        first=ALP.sources(chi,.01,[1,0,0],0,m0)
        second=ALP.sources(chi,.01,[0,1,0],1,m0)
        total=ALP.sources(chi,.01,[1,1,0],1,m0)
        scaled=ALP.sources(chi,.01,[1,1,0],1,m0,fa=2)
        for a,b,c,d in zip(first,second,total,scaled):
            np.testing.assert_allclose(a+b,c,atol=1e-13)
            np.testing.assert_allclose(c/2,d,atol=1e-13)
        U,Y,Z=total
        np.testing.assert_allclose(Y+.02*U-2*Z*np.eye(3),0,atol=1e-13)
    def test_neutral_pole_rejected_and_three_p_zero(self):
        chi=self.tree['pi0']['chi'];m0=self.ald['m0_squared_gev2']
        pole=np.linalg.eigvalsh(chi+(m0/3)*np.ones((3,3)))[0]
        with self.assertRaises(ValueError):ALP.sources(chi,pole,[1,0,0],0,m0)
        U,Y,_=ALP.sources(chi,.01,[1,0,0],0,m0)
        rows=ALP.target_rows(chi,U,Y,.01,self.sector,self.e)
        self.assertTrue(np.array_equal(rows['3P'],np.zeros(241)))
        self.assertEqual(len(rows),11)

    def test_physical_source_action_relations(self):
        chi=self.tree['pi0']['chi'];m0=self.ald['m0_squared_gev2']
        for w in (.01,.16,1.44):
            for uv in np.eye(4):
                U,Y,Z=ALP.sources(chi,w,uv[:3],uv[3],m0)
                certificate=ALP.relation_certificate(ALP.target_rows(chi,U,Y,w,self.sector,self.e))
                self.assertEqual(certificate['relation_count'],13)
                self.assertTrue(certificate['valid'],certificate)

    def test_null_witness_and_identifiable_target(self):
        measured=np.array([[1.,0.,0.],[0.,1.,0.]])
        invisible=AUDIT.null_diagnostic(measured,[[0.,0.,2.]],np.eye(3))
        self.assertEqual(invisible['null_image_rank'],1)
        self.assertLess(invisible['measured_null_relative_defect'],1e-12)
        self.assertGreater(invisible['target_witness_norm'],1.)
        visible=AUDIT.null_diagnostic(measured,[[1.,3.,0.]],np.eye(3))
        self.assertEqual(visible['null_image_rank'],0)


if __name__=='__main__':unittest.main()
