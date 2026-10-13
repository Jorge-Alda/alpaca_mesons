"""Spectral/source/action consistency and independent baseline reductions."""
import unittest
import tempfile
from pathlib import Path
import numpy as np
from .vector_sector import (VectorSector, physical_generators, source_residues,
                            neutral_mass_matrix)
from .vector_sector_amplitudes import VectorAmplitudeInputs, tff_components, radiative_components
from .vector_sector_amplitudes import q2_abs_tff_predictor, radiative_width, radiative_width_predictor
from .vector_sector_amplitudes import two_photon_width
from .odd_parity_rcht import OddParityCouplings, VJPCouplings, VVPCouplings, ideal_vector_generators
from .radiative_sm import Q, RadiativeSMInputs, radiative_sm_components
from .pseudoscalar_tff import TFFInputs, tff_components as old_tff

class VectorSectorTests(unittest.TestCase):
    def setUp(self):
        self.c=OddParityCouplings(VJPCouplings(.01,-.02,.03,.04,.05,-.06,.07),VVPCouplings(.08,-.09,.1,.11))
        self.F=.092;self.FV=.15;self.MV=.77
        self.chi=np.diag([.02,.025,.45]);self.tp=np.diag([1,-1,0])/np.sqrt(2)
        self.sector=VectorSector.universal(self.MV,self.FV)
        self.x=VectorAmplitudeInputs(self.F,self.MV,.135,self.tp,self.chi,self.sector)

    def test_baseline_tff_at_two_virtualities(self):
        fq=np.trace((Q@Q+Q@Q)@self.tp).real
        fh=np.trace((Q@Q+Q@Q)@(self.tp@self.chi+self.chi@self.tp)).real
        old=TFFInputs(self.F,self.FV,self.MV,.135,fq,fh)
        for s,t in [(0,0),(-.4,-1.3),(.1,-2)]:
            np.testing.assert_allclose(tff_components(self.c,self.x,s,t),old_tff(self.c,old,s,t),atol=1e-14)

    def test_radiative_baseline_all_generators_including_c4(self):
        g=ideal_vector_generators()
        for tv,tp in [(g['omega'],self.tp),(g['rho+'],g['rho-']), (g['K*+'],g['K*-'])]:
            old=RadiativeSMInputs(.3,self.F,self.FV,self.MV,.81,.2,tv,tp,self.chi)
            x=VectorAmplitudeInputs(self.F,self.MV,.2,tp,self.chi,self.sector)
            np.testing.assert_allclose(radiative_components(self.c,x,tv,.81,.3),radiative_sm_components(self.c,old),atol=1e-14)

    def test_split_spectral_response_independent_matrix_inverse(self):
        theta=.14;masses=(.775,.783,1.019);g=physical_generators(theta)
        residues=source_residues(g,2*self.chi,.15,.021,-.017)
        sector=VectorSector(masses,g,residues)
        ideal=physical_generators(0)
        rotation=np.array([[np.trace(a@b).real for a in ideal] for b in g])
        matrix=rotation.T@np.diag(np.array(masses)**2)@rotation
        source=rotation.T@np.array(residues)
        for s in [-3,0,.3]:
            coeff=np.linalg.solve(matrix-s*np.eye(3),source)
            expected=sum(c*t for c,t in zip(coeff,ideal))
            np.testing.assert_allclose(sector.response(s),expected,atol=1e-14)

    def test_state_phase_and_degenerate_basis_invariance(self):
        sector=VectorSector((.77,.8,1.02),physical_generators(.12),(.1,.03,-.04))
        phase=VectorSector(sector.masses,tuple(-t for t in sector.generators),tuple(-r for r in sector.residues))
        np.testing.assert_allclose(phase.response(-.2),sector.response(-.2))
        for theta in [0,.4,-.7]:
            np.testing.assert_allclose(VectorSector.universal(.77,.15,theta).response(-.2),.15*Q/(.77**2+.2),atol=1e-14)

    def test_source_slots_and_mass_diagonalization(self):
        xl,xs=.037,.944;g=physical_generators(.12);fv,z,y=.15,.02,-.01
        residues=source_residues(g,np.diag([xl,xl,xs]),fv,z,y)
        fl,fs=fv+2*z*xl,fv+2*z*xs
        ideal=np.array([fl/np.sqrt(2),fl/(3*np.sqrt(2))+y*np.sqrt(2)*(xl-xs)/3,-fs/3+y*(xl-xs)/3])
        rotation=np.array([[1,0,0],[0,np.cos(.12),-np.sin(.12)],[0,np.sin(.12),np.cos(.12)]])
        np.testing.assert_allclose(residues,rotation@ideal,atol=1e-15)
        matrix=neutral_mass_matrix(.8,.06,.03,.02,xl,xs)
        theta=.5*np.arctan2(2*matrix[0,1],matrix[1,1]-matrix[0,0])
        rot=np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
        np.testing.assert_allclose(rot@matrix@rot.T,np.diag(np.linalg.eigvalsh(matrix)),atol=1e-14)

    def test_bose_exchange_and_poles_and_neutral_selection(self):
        sector=VectorSector((.775,.783,1.019),physical_generators(.13),(.1,.03,-.04))
        x=VectorAmplitudeInputs(self.F,self.MV,.135,self.tp,self.chi,sector)
        np.testing.assert_allclose(tff_components(self.c,x,-.4,-1.2),tff_components(self.c,x,-1.2,-.4),atol=1e-14)
        for m in sector.masses:
            with self.assertRaises(ValueError): sector.response(m*m)
        charged=VectorAmplitudeInputs(self.F,self.MV,.135,ideal_vector_generators()['rho+'],self.chi,sector)
        with self.assertRaises(ValueError): tff_components(self.c,charged)
        with self.assertRaises(ValueError): VectorSector((.77,)*3,(self.tp,)*3,(.1,)*3)

    def test_stage3_scenarios_keep_both_branches(self):
        for source in ('universal','single_trace_positive','single_trace_negative','double_trace'):
            kw=dict(rho_input='rho_neutral',kstar_input='kstar_charged',phi_input='phi_evaluation',source_hypothesis=source)
            near=VectorSector.from_stage3(**kw)
            other=VectorSector.from_stage3(**kw,branch='other')
            self.assertNotEqual(near.hypothesis,other.hypothesis)
            self.assertLess(abs(near.generators[1][2,2]),abs(other.generators[1][2,2]))
            np.testing.assert_allclose(near.masses,other.masses)

    def test_source_complex_trace_and_tampered_input_rejected(self):
        with self.assertRaises(ValueError):
            source_residues((1j*self.tp,),self.chi,.15)
        path=Path(__file__).resolve().parents[2]/'data/vector_sector/pdg2024.json'
        with tempfile.TemporaryDirectory() as tmp:
            tampered=Path(tmp)/'input.json'
            tampered.write_bytes(path.read_bytes()+b' ')
            with self.assertRaisesRegex(ValueError,'checksum'):
                VectorSector.from_stage3(rho_input='rho_neutral',kstar_input='kstar_charged',
                    phi_input='phi_evaluation',source_hypothesis='double_trace',inputs_path=tampered)

    def test_fit_predictor_endpoints_and_crossed_width_spin_factors(self):
        import math
        parameters=dict(zip([f'c{i}' for i in range(1,8)]+[f'd{i}' for i in range(1,5)],self.c.as_array()))
        pred=q2_abs_tff_predictor(parameters,dict(inputs=self.x,q2=[0,1,3]))
        self.assertEqual(pred[0],0)
        self.assertFalse(pred.flags.writeable)
        np.testing.assert_allclose(pred,[q*abs(sum(tff_components(self.c,self.x,-q,0))) for q in (0,1,3)])
        for bad in ([],[float('nan')],[-1],[1j],[[1]]):
            with self.assertRaises(ValueError):q2_abs_tff_predictor(parameters,dict(inputs=self.x,q2=bad))
        tv=ideal_vector_generators()['omega']; mv=.8;e=.3
        for mp,orientation,denominator in ((.135,'V_to_P',96*math.pi*mv**3),(.95,'P_to_V',32*math.pi*.95**3)):
            x=VectorAmplitudeInputs(self.F,self.MV,mp,self.tp,self.chi,self.sector)
            g=sum(radiative_components(self.c,x,tv,mv,e))
            expected=abs(g)**2*abs(mv**2-mp**2)**3/denominator
            self.assertAlmostEqual(radiative_width(self.c,x,tv,mv,e,orientation),expected)
            pred=radiative_width_predictor(parameters,dict(inputs=x,TV=tv,Mvector=mv,e=e,orientation=orientation))
            self.assertEqual(pred.shape,(1,));self.assertFalse(pred.flags.writeable)
            self.assertAlmostEqual(pred[0],expected)
            with self.assertRaises(ValueError):radiative_width(self.c,x,tv,mv,e,'P_to_V' if orientation=='V_to_P' else 'V_to_P')
        endpoint=VectorAmplitudeInputs(self.F,self.MV,mv,self.tp,self.chi,self.sector)
        for orientation in ('V_to_P','P_to_V'):
            with self.assertRaises(ValueError):radiative_width(self.c,endpoint,tv,mv,e,orientation)

    def test_published_ope_pion_two_photon_width(self):
        import math
        f,fv,mv,mp=self.F,self.FV,self.MV,.135
        c1,c5,d2=.021,-.015,.008
        a=3*mv/(64*math.pi**2*np.sqrt(2)*fv)
        b=3*mv**2/(64*math.pi**2*fv**2)
        couplings=OddParityCouplings(VJPCouplings(c1,c1+c5,-c1/4,.04,c5,c5-a,.07),
                       VVPCouplings(-8*d2-b+f*f/(4*fv*fv),d2,-b+f*f/(8*fv*fv),.11))
        x=VectorAmplitudeInputs(f,mv,mp,self.tp,np.diag([mp*mp,mp*mp,.45]),self.sector)
        alpha=1/137.035999084
        published=alpha**2*mp**3/(64*math.pi**3*f*f)*(1-4*math.pi**2*f*f*mp*mp/(3*mv**4))**2
        self.assertAlmostEqual(two_photon_width(couplings,x,alpha)/published,1,places=13)
        chiral=VectorAmplitudeInputs(f,mv,0,self.tp,np.zeros((3,3)),self.sector)
        self.assertEqual(two_photon_width(couplings,chiral,alpha),0)
        for invalid in (0,-.1,float('nan'),float('inf'),1j):
            with self.assertRaises((ValueError,TypeError)):two_photon_width(couplings,x,invalid)
