"""Baseline normalization, trace selection, Bose/crossing and phase checks."""
import unittest
from dataclasses import replace
import numpy as np
from .odd_parity_rcht import OddParityCouplings,VJPCouplings,VVPCouplings
from .vector_sector import VectorSector,physical_generators,source_residues
from .vector_sector_amplitudes import VectorAmplitudeInputs
from . import vector_sector_amplitudes as old
from . import odd_trace_amplitudes as new

class TraceAmplitudes(unittest.TestCase):
    def setUp(self):
        self.gens=physical_generators(.08)
        residues=source_residues(self.gens,np.diag([.04,.04,.9]),.15,.003,.02)
        self.sector=VectorSector((.775,.783,1.019),self.gens,residues)
        self.x=VectorAmplitudeInputs(.093,.770,.135,self.gens[0],np.diag([.02,.03,.46]),self.sector)
        self.c=OddParityCouplings(VJPCouplings(.001,-.002,.003,.004,.005,-.006,.007),
                                  VVPCouplings(.008,-.009,.010,.011))

    def test_all_baseline_coefficients(self):
        for row in np.eye(11):
            c=OddParityCouplings(VJPCouplings(*row[:7]),VVPCouplings(*row[7:]))
            traced=new.OddTraceCouplings.from_baseline(c)
            for s,t in ((0,0),(-.1,-.3),(-1,0)):
                np.testing.assert_allclose(new.tff_components(traced,self.x,s,t),old.tff_components(c,self.x,s,t),atol=1e-13,rtol=1e-13)
            for TV in self.gens:
                np.testing.assert_allclose(new.radiative_components(traced,self.x,TV,.783,.3028),old.radiative_components(c,self.x,TV,.783,.3028),atol=1e-13,rtol=1e-13)
        # Charged c4 normalization as well as neutral limits.
        TV=np.zeros((3,3));TV[2,0]=1
        charged=replace(self.x,TP=TV.T,mp=.494)
        traced=new.OddTraceCouplings.from_baseline(self.c)
        np.testing.assert_allclose(new.radiative_components(traced,charged,TV,.892,.3028),
                                  old.radiative_components(self.c,charged,TV,.892,.3028),atol=1e-13)

    def test_bose_with_all_double_trace_terms(self):
        c=new.OddTraceCouplings((self.c.vjp,replace(self.c.vjp,c4=0),replace(self.c.vjp,c4=0),replace(self.c.vjp,c4=0)),
                               (self.c.vvp,self.c.vvp,self.c.vvp),.013)
        for s,t in ((-.1,-.3),(-1.,0.)):
            np.testing.assert_allclose(new.tff_components(c,self.x,s,t),new.tff_components(c,self.x,t,s),atol=1e-13)

    def test_singleton_photon_trace_is_exactly_blind(self):
        c=new.OddTraceCouplings((new.ZC,new.ZC,replace(self.c.vjp,c4=0),new.ZC),(new.ZD,)*3)
        self.assertAlmostEqual(abs(new.C(c,self.x,self.gens[1],new.Q,.6,0)),0,places=14)
        self.assertAlmostEqual(abs(new.tff_components(c,self.x,-.2,0)[1]),0,places=14)

    def test_double_trace_family_normalization(self):
        c=replace(new.ZC,c1=.03,c3=.01)
        traced=new.OddTraceCouplings((new.ZC,c,new.ZC,new.ZC),(new.ZD,)*3)
        A=self.gens[1];B=new.Q;Z=self.x.TP;mass=Z@self.x.chi+self.x.chi@Z
        expected=2*np.trace(A)*np.trace(B@Z)*c.c1*(self.x.mp**2-.5+.2)
        expected+=8*c.c3*np.trace(A)*np.trace(B@mass)
        self.assertAlmostEqual(abs(new.C(traced,self.x,A,B,.5,.2)-expected),0,places=14)

    def test_antisymmetric_family_sign_and_real_photon_zero(self):
        c=new.OddTraceCouplings((new.ZC,)*4,(new.ZD,)*3,.03)
        A=self.sector.response(-.1);B=self.sector.response(-.4)
        minus=np.trace(A)*np.trace(B@self.x.TP)-np.trace(B)*np.trace(A@self.x.TP)
        self.assertAlmostEqual(abs(new.H(c,self.x,A,B,-.1,-.4)+.03*minus*.3),0,places=14)
        self.assertEqual(new.H(c,self.x,A,B,0,0),0)

    def test_width_crossing_and_internal_phase(self):
        c=new.OddTraceCouplings.from_baseline(self.c)
        x=replace(self.x,mp=1.2)
        G=abs(sum(new.radiative_components(c,x,self.gens[0],.775,.3028)))**2
        self.assertAlmostEqual(new.radiative_width(c,x,self.gens[0],.775,.3028,'P_to_V'),
                              G*(1.2**2-.775**2)**3/(32*np.pi*1.2**3))
        flipped=VectorSector(self.sector.masses,tuple(-t for t in self.gens),tuple(-r for r in self.sector.residues))
        np.testing.assert_allclose(new.tff_components(c,self.x,-.2,0),new.tff_components(c,replace(self.x,sector=flipped),-.2,0),atol=1e-14)

    def test_invalid_and_poles(self):
        with self.assertRaises(ValueError):new.OddTraceCouplings((self.c.vjp,)*4,(self.c.vvp,)*3)
        with self.assertRaises(ValueError):new.OddTraceCouplings((new.ZC,)*4,(new.ZD,)*3,complex(1,1))
        with self.assertRaises(ValueError):new.tff_components(new.OddTraceCouplings.from_baseline(self.c),self.x,self.sector.masses[0]**2,0)

if __name__=='__main__':unittest.main()
