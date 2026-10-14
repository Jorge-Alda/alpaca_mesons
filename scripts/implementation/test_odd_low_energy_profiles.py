"""Finite profile algebra tests. Run .venv/bin/python scripts/implementation/test_odd_low_energy_profiles.py."""
import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import profile_odd_low_energy as P

class ProfileTests(unittest.TestCase):
    def test_identifiable_covariance_with_null_column(self):
        J=np.array([[1.,0.,0.],[0.,2.,0.]])
        H=np.array([[3.,4.,0.]])
        cov,rank,defect,_=P.identifiable_covariance(J,H)
        self.assertEqual(rank,2)
        self.assertLess(defect,1e-12)
        np.testing.assert_allclose(cov,[[13.]])
        H[0,2]=1.
        self.assertGreater(P.identifiable_covariance(J,H)[2],.1)

    def test_covariance_coordinate_scaling(self):
        J=np.array([[1.,2.,0.],[0.,1.,3.],[2.,0.,1.]])
        H=np.array([[2.,3.,4.]])
        scales=np.array([.01,100.,3.])
        np.testing.assert_allclose(P.identifiable_covariance(J,H)[0],P.identifiable_covariance(J*scales,H*scales)[0],rtol=1e-12)

    def test_degenerate_contrast_choice(self):
        choices=P.deterministic_eigenvectors(np.eye(3),2)
        np.testing.assert_allclose(choices[0][1],[1.,0.,0.])
        np.testing.assert_allclose(choices[1][1],[0.,1.,0.])

    def test_constraint_at_moving_nuisance(self):
        class Toy:
            gauge=np.eye(2)
            def blocks(self,v):
                return None,None,None,np.array([[1.+v[0],2.],[3.,4.-v[0]]]),np.array([v[0],-v[0]])
        toy=Toy();c=np.array([.3,.7]);reference=np.array([2.,1.])
        for nuisance in (-.2,.4):
            v=np.array([nuisance]);part,null=P.constraint_coordinates(toy,v,c,reference,1.25)
            coefficients=part+null@np.array([.9])
            rows,offset=toy.blocks(v)[3:5]
            self.assertAlmostEqual(float(c@(offset+rows@coefficients-reference)),1.25,places=12)

if __name__=='__main__':unittest.main()
