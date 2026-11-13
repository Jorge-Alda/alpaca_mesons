import unittest
import numpy as np
from calculate_alp_three_pseudoscalar import decay_width, prepare_amplitude

class ThreeBodyMasses(unittest.TestCase):
    def setUp(self):
        self.zero = np.zeros((3,3))
        self.scalar = [np.diag([.003,.005,.1]),self.zero,self.zero,self.zero,self.zero]
        self.pseudo = [self.zero,np.diag([.0001,.0002,0]),self.zero,self.zero]
        self.theta = [0.,0.,0.,0.]
        self.constants = dict(f_pi=.092,B_0=2.5,m_0_squared=.75)
    def amplitude(self, mesons):
        return prepare_amplitude(10.,.6,self.scalar,self.pseudo,self.theta,mesons,**self.constants)
    def width(self, mesons):
        return decay_width(10.,.6,self.scalar,self.pseudo,self.theta,mesons,**self.constants)
    def test_permutation_and_width(self):
        first=self.amplitude(('pi0','pi+','pi-'))
        second=self.amplitude(('pi+','pi0','pi-'))
        s12,s23=.09,.12
        s13=sum(first.masses[x]**2 for x in ('a','pi0','pi+','pi-'))-s12-s23
        self.assertAlmostEqual(abs(first(s12,s23)-second(s12,s13)),0.,places=13)
        self.assertGreater(self.width(('pi0','pi+','pi-')),0.)
        self.assertAlmostEqual(self.width(('pi0','pi+','pi-')),self.width(('pi+','pi0','pi-')),places=19)
    def test_closed_forbidden_and_zero(self):
        self.assertEqual(self.width(('pi+','pi+','pi0')),0.)
        self.assertEqual(self.width(('eta','eta','eta')),0.)
        self.assertEqual(decay_width(10.,.6,self.scalar,[self.zero]*4,self.theta,('pi0','pi+','pi-'),**self.constants),0.)
    def test_missing_higher_sources(self):
        with self.assertRaisesRegex(ValueError, 'requires'):
            prepare_amplitude(10.,.6,self.scalar[:3],self.pseudo,self.theta,('pi0','pi+','pi-'),**self.constants)

if __name__ == '__main__': unittest.main()
