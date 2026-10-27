"""Independent checks of the Born spin-one orientation integral.

Run: .venv/bin/python -m unittest discover -s scripts -p test_alp_three_pion_born_spin_one.py
"""

import unittest

import numpy as np
import sympy as sp
from scipy.integrate import quad
from derive_alp_three_pion_born_spin_one import derive


class BornSpinOneTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.results = derive()

    def test_polar_primitive(self):
        cosine = sp.Symbol('z', real=True)
        self.assertEqual(sp.simplify(sp.diff(self.results['primitive'], cosine)
                                     - self.results['polar_integrand']), 0)

    def test_regular_zero_momentum_limit(self):
        constant = sp.Symbol('a', positive=True)
        self.assertEqual(self.results['zero_momentum_limit'], 4/(3*constant))

    def test_cyclic_orientation_and_azimuth(self):
        energy = sp.Symbol('sqrt_x', positive=True)
        self.assertEqual(self.results['cyclic_ratios'], (2/energy,)*3)
        cosine = sp.Symbol('z', real=True)
        photon_momentum = sp.Symbol('abs_q', positive=True)
        self.assertEqual(sp.factor(self.results['azimuth_ratio']
                                  - photon_momentum**2*(1-cosine**2)/2), 0)

    def test_covariant_epsilon_current_sign(self):
        self.assertEqual(self.results['epsilon_difference'], sp.zeros(4, 1))

    def test_angular_overlap_against_direct_current(self):
        # An asymmetric triangle and energies avoid accidental exchange zeros.
        pions = np.array([[0.7, 0., 0.], [-0.2, 0.6, 0.], [-0.5, -0.6, 0.]])
        pion_mass = 0.14
        energies = np.sqrt(np.sum(pions*pions, axis=1)+pion_mass**2)
        total_energy = np.sum(energies)
        photon_energy, photon_momentum = 0., 2.1
        photon_squared = photon_energy**2-photon_momentum**2
        normal = -total_energy*np.cross(pions[0], pions[1])
        nodes, weights = np.polynomial.legendre.leggauss(40)
        average = 0.
        for cosine, weight in zip(nodes, weights):
            for azimuth in 2*np.pi*np.arange(80)/80:
                photon = photon_momentum*np.array([np.sqrt(1-cosine**2)*np.cos(azimuth),
                                                   np.sqrt(1-cosine**2)*np.sin(azimuth), cosine])
                current = np.zeros(3)
                for r, s, t in ((0, 2, 1), (1, 0, 2), (2, 1, 0)):
                    difference = ((total_energy-photon_energy)*(energies[s]-energies[t])
                                  + np.dot(photon, pions[s]-pions[t]))
                    denominator = photon_squared-2*photon_energy*energies[r]+2*np.dot(photon, pions[r])
                    current += difference*(2*pions[r]-photon)/denominator
                average += weight*np.dot(normal, np.cross(current, photon))/(2*80)
        integrals = []
        for pion, energy in zip(pions, energies):
            constant = photon_squared-2*photon_energy*energy
            linear = 2*photon_momentum*np.linalg.norm(pion)
            integrals.append(quad(lambda cosine: (1-cosine**2)/(constant+linear*cosine), -1, 1)[0])
        predicted = photon_momentum**2*np.dot(normal, normal)/total_energy*sum(integrals)
        self.assertAlmostEqual(average, predicted, places=9)


if __name__ == '__main__':
    unittest.main()
