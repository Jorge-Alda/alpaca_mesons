"""Checks for the ALP mass axis in decay plots."""
import unittest
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scripts.implementation.murayama import alp_mass_axis


class AlpMassAxisTest(unittest.TestCase):
    def test_ticks_stay_within_mass_range_and_map_to_their_values(self):
        kappa = np.logspace(-3, 0, 500)
        masses = 2.34 * kappa ** 1.03
        figure, axis = plt.subplots()
        axis.set_xscale('log')
        axis.set_xlim(kappa[0], kappa[-1])
        mass_axis = alp_mass_axis(axis, kappa, masses)
        figure.canvas.draw()
        major = mass_axis.xaxis.get_majorticklocs()
        minor = mass_axis.xaxis.get_minorticklocs()
        self.assertTrue(np.all((major >= masses[0]) & (major <= masses[-1])))
        self.assertTrue(np.all((minor >= masses[0]) & (minor <= masses[-1])))
        self.assertNotIn(10., major)
        positions = mass_axis._functions[1](major)
        np.testing.assert_allclose(mass_axis._functions[0](positions), major)
        plt.close(figure)

if __name__ == '__main__':
    unittest.main()
