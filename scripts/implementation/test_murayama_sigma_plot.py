import unittest
from unittest.mock import patch

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from scripts.implementation import murayama


class SigmaPlotTest(unittest.TestCase):
    def test_mass_labels_keep_both_branches(self):
        kappa = np.logspace(-3, 0, 500)
        masses = (np.log10(kappa) + 1.5)**2 + 1
        with patch.object(murayama, 'mesons', ['pi0']), \
             patch.object(murayama, 'sigma_two_pseudoscalar_width', return_value=1.), \
             patch.object(plt, 'tight_layout'), patch.object(plt, 'savefig'):
            murayama.sigma_PP(masses)
        axis = plt.gcf().axes[1]
        positions = axis.get_xticks()
        labels = np.array([float(label.get_text()) for label in axis.get_xticklabels()])
        np.testing.assert_allclose(np.interp(positions, kappa, masses), labels, rtol=1e-4)
        self.assertTrue(all(value in (1., 2., 5.) for value in labels))
        self.assertTrue(np.any(positions < 10**-1.5))
        self.assertTrue(np.any(positions > 10**-1.5))
        self.assertTrue(np.all(np.diff(np.log10(positions)) >= 0.16 - 1e-10))
        plt.close('all')


if __name__ == '__main__':
    unittest.main()
