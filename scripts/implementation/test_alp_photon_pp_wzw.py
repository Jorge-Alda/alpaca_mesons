"""Physical-channel tests of the no-mixing WZW coefficient."""
import math
import unittest
from unittest.mock import patch

import numpy as np

from alpaca_chiral.pseudoscalar_diag import ALPcouplings
from alpaca_chiral.three_pseudoscalar import F_PI, _pseudoscalars
from alpaca_chiral.vector_inputs import _read_vector_inputs
from .alp_photon_pp_wzw import photon_two_pseudoscalar_wzw_amplitude


class WZWAmplitudeTests(unittest.TestCase):
    """Check anomaly normalization, physical projections and basis input."""

    def setUp(self):
        """Set nonzero diagonal and flavour-changing no-mixing couplings."""
        self.couplings = ALPcouplings(
            0., np.diag([2., 0., 0.]), np.zeros((3, 3)),
            np.array([[4., .3+.7j, 0.], [.3-.7j, 6., 0.], [0., 0., 0.]]),
            np.zeros((3, 3)))
        self.prefactor = 3*math.sqrt(4*math.pi*_read_vector_inputs()['alpha_em'])/(12*math.pi**2*1000*F_PI**2)

    def amplitude(self, first, second, **kwargs):
        """Evaluate safely above all two-meson thresholds."""
        return photon_two_pseudoscalar_wzw_amplitude(
            3., 1000., self.couplings, first, second, **kwargs)

    def test_charged_anomaly_normalization(self):
        """Use independent analytic coefficients from the channel table."""
        for first, second in [('pi+', 'pi-'), ('K+', 'K-')]:
            self.assertAlmostEqual(self.amplitude(first, second), -2*self.prefactor)
        self.assertAlmostEqual(self.amplitude('K0', 'K0bar'), -2*self.prefactor)

    def test_physical_neutral_generators(self):
        """Evaluate the diagonal-generator identity without the six traces."""
        for neutral in ('pi0', 'eta', 'eta_prime'):
            diagonal = np.diag(_pseudoscalars()[neutral][1])
            expected = (diagonal[2]-diagonal[1])*math.sqrt(2)*(.3+.7j)/2
            self.assertAlmostEqual(self.amplitude('K0', neutral), self.prefactor*expected)

    def test_kaon_projection(self):
        """Check coherent outgoing p/q projection and antisymmetry."""
        epsilon_k = 2.228e-3*np.exp(1j*np.deg2rad(43.52))
        normalization = math.sqrt(2*(1+abs(epsilon_k)**2))
        kaon_p, kaon_q = (1+epsilon_k)/normalization, (1-epsilon_k)/normalization
        for state, sign in [('KS', 1), ('KL', -1)]:
            expected = np.conjugate(kaon_p)*self.amplitude('K0', 'eta') + sign*np.conjugate(kaon_q)*self.amplitude('K0bar', 'eta')
            self.assertAlmostEqual(self.amplitude(state, 'eta'), expected)
        self.assertAlmostEqual(self.amplitude('KS', 'KL'), -2*np.conjugate(kaon_p*kaon_q)*self.amplitude('K0', 'K0bar'))

    def test_no_basis_conversion(self):
        """Fail if a coupling translation is attempted."""
        with patch.object(self.couplings, 'translate', side_effect=AssertionError):
            self.assertAlmostEqual(self.amplitude('pi+', 'pi-'), -2*self.prefactor)

    def test_order_threshold_and_scale(self):
        """Verify Bose antisymmetry, threshold, orientation and fa scaling."""
        self.assertAlmostEqual(self.amplitude('eta', 'K0'), -self.amplitude('K0', 'eta'))
        self.assertEqual(self.amplitude('KS', 'KS'), 0j)
        self.assertEqual(self.amplitude('pi+', 'K-'), 0j)
        self.assertAlmostEqual(self.amplitude('pi+', 'pi-', orientation=-1), 2*self.prefactor)
        self.assertEqual(photon_two_pseudoscalar_wzw_amplitude(.2, 1000., self.couplings, 'pi+', 'pi-'), 0j)
        self.assertAlmostEqual(photon_two_pseudoscalar_wzw_amplitude(3., 2000., self.couplings, 'pi+', 'pi-'), -self.prefactor)


if __name__ == '__main__':
    unittest.main()
