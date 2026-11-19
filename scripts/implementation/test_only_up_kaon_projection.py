"""Check outgoing physical-kaon projection in the only-up amplitude."""
import unittest
from unittest.mock import patch
import numpy as np
from .only_up_decays import only_up_inputs, three_pseudoscalar


class KaonProjectionTest(unittest.TestCase):
    def test_complex_projection_for_each_leg(self):
        values = {('K0_1', 'K0_1', 'pi+'): 1+2j,
                  ('K0_1', 'K0_2', 'pi+'): 3-1j,
                  ('K0_2', 'K0_1', 'pi+'): -2+4j,
                  ('K0_2', 'K0_2', 'pi+'): 5+3j}
        def prepare(*args, **kwargs):
            value = values[tuple(args[5])]
            return lambda s12, s23: value
        epsilon = 2.228e-3*np.exp(1j*np.deg2rad(43.52))
        p = np.conjugate((1+epsilon)/np.sqrt(2*(1+abs(epsilon)**2)))
        q = np.conjugate((1-epsilon)/np.sqrt(2*(1+abs(epsilon)**2)))
        short = ((p+q)/np.sqrt(2), 1j*(-p+q)/np.sqrt(2))
        long = ((p-q)/np.sqrt(2), 1j*(-p-q)/np.sqrt(2))
        expected = sum(short[i]*long[j]*values[(f'K0_{i+1}', f'K0_{j+1}', 'pi+')]
                       for i in range(2) for j in range(2))
        with patch('calculate_alp_three_pseudoscalar_finite.prepare_amplitude', side_effect=prepare):
            actual = three_pseudoscalar(only_up_inputs(.01), ('K_S', 'K_L', 'pi+'))(.3, .4)
        self.assertAlmostEqual(actual, expected)

    def test_real_calculator_accepts_physical_kaons(self):
        amplitude = three_pseudoscalar(only_up_inputs(.01), ('pi0', 'K_S', 'K_L'))
        self.assertTrue(np.isfinite(amplitude(.3, .4)))

class KaonWidthTest(unittest.TestCase):
    def test_closed_channels(self):
        from .only_up_decays import three_pseudoscalar_width
        for radial_mass_squared in (None, 100.):
            self.assertEqual(three_pseudoscalar_width(only_up_inputs(.01),
                             ('pi0', 'K_S', 'K_L'), radial_mass_squared), 0.)

    def test_open_channels_and_identical_factor(self):
        from .only_up_decays import three_pseudoscalar_width
        inputs = only_up_inputs(1.)
        def integral(integrand, limits, lower, upper, mass_sum, poles, epsrel):
            s = (lower + upper)/2
            minimum, maximum = limits(s)
            value = integrand(s)((minimum + maximum)/2)
            self.assertTrue(np.isfinite(value))
            return 1.
        with patch('alpaca_chiral._three_body_quadrature.integrate_three_body', side_effect=integral):
            for radial_mass_squared in (None, 100.):
                mixed = three_pseudoscalar_width(inputs, ('pi0', 'K_S', 'K_L'), radial_mass_squared)
                identical = three_pseudoscalar_width(inputs, ('pi0', 'K_S', 'K_S'), radial_mass_squared)
                self.assertEqual(mixed, 0.)
                self.assertGreater(identical, 0.)

    def test_symmetry_zero_before_integration(self):
        from .only_up_decays import three_pseudoscalar_width
        inputs = only_up_inputs(1.)
        with patch('alpaca_chiral._three_body_quadrature.integrate_three_body', side_effect=AssertionError('integration reached')):
            for spectator in ('pi0', 'eta', 'eta_prime'):
                for radial in (None, 100.):
                    labels = ('K_L', spectator, 'K_S')
                    self.assertEqual(three_pseudoscalar(inputs, labels, radial)(.4, .5), 0j)
                    self.assertEqual(three_pseudoscalar_width(inputs, labels, radial), 0.)
