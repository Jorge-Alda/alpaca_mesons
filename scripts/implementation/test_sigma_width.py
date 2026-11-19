import math
import unittest

from scripts.implementation.only_up_decays import sigma_two_pseudoscalar_width
from scripts.implementation.only_up_inputs import only_up_inputs
from scripts.implementation.only_up_decays import spectrum
from calculate_alp_three_pseudoscalar_finite import prepare_amplitude


class SigmaWidthTest(unittest.TestCase):
    def test_all_physical_pairs_and_selection_rules(self):
        from itertools import combinations_with_replacement
        labels = ('a', 'pi0', 'eta', 'eta_prime', 'pi+', 'pi-', 'K+', 'K-',
                  'K0', 'K0bar', 'K_S', 'K_L')
        neutral = {'a', 'pi0', 'eta', 'eta_prime'}
        for pair in combinations_with_replacement(labels, 2):
            width = sigma_two_pseudoscalar_width(0.01, 10., pair)
            self.assertGreaterEqual(width, 0.)
            if not (set(pair) <= neutral or set(pair) in ({'pi+', 'pi-'}, {'K+', 'K-'})):
                self.assertEqual(width, 0., pair)
        self.assertGreater(sigma_two_pseudoscalar_width(0.01, 10., ('K+', 'K-')), 0.)

    def test_pole_masses_match_finite_mixing_calculator(self):
        inputs = only_up_inputs(0.3)
        amplitude = prepare_amplitude(inputs['f_a'], inputs['m_a0'], inputs['hat_S'],
                                      inputs['hat_P'], inputs['hat_Theta'],
                                      ['pi0', 'eta', 'eta_prime'], f_pi=inputs['f_pi'],
                                      B_0=inputs['B_0'], m_0_squared=inputs['m_0_squared'])
        for name in ('a', 'pi0', 'eta', 'eta_prime'):
            self.assertAlmostEqual(amplitude.masses[name], spectrum(inputs)[name])

    def test_charged_vertex_and_threshold(self):
        kappa, coupling = 0.01, 10.
        inputs = only_up_inputs(kappa)
        mass = math.sqrt(inputs['B_0']*inputs['hat_S'][0][0, 0]*inputs['f_pi']**2/inputs['f_a']**2 + 2*coupling*inputs['f_a']**2)
        pion = spectrum(inputs)['pi+']
        vertex = inputs['B_0']*inputs['hat_S'][0][0, 0]/inputs['f_a']
        expected = vertex**2*math.sqrt(1-4*pion**2/mass**2)/(16*math.pi*mass)
        self.assertAlmostEqual(sigma_two_pseudoscalar_width(kappa, coupling, ['pi+', 'pi-']), expected)
        self.assertEqual(sigma_two_pseudoscalar_width(kappa, 0.001, ['pi+', 'pi-']), 0.)

    def test_neutral_pair_order(self):
        first = sigma_two_pseudoscalar_width(0.01, 10., ['a', 'pi0'])
        second = sigma_two_pseudoscalar_width(0.01, 10., ['pi0', 'a'])
        self.assertAlmostEqual(first, second)
        self.assertGreater(sigma_two_pseudoscalar_width(0.01, 10., ['a', 'a']), 0.)
