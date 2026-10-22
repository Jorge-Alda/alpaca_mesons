"""Physical invariants and pre-integration zeros of the photon PP width."""
import unittest
from dataclasses import replace
from unittest.mock import patch

import numpy as np

from .alp_photon_pp_width import (
    PhotonPPAmplitude, alp_photon_pp_width, build_photon_pp_amplitude,
    photon_pp_t_limits,
)
from .alp_vv_width import read_preliminary_fit_parameters
from .three_meson_width import _kaon_mixing_coefficients


class PhotonPPTests(unittest.TestCase):
    """Checks independent kinematics, Bose ordering and coherent kaon states."""

    @classmethod
    def setUpClass(cls):
        cls.fit = read_preliminary_fit_parameters()
        cls.axial = np.diag([0.8, -0.3, 0.15])
        cls.pseudoscalar = np.diag([0.1, 0.2, -0.1])

    def build(self, pair, **kwargs):
        return build_photon_pp_amplitude(1.5, self.axial, self.pseudoscalar,
                                         pair, 100., fit_parameters=self.fit, **kwargs)

    def test_gram_determinant(self):
        amplitude = PhotonPPAmplitude(1.5, (0.14, 0.49), ('pi+', 'K-'),
                                      ((0, 0.77, 0.15, 1., 0.),))
        s = 0.9
        lower, upper = photon_pp_t_limits(1.5, 0.14, 0.49, s)
        t = lower+0.37*(upper-lower)
        u = 1.5**2+0.14**2+0.49**2-s-t
        gram = np.array([[0., (u-0.14**2)/2, (t-0.49**2)/2],
                         [(u-0.14**2)/2, 0.14**2, (s-0.14**2-0.49**2)/2],
                         [(t-0.49**2)/2, (s-0.14**2-0.49**2)/2, 0.49**2]])
        self.assertAlmostEqual(amplitude.squared(s, t), np.linalg.det(gram)*abs(amplitude(s,t))**2)
        self.assertEqual(amplitude.squared(s, lower), 0.)
        self.assertEqual(amplitude.squared(s, upper), 0.)

    def test_width_normalization_massless_constant(self):
        # A constant scalar coefficient has integral m_a^10/480.
        class ConstantAmplitude(PhotonPPAmplitude):
            def __call__(self, s, t):
                return 1.+0j
        amplitude = ConstantAmplitude(1.5, (0., 0.), ('pi+', 'pi-'),
                                      ((0, 3., 0., 1., 0.),))
        with patch('scripts.implementation.alp_photon_pp_width.build_photon_pp_amplitude',
                   return_value=amplitude):
            width = alp_photon_pp_width(1.5, self.axial, self.pseudoscalar,
                                       ('pi+', 'pi-'), 100.)
        expected = 1.5**7/(122880*np.pi**3)
        self.assertAlmostEqual(width/expected, 1., places=12)

    def test_ordering(self):
        first, reversed_pair = self.build(('K+', 'K-')), self.build(('K-', 'K+'))
        s, t = 1.2, 0.6
        u = 1.5**2+sum(m*m for m in first.masses)-s-t
        self.assertAlmostEqual(abs(first(s,t)+reversed_pair(s,u)), 0., places=13)

    def test_kaon_bra_projection(self):
        epsilon = 0.12+0.08j  # Resolve conjugation unambiguously.
        flavor = self.build(('K0', 'K0bar'))
        reverse = self.build(('K0bar', 'K0'))
        physical = self.build(('KL', 'KS'), epsilon_k=epsilon)
        p_mixing, q_mixing = _kaon_mixing_coefficients(epsilon)
        expected = (p_mixing*q_mixing).conjugate()*(flavor(1.2,0.6)-reverse(1.2,0.6))
        self.assertAlmostEqual(abs(physical(1.2,0.6)-expected), 0., places=13)
        self.assertTrue(physical.terms)

    def test_forbidden_before_integration(self):
        with patch('scripts.implementation.alp_photon_pp_width.quad', side_effect=AssertionError('quadrature reached')):
            for pair, mass in ((('pi+', 'pi-'), 0.2), (('pi+', 'pi0'), 1.5),
                               (('pi0', 'eta'), 1.5), (('KS', 'KS'), 1.5),
                               (('KL', 'KL'), 1.5), (('K0', 'K0'), 1.5)):
                with self.subTest(pair=pair):
                    self.assertEqual(alp_photon_pp_width(mass, self.axial, self.pseudoscalar,
                        pair, 100., fit_parameters=self.fit), 0.)
        empty = self.build(('KS', 'KS'))
        self.assertEqual(empty(object(), object()), 0j)

    def test_conditional_gv_and_source_scaling(self):
        amplitude = self.build(('pi+', 'pi-'))
        scaled_gv = self.build(('pi+', 'pi-'), g_v=2*self.fit.even_parameters['G_V_GeV'])
        scaled_source = build_photon_pp_amplitude(1.5, 2*self.axial, 2*self.pseudoscalar,
                     ('pi+', 'pi-'), 100., fit_parameters=self.fit)
        self.assertAlmostEqual(abs(scaled_gv(.6,.7)-2*amplitude(.6,.7)), 0., places=13)
        self.assertAlmostEqual(abs(scaled_source(.6,.7)-2*amplitude(.6,.7)), 0., places=13)

    def test_charged_c4_contributes(self):
        coefficients = dict(self.fit.odd_coefficients)
        coefficients['c4'] += 1.
        altered = replace(self.fit, odd_coefficients=coefficients)
        baseline = self.build(('K+', 'K-'))
        changed = build_photon_pp_amplitude(1.5, self.axial, self.pseudoscalar,
                     ('K+', 'K-'), 100., fit_parameters=altered)
        self.assertGreater(abs(changed(1.2,.6)-baseline(1.2,.6)), 1e-8)
        neutral = self.build(('K0', 'K0bar'))
        neutral_changed = build_photon_pp_amplitude(1.5, self.axial, self.pseudoscalar,
                     ('K0', 'K0bar'), 100., fit_parameters=altered)
        self.assertEqual(neutral.terms, neutral_changed.terms)


if __name__ == '__main__':
    unittest.main()
