"""Alda normalization and mass-source checks for the strong VPP calculation."""
from dataclasses import replace
import math
import unittest
from unittest.mock import patch

import numpy as np

from . import alp_vpp_width as implementation
from .alp_vv_width import read_preliminary_fit_parameters
from alpaca_chiral import vector_two_pseudoscalar as library


class VPPConventionTests(unittest.TestCase):
    """Check generator norms and a pure d2 amplitude against the trace formula."""

    def test_alda_generator_norms(self):
        for generators, _ in (implementation._meson_inputs(implementation.EPSILON_K),
                              library._meson_inputs()):
            for name, generator in generators.items():
                with self.subTest(name=name):
                    self.assertAlmostEqual(np.trace(generator.conj().T @ generator).real, 2.)

    def test_crossed_d2_mass_source(self):
        fitted = read_preliminary_fit_parameters()
        fitted = replace(fitted, odd_coefficients=dict(d1=0., d2=1., d3=0., d4=0.),
                         shared_total_widths={**fitted.shared_total_widths, 'rhoplus': .1491})
        axial = np.diag([1., 0., 0.])
        zero = np.zeros((3, 3))
        amplitude = implementation.build_vpp_amplitude(
            1.5, axial, zero, ('omega', 'pi+', 'pi-'), 1000.,
            g_v=.065, fit_parameters=fitted)
        generators, _ = implementation._meson_inputs(implementation.EPSILON_K)
        first, second = (generators[name].conj().T for name in ('pi+', 'pi-'))
        external = fitted.vector_generators['omega']
        s, t = .25, .95
        u = 1.5**2+sum(m*m for m in amplitude.masses)-s-t
        expected = 0j
        mass = fitted.even_parameters['m_rho_GeV']
        width = fitted.shared_total_widths['rhoplus']
        for internal in implementation._real_components(0, 1):
            for emitted, source, virtuality, sign in ((first, second, t, 1),
                                                      (second, first, u, -1)):
                # i chi_- = +{lambda,chi}/f_pi; N(d2)=-8 Tr({TV,TR}Y).
                mass_trace = np.trace((external@internal+internal@external)
                                     @ (source@fitted.quark_mass_matrix
                                        + fitted.quark_mass_matrix@source))/.093
                expected += (sign*math.sqrt(2)*.065/(.093*1000.*.78266)
                             *np.trace(internal@(emitted@axial-axial@emitted))
                             *(-8*mass_trace)/(virtuality-mass**2+1j*mass*width))
        self.assertGreater(abs(expected), 0.)
        np.testing.assert_allclose(amplitude(s, t)[0], expected, rtol=1e-13)
        library_fit = dict(vector_generators=fitted.vector_generators,
                           vector_masses=fitted.vector_masses,
                           quark_mass_matrix=fitted.quark_mass_matrix,
                           odd_coefficients=fitted.odd_coefficients,
                           even_parameters={**fitted.even_parameters, 'G_V_GeV': .065})
        with patch.object(library, '_read_vector_inputs', return_value=library_fit), \
             patch.object(library, '_no_mixing_couplings', return_value=(axial, zero, zero)), \
             patch.object(library, 'cqV_no_mixing', return_value=zero):
            library_amplitude = library.vector_two_pseudoscalar_amplitude(
                1.5, 1000., None, 'omega', 'pi+', 'pi-')
        np.testing.assert_allclose(library_amplitude(s, t)[0], expected, rtol=1e-13)


class VPPNeutralKaonTests(unittest.TestCase):
    """Compare physical-state coefficients with explicit outgoing flavor bras."""

    @classmethod
    def setUpClass(cls):
        cls.fitted = read_preliminary_fit_parameters()
        cls.epsilon_k = .13+.21j
        normalization = math.sqrt(2*(1+abs(cls.epsilon_k)**2))
        cls.p_mixing = (1+cls.epsilon_k)/normalization
        cls.q_mixing = (1-cls.epsilon_k)/normalization
        cls.axial = np.diag([.2, .4, .8]).astype(complex)
        cls.axial[1, 2] = .31+.17j
        cls.axial[2, 1] = cls.axial[1, 2].conjugate()
        cls.pseudoscalar = np.diag([.9, -.3, .5]).astype(complex)
        cls.pseudoscalar[1, 2] = -.12+.28j
        cls.pseudoscalar[2, 1] = cls.pseudoscalar[1, 2].conjugate()
        cls.vector_source = np.zeros((3, 3), complex)
        cls.vector_source[1, 2] = .7+.4j
        cls.vector_source[2, 1] = .7-.4j

    def amplitude(self, vector, first, second, axial=None, pseudoscalar=None):
        return implementation.build_vpp_amplitude(
            2.5, self.axial if axial is None else axial,
            self.pseudoscalar if pseudoscalar is None else pseudoscalar,
            (vector, first, second), 1000., g_v=.065,
            c_q_VS=self.vector_source, fit_parameters=self.fitted,
            epsilon_k=self.epsilon_k)

    def flavor_weights(self, name):
        if name not in ('KL', 'KS'):
            return ((name, 1.),)
        return (('K0', self.p_mixing.conjugate()),
                ('K0bar', (-1 if name == 'KL' else 1)*self.q_mixing.conjugate()))

    def test_complex_mixing_projects_all_coefficients(self):
        for vector in ('omega', 'phi'):
            for pair in (('KL', 'KL'), ('KL', 'KS'), ('KS', 'KS'),
                         ('pi0', 'KL'), ('KS', 'eta')):
                physical = self.amplitude(vector, *pair)
                flavor_amplitudes = [
                    (first_weight*second_weight,
                     self.amplitude(vector, first_flavor, second_flavor))
                    for first_flavor, first_weight in self.flavor_weights(pair[0])
                    for second_flavor, second_weight in self.flavor_weights(pair[1])]
                for s, t in ((1.2, 2.0), (1.4, 2.2), (1.7, 1.9)):
                    with self.subTest(vector=vector, pair=pair, s=s, t=t):
                        expected = sum((weight*np.asarray(amplitude(s, t))
                                        for weight, amplitude in flavor_amplitudes),
                                       np.zeros(3, complex))
                        np.testing.assert_allclose(physical(s, t), expected,
                                                   rtol=3e-13, atol=1e-16)

    def test_odd_equal_kaons_vanish_and_mixed_kaons_survive(self):
        for vector in ('omega', 'phi'):
            for pair in (('KL', 'KL'), ('KS', 'KS'), ('K0', 'K0'), ('K0bar', 'K0bar')):
                with self.subTest(vector=vector, pair=pair):
                    self.assertEqual(self.amplitude(vector, *pair).terms, ())
            mixed = self.amplitude(vector, 'KL', 'KS')
            flavor = self.amplitude(vector, 'K0', 'K0bar')
            for s, t in ((1.2, 2.), (1.4, 2.2)):
                self.assertGreater(abs(mixed(s, t)[0]), 1e-8)
                np.testing.assert_allclose(mixed(s, t)[0],
                    2*self.p_mixing.conjugate()*self.q_mixing.conjugate()*flavor(s, t)[0],
                    rtol=3e-13)
                reverse = self.amplitude(vector, 'KS', 'KL')
                u = mixed.alp_mass**2+sum(m*m for m in mixed.masses)-s-t
                np.testing.assert_allclose(reverse(s, u)[0], -mixed(s, t)[0], rtol=3e-13)
                self.assertAlmostEqual(reverse.squared(s, u), mixed.squared(s, t), places=14)

    def test_forbidden_equal_kaon_width_skips_quadrature(self):
        with patch.object(implementation, 'quad', side_effect=AssertionError('quadrature called')):
            for vector in ('omega', 'phi'):
                width = implementation.alp_vpp_width(
                    2.5, self.axial, self.pseudoscalar, (vector, 'KL', 'KL'),
                    1000., g_v=.065, c_q_VS=self.vector_source,
                    fit_parameters=self.fitted, epsilon_k=self.epsilon_k)
                self.assertEqual(width, 0.)

    def test_each_odd_source_has_the_neutral_pair_selection_rule(self):
        zero = np.zeros((3, 3))
        for vector in ('omega', 'phi'):
            for axial, pseudoscalar in ((self.axial, zero), (zero, self.pseudoscalar)):
                with self.subTest(vector=vector, axial=bool(np.any(axial))):
                    self.assertEqual(self.amplitude(vector, 'KL', 'KL', axial, pseudoscalar).terms, ())
                    self.assertGreater(abs(self.amplitude(vector, 'KL', 'KS', axial, pseudoscalar)(1.2, 2.)[0]), 1e-8)


if __name__ == '__main__':
    unittest.main()
