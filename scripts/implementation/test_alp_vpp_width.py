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


if __name__ == '__main__':
    unittest.main()
