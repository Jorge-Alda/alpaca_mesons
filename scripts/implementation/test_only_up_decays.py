"""Cross-check only-up inputs and finite-mixing decay adapters."""
import contextlib
import io
import math
import sys
import unittest
from pathlib import Path

import numpy as np
from scipy.integrate import quad

from .only_up_inputs import (B_ZERO, DOWN_MASS, M1_SQUARED, PION_DECAY_CONSTANT,
                             PION_MASS, STRANGE_MASS, UP_MASS, only_up_inputs)
from .only_up_decays import (photon_two_pseudoscalar, photon_two_pseudoscalar_width, spectrum,
                             three_pseudoscalar, three_pseudoscalar_width,
                             two_photon_coefficient, two_photon_width)
from .alp_photon_pp_wzw_finite_mixing import alternating_trace, external_legs, wzw_coefficient

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from derive_wzw_external_currents import derive_neutral_pion_two_photons


class OnlyUpTests(unittest.TestCase):
    def test_phase_derivatives_and_inputs(self):
        data = only_up_inputs(.01)
        self.assertAlmostEqual(data['f_a'], math.sqrt(2)*UP_MASS/.01)
        self.assertAlmostEqual(B_ZERO*(UP_MASS+DOWN_MASS), PION_MASS**2)
        self.assertAlmostEqual(M1_SQUARED+2*B_ZERO*(UP_MASS+DOWN_MASS+STRANGE_MASS)/3,
                               .95778**2)
        self.assertTrue(np.array_equal(data['hat_P'][1], np.diag([-UP_MASS, 0., 0.])))
        self.assertTrue(np.array_equal(data['hat_P'][3], -data['hat_P'][1]))
        self.assertTrue(np.array_equal(data['hat_S'][4], -data['hat_S'][2]))

    def test_poles_against_paper_matrix(self):
        data = only_up_inputs(.01)
        ratio = PION_DECAY_CONSTANT/data['f_a']
        u, d, s = UP_MASS, DOWN_MASS, STRANGE_MASS
        root3 = math.sqrt(3)
        matrix = B_ZERO*np.array([
            [ratio**2*u, -ratio*u, -ratio*u/root3, -ratio*u*math.sqrt(2/3)],
            [-ratio*u, u+d, (u-d)/root3, math.sqrt(2/3)*(u-d)],
            [-ratio*u/root3, (u-d)/root3, (u+d+4*s)/3, math.sqrt(2)*(u+d-2*s)/3],
            [-ratio*u*math.sqrt(2/3), math.sqrt(2/3)*(u-d), math.sqrt(2)*(u+d-2*s)/3,
             M1_SQUARED/B_ZERO+2*(u+d+s)/3]])
        masses = spectrum(data)
        found = sorted(masses[key]**2 for key in ('a', 'pi0', 'eta', 'eta_prime'))
        eigenvalues, eigenvectors = np.linalg.eigh(matrix)
        np.testing.assert_allclose(found, eigenvalues, rtol=1e-12, atol=1e-14)
        column = eigenvectors[:, np.argmin(abs(eigenvalues-masses['a']**2))]
        if column[0] < 0:
            column = -column
        generators = [np.diag([1., -1., 0.]), np.diag([1., 1., -2.])/root3,
                      np.eye(3)*math.sqrt(2/3)]
        expected_leg = -sum((component*generator for component, generator in
                             zip(column[1:], generators)), np.zeros((3, 3)))
        legs, _ = external_legs('pi+', 'pi-', data['m_a0'], data['f_a'],
                                 data['hat_S'], data['hat_P'], data['hat_Theta'],
                                 f_pi=data['f_pi'], B_0=B_ZERO, m_0_squared=M1_SQUARED)
        np.testing.assert_allclose(legs[0], expected_leg, rtol=1e-11, atol=1e-13)

    def test_common_wzw_orientation(self):
        with contextlib.redirect_stdout(io.StringIO()):
            orientation = int(derive_neutral_pion_two_photons())
        self.assertEqual(orientation, -1)
        data = only_up_inputs(.01)
        args = (data['m_a0'], data['f_a'], data['hat_S'], data['hat_P'], data['hat_Theta'])
        kwargs = {'f_pi': data['f_pi'], 'B_0': data['B_0'], 'm_0_squared': data['m_0_squared']}
        legs, _ = external_legs('pi+', 'pi-', *args, **kwargs)
        expected = -orientation*3*math.sqrt(4*math.pi/137)/(24*math.pi**2*data['f_pi']**3)*alternating_trace(*legs)
        self.assertAlmostEqual(photon_two_pseudoscalar(data, ('pi+', 'pi-')), expected)

    def test_photon_widths(self):
        data = only_up_inputs(.3)
        masses = spectrum(data)
        alpha = 1/137
        expected_two_photon = (alpha**2*masses['a']**3*two_photon_coefficient(data)**2
                               /(64*math.pi**3*data['f_a']**2))
        self.assertAlmostEqual(two_photon_width(data), expected_two_photon, places=20)
        self.assertEqual(photon_two_pseudoscalar_width(only_up_inputs(.01), ('pi+', 'pi-')), 0.)
        self.assertEqual(photon_two_pseudoscalar_width(data, ('pi0', 'pi0')), 0.)
        first_mass, second_mass, parent = masses['pi+'], masses['pi-'], masses['a']
        coefficient = photon_two_pseudoscalar(data, ('pi+', 'pi-'))
        def integrate_t(s):
            kallen = (s-(first_mass+second_mass)**2)*(s-(first_mass-second_mass)**2)
            center = second_mass**2+(parent**2-s)*(s+second_mass**2-first_mass**2)/(2*s)
            half = (parent**2-s)*math.sqrt(max(0., kallen))/(2*s)
            return quad(lambda t: s*(t-(center-half))*(center+half-t)/4,
                        center-half, center+half, epsabs=0., epsrel=1e-10)[0]
        integral = quad(integrate_t, (first_mass+second_mass)**2, parent**2,
                        epsabs=0., epsrel=1e-10)[0]
        expected = abs(coefficient)**2*integral/(256*math.pi**3*parent**3)
        actual = photon_two_pseudoscalar_width(data, ('pi+', 'pi-'))
        self.assertGreater(actual, 0.)
        self.assertAlmostEqual(actual/expected, 1., places=8)
        self.assertAlmostEqual(photon_two_pseudoscalar_width(data, ('pi-', 'pi+'))/actual, 1.)

    def test_open_three_pion_width(self):
        data = only_up_inputs(.3)
        mesons = ('pi0', 'pi+', 'pi-')
        masses = spectrum(data)
        self.assertGreater(masses['a'], sum(masses[meson] for meson in mesons))
        from calculate_alp_three_pseudoscalar_finite import decay_width
        direct = decay_width(data['f_a'], data['m_a0'], data['hat_S'], data['hat_P'],
                             data['hat_Theta'], mesons, f_pi=data['f_pi'],
                             B_0=data['B_0'], m_0_squared=data['m_0_squared'])
        self.assertGreater(direct, 0.)
        self.assertAlmostEqual(three_pseudoscalar_width(data, mesons), direct, places=15)

    def test_photon_contractions_and_width_threshold(self):
        data = only_up_inputs(.01)
        args = (data['m_a0'], data['f_a'], data['hat_S'], data['hat_P'], data['hat_Theta'])
        kwargs = {'f_pi': data['f_pi'], 'B_0': data['B_0'], 'm_0_squared': data['m_0_squared']}
        legs, _ = external_legs('pi+', 'pi-', *args, **kwargs)
        expected = 3*data['f_a']/data['f_pi']*np.trace(np.diag([4/9, 1/9, 1/9])@legs[0])
        self.assertAlmostEqual(two_photon_coefficient(data), expected.real)
        self.assertAlmostEqual(photon_two_pseudoscalar(data, ('pi+', 'pi-')),
                               wzw_coefficient('pi+', 'pi-', *args, **kwargs, e=math.sqrt(4*math.pi/137)))
        self.assertTrue(np.isfinite(three_pseudoscalar(data, ('pi0', 'pi+', 'pi-'))(.1, .1)))
        self.assertEqual(three_pseudoscalar_width(data, ('pi0', 'pi+', 'pi-')), 0.)


if __name__ == '__main__':
    unittest.main()
