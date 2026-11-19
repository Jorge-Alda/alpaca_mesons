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
                             two_photon_coefficient, two_photon_width,
                             BARYON_MASSES, TOPOLOGICAL_SINGLET_BARYON,
                             TOPOLOGICAL_SINGLET_BARYON_NAIVE_UNCERTAINTY, baryon_pair_couplings,
                             baryon_pair_width, baryon_pair_width_from_couplings, baryon_pair_widths,
                             neutral_flavour_conserving_channels)
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


class RadialExchangeTests(unittest.TestCase):
    """sigma exchange in a -> 3P against the direct polar tree amplitude of
    scripts/derive_murayama_only_up_three_pseudoscalar_amplitude.py (same Lagrangian, same legs)."""

    KAPPA = .3
    RADIAL_MASS_SQUARED = 2.5  # GeV^2

    @classmethod
    def setUpClass(cls):
        from derive_murayama_only_up_three_pseudoscalar_amplitude import charge_model_parameters, spectrum
        cls.inputs = only_up_inputs(cls.KAPPA)
        f_a = cls.inputs['f_a']
        quartic = (cls.RADIAL_MASS_SQUARED - B_ZERO*UP_MASS*PION_DECAY_CONSTANT**2/f_a**2)/(2*f_a**2)
        cls.parameters = charge_model_parameters((1, 0, 0), (UP_MASS, DOWN_MASS, STRANGE_MASS), f_a, quartic,
                                                 B_ZERO, PION_DECAY_CONSTANT, M1_SQUARED)
        cls.data = spectrum(cls.parameters)

    def test_radial_mass(self):
        from derive_alp_three_pseudoscalar_radial_exchange import radial_mass_squared_trace_form
        self.assertAlmostEqual(radial_mass_squared_trace_form(self.parameters), self.RADIAL_MASS_SQUARED, places=12)

    def test_amplitude_against_direct_polar_tree(self):
        from derive_murayama_only_up_three_pseudoscalar_amplitude import direct_pieces, flavour_direction
        for mesons in (('pi0', 'pi+', 'pi-'), ('pi0', 'eta', 'eta_prime'), ('eta', 'pi0', 'pi0'),
                       ('pi0', 'K+', 'K-'), ('pi0', 'K0_1', 'K0_1')):
            total = three_pseudoscalar(self.inputs, mesons, self.RADIAL_MASS_SQUARED)
            amplitude = total.without_radial
            masses_squared = [amplitude.masses['a']**2] + [amplitude.masses[name]**2 for name in mesons]
            for mass_squared in masses_squared:  # same pole masses as the eleven-field polar Hessian
                self.assertLess(np.min(np.abs(self.data['masses_squared'] - mass_squared)), 1e-12*mass_squared)
            directions = [flavour_direction(matrix, alp) for matrix, alp in amplitude.external]
            s12, s23 = .31, .27
            s13 = sum(masses_squared) - s12 - s23
            direct = direct_pieces(self.parameters, self.data, directions, masses_squared, s12, s13)
            with self.subTest(mesons=mesons):
                scale = abs(direct['total'])
                self.assertLess(abs(total(s12, s23) - direct['total']), 1e-9*scale)
                self.assertLess(abs(total.radial(s12, s23) - direct['radial_exchange']), 1e-9*scale)
                self.assertLess(abs(three_pseudoscalar(self.inputs, mesons)(s12, s23) - amplitude(s12, s23)), 1e-15*scale)

    def test_width_defaults_and_limits(self):
        mesons = ('pi0', 'pi+', 'pi-')
        default = three_pseudoscalar_width(self.inputs, mesons)
        heavy = three_pseudoscalar_width(self.inputs, mesons, 1e8)
        self.assertAlmostEqual(heavy/default, 1., places=6)
        with_radial = three_pseudoscalar_width(self.inputs, mesons, self.RADIAL_MASS_SQUARED, epsrel=1e-8)
        self.assertGreater(abs(with_radial/default - 1), 1e-4)
        # Independent integration of the same coherent amplitude.
        from scipy.integrate import dblquad
        total = three_pseudoscalar(self.inputs, mesons, self.RADIAL_MASS_SQUARED)
        masses = spectrum(self.inputs)
        parent, (first, second, third) = masses['a'], [masses[name] for name in mesons]
        def limits(s, index):
            kallen = lambda x, y, z: (x - y - z)**2 - 4*y*z
            root = math.sqrt(max(0., kallen(s, first**2, second**2)*kallen(parent**2, s, third**2)))
            midpoint = second**2 + third**2 + (s + second**2 - first**2)*(parent**2 - s - third**2)/(2*s)
            return midpoint + (2*index - 1)*root/(2*s)
        integral = dblquad(lambda t, s: abs(total(s, t))**2, (first + second)**2, (parent - third)**2,
                           lambda s: limits(s, 0), lambda s: limits(s, 1), epsabs=0, epsrel=1e-9)[0]
        self.assertAlmostEqual(with_radial/(integral/(256*math.pi**3*parent**3)), 1., places=6)

    def test_radial_pole_in_dalitz_region_raises(self):
        with self.assertRaises(ValueError):
            three_pseudoscalar_width(self.inputs, ('pi0', 'pi+', 'pi-'), .1)


class BaryonPairWidthTests(unittest.TestCase):
    """a -> B_1 bar(B_2); reference widths in GeV from scripts/derive_murayama_baryon_decays.py
    before its numerical part was moved here (memory/alp/baryon/murayama-only-up-baryon-pair-widths.md)."""

    REFERENCE = {
        .9: {('p', 'p'): 0.01739913355664702, ('n', 'n'): 0.0028891438799978136},
        1.: {('p', 'p'): 0.02024736579170888, ('n', 'n'): 0.0035501019230788477,
             ('Lambda', 'Lambda'): 0.00016339717041817878,
             ('Lambda', 'Sigma0'): 0.0024840637080593437,
             ('Sigma0', 'Lambda'): 0.0024840637080593446},
    }

    def test_widths_against_derivation_values(self):
        for kappa, reference in self.REFERENCE.items():
            inputs = only_up_inputs(kappa)
            widths = baryon_pair_widths(inputs)
            self.assertEqual(set(widths), set(neutral_flavour_conserving_channels()))
            for pair, width in widths.items():
                with self.subTest(kappa=kappa, pair=pair):
                    expected = reference.get(pair, 0.)
                    self.assertAlmostEqual(width, expected, delta=1e-12*expected)
                    self.assertEqual(baryon_pair_width(inputs, pair), width)

    def test_below_threshold_is_zero(self):
        inputs = only_up_inputs(.3)
        alp_mass = spectrum(inputs)['a']
        self.assertLess(alp_mass, 2*BARYON_MASSES['p'])
        self.assertTrue(all(width == 0. for width in baryon_pair_widths(inputs).values()))
        self.assertEqual(baryon_pair_width_from_couplings(.1, .1, 1., 1., 2., .003), 0.)

    def test_input_dictionary_and_singlet_argument(self):
        inputs = only_up_inputs(1.)
        # kappa_u enters only through inputs['f_a']: overriding f_a reproduces kappa_u = 0.9.
        modified = dict(inputs, f_a=math.sqrt(2)*UP_MASS/.9)
        for pair, expected in self.REFERENCE[.9].items():
            self.assertAlmostEqual(baryon_pair_width(modified, pair), expected, delta=1e-12*expected)
        central = baryon_pair_width(inputs, ('Lambda', 'Lambda'))
        shifted = baryon_pair_width(inputs, ('Lambda', 'Lambda'), TOPOLOGICAL_SINGLET_BARYON
                                    - TOPOLOGICAL_SINGLET_BARYON_NAIVE_UNCERTAINTY)
        self.assertNotAlmostEqual(central/shifted, 1., places=2)
        # The Lambda Sigma0bar coupling has no singlet part.
        self.assertEqual(baryon_pair_width(inputs, ('Lambda', 'Sigma0'), .3),
                         baryon_pair_width(inputs, ('Lambda', 'Sigma0')))
        vector, axial = baryon_pair_couplings(inputs, ('p', 'p'))
        self.assertEqual(vector, 0)
        self.assertAlmostEqual(axial.real, -1.00e-3, delta=1e-5)


if __name__ == '__main__':
    unittest.main()
