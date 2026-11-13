"""Independent numerical checks of the finite-mixing mass interface."""
import unittest

import numpy as np
import sympy as sp

from calculate_alp_meson_masses import calculate_masses
from derive_alp_lo_neutral_blocks import closed_form_blocks
from solve_alp_meson_spectrum import solve_spectrum_numerical


class TestMasses(unittest.TestCase):
    def setUp(self):
        self.zero = np.zeros((3, 3))
        self.scalar_zero = np.diag([0.003, 0.005, 0.095])
        self.scalar_second = np.diag([-0.001, 0.002, -0.003])
        self.pseudoscalar_first = np.array([[0.001, 0, 0],
                                            [0, -0.002, 0.002 + 0.001j],
                                            [0, 0.002 - 0.001j, 0.003]])
        self.kwargs = dict(f_pi=0.092, B_0=2.5, m_0_squared=0.75)

    def masses(self, f_a=0.2, m_a0=0.8, theta=0.3, pseudoscalar=None):
        if pseudoscalar is None:
            pseudoscalar = self.pseudoscalar_first
        return calculate_masses(f_a, m_a0, [self.scalar_zero, self.zero, self.scalar_second],
                                [self.zero, pseudoscalar], [0., theta, 2.], **self.kwargs)

    def test_decoupled_and_scalar_second_derivative(self):
        masses = self.masses(theta=0., pseudoscalar=self.zero)
        correction = -(self.kwargs['f_pi'] / 0.2)**2 * self.kwargs['B_0'] * np.trace(self.scalar_second)
        self.assertAlmostEqual(masses['a']**2, 0.8**2 + correction)
        self.assertAlmostEqual(masses['pi+']**2, 2.5 * (0.003 + 0.005))
        self.assertAlmostEqual(masses['K0_1']**2, 2.5 * (0.005 + 0.095))
        self.assertAlmostEqual(masses['K0_2'], masses['K0_1'])
        self.assertEqual(masses['K+'], masses['K-'])
        self.assertEqual(masses['pi+'], masses['pi-'])

    def test_large_two_level_mixing(self):
        """Compare to the exact two-field eigenvalues in the isospin limit."""
        scalar_zero = np.diag([0.004, 0.004, 0.095])
        pseudoscalar = np.diag([0.003, -0.003, 0.])
        f_a, bare_alp_mass = 0.1, 0.15
        result = calculate_masses(f_a, bare_alp_mass,
                                  [scalar_zero, self.zero, self.zero],
                                  [self.zero, pseudoscalar], [0., 0.], **self.kwargs)
        pion_squared = 2*self.kwargs['B_0']*0.004
        mixing = -self.kwargs['f_pi']/f_a * self.kwargs['B_0']*0.006
        expected = [(pion_squared + bare_alp_mass**2 + sign *
                     np.sqrt((pion_squared - bare_alp_mass**2)**2 + 4*mixing**2))/2
                    for sign in (-1, 1)]
        np.testing.assert_allclose(sorted([result['pi0']**2, result['a']**2]), expected)

    def test_complex_source_against_generalized_pole_solver(self):
        """Use independent legacy blocks and the generalized eigenvalue solver."""
        b0 = self.kwargs['B_0']
        scalar = sp.Matrix(self.scalar_zero)
        pseudoscalar = sp.Matrix(self.pseudoscalar_first)
        # P_1 = {c^P,S_0}; solve elementwise for c^P.
        coupling = sp.Matrix(3, 3, lambda i, j: pseudoscalar[i, j] / (scalar[i, i] + scalar[j, j]))
        kinetic, mass_squared = closed_form_blocks(
            sp.zeros(3), coupling, 2*b0*scalar, sp.Float(0.3),
            sp.Float(0.092), sp.Float(0.2), sp.Float(0.75), sp.Float(0.64))
        mass_squared[5, 5] -= (0.092/0.2)**2 * b0 * np.trace(self.scalar_second)
        spectrum = solve_spectrum_numerical(np.array(kinetic, dtype=float),
                                            np.array(mass_squared, dtype=float))
        actual = self.masses()
        neutral = [actual[label]**2 for label in ('pi0', 'eta', 'eta_prime', 'K0_1', 'K0_2', 'a')]
        np.testing.assert_allclose(sorted(neutral), spectrum.pole_squared_masses, rtol=1e-12, atol=1e-13)
        self.assertAlmostEqual(actual['K0_2']**2, b0*(0.005+0.095))

    def test_invalid_vacuum_charge_and_negative_pole(self):
        violating = self.pseudoscalar_first.copy()
        violating[0, 1] = violating[1, 0] = 0.001
        with self.assertRaisesRegex(ValueError, 'electric charge'):
            self.masses(pseudoscalar=violating)
        with self.assertRaisesRegex(ValueError, 'negative squared mass'):
            self.masses(f_a=0.01, m_a0=0.01, theta=0., pseudoscalar=100*self.pseudoscalar_first)
        with self.assertRaisesRegex(ValueError, 'aligned vacuum'):
            calculate_masses(1., 1., [self.scalar_zero, self.zero, self.zero],
                             [np.eye(3), self.zero], [0., 0.], **self.kwargs)

    def test_irrelevant_higher_derivatives(self):
        reference = self.masses()
        changed = calculate_masses(0.2, 0.8, [self.scalar_zero, 9*self.zero, self.scalar_second, np.eye(3)],
                                   [self.zero, self.pseudoscalar_first, np.eye(3)],
                                   [0., 0.3, -100.], **self.kwargs)
        self.assertEqual(reference, changed)


if __name__ == '__main__':
    unittest.main()
