"""Tests of the double-trace singlet-source option y_V of the joint vector likelihood.

Run: .venv/bin/python -m unittest scripts.test_joint_vector_double_trace_singlet_source_option
Cheap quadrature orders (1, 3, 4); points are the best accepted endpoint
basin_probes/even_shared_widths extended by y_V in GeV^-1.
"""
import sys
import unittest
from pathlib import Path

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import derive_double_trace_vector_mass_and_kkpi as derivation  # noqa: E402
import fit_babar_kkpi_global_vector_sector as joint_fit  # noqa: E402
import joint_vector_double_trace_singlet_source_option as option  # noqa: E402
import refit_joint_vector_sensitivity as sensitivity  # noqa: E402
from fit_babar_kkpi_g_coefficients import CHANNELS  # noqa: E402
from implementation import vector_sector  # noqa: E402

ORDERS = (1, 3, 4)
TEST_COEFFICIENT = .013
BEST = np.asarray(sensitivity.accepted_context()['accepted'][0]['endpoint'], dtype=float)


def ideal_residues(model, values):
    """(R_rho, R_omega_I, R_phi_I) in GeV from vector_sector.source_residues with the ideal generators."""
    point = model.point+joint_fit.EVEN_SCALES*np.asarray(values)[:6]
    generators = [np.array(matrix, dtype=float) for matrix in derivation.ideal_generators().values()]
    spurion = np.diag([model.assumptions['xl_gev2']]*2+[model.assumptions['xs_gev2']])
    return vector_sector.source_residues(generators, spurion, point[4], point[5], values[option.COORDINATE_INDEX])


def ideal_propagators(model, state, channel, values):
    """rho and ideal-basis omega/phi propagators at q, assembled independently of the option code."""
    point, sector, mk, _, totals = state
    reference = dict(model.preliminary, m_rho_GeV=point[0], m_Kstar_GeV=mk,
                     Gamma_total_omega_GeV=totals['omega'], Gamma_total_phi_GeV=totals['phi'])
    index = CHANNELS.index(channel)
    beta = (values[26+index], values[28+index])
    q = model.phase[channel]['q']
    light_basis = np.diag([1, 1, 0])/np.sqrt(2)
    strange_basis = np.diag([0, 0, 1])
    rotation = np.array([[np.trace(generator@basis).real for generator in sector.generators[1:]]
                         for basis in (light_basis, strange_basis)])
    poles = [1/(sector.masses[n]**2-q-1j*sector.masses[n]*totals[name]) for n, name in ((1, 'omega'), (2, 'phi'))]

    def excited(name):
        ground, with_excitations = joint_fit.reference_propagator_pair(name, q, beta, reference)
        return with_excitations-ground
    light_light = rotation[0, 0]**2*poles[0]+rotation[0, 1]**2*poles[1]+excited('omega')
    strange_strange = rotation[1, 0]**2*poles[0]+rotation[1, 1]**2*poles[1]+excited('phi')
    light_strange = rotation[0, 0]*rotation[1, 0]*poles[0]+rotation[0, 1]*rotation[1, 1]*poles[1]
    rho = joint_fit.running_propagator('rho', q, beta, point[0], totals['rho0'], reference)
    return rho, light_light, strange_strange, light_strange


class DefaultReproducesBaseline(unittest.TestCase):
    """Option off, and option on at y_V = 0, give the baseline residual and Jacobian bitwise."""

    @classmethod
    def setUpClass(cls):
        cls.baseline = joint_fit.GlobalLikelihood(*ORDERS, bin_block_size=4)
        cls.baseline_residual = cls.baseline.residual(BEST)
        cls.baseline_jacobian = cls.baseline.jacobian(BEST)

    def test_option_off_bitwise(self):
        model = option.DoubleTraceSingletSourceLikelihood(*ORDERS, bin_block_size=4)
        self.assertEqual(len(model.parameter_names), 32)
        self.assertTrue(np.array_equal(model.residual(BEST), self.baseline_residual))
        self.assertTrue(np.array_equal(model.jacobian(BEST), self.baseline_jacobian))

    def test_option_on_at_zero_bitwise(self):
        model = option.DoubleTraceSingletSourceLikelihood(*ORDERS, bin_block_size=4, double_trace_singlet_source=True)
        values = model.extend(BEST, 0.)
        self.assertTrue(np.array_equal(model.residual(values), self.baseline_residual))
        jacobian = model.jacobian(values)
        self.assertEqual(jacobian.shape, (75, 33))
        self.assertTrue(np.array_equal(jacobian[:, :32], self.baseline_jacobian))

    def test_option_on_at_zero_full_grids_bitwise(self):
        baseline = joint_fit.GlobalLikelihood(*ORDERS)
        model = option.DoubleTraceSingletSourceLikelihood(*ORDERS, double_trace_singlet_source=True)
        values = model.extend(BEST, 0.)
        self.assertTrue(np.array_equal(model.residual(values), baseline.residual(BEST)))
        state = model.unpack(values)
        self.assertTrue(np.array_equal(model.babar_blocks(values, state, derivatives=True),
                                       baseline.babar_blocks(BEST, baseline.unpack(BEST), derivatives=True)))


class OnPathMatchesDerivation(unittest.TestCase):
    """At y_V != 0 the residues and the BaBar insertions follow the 4.1 derivation."""

    @classmethod
    def setUpClass(cls):
        cls.model = option.DoubleTraceSingletSourceLikelihood(*ORDERS, double_trace_singlet_source=True)
        cls.values = cls.model.extend(BEST, TEST_COEFFICIENT)
        cls.residual = cls.model.residual(cls.values)
        cls.state = cls.model.unpack(cls.values)
        cls.cached = dict(cls.model.babar_cached)

    def test_isoscalar_factors_equal_symbolic_replacements(self):
        factors = derivation.kkpi_source_factors()
        symbols = sp.symbols('F_l F_s', real=True)+(derivation.y_V, derivation.x_l, derivation.x_s)
        light_expression = sp.expand(factors['shifted_light'].subs(
            derivation.F_V, symbols[0]-2*derivation.z_V*derivation.x_l))
        strange_expression = sp.expand(factors['shifted_strange'].subs(
            derivation.F_V, symbols[1]-2*derivation.z_V*derivation.x_s))
        self.assertNotIn(derivation.z_V, light_expression.free_symbols | strange_expression.free_symbols)
        light_function = sp.lambdify(symbols, light_expression)
        strange_function = sp.lambdify(symbols, strange_expression)
        for arguments in ((.156, -.154, .013, .0366, .9437), (.2, .1, -.3, .02, .5)):
            light, strange = option.isoscalar_source_factors(*arguments)
            self.assertAlmostEqual(light, light_function(*arguments), places=14)
            self.assertAlmostEqual(strange, strange_function(*arguments), places=14)

    def test_symbolic_replacements_only_shift_by_y(self):
        factors = derivation.kkpi_source_factors()
        self.assertTrue(all(value == 0 for value in factors['differences'].values()))

    def test_physical_residues_are_rotated_ideal_residues(self):
        rho, omega, phi = ideal_residues(self.model, self.values)
        _, sector, _, _, _ = self.state
        basis = [np.array(matrix, dtype=float) for matrix in derivation.ideal_generators().values()]
        for generator, residue in zip(sector.generators, sector.residues):
            expected = sum(np.trace(generator@matrix).real*value for matrix, value in zip(basis, (rho, omega, phi)))
            self.assertAlmostEqual(residue, expected, places=13)

    def test_ideal_residues_equal_symbolic_residues(self):
        symbolic = derivation.photon_residues()
        point = self.model.point+joint_fit.EVEN_SCALES*self.values[:6]
        substitutions = {derivation.F_V: point[4], derivation.z_V: point[5],
                         derivation.y_V: TEST_COEFFICIENT, derivation.x_l: self.model.assumptions['xl_gev2'],
                         derivation.x_s: self.model.assumptions['xs_gev2']}
        for name, value in zip(('rho', 'omega_I', 'phi_I'), ideal_residues(self.model, self.values)):
            self.assertAlmostEqual(float(symbolic[name].subs(substitutions)), value, places=13)

    def test_babar_insertions_are_propagator_times_residues(self):
        rho_residue, omega_residue, phi_residue = ideal_residues(self.model, self.values)
        for channel in CHANNELS:
            rho, light_light, strange_strange, light_strange = ideal_propagators(
                self.model, self.state, channel, self.values)
            _, _, _, _, source_rho, source_light, source_strange = self.cached[channel]
            np.testing.assert_allclose(source_rho/np.sqrt(2), rho*rho_residue, rtol=1e-12, atol=0)
            np.testing.assert_allclose(source_light/(3*np.sqrt(2)),
                                       light_light*omega_residue+light_strange*phi_residue, rtol=1e-11, atol=1e-14)
            np.testing.assert_allclose(-source_strange/3,
                                       light_strange*omega_residue+strange_strange*phi_residue,
                                       rtol=1e-11, atol=1e-14)

    def test_isovector_insertion_unchanged(self):
        zero = option.DoubleTraceSingletSourceLikelihood(*ORDERS, double_trace_singlet_source=True)
        zero.residual(zero.extend(BEST, 0.))
        for channel in CHANNELS:
            for position in (0, 1, 2, 3, 4):
                self.assertTrue(np.array_equal(zero.babar_cached[channel][position], self.cached[channel][position]))
            self.assertFalse(np.array_equal(zero.babar_cached[channel][5], self.cached[channel][5]))

    def test_only_expected_rows_change(self):
        baseline = joint_fit.GlobalLikelihood(*ORDERS).residual(BEST)
        exactly_unchanged = [*range(12, 16), *range(19, 28), 73, 74]
        self.assertTrue(np.array_equal(self.residual[exactly_unchanged], baseline[exactly_unchanged]))
        # rho leptonic width: <T_rho> = 0, so y_V drops out up to rounding of the trace.
        self.assertLess(abs(self.residual[16]-baseline[16]), 1e-12)
        changed = [*range(0, 12), 17, 18, *range(28, 73)]
        self.assertTrue(np.all(self.residual[changed] != baseline[changed]))

    def test_cache_follows_y(self):
        model = option.DoubleTraceSingletSourceLikelihood(*ORDERS, double_trace_singlet_source=True)
        model.residual(model.extend(BEST, 0.))
        self.assertTrue(np.array_equal(model.residual(self.values), self.residual))


class CoordinateAndDomain(unittest.TestCase):
    """Names, bounds, start, length checks, penalty and refused combinations."""

    @classmethod
    def setUpClass(cls):
        cls.model = option.DoubleTraceSingletSourceLikelihood(*ORDERS, bin_block_size=4,
                                                              double_trace_singlet_source=True)

    def test_names_bounds_start(self):
        self.assertEqual(self.model.parameter_names[32:], option.ADDED_COORDINATE_NAMES)
        lower, upper = self.model.parameter_bounds()
        baseline_lower, baseline_upper = joint_fit.parameter_bounds(self.model)
        self.assertTrue(np.array_equal(lower[:32], baseline_lower) and np.array_equal(upper[:32], baseline_upper))
        self.assertEqual((lower[32], upper[32]), (-.5, .5))
        self.assertEqual(self.model.initial()[32], 0.)

    def test_wrong_length_raises(self):
        with self.assertRaises(ValueError):
            self.model.residual(BEST)

    def test_outside_bounds_gives_penalty(self):
        values = self.model.extend(BEST, .6)
        residual = self.model.residual(values)
        self.assertTrue(np.all(residual == 1e8+np.linalg.norm(values)))

    def test_extend_needs_option(self):
        with self.assertRaises(ValueError):
            option.DoubleTraceSingletSourceLikelihood(*ORDERS, bin_block_size=4).extend(BEST)

    def test_incompatible_options_refused(self):
        for options in (dict(regulator_gev=2.), dict(source_hypothesis='universal'),
                        dict(pole_shift=dict(state='rho(1450)', quantity='mass', shift_gev=.01))):
            with self.assertRaises(ValueError):
                option.DoubleTraceSingletSourceLikelihood(*ORDERS, bin_block_size=4,
                                                          double_trace_singlet_source=True, **options)

    def test_jacobian_column_and_steps(self):
        values = self.model.extend(BEST, TEST_COEFFICIENT)
        jacobian = self.model.jacobian(values)
        self.assertTrue(np.all(np.isfinite(jacobian)))
        step = 1e-7
        steps = np.sqrt(np.finfo(float).eps)*np.maximum(1., abs(values))
        steps[32] = step
        shifted = values.copy()
        shifted[32] += step
        expected = (self.model.residual(shifted)-self.model.residual(values))/step
        self.assertTrue(np.array_equal(self.model.jacobian(values, steps)[:, 32], expected))
        self.assertTrue(np.all(jacobian[19:28, 32] == 0))


if __name__ == '__main__':
    unittest.main()
