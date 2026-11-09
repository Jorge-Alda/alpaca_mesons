"""Tests of the switchable isospin-breaking VPP vertex of the joint likelihood.

Run: .venv/bin/python -m unittest scripts.test_joint_vector_isospin_breaking_vpp_option

Expected factors are computed with the functions of
scripts/scan_isospin_breaking_vpp_simultaneous.py and
scripts/derive_isospin_breaking_vpp_operators.py (SymPy traces), and the BaBar
substitution is compared with babar_amplitudes of
scripts/derive_kstar_kpi_saturation_freedom.py. Reduced quadrature orders
(1, 3, 4) except for the bitwise default test, which uses the orders of the
earlier K* option test.
"""
import math
import sys
import unittest
from pathlib import Path

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import derive_isospin_breaking_vpp_operators as derivation  # noqa: E402
import derive_kstar_kpi_saturation_freedom as saturation  # noqa: E402
import fit_babar_kkpi_global_vector_sector as joint_fit  # noqa: E402
import joint_vector_isospin_breaking_vpp_option as vpp_option  # noqa: E402
import scan_isospin_breaking_vpp_simultaneous as scan  # noqa: E402
from test_joint_vector_likelihood_options import accepted_endpoint  # noqa: E402

SCAN_POINT = np.array(vpp_option.SCAN_POINT_GEV_MINUS2)
PHI_ROWS = {26, 27}
RADIATIVE_KSTAR_ROWS = {6, 7}  # K*0 -> K0 gamma, K*+ -> K+ gamma


def width_row(model, name):
    """Residual index of a total-width row."""
    return 19+model.total_order.index(name)


def scan_width_ratio(modes, coefficients, point):
    """Gamma(r)/Gamma(0) of one K* charge state from the scan script's traces and phase space."""
    traces = scan.channel_traces(coefficients, *point)
    return (sum(scan.PHASE_SPACE[m]*traces[m]**2 for m in modes)
            / sum(scan.PHASE_SPACE[m]*coefficients[m][0]**2 for m in modes))


def symbolic_phi_factors(mixing_angle, point):
    """T(r)/T(0) of phi -> K+K-, K0K0bar from the SymPy traces of the derivation script."""
    published = derivation.alda_inputs()['published_lambda_diagonals']['pi0']
    traces = derivation.effective_traces([sp.Float(x, 30) for x in published], sp.Float(mixing_angle, 30))
    numbers = dict(zip((scan.XU, scan.XD, scan.XS), derivation.vacuum_spurion_values()))
    values = dict(zip((scan.R1, scan.R2, scan.R3), point))
    zero = dict.fromkeys((scan.R1, scan.R2, scan.R3), 0)
    return {mode: float(traces[mode].subs(numbers).subs(values)/traces[mode].subs(numbers).subs(zero))
            for mode in vpp_option.PHI_MODES}


class DefaultsReproduceBaseline(unittest.TestCase):
    def test_option_off_bitwise_equal_residual_and_jacobian(self):
        values = accepted_endpoint()
        baseline = joint_fit.GlobalLikelihood(2, 16, 16, bin_block_size=4)
        for model in (vpp_option.IsospinBreakingVppLikelihood(2, 16, 16, bin_block_size=4),
                      vpp_option.IsospinBreakingVppLikelihood(
                          2, 16, 16, bin_block_size=4, isospin_breaking_vpp=False,
                          vpp_kstar_excitation_scaling='ground_only')):
            self.assertEqual(model.parameter_names, joint_fit.PARAMETER_NAMES)
            np.testing.assert_array_equal(model.initial(), baseline.initial())
            self.assertTrue(np.array_equal(baseline.residual(values), model.residual(values)))
            self.assertTrue(np.array_equal(baseline.jacobian(values), model.jacobian(values)))
            for expected, actual in zip(joint_fit.parameter_bounds(baseline), model.parameter_bounds()):
                np.testing.assert_array_equal(expected, actual)
        self.assertIs(joint_fit.running_propagator, vpp_option.RUNNING_PROPAGATOR)

    def test_option_on_at_zero_coefficients_bitwise_equal(self):
        values = accepted_endpoint()
        baseline = joint_fit.GlobalLikelihood(1, 3, 4, bin_block_size=4)
        for scaling in vpp_option.EXCITATION_SCALINGS:
            model = vpp_option.IsospinBreakingVppLikelihood(
                1, 3, 4, bin_block_size=4, isospin_breaking_vpp=True, vpp_kstar_excitation_scaling=scaling)
            point = model.extend(values)
            np.testing.assert_array_equal(model.initial()[:32], baseline.initial())
            np.testing.assert_array_equal(model.initial()[32:], np.zeros(3))
            self.assertTrue(np.array_equal(baseline.residual(values), model.residual(point)))
            jacobian = model.jacobian(point)
            self.assertEqual(jacobian.shape, (75, 35))
            self.assertTrue(np.array_equal(baseline.jacobian(values), jacobian[:, :32]))

    def test_coordinates_and_bounds(self):
        model = vpp_option.IsospinBreakingVppLikelihood(1, 3, 4, isospin_breaking_vpp=True)
        self.assertEqual(model.parameter_names[32:], vpp_option.ADDED_COORDINATE_NAMES)
        self.assertEqual(75-len(model.parameter_names), 40)
        lower, upper = model.parameter_bounds()
        np.testing.assert_array_equal(lower[32:], [-5., -5., -5.])
        np.testing.assert_array_equal(upper[32:], [5., 5., 5.])


class FactorsAgreeWithScanScript(unittest.TestCase):
    def test_kstar_and_rho_factors_equal_scan_trace_coefficients(self):
        coefficients = scan.trace_coefficients(derivation.vacuum_spurion_values())
        factors = vpp_option.vertex_factors(SCAN_POINT, derivation.vacuum_spurion_values(),
                                            derivation.alda_inputs()['published_lambda_diagonals']['pi0'])
        traces = scan.channel_traces(coefficients, *SCAN_POINT)
        for mode, factor in factors.items():
            self.assertAlmostEqual(factor, traces[mode]/coefficients[mode][0], places=12)

    def test_phi_factors_at_fixed_angle_equal_scan_trace_coefficients(self):
        angle = scan.gv_fit.fixed_parameters()[2]
        generator = (math.sin(angle)*np.diag([1., 1., 0.])/math.sqrt(2)
                     + math.cos(angle)*np.diag([0., 0., 1.]))
        coefficients = scan.trace_coefficients(derivation.vacuum_spurion_values())
        traces = scan.channel_traces(coefficients, *SCAN_POINT)
        factors = vpp_option.vertex_factors(SCAN_POINT, derivation.vacuum_spurion_values(),
                                            derivation.alda_inputs()['published_lambda_diagonals']['pi0'],
                                            generator)
        for mode in vpp_option.PHI_MODES:
            self.assertAlmostEqual(factors[mode], traces[mode]/coefficients[mode][0], places=12)

    def test_charge_conjugate_kstar_modes_have_equal_factors(self):
        model = vpp_option.IsospinBreakingVppLikelihood(1, 3, 4)
        spurion = np.diag(model.spurion_diagonal)
        unit = vpp_option.matrix_unit
        pion_zero = np.diag(model.pion_zero_diagonal)/math.sqrt(2)
        conjugates = {'kstar_charged_K+_pi0': (unit(0, 2), unit(2, 0), pion_zero),
                      'kstar_charged_K0_pi+': (unit(0, 2), unit(2, 1), unit(1, 0)),
                      'kstar_neutral_K+_pi-': (unit(1, 2), unit(2, 0), unit(0, 1))}
        factors = model.factors(SCAN_POINT)
        for mode, matrices in conjugates.items():
            ratio = (vpp_option.effective_trace(*matrices, spurion, SCAN_POINT)
                     / vpp_option.effective_trace(*matrices, spurion, np.zeros(3)))
            self.assertAlmostEqual(ratio, factors[mode], places=13)


class ScanPointChangesPredictedRows(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.values = accepted_endpoint()
        cls.reference = joint_fit.GlobalLikelihood(1, 3, 4)
        cls.baseline = cls.reference.residual(cls.values)
        cls.reference_state = cls.reference.unpack(cls.values)
        _, _, cls.reference_details = cls.reference.old_blocks(cls.values, cls.reference_state)
        cls.model = vpp_option.IsospinBreakingVppLikelihood(1, 3, 4, isospin_breaking_vpp=True)
        cls.point = cls.model.extend(cls.values, SCAN_POINT)
        cls.residual = cls.model.residual(cls.point)
        cls.state = cls.model.unpack(cls.point)
        _, _, cls.details = cls.model.old_blocks(cls.point, cls.state)
        cls.scan_coefficients = scan.trace_coefficients(derivation.vacuum_spurion_values())
        cls.scan_values = scan.observables(cls.scan_coefficients, *SCAN_POINT)

    def test_changed_rows(self):
        changed = set(np.flatnonzero(self.residual != self.baseline).tolist())
        expected = (RADIATIVE_KSTAR_ROWS | PHI_ROWS
                    | {width_row(self.model, name) for name in vpp_option.KSTAR_NAMES})
        self.assertTrue(expected <= changed)
        self.assertTrue(changed & set(range(28, 28+len(self.model.babar['kkpi0']['observed']))))
        self.assertTrue(changed & set(range(28+len(self.model.babar['kkpi0']['observed']), 73)))
        self.assertTrue(changed <= expected | set(range(28, 73)))

    def test_kstar_widths_scale_as_scan_prediction(self):
        totals, reference_totals = self.state[4], self.reference_state[4]
        ratios = {}
        for name in vpp_option.KSTAR_NAMES:
            modes = [mode for _, _, mode in vpp_option.KSTAR_MODES[name]]
            expected = scan_width_ratio(modes, self.scan_coefficients, SCAN_POINT)
            ratios[name] = totals[name]/reference_totals[name]
            self.assertAlmostEqual(ratios[name], expected, places=12)
            row = width_row(self.model, name)
            sigma = self.model.total_sigma[self.model.total_order.index(name)]
            self.assertAlmostEqual(self.residual[row]-self.baseline[row],
                                   reference_totals[name]*(expected-1)/sigma, places=9)
        # Direction: the charged/neutral width ratio rises by the scan's R.
        self.assertAlmostEqual(ratios['Kstarplus']/ratios['Kstar0'], self.scan_values['kstar_ratio'], places=12)
        self.assertGreater(ratios['Kstarplus']/ratios['Kstar0'], 1.)
        # Both K* couplings are reduced at this point (common-mode factor < 1).
        self.assertLess(ratios['Kstarplus'], 1.)
        self.assertLess(ratios['Kstar0'], 1.)

    def test_radiative_kstar_rows_follow_total_widths(self):
        for row, name in zip((6, 7), ('Kstar0', 'Kstarplus')):
            self.assertEqual(self.model.observations[row]['total_width_key'], name)
            expected = (self.reference_details['odd_predictions'][row]
                        * self.reference_state[4][name]/self.state[4][name])
            self.assertAlmostEqual(self.details['odd_predictions'][row]/expected, 1., places=12)

    def test_phi_rows_scale_as_derivation_traces(self):
        generator = np.real(self.state[1].generators[2])
        angle = math.atan2(math.sqrt(2)*generator[0, 0], generator[2, 2])
        np.testing.assert_allclose(generator, math.sin(angle)*np.diag([1., 1., 0.])/math.sqrt(2)
                                   + math.cos(angle)*np.diag([0., 0., 1.]), atol=1e-14)
        expected = symbolic_phi_factors(angle, SCAN_POINT)
        predictions = self.details['phi_branching_predictions']
        reference = self.reference_details['phi_branching_predictions']
        for mode, new, old in zip(vpp_option.PHI_MODES, predictions, reference):
            self.assertAlmostEqual(new/old, expected[mode]**2, places=12)
        change = (predictions[0]/predictions[1])/(reference[0]/reference[1])-1
        # Direction predicted by the scan (fixed provisional angle): -4.43 %.
        self.assertLess(change, 0.)
        self.assertLess(self.scan_values['phi_change'], 0.)
        self.assertAlmostEqual(change, self.scan_values['phi_change'], delta=.01)
        expected_residual = np.linalg.solve(self.model.phi_cholesky, np.array(predictions)-self.model.phi_target)
        np.testing.assert_allclose(self.residual[26:28], expected_residual, rtol=1e-13)

    def test_babar_terms_equal_per_mode_saturation_substitution(self):
        factors = self.model.factors(SCAN_POINT)
        widths = (self.state[4]['Kstarplus'], self.state[4]['Kstar0'])
        count = len(self.model.babar['kkpi0']['observed'])
        rows = {'kkpi0': slice(0, count), 'kskpi': slice(count, None)}
        kappas = {'kkpi0': (factors['kstar_charged_K+_pi0']**2, 1.),
                  'kskpi': (factors['kstar_charged_K0_pi+']**2, factors['kstar_neutral_K+_pi-']**2)}
        for scaling, record_scaling in (('ground_only', 'ground_only'), ('whole_propagator', 'whole')):
            model = vpp_option.IsospinBreakingVppLikelihood(1, 3, 4, isospin_breaking_vpp=True,
                                                            vpp_kstar_excitation_scaling=scaling)
            babar = model.residual(self.point)[28:73]
            for channel, (charged, neutral) in kappas.items():
                _, expected = saturation.babar_amplitudes(
                    self.reference, self.values, self.reference_state, charged, neutral,
                    excited_scaling=record_scaling, kstar_widths=widths)
                np.testing.assert_allclose(babar[rows[channel]], expected[rows[channel]], rtol=1e-11, atol=1e-11)
        self.assertIs(joint_fit.running_propagator, vpp_option.RUNNING_PROPAGATOR)

    def test_r3_alone_changes_only_phi_rows(self):
        point = self.model.extend(self.values, (0., 0., SCAN_POINT[2]))
        changed = set(np.flatnonzero(self.model.residual(point) != self.baseline).tolist())
        self.assertEqual(changed, PHI_ROWS)


class IsospinLimitRho(unittest.TestCase):
    def test_rho_plus_over_rho_zero_unchanged_in_isospin_limit(self):
        spurions = derivation.vacuum_spurion_values()
        for point in (SCAN_POINT, (1., -3., 2.), (-5., 5., -5.)):
            factors = vpp_option.vertex_factors(point, spurions, (1., -1., 0.))
            self.assertAlmostEqual(factors['rho+_pi+_pi0']/factors['rho0_pi+_pi-'], 1., places=14)
        # Same statement from the SymPy traces of the derivation script (exact, all r_i).
        traces = derivation.effective_traces((1, -1, 0), 0)
        self.assertEqual(sp.simplify(traces['rho+_pi+_pi0']+traces['rho0_pi+_pi-']), 0)

    def test_published_pi0_gives_small_rho_split(self):
        model = vpp_option.IsospinBreakingVppLikelihood(1, 3, 4)
        factors = model.factors(SCAN_POINT)
        coefficients = scan.trace_coefficients(derivation.vacuum_spurion_values())
        expected = scan.observables(coefficients, *SCAN_POINT)['rho_split']
        self.assertAlmostEqual(factors['rho+_pi+_pi0']/factors['rho0_pi+_pi-'], expected, places=12)
        self.assertLess(abs(expected-1), 2e-3)


class JacobianAtScanPoint(unittest.TestCase):
    def test_finite_jacobian_and_added_columns(self):
        values = accepted_endpoint()
        model = vpp_option.IsospinBreakingVppLikelihood(1, 3, 4, bin_block_size=4, isospin_breaking_vpp=True)
        point = model.extend(values, SCAN_POINT)
        jacobian = model.jacobian(point)
        self.assertEqual(jacobian.shape, (75, 35))
        self.assertTrue(np.all(np.isfinite(jacobian)))
        for column in (32, 33):
            self.assertTrue(np.any(jacobian[28:73, column] != 0))
            self.assertTrue(np.all(jacobian[[width_row(model, n) for n in vpp_option.KSTAR_NAMES], column] != 0))
        # r3 enters only the phi rows among the implemented ones.
        self.assertTrue(np.all(jacobian[26:28, 34] != 0))
        np.testing.assert_array_equal(np.delete(jacobian[:, 34], [26, 27]), 0.)


class InvalidInputs(unittest.TestCase):
    def test_constructor_options(self):
        with self.assertRaises(ValueError):
            vpp_option.IsospinBreakingVppLikelihood(1, 3, 4, isospin_breaking_vpp=True,
                                                    vpp_kstar_excitation_scaling='excited_only')
        with self.assertRaises(ValueError):
            vpp_option.IsospinBreakingVppLikelihood(1, 3, 4, vpp_kstar_excitation_scaling='whole_propagator')

    def test_coordinates(self):
        values = accepted_endpoint()
        off = vpp_option.IsospinBreakingVppLikelihood(1, 3, 4)
        on = vpp_option.IsospinBreakingVppLikelihood(1, 3, 4, isospin_breaking_vpp=True)
        with self.assertRaises(ValueError):
            off.extend(values)
        with self.assertRaises(ValueError):
            on.extend(values, (1., 2.))
        with self.assertRaises(ValueError):
            on.residual(values)
        with self.assertRaises(ValueError):
            on.jacobian(values)
        with self.assertRaises(ValueError):
            off.residual(np.r_[values, 0., 0., 0.])
        with self.assertRaises(ValueError):
            vpp_option.vertex_factors((1., 2.), derivation.vacuum_spurion_values(), (1., -1., 0.))

    def test_coefficient_domain(self):
        model = vpp_option.IsospinBreakingVppLikelihood(1, 3, 4, isospin_breaking_vpp=True)
        point = model.extend(accepted_endpoint(), (0., 5.5, 0.))
        with self.assertRaises(ValueError):
            model.unpack(point)
        np.testing.assert_array_equal(model.residual(point), np.full(75, 1e8+np.linalg.norm(point)))


if __name__ == '__main__':
    unittest.main()
