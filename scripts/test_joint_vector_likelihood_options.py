"""Tests of the switchable joint-likelihood variations.

Run: .venv/bin/python -m unittest scripts.test_joint_vector_likelihood_options
"""
import json
import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fit_babar_kkpi_global_vector_sector as joint_fit
import joint_vector_likelihood_options as options_module

ROOT = joint_fit.ROOT


def accepted_endpoint():
    comparison = json.loads((ROOT/'data/vector_sector/joint_vector_endpoint_comparison.json').read_text())
    candidate = next(c for c in comparison['candidates'] if c['id'] == 'basin_probes/even_shared_widths')
    assert candidate['accepted']
    return np.asarray(candidate['endpoint'])


class DefaultReproducesBaseline(unittest.TestCase):
    def test_bitwise_equal_residual_and_jacobian_at_accepted_endpoint(self):
        values = accepted_endpoint()
        baseline = joint_fit.GlobalLikelihood(2, 16, 16, bin_block_size=4)
        varied = options_module.VariedLikelihood(2, 16, 16, bin_block_size=4)
        self.assertTrue(np.array_equal(baseline.residual(values), varied.residual(values)))
        self.assertTrue(np.array_equal(baseline.jacobian(values), varied.jacobian(values)))
        self.assertTrue(np.array_equal(baseline.residual(values), varied.full_residual(values)))
        explicit = options_module.VariedLikelihood(2, 16, 16, bin_block_size=4, pole_shift=None)
        self.assertTrue(np.array_equal(baseline.residual(values), explicit.residual(values)))
        self.assertTrue(np.array_equal(baseline.jacobian(values), explicit.jacobian(values)))
        self.assertIs(joint_fit.reference_propagator, options_module.REFERENCE_PROPAGATOR)


class VariationsChangeOnlyIntendedRows(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.values = accepted_endpoint()
        cls.baseline = joint_fit.GlobalLikelihood(1, 3, 4, bin_block_size=4).residual(cls.values)

    def residual(self, **options):
        return options_module.VariedLikelihood(1, 3, 4, bin_block_size=4, **options).full_residual(self.values)

    def changed(self, residual):
        return set(np.flatnonzero(residual != self.baseline).tolist())

    def test_phi_covariance_options(self):
        self.assertEqual(self.changed(self.residual(phi_kk_correlation=False)), {27})
        self.assertEqual(self.changed(self.residual(phi_electronic_width='average')), {18})

    def test_babar_error_and_regulator_options(self):
        self.assertTrue(self.changed(self.residual(babar_fractional_systematic=.05)) <= set(range(28, 73)))
        self.assertTrue(self.changed(self.residual(regulator_gev=3.0)) <= set(range(28, 73)))
        self.assertTrue(self.changed(self.residual(regulator_gev=3.0)))

    def test_window_and_exclusion_masks(self):
        model = options_module.VariedLikelihood(1, 3, 4, bin_block_size=4, babar_upper_energy_gev=1.70)
        dropped = sum(int(np.sum(model.babar[c]['edges'][:, 1] > 1.70+1e-9)) for c in joint_fit.CHANNELS)
        self.assertEqual(dropped, 8)
        residual = model.residual(self.values)
        self.assertEqual(residual.shape, (75-dropped,))
        np.testing.assert_array_equal(residual, self.baseline[model.row_mask])
        self.assertEqual(model.jacobian(self.values).shape, (75-dropped, 32))
        excluded = options_module.VariedLikelihood(1, 3, 4, bin_block_size=4, excluded_blocks=('phi_KK',))
        np.testing.assert_array_equal(excluded.residual(self.values), np.delete(self.baseline, [26, 27]))

    def test_source_and_branch_starts_are_in_domain(self):
        for options, hypothesis, branch in (
                (dict(source_hypothesis='universal'), 'universal', 'near_ideal'),
                (dict(source_hypothesis='single_trace_positive'), 'single_trace_positive', 'near_ideal'),
                (dict(inverse_mass_branch='second'), 'single_trace_negative', 'second')):
            model = options_module.VariedLikelihood(1, 3, 4, bin_block_size=4, **options)
            start = self.values.copy()
            if branch == 'second':
                start[2:4] = 0.
            coupling, slope = model.stage3_source(hypothesis, branch)
            start[4] = (coupling-model.point[4])/joint_fit.EVEN_SCALES[4]
            start[5] = (slope-model.point[5])/joint_fit.EVEN_SCALES[5]
            point = model.unpack(start)[0]
            self.assertTrue(np.all(np.isfinite(model.residual(start))))
            if hypothesis == 'universal':
                self.assertEqual(point[5], 0.)
            if hypothesis == 'single_trace_positive':
                # The baseline endpoint has a negative strange residue.
                with self.assertRaises(ValueError):
                    model.unpack(self.values)

    def test_unknown_source_hypothesis_is_rejected(self):
        with self.assertRaises(ValueError):
            options_module.VariedLikelihood(1, 3, 4, bin_block_size=4, source_hypothesis='double_trace')


class PoleShifts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.values = accepted_endpoint()
        cls.baseline = joint_fit.GlobalLikelihood(1, 3, 4, bin_block_size=4).residual(cls.values)
        cls.kkpi0_rows = len(options_module.VariedLikelihood(1, 3, 4).babar['kkpi0']['edges'])

    def residual(self, state, quantity, shift, values=None):
        model = options_module.VariedLikelihood(
            1, 3, 4, bin_block_size=4, pole_shift=dict(state=state, quantity=quantity, shift_gev=shift))
        return model.residual(self.values if values is None else values)

    def test_single_pole_propagator_matches_reference_propagator(self):
        invariant = np.linspace(.3, 3.2, 41)
        reference = dict(m_rho_GeV=.775, m_Kstar_GeV=.892, Gamma_total_omega_GeV=.0087,
                         Gamma_total_phi_GeV=.0042)
        for state, (species, slot, mass, width) in options_module.EXCITED_POLE_BASELINES.items():
            beta = (1., 0.) if slot == 0 else (0., 1.)
            expected = (options_module.REFERENCE_PROPAGATOR(species, invariant, reference, beta)
                        - options_module.REFERENCE_PROPAGATOR(species, invariant, reference, 0.))
            np.testing.assert_allclose(options_module.excited_pole_propagator(species, invariant, mass, width),
                                       expected, rtol=1e-9, atol=1e-12, err_msg=state)

    def test_zero_shift_reproduces_baseline(self):
        for state in (*options_module.EXCITED_POLE_BASELINES, 'a2(1320)'):
            np.testing.assert_allclose(self.residual(state, 'mass', 0.), self.baseline, rtol=0, atol=1e-10,
                                       err_msg=state)
        self.assertIs(joint_fit.reference_propagator, options_module.REFERENCE_PROPAGATOR)

    def test_shifts_change_only_babar_rows(self):
        for state, quantity in (('K*(1410)', 'mass'), ('phi(2170)', 'width'), ('rho(1450)', 'width')):
            changed = set(np.flatnonzero(self.residual(state, quantity, .02) != self.baseline).tolist())
            self.assertTrue(changed and changed <= set(range(28, 73)), state)
        changed = set(np.flatnonzero(self.residual('a2(1320)', 'width', .005) != self.baseline).tolist())
        self.assertTrue(changed and changed <= set(range(28+self.kkpi0_rows, 73)))

    def test_tensor_shift_vanishes_without_tensor_products(self):
        values = self.values.copy()
        values[24:26] = 0.
        baseline = joint_fit.GlobalLikelihood(1, 3, 4, bin_block_size=4).residual(values)
        np.testing.assert_allclose(self.residual('a2(1320)', 'mass', .01, values), baseline, rtol=0, atol=1e-12)

    def test_shifted_jacobian_is_finite(self):
        model = options_module.VariedLikelihood(
            1, 3, 4, bin_block_size=4, pole_shift=dict(state='a2(1320)', quantity='width', shift_gev=.005))
        jacobian = model.jacobian(self.values)
        self.assertEqual(jacobian.shape, (75, 32))
        self.assertTrue(np.all(np.isfinite(jacobian)))

    def test_declared_variations_read_pdg_uncertainties(self):
        variations = {item['name']: item for item in options_module.pole_shift_variations()}
        self.assertEqual(len(variations), 36)
        self.assertAlmostEqual(variations['pole_phi_2170_width_plus_one_sigma']['options']['pole_shift']['shift_gev'], .024)
        self.assertAlmostEqual(variations['pole_phi_2170_width_minus_one_sigma']['options']['pole_shift']['shift_gev'], -.018)
        self.assertAlmostEqual(variations['pole_a2_1320_width_minus_one_sigma']['options']['pole_shift']['shift_gev'], -.005)
        self.assertEqual(variations['pole_phi_2170_mass_plus_one_sigma']['baseline_value_gev'], 2.162)


if __name__ == '__main__':
    unittest.main()
