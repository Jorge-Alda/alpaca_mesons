"""Tests of the switchable charged/neutral K* treatment of the joint likelihood.

Run: .venv/bin/python -m unittest scripts.test_joint_vector_kstar_charge_options
"""
import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fit_babar_kkpi_global_vector_sector as joint_fit
import joint_vector_kstar_charge_options as kstar_options
from test_joint_vector_likelihood_options import accepted_endpoint

KKPI0_ROWS = 28+np.arange(len(joint_fit.GlobalLikelihood(1, 3, 4).babar['kkpi0']['edges']))


def width_row(model, name):
    """Residual index of a total-width row."""
    return 19+model.total_order.index(name)


class DefaultsReproduceBaseline(unittest.TestCase):
    def test_bitwise_equal_residual_and_jacobian_at_accepted_endpoint(self):
        values = accepted_endpoint()
        baseline = joint_fit.GlobalLikelihood(2, 16, 16, bin_block_size=4)
        for model in (kstar_options.KstarChargeLikelihood(2, 16, 16, bin_block_size=4),
                      kstar_options.KstarChargeLikelihood(
                          2, 16, 16, bin_block_size=4, separate_kstar_charge_states=False,
                          free_kstar_mass_difference=False, kstar_width_treatment='kpi_saturation')):
            self.assertEqual(model.parameter_names, joint_fit.PARAMETER_NAMES)
            np.testing.assert_array_equal(model.initial(), baseline.initial())
            self.assertTrue(np.array_equal(baseline.residual(values), model.residual(values)))
            self.assertTrue(np.array_equal(baseline.jacobian(values), model.jacobian(values)))
        for expected, actual in zip(joint_fit.parameter_bounds(baseline), model.parameter_bounds()):
            np.testing.assert_array_equal(expected, actual)
        self.assertIs(joint_fit.running_propagator, kstar_options.RUNNING_PROPAGATOR)

    def test_bitwise_equal_without_energy_blocks(self):
        values = accepted_endpoint()
        baseline = joint_fit.GlobalLikelihood(1, 3, 4)
        model = kstar_options.KstarChargeLikelihood(1, 3, 4)
        self.assertTrue(np.array_equal(baseline.residual(values), model.residual(values)))


class OptionsChangeIntendedRows(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.values = accepted_endpoint()
        cls.reference = joint_fit.GlobalLikelihood(1, 3, 4, bin_block_size=4)
        cls.baseline = cls.reference.residual(cls.values)

    def model(self, **options):
        return kstar_options.KstarChargeLikelihood(1, 3, 4, bin_block_size=4, **options)

    def changed(self, residual):
        return set(np.flatnonzero(residual != self.baseline).tolist())

    def test_separate_charge_states(self):
        model = self.model(separate_kstar_charge_states=True)
        self.assertEqual(len(model.parameter_names), 32)
        self.assertAlmostEqual(model.pdg_kstar_mass_difference_gev, .00388, places=12)
        changed = self.changed(model.residual(self.values))
        neutral_rows = {width_row(model, 'Kstar0'), 6}  # width row and K*0 -> K0 gamma
        self.assertTrue(neutral_rows <= changed)
        self.assertTrue(changed & set(KKPI0_ROWS.tolist()))
        self.assertTrue(changed & set(range(28+len(KKPI0_ROWS), 73)))
        self.assertTrue(changed <= neutral_rows | set(range(28, 73)))
        self.assertIs(joint_fit.running_propagator, kstar_options.RUNNING_PROPAGATOR)

    def test_free_mass_difference_at_start_equals_fixed_difference(self):
        fixed = self.model(separate_kstar_charge_states=True)
        free = self.model(separate_kstar_charge_states=True, free_kstar_mass_difference=True)
        self.assertEqual(free.parameter_names[-1], 'kstar_neutral_minus_charged_mass_GeV')
        start = free.extend(self.values)
        self.assertEqual(start.shape, (33,))
        np.testing.assert_array_equal(fixed.residual(self.values), free.residual(start))
        lower, upper = free.parameter_bounds()
        self.assertEqual((lower[32], upper[32]), kstar_options.MASS_DIFFERENCE_BOUNDS_GEV)
        moved = start.copy()
        moved[32] += .002
        changed = set(np.flatnonzero(free.residual(moved) != free.residual(start)).tolist())
        # Only the K*0 width row, the K*0 -> K0 gamma row and the K_S K pi bins.
        self.assertTrue(changed and changed <= {width_row(free, 'Kstar0'), 6}
                        | set(range(28+len(KKPI0_ROWS), 73)))

    def test_free_total_widths(self):
        model = self.model(kstar_width_treatment='free_total_widths')
        self.assertEqual(model.parameter_names[32:], ('total_width_pull_Kstarplus', 'total_width_pull_Kstar0'))
        start = model.extend(self.values)
        residual = model.residual(start)
        for name in kstar_options.KSTAR_NAMES:
            self.assertEqual(residual[width_row(model, name)], 0.)
        changed = self.changed(residual)
        self.assertTrue(changed <= {6, 7, width_row(model, 'Kstar0'), width_row(model, 'Kstarplus')}
                        | set(range(28, 73)))
        self.assertTrue(changed & set(range(28, 73)))
        lower, upper = model.parameter_bounds()
        for index, name in zip((32, 33), kstar_options.KSTAR_NAMES):
            width = model.odd_data['total_widths'][name]
            self.assertEqual(lower[index], -.99*width['value_gev']/width['sigma_gev'])
            self.assertEqual(upper[index], np.inf)

    def test_split_reduces_to_shared_for_equal_widths_and_zero_difference(self):
        options = dict(free_kstar_mass_difference=True, kstar_width_treatment='free_total_widths')
        split = self.model(separate_kstar_charge_states=True, **options)
        shared = self.model(kstar_width_treatment='free_total_widths')
        widths = shared.odd_data['total_widths']
        point = split.extend(self.values)
        point[33] = (widths['Kstarplus']['value_gev']-widths['Kstar0']['value_gev'])/widths['Kstar0']['sigma_gev']
        point[34] = 0.
        np.testing.assert_allclose(split.residual(point), shared.residual(point[:34]), rtol=1e-12, atol=1e-12)

    def test_added_jacobian_columns(self):
        model = self.model(separate_kstar_charge_states=True, free_kstar_mass_difference=True,
                           kstar_width_treatment='free_total_widths')
        point = model.extend(self.values)
        jacobian = model.jacobian(point)
        self.assertEqual(jacobian.shape, (75, 35))
        self.assertTrue(np.all(np.isfinite(jacobian)))
        for column, name in zip((32, 33), kstar_options.KSTAR_NAMES):
            self.assertAlmostEqual(jacobian[width_row(model, name), column], 1., places=6)
        self.assertTrue(np.any(jacobian[28:73, 34] != 0))

    def test_invalid_combinations_are_rejected(self):
        with self.assertRaises(ValueError):
            self.model(free_kstar_mass_difference=True)
        with self.assertRaises(ValueError):
            self.model(kstar_width_treatment='shared')


if __name__ == '__main__':
    unittest.main()
