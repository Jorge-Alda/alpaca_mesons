"""Tests of the block-exclusion and variation refit driver on cheap orders.

Run: .venv/bin/python -m unittest scripts.test_refit_joint_vector_sensitivity
"""
import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import assess_joint_vector_fit_adequacy as adequacy
import joint_vector_likelihood_options as options_module
import refit_joint_vector_sensitivity as driver

SETTINGS = dict(driver.OPTIMIZER_SETTINGS, max_nfev=2)


def cheap_context():
    comparison = json.loads((adequacy.ROOT/adequacy.COMPARISON).read_text())
    accepted = adequacy.load_accepted_endpoints(comparison)
    return dict(accepted=[dict(id=c['id'], endpoint=c['endpoint']) for c in accepted],
                scales=comparison['optimizer_coordinate_scales'], local_errors=[1.]*32,
                baseline_blocks={}, baseline_chi2=0.)


class ExclusionTests(unittest.TestCase):
    def test_reduced_block_exclusion_record(self):
        item = next(item for item in driver.EXCLUSIONS if item['name'] == 'without_phi_KK')
        models = [options_module.VariedLikelihood(1, 3, 4, bin_block_size=4, excluded_blocks=('phi_KK',))
                  for _ in range(2)]
        record = driver.run_exclusion(item, cheap_context(), models, SETTINGS)
        self.assertEqual(record['refit']['state'], 'finished')
        self.assertIn('status', record['refit'])
        high = record['high_order']
        self.assertEqual(high['residual_count'], 73)
        self.assertEqual(high['degrees_of_freedom'], 41)
        self.assertNotIn('phi_KK', [name for name, block in high['blocks'].items() if block['count']])
        self.assertEqual([row['index'] for row in high['removed_row_pulls']], [26, 27])
        full = np.asarray(high['full_residuals'])
        self.assertAlmostEqual(high['chi2'], float(np.delete(full, [26, 27]) @ np.delete(full, [26, 27])))
        self.assertIn('conditional', record)

    def test_singular_column_detection(self):
        jacobian = np.zeros((40, 32))
        jacobian[:32, :] = np.eye(32)
        jacobian[:, 5] = 0.
        result = driver.near_singular_columns(jacobian, np.ones(32), np.arange(32))
        self.assertFalse(driver.near_singular_columns(np.delete(jacobian, 5, axis=1), np.ones(32),
                                                      np.delete(np.arange(32), 5))['conditional'])
        self.assertTrue(result['conditional'])
        self.assertEqual(result['zero_columns'], ['scaled_z_V_GeV_inverse'])

    def test_babar_exclusion_has_no_nominal_p_value_but_constrained_count(self):
        item = next(item for item in driver.EXCLUSIONS if item['name'] == 'without_BaBar')
        models = [options_module.VariedLikelihood(1, 3, 4, bin_block_size=4, excluded_blocks=('BaBar',))
                  for _ in range(2)]
        record = driver.run_exclusion(item, cheap_context(), models, SETTINGS)
        high = record['high_order']
        self.assertEqual(high['residual_count'], 30)
        self.assertEqual(high['degrees_of_freedom'], -2)
        self.assertIsNone(high['nominal_p_value'])
        zero = len(record['refit']['jacobian_conditioning']['zero_columns'])
        self.assertGreater(zero, 0)
        self.assertEqual(high['degrees_of_freedom_without_zero_columns'], 30-(32-zero))
        self.assertTrue(record['conditional'])
        json.dumps(record, allow_nan=False)


class CheckpointTests(unittest.TestCase):
    def test_refit_resumes_from_checkpoint_with_remaining_budget(self):
        context = cheap_context()
        model = options_module.VariedLikelihood(1, 3, 4, bin_block_size=4)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'item.progress.json'
            first = driver.refit(model, context['accepted'][0]['endpoint'], {}, context['scales'], SETTINGS, path)
            self.assertTrue(path.exists())
            second = driver.refit(model, context['accepted'][0]['endpoint'], {}, context['scales'], SETTINGS, path)
        self.assertNotIn('resumed_from_checkpoint', first)
        self.assertTrue(second['resumed_from_checkpoint'])
        self.assertGreater(second['prior_nfev'], 0)


class VariationTests(unittest.TestCase):
    def test_variation_runs_without_audit_and_fixes_tensor(self):
        variation = next(item for item in driver.VARIATIONS if item['name'] == 'tensor_products_zero')
        models = [options_module.VariedLikelihood(1, 3, 4, bin_block_size=4) for _ in range(3)]
        record = driver.run_variation(variation, cheap_context(), models, SETTINGS, with_audit=False)
        run = record['runs'][0]
        self.assertEqual(run['free_coordinate_count'], 30)
        self.assertEqual(run['endpoint'][24:26], [0., 0.])
        self.assertEqual(run['high_order']['degrees_of_freedom'], 45)

    def test_second_branch_shift_uses_baseline_reference(self):
        model = options_module.VariedLikelihood(1, 3, 4, bin_block_size=4, inverse_mass_branch='second')
        values = np.zeros(32)
        converted = driver.baseline_coordinates(model, values)
        self.assertAlmostEqual(converted[2]*0.01+model.preliminary['Delta_V_GeV2'], model.point[2])

    def test_not_implemented_entries_state_missing_ingredient(self):
        self.assertEqual([entry['name'] for entry in driver.NOT_IMPLEMENTED], ['source_single_plus_double_trace'])
        self.assertIn('deferred by user decision 2026-11-07', driver.NOT_IMPLEMENTED[0]['reason'])

    def test_pole_shift_variations_are_declared_and_refit(self):
        names = [item['name'] for item in driver.VARIATIONS]
        self.assertEqual(len(names), 11+36)
        self.assertEqual(len(set(names)), len(names))
        variation = next(item for item in driver.VARIATIONS if item['name'] == 'pole_a2_1320_width_plus_one_sigma')
        models = [options_module.VariedLikelihood(1, 3, 4, bin_block_size=4, **variation['options']) for _ in range(3)]
        record = driver.run_variation(variation, cheap_context(), models, SETTINGS, with_audit=False)
        self.assertEqual(record['runs'][0]['state'], 'finished')
        self.assertEqual(record['runs'][0]['high_order']['degrees_of_freedom'], 43)
        json.dumps(record, allow_nan=False)

    def test_selected_run_prefers_accepted_then_lowest_chi2(self):
        runs = [dict(accepted=False, high_order=dict(valid=True, chi2=1.)),
                dict(accepted=True, high_order=dict(valid=True, chi2=3.)),
                dict(accepted=True, high_order=dict(valid=True, chi2=2.)), dict(state='error', accepted=False)]
        record = driver.variation_record(dict(name='x'), runs, 5.)
        self.assertEqual((record['selected_run'], record['selected_accepted'], record['basin_switch_reruns']),
                         (2, True, True))


if __name__ == '__main__':
    unittest.main()
