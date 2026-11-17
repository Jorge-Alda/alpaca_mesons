"""Tests of the driver of the joint vector fit with the double-trace singlet-source coordinate y_V.

Run: .venv/bin/python -m unittest scripts.test_fit_joint_vector_double_trace_singlet_source
Cheap quadrature orders (1, 3, 4); the optimizer is mocked or forbidden.
"""
import json
import sys
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fit_babar_kkpi_global_vector_sector as joint_fit  # noqa: E402
import fit_joint_vector_double_trace_singlet_source as driver  # noqa: E402

ORDERS = (1, 3, 4)


def forbidden(*arguments, **keywords):
    raise AssertionError('the optimizer must not be called')


class Starts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.context = driver.sensitivity.accepted_context()
        cls.model = driver.make_models((ORDERS,))[0]
        cls.definitions = driver.start_definitions(cls.context)

    def test_start_definitions(self):
        names = [item['name'] for item in self.definitions]
        self.assertEqual(len(names), 7)
        self.assertEqual(len(set(names)), 7)
        zero_free = [item for item in self.definitions
                     if item['coefficient_gev_inverse'] == 0. and not item['coefficient_fixed']]
        self.assertEqual([item['start_endpoint'] for item in zero_free],
                         [item['id'] for item in self.context['accepted']])
        offsets = sorted(item['coefficient_gev_inverse'] for item in self.definitions
                         if item['coefficient_gev_inverse'] != 0.)
        self.assertEqual(offsets, [-.022, .022])
        fixed = [item['name'] for item in self.definitions if item['coefficient_fixed']]
        self.assertEqual(fixed, [driver.CONTROL])

    def test_check_model(self):
        summary = driver.check_model(self.model)
        self.assertEqual(summary['coordinate_count'], 33)
        self.assertEqual(summary['nominal_degrees_of_freedom'], 42)
        with self.assertRaises(ValueError):
            driver.check_model(driver.option.DoubleTraceSingletSourceLikelihood(*ORDERS, bin_block_size=4))

    def test_zero_starts_reproduce_baseline(self):
        baseline = joint_fit.GlobalLikelihood(*ORDERS, bin_block_size=4)
        for definition in self.definitions:
            if definition['coefficient_gev_inverse'] != 0.:
                continue
            values, fixed = driver.start_from(self.model, self.context, definition)
            endpoint = np.asarray(driver.endpoint_by_identifier(self.context, definition['start_endpoint']))
            self.assertTrue(np.array_equal(self.model.residual(values), baseline.residual(endpoint)))
            self.assertEqual(fixed, {32: 0.} if definition['coefficient_fixed'] else {})

    def test_validate_point(self):
        lower, upper = self.model.parameter_bounds()
        values, _ = driver.start_from(self.model, self.context, self.definitions[0])
        with self.assertRaises(ValueError):
            driver.validate_point(self.model, values[:32], lower, upper)
        outside = values.copy()
        outside[32] = .6
        with self.assertRaises(ValueError):
            driver.validate_point(self.model, outside, lower, upper)

    def test_optimizer_scales(self):
        jacobian = np.zeros((75, 33))
        jacobian[0, 32] = 4.
        scales = driver.optimizer_scales(self.context, jacobian)
        self.assertEqual(len(scales), 33)
        self.assertEqual(scales[32], .25)
        self.assertTrue(np.array_equal(scales[:32], self.context['scales']))
        with self.assertRaises(ValueError):
            driver.optimizer_scales(self.context, np.zeros((75, 33)))


class Refit(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.context = driver.sensitivity.accepted_context()
        cls.model = driver.make_models((ORDERS,))[0]
        cls.definitions = driver.start_definitions(cls.context)
        cls.scales = np.ones(33)

    def fake_result(self, values, free):
        residual = self.model.residual(values)
        return SimpleNamespace(x=values[free]+1e-6, fun=residual, success=True, status=2, message='mock',
                               nfev=3, njev=2, optimality=.1, active_mask=np.zeros(len(free), dtype=int),
                               jac=self.model.jacobian(values)[:, free])

    def test_refit_record_and_fixed_coordinate(self):
        definition = next(item for item in self.definitions if item['name'] == driver.CONTROL)
        values, fixed = driver.start_from(self.model, self.context, definition)
        free = np.arange(32)
        with mock.patch.object(driver, 'least_squares', return_value=self.fake_result(values, free)) as optimizer:
            record = driver.refit(self.model, values, fixed, self.scales, driver.OPTIMIZER_SETTINGS)
        self.assertEqual(optimizer.call_args.args[1].shape, (32,))
        self.assertEqual(record['state'], 'finished')
        self.assertEqual(record['endpoint'][32], 0.)
        self.assertEqual(record['free_coordinate_count'], 32)
        for key in ('status', 'message', 'nfev', 'njev', 'optimality', 'active_bounds', 'jacobian_conditioning'):
            self.assertIn(key, record)
        json.dumps(record)

    def test_refit_error_is_kept(self):
        values, fixed = driver.start_from(self.model, self.context, self.definitions[0])
        with mock.patch.object(driver, 'least_squares', side_effect=ValueError('mock failure')):
            record = driver.refit(self.model, values, fixed, self.scales, driver.OPTIMIZER_SETTINGS)
        self.assertEqual((record['state'], record['success'], record['message']), ('error', False, 'mock failure'))

    def test_selected_run(self):
        runs = [dict(coefficient_fixed=False, accepted=True, high_order=dict(valid=True, chi2=200.)),
                dict(coefficient_fixed=False, accepted=False, high_order=dict(valid=True, chi2=100.)),
                dict(coefficient_fixed=True, accepted=True, high_order=dict(valid=True, chi2=150.)),
                dict(coefficient_fixed=False, accepted=True, high_order=dict(valid=True, chi2=190.)),
                dict(coefficient_fixed=False, state='error', accepted=False)]
        self.assertEqual(driver.selected_run(runs), 3)
        self.assertIsNone(driver.selected_run(runs[1:3]))

    def test_prepared_ledger(self):
        ledger = driver.prepare_ledger(self.model, self.context, self.scales, self.definitions, dict(passed=True))
        self.assertEqual(ledger['state'], 'prepared')
        self.assertEqual(ledger['nominal_degrees_of_freedom'], 42)
        self.assertEqual(ledger['bounds']['lower'][32], -.5)
        self.assertEqual(ledger['bounds']['upper'][32], .5)
        self.assertEqual(len(ledger['starts']), 7)
        self.assertIn('scripts/joint_vector_double_trace_singlet_source_option.py', ledger['source_sha256'])
        self.assertEqual(json.loads(json.dumps(ledger)), ledger)


class DryRun(unittest.TestCase):
    def test_dry_run_never_calls_optimizer(self):
        with mock.patch.object(driver, 'least_squares', side_effect=forbidden):
            summary = driver.dry_run(ORDERS)
        self.assertFalse(summary['optimizer_called'])
        self.assertEqual(len(summary['starts']), 7)
        self.assertTrue(summary['source_hash_gate']['passed'])
        for name, start in summary['starts'].items():
            self.assertEqual(start['columns'], 33)
            if start['coefficient_gev_inverse'] == 0.:
                self.assertTrue(all(start['equals_baseline_bitwise'].values()))
        self.assertEqual(summary['starts'][driver.CONTROL]['degrees_of_freedom'], 43)
        self.assertGreater(summary['added_optimizer_scale_gev_inverse'], 0.)


if __name__ == '__main__':
    unittest.main()
