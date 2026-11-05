"""Checks for deterministic joint-fit starts and attempt records.

Run: .venv/bin/python -m unittest scripts.test_fit_joint_vector_multistart
"""
import json
import hashlib
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fit_babar_kkpi_global_vector_sector as joint_fit
import fit_joint_vector_multistart as multi
from scripts import test_global_kkpi_running_propagator as original_path


class JointVectorStarts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = joint_fit.GlobalLikelihood(1, 3, 4, bin_block_size=4)
        cls.reference = np.array(json.loads((multi.ROOT / multi.REFERENCE).read_text())['parameter_values'])
        cls.lower, cls.upper = joint_fit.parameter_bounds(cls.model)

    def test_bounds_match_original_fit(self):
        np.testing.assert_array_equal(self.lower[:7], [-100, -100, -20, -10, -15, -10, .005])
        np.testing.assert_array_equal(self.upper[:7], [100, 100, 20, 10, 15, 10, .2])
        for index, name in enumerate(joint_fit.FREE_TOTALS):
            row = self.model.odd_data['total_widths'][name]
            self.assertEqual(self.lower[7+index], -.99 * row['value_gev'] / row['sigma_gev'])
        np.testing.assert_array_equal(self.lower[26:30], -np.ones(4)*5)
        np.testing.assert_array_equal(self.upper[26:30], np.ones(4)*5)
        self.assertTrue(np.isneginf(self.lower[12:26]).all())
        self.assertTrue(np.isposinf(self.upper[7:26]).all())

    def test_starts_are_deterministic_and_physically_valid(self):
        first = multi.named_starts(self.model, self.reference)
        second = multi.named_starts(self.model, self.reference)
        self.assertEqual(list(first), ['stored', 'preliminary', 'stored_negative_c4', 'stored_negative_odd'])
        for name, values in first.items():
            np.testing.assert_array_equal(values, second[name])
            np.testing.assert_array_equal(multi.validate_start(self.model, values, self.lower, self.upper), values)
        np.testing.assert_array_equal(first['stored'], self.reference)
        self.assertFalse(np.array_equal(first['stored_negative_odd'], first['stored']))

    def test_invalid_starts_are_rejected(self):
        for invalid in (self.reference[:31], np.full(32, np.nan)):
            with self.assertRaises(ValueError):
                multi.validate_start(self.model, invalid, self.lower, self.upper)
        invalid = self.reference.copy(); invalid[26] = 5.1
        with self.assertRaisesRegex(ValueError, 'bounds'):
            multi.validate_start(self.model, invalid, self.lower, self.upper)
        # A sign-reversed source is within the rectangular bounds but crosses
        # the source branch enforced by the likelihood.
        invalid = self.reference.copy(); invalid[5] = 10
        self.assertTrue(np.all(invalid >= self.lower) and np.all(invalid <= self.upper))
        with self.assertRaises(ValueError):
            multi.validate_start(self.model, invalid, self.lower, self.upper)

    def test_fixed_scales_and_attempt_serialization(self):
        first = multi.optimizer_scales(self.model, self.reference)
        second = multi.optimizer_scales(self.model, self.reference)
        np.testing.assert_array_equal(first, second)
        self.assertTrue(np.all(np.isfinite(first)) and np.all(first > 0))
        settings = {'method': 'trf', 'max_nfev': 2, 'ftol': 1e-8, 'xtol': 1e-9, 'gtol': 1e-7}
        with patch.object(multi, 'optimizer_scales', return_value=first):
            ledger = multi.prepare_ledger(self.model, (1, 3, 4), 4, settings)
        self.assertEqual(len(ledger['parameter_names']), 32)
        self.assertEqual(ledger['bounds']['lower'][12], None)
        self.assertEqual(ledger['quadrature']['bin_block_size'], 4)
        outcomes = []
        for index in range(2):
            endpoint = self.reference.copy(); endpoint[16] += .01 * index
            result = SimpleNamespace(x=endpoint, fun=np.full(75, index + 1.),
                success=bool(index), status=index, message=f'outcome {index}', nfev=2,
                njev=1, optimality=3., grad=np.ones(32), active_mask=np.zeros(32, dtype=int))
            with patch.object(multi, 'least_squares', return_value=result) as optimizer:
                multi.run_attempt(self.model, ledger['attempts'][index], self.lower, self.upper, first, settings)
                np.testing.assert_array_equal(optimizer.call_args.kwargs['x_scale'], first)
            outcomes.append(endpoint.tolist())
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'ledger.json'
            multi.write_ledger(path, ledger)
            restored = json.loads(path.read_text())
        self.assertEqual(restored['attempts'][0]['endpoint'], outcomes[0])
        self.assertEqual(restored['attempts'][1]['endpoint'], outcomes[1])
        self.assertFalse(restored['attempts'][0]['success'])
        self.assertTrue(restored['attempts'][1]['success'])
        self.assertFalse(restored['attempts'][1]['accepted'])
        self.assertEqual(restored['attempts'][0]['message'], 'outcome 0')
        self.assertTrue(restored['attempts'][0]['endpoint_diagnostics']['valid'])
        self.assertEqual(sum(item['count'] for item in restored['attempts'][0]['endpoint_diagnostics']['blocks'].values()), 75)
        self.assertEqual(len(restored['attempts'][0]['endpoint_diagnostics']['babar_predictions']['kkpi0']), 21)

    def test_invalid_endpoint_and_optimizer_error_are_recorded(self):
        invalid = self.reference.copy(); invalid[5] = 10
        result = multi.endpoint_diagnostics(self.model, invalid, self.lower, self.upper)
        self.assertFalse(result['valid'])
        self.assertIn('branch', result['message'])
        attempt = {'name': 'stored', 'start': self.reference.tolist(), 'accepted': False}
        with patch.object(multi, 'least_squares', side_effect=ValueError('optimizer failure')):
            multi.run_attempt(self.model, attempt, self.lower, self.upper, np.ones(32), {})
        self.assertEqual(attempt['state'], 'error')
        self.assertFalse(attempt['success'])
        self.assertEqual(attempt['message'], 'optimizer failure')
        self.assertFalse(attempt['accepted'])

    def test_refinement_uses_certified_orders_and_guards_failed_endpoints(self):
        attempt = {'success': False, 'endpoint': self.reference.tolist(),
                   'endpoint_diagnostics': {'valid': True}}
        with patch.object(multi.joint_fit, 'audit_quadrature') as audit:
            result = multi.refine_attempt(self.model, attempt, self.lower, self.upper)
            self.assertFalse(result['passed'])
            audit.assert_not_called()
        attempt['success'] = True
        with patch.object(multi.joint_fit, 'audit_quadrature',
                          return_value={'status': 'certified', 'passed': True}) as audit:
            result = multi.refine_attempt(self.model, attempt, self.lower, self.upper)
        self.assertTrue(result['passed'])
        self.assertEqual(audit.call_args.args[2], {'energy': [2, 4], 'pair': [64, 128], 'angle': [64, 128]})
        np.testing.assert_array_equal(audit.call_args.kwargs['parameter_scales'], np.ones(32))
        self.assertEqual(result['bin_block_size'], 4)

    def test_all_starts_match_original_full_grid_calculation(self):
        full = joint_fit.GlobalLikelihood(1, 3, 4)
        blocked = joint_fit.GlobalLikelihood(1, 3, 4, bin_block_size=4)
        starts = multi.named_starts(blocked, self.reference)
        for name, values in starts.items():
            with self.subTest(start=name):
                # The earlier code evaluates each reference propagator
                # separately and finite-differences complete residuals.
                with (patch.object(original_path, 'fit', joint_fit),
                      patch.object(joint_fit, 'running_propagator',
                                   side_effect=original_path.original_running_propagator),
                      patch.object(joint_fit, 'reference_propagator_pair',
                                   side_effect=original_path.original_reference_propagator_pair)):
                    expected_residual = full.residual(values)
                    _, expected_predictions = full.babar_blocks(values, full.unpack(values))
                    expected_jacobian = original_path.original_jacobian(full, values)
                actual_residual = blocked.residual(values)
                _, actual_predictions = blocked.babar_blocks(values, blocked.unpack(values))
                actual_jacobian = blocked.jacobian(values)
                np.testing.assert_array_equal(actual_residual.view(np.uint64),
                                              expected_residual.view(np.uint64))
                for channel in joint_fit.CHANNELS:
                    np.testing.assert_array_equal(actual_predictions[channel].view(np.uint64),
                                                  expected_predictions[channel].view(np.uint64))
                np.testing.assert_array_equal(actual_jacobian.view(np.uint64),
                                              expected_jacobian.view(np.uint64))

    def test_continuation_preserves_saved_endpoint_and_scales(self):
        source_path = multi.ROOT / 'data/vector_sector/joint_vector_multistart.json'
        source_bytes = source_path.read_bytes()
        source = json.loads(source_bytes)
        settings = {'method': 'trf', 'max_nfev': 5000, 'ftol': 1e-8,
                    'xtol': 1e-9, 'gtol': 1e-7}
        ledger = multi.prepare_continuation(self.model, source_path, 'stored_negative_c4',
                                            (1, 3, 4), 4, settings)
        prior = next(a for a in source['attempts'] if a['name'] == 'stored_negative_c4')
        self.assertEqual(ledger['attempts'][0]['start'], prior['endpoint'])
        self.assertEqual(ledger['optimizer_coordinate_scales'], source['optimizer_coordinate_scales'])
        self.assertEqual(ledger['attempts'][0]['prior_nfev'], prior['nfev'])
        self.assertEqual(ledger['resume_source']['sha256'], hashlib.sha256(source_bytes).hexdigest())
        self.assertFalse(ledger['resume_source']['source_success'])
        self.assertEqual(source_path.read_bytes(), source_bytes)

    def test_continuation_rejects_bounds_and_branch_crossings(self):
        source = json.loads((multi.ROOT / 'data/vector_sector/joint_vector_multistart.json').read_text())
        settings = {'method': 'trf', 'max_nfev': 5, 'ftol': 1e-8, 'xtol': 1e-9, 'gtol': 1e-7}
        for coordinate, value in ((26, 5.1), (5, 10.0)):
            altered = json.loads(json.dumps(source))
            attempt = next(a for a in altered['attempts'] if a['name'] == 'stored_negative_c4')
            attempt['endpoint'][coordinate] = value
            with tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / 'source.json'
                path.write_text(json.dumps(altered))
                with self.assertRaises(ValueError):
                    multi.prepare_continuation(self.model, path, 'stored_negative_c4', (1, 3, 4), 4, settings)

    def test_checkpoint_contains_iteration_and_optional_optimality(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'continuation.progress.json'
            result = SimpleNamespace(x=self.reference, fun=np.array([3., 4.]), nfev=7, nit=4)
            multi.write_checkpoint(path, result, 1500)
            checkpoint = json.loads(path.read_text())
            self.assertEqual(checkpoint['parameter_values'], self.reference.tolist())
            self.assertEqual(checkpoint['chi2'], 25.)
            self.assertEqual((checkpoint['nfev'], checkpoint['cumulative_nfev']), (7, 1507))
            self.assertIsNone(checkpoint['optimality'])
            result.optimality = .125
            multi.write_checkpoint(path, result, 1500)
            self.assertEqual(json.loads(path.read_text())['optimality'], .125)
            self.assertFalse(path.with_name(path.name + '.tmp').exists())

    def test_run_attempt_passes_checkpoint_callback(self):
        attempt = {'name': 'stored', 'start': self.reference.tolist(), 'accepted': False}
        result = SimpleNamespace(x=self.reference, fun=np.ones(75), success=False, status=0,
                                 message='limit', nfev=2, njev=1, optimality=3.,
                                 grad=np.ones(32), active_mask=np.zeros(32, dtype=int))
        observed = []
        def callback(intermediate_result):
            observed.append(intermediate_result.nfev)
        with patch.object(multi, 'least_squares', return_value=result) as optimizer:
            multi.run_attempt(self.model, attempt, self.lower, self.upper, np.ones(32), {}, callback)
            optimizer.call_args.kwargs['callback'](SimpleNamespace(nfev=2))
        self.assertEqual(observed, [2])


if __name__ == '__main__':
    unittest.main()
