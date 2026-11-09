"""Fast tests of the isospin-breaking VPP fit driver (no optimizer is run).

Run: .venv/bin/python -m unittest scripts.test_fit_joint_vector_isospin_breaking_vpp
Cheap quadrature (1,3,4) with four-bin blocks; least_squares is always mocked.
"""
import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fit_babar_kkpi_global_vector_sector as joint_fit  # noqa: E402
import fit_joint_vector_isospin_breaking_vpp as driver  # noqa: E402
import fit_joint_vector_multistart as multi  # noqa: E402
import refit_joint_vector_sensitivity as sensitivity  # noqa: E402

CHEAP_ORDERS = (1, 3, 4)


class ArgumentTests(unittest.TestCase):
    def test_defaults(self):
        args = driver.parse_arguments([])
        self.assertFalse(args.dry_run)
        self.assertEqual(args.output, 'data/vector_sector/joint_vector_isospin_breaking_vpp_fit.json')
        self.assertEqual(args.workers, 2)

    def test_dry_run_and_options(self):
        args = driver.parse_arguments(['--dry-run', '--workers', '3', '--output', 'x.json'])
        self.assertTrue(args.dry_run)
        self.assertEqual((args.workers, args.output), (3, 'x.json'))

    def test_settings_follow_baseline_refits(self):
        self.assertEqual(driver.FIT_ORDERS, (2, 16, 16))
        self.assertEqual(driver.HIGH_ORDERS, (4, 128, 128))
        self.assertEqual(driver.REFINEMENT_ORDERS, (2, 64, 64))
        self.assertEqual(driver.OPTIMIZER_SETTINGS,
                         dict(method='trf', max_nfev=5000, ftol=1e-8, xtol=1e-9, gtol=1e-7))


class ModelAndStartTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.context = sensitivity.accepted_context()
        cls.model = driver.FitLikelihood(*CHEAP_ORDERS, bin_block_size=4, **driver.MODEL_OPTIONS)
        cls.baseline = joint_fit.GlobalLikelihood(*CHEAP_ORDERS, bin_block_size=4)
        cls.starts = driver.start_points(cls.model, cls.context)
        cls.best = np.asarray(cls.context['accepted'][0]['endpoint'])

    def test_model_check_and_bounds(self):
        summary = driver.check_model(self.model)
        self.assertEqual(summary['coordinate_count'], 35)
        self.assertEqual(summary['nominal_degrees_of_freedom'], 40)
        lower, upper = self.model.parameter_bounds()
        baseline_lower, baseline_upper = joint_fit.parameter_bounds(self.model)
        np.testing.assert_array_equal(lower[:32], baseline_lower)
        np.testing.assert_array_equal(upper[:32], baseline_upper)
        np.testing.assert_array_equal(lower[32:], [-5.]*3)
        np.testing.assert_array_equal(upper[32:], [5.]*3)

    def test_model_check_rejects_switches(self):
        for options in (dict(isospin_breaking_vpp=False),
                        dict(isospin_breaking_vpp=True, vpp_s_channel_factors=False),
                        dict(isospin_breaking_vpp=True, vpp_kstar_excitation_scaling='whole_propagator')):
            with self.assertRaises(ValueError):
                driver.check_model(driver.FitLikelihood(*CHEAP_ORDERS, **options))

    def test_start_points(self):
        self.assertEqual([item['name'] for item in self.starts], ['scan_point', 'zero_control'])
        scan, control = self.starts
        for item in self.starts:
            self.assertEqual(item['start_endpoint'], 'basin_probes/even_shared_widths')
            np.testing.assert_array_equal(item['start'][:32], self.best)
        np.testing.assert_array_equal(scan['start'][32:], [-.350, 1.855, .729])
        self.assertEqual(scan['fixed_coordinates'], {})
        self.assertEqual(scan['free_coordinate_count'], 35)
        np.testing.assert_array_equal(control['start'][32:], [0., 0., 0.])
        self.assertEqual(control['fixed_coordinates'],
                         {name: 0. for name in driver.vpp_option.ADDED_COORDINATE_NAMES})
        self.assertEqual(control['free_coordinate_count'], 32)

    def test_start_outside_bounds_is_rejected(self):
        definition = dict(driver.STARTS[0], coefficients_gev_minus2=(0., 6., 0.))
        with self.assertRaises(ValueError):
            driver.start_from(self.model, self.best, definition)

    def test_zero_control_reproduces_baseline_residual(self):
        control = np.asarray(self.starts[1]['start'])
        np.testing.assert_array_equal(self.model.residual(control), self.baseline.residual(self.best))

    def test_parameter_names_are_fixed_during_audit_renaming(self):
        names = self.model.parameter_names
        with driver.extended_parameter_names(names):
            self.assertEqual(joint_fit.PARAMETER_NAMES, names)
            self.assertEqual(self.model.parameter_names, names)
        self.assertEqual(len(joint_fit.PARAMETER_NAMES), 32)

    def test_optimizer_scales(self):
        jacobian = np.ones((75, 35))
        jacobian[:, 32:] *= [1., 2., 4.]
        scales = driver.optimizer_scales(self.context, jacobian)
        np.testing.assert_array_equal(scales[:32], self.context['scales'])
        np.testing.assert_allclose(scales[32:], 1/(np.sqrt(75)*np.array([1., 2., 4.])))
        jacobian[:, 33] = 0.
        with self.assertRaises(ValueError):
            driver.optimizer_scales(self.context, jacobian)

    def test_prepared_ledger_schema(self):
        scales = np.r_[self.context['scales'], 1., 1., 1.]
        gate = dict(passed=True)
        ledger = driver.prepare_ledger(self.model, self.context, scales, self.starts, gate)
        for key in ('study', 'state', 'model_options', 'parameter_names', 'bounds', 'starts',
                    'nominal_degrees_of_freedom', 'fit_quadrature', 'evaluation_quadrature',
                    'refinement_quadrature', 'optimizer_settings', 'optimizer_coordinate_scales',
                    'acceptance_rule', 'basin_switch_rule', 'adequacy_criterion', 'control_rule',
                    'source_hash_gate', 'entries', 'input_sha256', 'source_sha256',
                    'baseline_chi2_high_order', 'accepted_endpoints'):
            self.assertIn(key, ledger)
        self.assertEqual(ledger['state'], 'prepared')
        self.assertEqual(ledger['nominal_degrees_of_freedom'], 40)
        self.assertEqual(len(ledger['parameter_names']), 35)
        self.assertEqual(len(ledger['bounds']['lower']), 35)
        self.assertEqual(ledger['bounds']['upper'][32:], [5., 5., 5.])
        self.assertEqual(len(ledger['optimizer_coordinate_scales']), 35)
        self.assertEqual((ledger['fit_quadrature']['energy'], ledger['fit_quadrature']['pair_per_segment'],
                          ledger['fit_quadrature']['angle_per_segment']), (2, 16, 16))
        self.assertEqual(ledger['evaluation_quadrature']['pair_per_segment'], 128)
        self.assertIn('scripts/fit_joint_vector_isospin_breaking_vpp.py', ledger['source_sha256'])
        with tempfile.TemporaryDirectory() as directory:
            multi.write_ledger(Path(directory)/'ledger.json', ledger)
            reread = json.loads((Path(directory)/'ledger.json').read_text())
        self.assertEqual(reread['starts'][0]['start'], self.starts[0]['start'])

    def test_refit_record_with_mocked_optimizer(self):
        control = self.starts[1]
        fixed = {32: 0., 33: 0., 34: 0.}
        free = np.arange(32)
        fake = SimpleNamespace(x=np.asarray(control['start'])[free], fun=np.ones(75), success=True, status=2,
                               message='mock', nfev=1, njev=1, optimality=0., active_mask=np.zeros(32, int),
                               jac=np.eye(75, 32))
        with mock.patch.object(driver, 'least_squares', return_value=fake) as optimizer:
            record = driver.refit(self.model, control['start'], fixed, np.ones(35), driver.OPTIMIZER_SETTINGS)
        optimizer.assert_called_once()
        for key in ('start', 'fixed_coordinates', 'free_coordinate_count', 'state', 'endpoint', 'fit_order_chi2',
                    'success', 'status', 'message', 'nfev', 'njev', 'optimality', 'active_bounds',
                    'penalty_endpoint', 'jacobian_conditioning', 'refit_seconds'):
            self.assertIn(key, record)
        self.assertEqual(record['free_coordinate_count'], 32)
        self.assertEqual(record['endpoint'], control['start'])
        self.assertEqual(len(joint_fit.PARAMETER_NAMES), 32)

    def test_failed_refit_is_kept(self):
        with mock.patch.object(driver, 'least_squares') as optimizer:
            record = driver.refit(self.model, np.zeros(32), {}, np.ones(35), driver.OPTIMIZER_SETTINGS)
        optimizer.assert_not_called()
        self.assertEqual(record['state'], 'error')
        self.assertFalse(record['success'])


class DryRunTests(unittest.TestCase):
    def test_dry_run_never_calls_optimizer(self):
        with mock.patch.object(driver, 'least_squares') as optimizer, \
                mock.patch.object(sensitivity, 'least_squares') as sensitivity_optimizer, \
                mock.patch.object(multi, 'least_squares') as multi_optimizer, \
                mock.patch.object(driver, 'run_start') as run_start, \
                mock.patch.object(driver, 'ProcessPoolExecutor') as executor:
            summary = driver.dry_run(CHEAP_ORDERS)
        for patched in (optimizer, sensitivity_optimizer, multi_optimizer, run_start, executor):
            patched.assert_not_called()
        self.assertTrue(summary['dry_run'])
        self.assertTrue(summary['source_hash_gate']['passed'])
        self.assertTrue(all(summary['zero_control_equals_baseline'].values()))
        scan, control = summary['starts']['scan_point'], summary['starts']['zero_control']
        self.assertEqual((scan['rows'], scan['columns'], scan['degrees_of_freedom']), (75, 35, 40))
        self.assertEqual((control['rows'], control['columns'], control['degrees_of_freedom']), (75, 35, 43))
        self.assertAlmostEqual(sum(block['chi2'] for block in scan['blocks'].values()), scan['chi2'])

    def test_main_dry_run_writes_nothing(self):
        with mock.patch.object(driver, 'dry_run', return_value=dict(dry_run=True)) as dry, \
                mock.patch.object(multi, 'write_ledger') as writer, \
                mock.patch.object(driver, 'least_squares') as optimizer, \
                mock.patch('builtins.print'):
            driver.main(['--dry-run'])
        dry.assert_called_once()
        writer.assert_not_called()
        optimizer.assert_not_called()


class GateTests(unittest.TestCase):
    def test_adequacy_verdict(self):
        def high(total, block):
            return dict(valid=True, nominal_p_value=total,
                        blocks={'a': dict(nominal_p_value_row_count_dof=block),
                                'b': dict(nominal_p_value_row_count_dof=.5)})
        self.assertTrue(driver.adequacy_verdict(high(.02, 2e-3))['meets_criterion'])
        self.assertFalse(driver.adequacy_verdict(high(.005, .5))['meets_criterion'])
        self.assertEqual(driver.adequacy_verdict(high(.5, 5e-4))['failing_blocks'], ['a'])
        self.assertFalse(driver.adequacy_verdict(dict(valid=False))['meets_criterion'])

    def test_control_reproduction_rule(self):
        context = dict(baseline_chi2=225.)
        good = dict(high_order=dict(valid=True, chi2=225.005), scaled_distance_to_best=.005)
        self.assertTrue(driver.control_reproduces_baseline(good, context))
        self.assertFalse(driver.control_reproduces_baseline(dict(good, scaled_distance_to_best=.02), context))
        self.assertFalse(driver.control_reproduces_baseline(
            dict(good, high_order=dict(valid=True, chi2=225.02)), context))

    def test_source_hash_gate_detects_change(self):
        baseline = joint_fit.GlobalLikelihood(*CHEAP_ORDERS, bin_block_size=4)
        self.assertTrue(driver.source_hash_gate(baseline)['passed'])
        reference = json.loads((driver.ROOT/driver.REFERENCE_LEDGER).read_text())
        reference['source_sha256']['scripts/fit_vector_gv.py'] = '0'*64
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'reference.json'  # absolute: ROOT/path is path itself
            path.write_text(json.dumps(dict(input_sha256=reference['input_sha256'],
                                            source_sha256=reference['source_sha256'])))
            with self.assertRaisesRegex(ValueError, 'fit_vector_gv'):
                driver.source_hash_gate(baseline, str(path))

    def test_unsuccessful_run_is_not_audited(self):
        result, accepted = driver.audit(None, dict(success=False, state='finished'))
        self.assertEqual(result['status'], 'not_evaluated')
        self.assertFalse(accepted)


if __name__ == '__main__':
    unittest.main()
