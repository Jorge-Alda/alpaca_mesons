"""Tests of the counting wrappers of benchmark_isospin_breaking_vpp_cost.py (cheap quadrature).

Run: .venv/bin/python -m unittest scripts.test_benchmark_isospin_breaking_vpp_cost
"""
import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import benchmark_isospin_breaking_vpp_cost as benchmark  # noqa: E402
import fit_babar_kkpi_global_vector_sector as joint_fit  # noqa: E402
import joint_vector_isospin_breaking_vpp_option as vpp_option  # noqa: E402

CHEAP_ORDERS = (1, 3, 4)


class InstrumentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.index, endpoint = benchmark.restart_point()
        cls.cases = benchmark.configurations(endpoint, CHEAP_ORDERS)

    def test_restart_run_index(self):
        self.assertEqual(self.index, 2)

    def test_full_residual_counts_per_jacobian(self):
        expected = dict(baseline=(11, 0), option_r_zero=(11, 4), option_restart=(11, 4),
                        reuse_r_zero=(11, 3), reuse_restart=(11, 3))
        for name, (model, point) in self.cases.items():
            calls, _ = benchmark.count_calls(model, point, 'jacobian')
            with self.subTest(configuration=name):
                self.assertEqual((calls.get('full_residual_inside_baseline_jacobian', 0),
                                  calls.get('full_residual_outside', 0)), expected[name])
                self.assertEqual(calls['babar_blocks_channel_restricted'], 6)
                self.assertEqual(calls['babar_blocks_derivative_pass'], 1)

    def test_zero_excitation_calls_only_with_nonzero_r(self):
        model, point = self.cases['option_restart']
        calls, _ = benchmark.count_calls(model, point, 'residual')
        self.assertEqual(calls['running_propagator_kstar_zero_excitation'], calls['running_propagator_kstar'])
        model, point = self.cases['option_r_zero']
        calls, _ = benchmark.count_calls(model, point, 'residual')
        self.assertNotIn('running_propagator_kstar_zero_excitation', calls)

    def test_wrappers_preserve_outputs_and_restore_globals(self):
        running, pair = joint_fit.running_propagator, joint_fit.reference_propagator_pair
        residual = joint_fit.GlobalLikelihood.residual
        self.assertTrue(all(benchmark.outputs_unchanged_by_instrument(self.cases).values()))
        self.assertIs(joint_fit.running_propagator, running)
        self.assertIs(vpp_option.RUNNING_PROPAGATOR, running)
        self.assertIs(joint_fit.reference_propagator_pair, pair)
        self.assertIs(joint_fit.GlobalLikelihood.residual, residual)

    def test_reuse_bitwise_equal(self):
        result = benchmark.reuse_bitwise_equal(self.cases)
        self.assertTrue(all(all(item.values()) for item in result.values()))

    def test_wall_time_estimate(self):
        timings = dict(residual=dict(option_restart=dict(median_seconds=.1)),
                       jacobian=dict(option_restart=dict(median_seconds=1.)))
        self.assertTrue(np.isclose(benchmark.wall_time_estimate(timings, 10, 8)['seconds'], 9.))


if __name__ == '__main__':
    unittest.main()
