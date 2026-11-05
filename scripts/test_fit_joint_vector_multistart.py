"""Checks for deterministic joint-fit starts and attempt records.

Run: .venv/bin/python -m unittest scripts.test_fit_joint_vector_multistart
"""
import json
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


if __name__ == '__main__':
    unittest.main()
