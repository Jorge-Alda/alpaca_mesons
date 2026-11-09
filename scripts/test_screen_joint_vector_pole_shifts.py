"""Tests of the linear pole-shift screen.

Run: .venv/bin/python -m unittest scripts.test_screen_joint_vector_pole_shifts
"""
import json
import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import screen_joint_vector_pole_shifts as screen_module


class LinearMinimumTests(unittest.TestCase):
    def test_linear_minimum_is_projection_residual(self):
        generator = np.random.default_rng(1)
        jacobian = generator.normal(size=(20, 4))
        residual = generator.normal(size=20)
        step, value = screen_module.linear_minimum(residual, jacobian)
        projector = jacobian@np.linalg.pinv(jacobian)
        remainder = residual-projector@residual
        self.assertAlmostEqual(value, float(remainder @ remainder))
        np.testing.assert_allclose(jacobian.T@(residual+jacobian@step), 0., atol=1e-12)

    def test_stored_screen_covers_all_shifts_and_orders_by_size(self):
        ledger = json.loads((screen_module.ROOT/'data/vector_sector/joint_vector_pole_shift_screen.json').read_text())
        self.assertEqual(len(ledger['entries']), 36)
        sizes = {entry['name']: abs(entry['linear_refit_delta_chi2']) for entry in ledger['entries']}
        ordered = [sizes[name] for name in ledger['refit_order']]
        self.assertEqual(ordered, sorted(ordered, reverse=True))


if __name__ == '__main__':
    unittest.main()
