"""Tests of the deterministic basin starts; no optimizer fit is run.

Run: .venv/bin/python -m unittest scripts.test_fit_joint_vector_basin_probes
"""
import hashlib
import json
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fit_babar_kkpi_global_vector_sector as joint_fit
import fit_joint_vector_basin_probes as probes
import fit_joint_vector_multistart as multi


class BasinProbeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = multi.ROOT / 'data/vector_sector/joint_vector_continuation_e4p128a128.json'
        cls.source_bytes = cls.source.read_bytes()
        cls.saved = json.loads(cls.source_bytes)
        cls.seed = np.asarray(cls.saved['attempts'][0]['endpoint'])
        cls.scales = np.asarray(cls.saved['optimizer_coordinate_scales'])
        cls.model = joint_fit.GlobalLikelihood(1, 3, 4, bin_block_size=4)
        cls.lower, cls.upper = joint_fit.parameter_bounds(cls.model)

    def test_starts_are_deterministic_distinct_and_confined_to_separate_blocks(self):
        first = probes.basin_starts(self.model, self.seed, self.scales)
        self.assertEqual(first, probes.basin_starts(self.model, self.seed, self.scales))
        self.assertEqual([p['name'] for p in first], list(probes.BLOCKS))
        for attempt in first:
            values = np.asarray(attempt['start'])
            multi.validate_start(self.model, values, self.lower, self.upper)
            changed = set(np.flatnonzero(values != self.seed))
            self.assertEqual(changed, set(probes.BLOCKS[attempt['name']]))
            self.assertFalse(attempt['accepted'])
            expected_radius = probes.SCALED_RADIUS * 0.5 ** attempt['halvings']
            self.assertAlmostEqual(attempt['realized_scaled_distance'], expected_radius, places=11)
        self.assertEqual(len({tuple(attempt['start']) for attempt in first}), 3)

    def test_invalid_displacement_is_halved_deterministically(self):
        original = multi.validate_start
        def reject_large(model, values, lower, upper):
            if np.linalg.norm((np.asarray(values)-self.seed)/self.scales) > .2:
                raise ValueError('synthetic physical-domain boundary')
            return original(model, values, lower, upper)
        with patch.object(multi, 'validate_start', side_effect=reject_large):
            starts = probes.basin_starts(self.model, self.seed, self.scales)
        self.assertEqual([attempt['halvings'] for attempt in starts], [1, 1, 1])
        self.assertTrue(all(attempt['realized_scaled_distance'] < .2 for attempt in starts))

    def test_preparation_preserves_source_scales_and_provenance(self):
        ledger = probes.prepare_probes(self.model, self.source, 'stored_negative_c4', 1500)
        self.assertEqual(ledger['probe_seed']['endpoint'], self.seed.tolist())
        self.assertEqual(ledger['optimizer_coordinate_scales'], self.scales.tolist())
        self.assertEqual(ledger['resume_source']['sha256'], hashlib.sha256(self.source_bytes).hexdigest())
        self.assertEqual(ledger['optimizer_settings']['max_nfev'], 1500)
        self.assertEqual(ledger['quadrature'], {'energy': 2, 'pair_per_segment': 16,
                                               'angle_per_segment': 16, 'bin_block_size': 4})
        self.assertEqual(len(ledger['attempts']), 3)
        self.assertIn('input_sha256', ledger)
        self.assertIn('source_sha256', ledger)
        self.assertEqual(self.source.read_bytes(), self.source_bytes)


if __name__ == '__main__':
    unittest.main()
