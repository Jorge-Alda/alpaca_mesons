"""Tests of the goodness-of-fit assessment at accepted minima (low order).

Run: .venv/bin/python -m unittest scripts.test_assess_joint_vector_fit_adequacy
"""
import json
import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import assess_joint_vector_fit_adequacy as adequacy
import fit_babar_kkpi_global_vector_sector as joint_fit


class AdequacyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.comparison = json.loads((adequacy.ROOT/adequacy.COMPARISON).read_text())
        cls.model = joint_fit.GlobalLikelihood(1, 3, 4, bin_block_size=4)

    def test_accepted_endpoints_are_the_four_ranked_ones(self):
        accepted = adequacy.load_accepted_endpoints(self.comparison)
        self.assertEqual(len(accepted), 4)
        self.assertEqual(accepted[0]['id'], 'basin_probes/even_shared_widths')
        self.assertTrue(all(candidate['accepted'] for candidate in accepted))

    def test_refuses_non_accepted_and_unknown_endpoints(self):
        with self.assertRaisesRegex(ValueError, 'not accepted'):
            adequacy.load_accepted_endpoints(self.comparison, ['continuation/stored_negative_c4'])
        with self.assertRaisesRegex(ValueError, 'not in the endpoint'):
            adequacy.load_accepted_endpoints(self.comparison, ['missing/endpoint'])

    def test_assess_requires_common_order(self):
        with self.assertRaises(ValueError):
            adequacy.assess(self.model, self.comparison)

    def test_record_contents_and_block_sums(self):
        labels = adequacy.row_labels(self.model)
        candidate = dict(adequacy.load_accepted_endpoints(self.comparison)[0])
        residual = self.model.residual(np.asarray(candidate['endpoint']))
        candidate['common_chi2'] = float(residual @ residual)
        record = adequacy.assess_endpoint(self.model, candidate, labels, with_local_errors=False)
        self.assertTrue(record['chi2_matches_comparison'])
        self.assertEqual(record['degrees_of_freedom'], 43)
        self.assertAlmostEqual(sum(block['chi2'] for block in record['blocks'].values()), record['chi2'])
        self.assertEqual([block['count'] for block in record['blocks'].values()], [12, 7, 7, 2, 45, 2])
        self.assertEqual(len(record['standardized_residuals']), 75)
        babar = [row for row in record['standardized_residuals'] if row['block'] == 'BaBar']
        self.assertEqual(len(babar), 45)
        self.assertTrue(all('channel' in row and 'energy_center_gev' in row for row in babar))
        pulls = [abs(row['pull']) for row in record['largest_pulls']]
        self.assertEqual(pulls, sorted(pulls, reverse=True))
        pattern = record['babar_energy_pattern']
        self.assertEqual(sum(part['count'] for part in pattern['kkpi0']['energy_ranges']), 21)
        self.assertEqual(sum(part['count'] for part in pattern['kskpi']['energy_ranges']), 24)

    def test_runs_statistics(self):
        result = adequacy.runs_statistics([1, 1, -1, -1, 1, -1])
        self.assertEqual(result['runs'], 4)
        self.assertAlmostEqual(result['expected_runs'], 1+2*3*3/6)
        self.assertIsNone(adequacy.runs_statistics([1, 1, 1])['expected_runs'])


if __name__ == '__main__':
    unittest.main()
