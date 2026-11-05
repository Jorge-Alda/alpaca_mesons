"""Check joint vector continuation selection without running an optimizer.

Run: .venv/bin/python -m unittest scripts.test_select_joint_vector_continuation
"""
import copy
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

from scripts import select_joint_vector_continuation as select


ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / 'data/vector_sector/joint_vector_multistart.json'


class ContinuationSelectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ledger_bytes = LEDGER.read_bytes()
        cls.ledger = json.loads(cls.ledger_bytes)
        cls.ledger_hash = hashlib.sha256(cls.ledger_bytes).hexdigest()

    def test_selection_ranking_distances_and_status(self):
        result = select.compare_endpoints(self.ledger, str(LEDGER), self.ledger_hash)
        self.assertEqual(result['source_ledger']['sha256'], self.ledger_hash)
        self.assertEqual(result['continuation_seed']['name'], 'stored_negative_c4')
        selected = next(a for a in self.ledger['attempts'] if a['name'] == 'stored_negative_c4')
        self.assertEqual(result['continuation_seed']['endpoint'], selected['endpoint'])
        self.assertIs(result['continuation_seed']['success'], False)
        self.assertIs(result['continuation_seed']['accepted'], False)
        self.assertEqual([entry['name'] for entry in result['objective_ranking']],
                         ['stored_negative_c4', 'stored', 'stored_negative_odd', 'preliminary'])
        self.assertEqual(len(result['pairwise_distances']), 6)
        first = result['pairwise_distances'][0]
        by_name = {a['name']: np.asarray(a['endpoint']) for a in self.ledger['attempts']}
        expected = np.linalg.norm((by_name[first['first']] - by_name[first['second']]) /
                                  np.asarray(self.ledger['optimizer_coordinate_scales']))
        self.assertEqual(first['scaled_distance'], expected)

    def test_rejects_invalid_ambiguous_and_changed_selection(self):
        invalid = copy.deepcopy(self.ledger)
        invalid['attempts'][0]['endpoint_diagnostics']['valid'] = False
        with self.assertRaisesRegex(ValueError, 'not physically valid'):
            select.compare_endpoints(invalid, 'ledger', 'hash')
        invalid = copy.deepcopy(self.ledger)
        invalid['optimizer_coordinate_scales'][0] = 0
        with self.assertRaisesRegex(ValueError, 'positive scales'):
            select.compare_endpoints(invalid, 'ledger', 'hash')
        invalid = copy.deepcopy(self.ledger)
        invalid['attempts'][0]['chi2'] = 0
        with self.assertRaisesRegex(ValueError, 'not the lowest'):
            select.compare_endpoints(invalid, 'ledger', 'hash')

    def test_command_reproducible_and_source_unchanged(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / 'selection.json'
            command = [sys.executable, str(ROOT / 'scripts/select_joint_vector_continuation.py'),
                       '--ledger', str(LEDGER), '--output', str(output)]
            subprocess.run(command, check=True, capture_output=True, text=True)
            first = output.read_bytes()
            subprocess.run(command, check=True, capture_output=True, text=True)
            self.assertEqual(output.read_bytes(), first)
            self.assertEqual(json.loads(first)['source_ledger']['sha256'], self.ledger_hash)
        self.assertEqual(LEDGER.read_bytes(), self.ledger_bytes)


if __name__ == '__main__':
    unittest.main()
