"""Tests of the Markdown summary of the joint vector-fit adequacy ledgers.

Run: .venv/bin/python -m unittest scripts.test_summarize_joint_vector_fit_adequacy
"""
import json
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import summarize_joint_vector_fit_adequacy as summary_module


class SummaryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ledgers = [json.loads((summary_module.ROOT/path).read_text())
                       for path in (summary_module.ADEQUACY, summary_module.EXCLUSION, summary_module.VARIATIONS)]

    def test_p_value_format(self):
        self.assertEqual(summary_module.format_p_value(1.2670481319e-26), '1.267e-26')
        self.assertEqual(summary_module.format_p_value(None), 'n/a')

    def test_function_evaluation_limit_flag(self):
        self.assertTrue(summary_module.hit_function_evaluation_limit(dict(cumulative_nfev=5000, message='')))
        self.assertTrue(summary_module.hit_function_evaluation_limit(
            dict(nfev=12, message='The maximum number of function evaluations is exceeded.')))
        self.assertFalse(summary_module.hit_function_evaluation_limit(
            dict(nfev=1282, message='`ftol` termination condition is satisfied.')))

    def test_worst_block_is_smallest_p_value(self):
        blocks = dict(a=dict(chi2=1., count=1, nominal_p_value_row_count_dof=.3),
                      b=dict(chi2=9., count=2, nominal_p_value_row_count_dof=.01))
        self.assertEqual(summary_module.worst_block(blocks)[0], 'b')

    def test_sign_subsequences_match_ledger_counts(self):
        for endpoint in self.ledgers[0]['endpoints']:
            for pattern in endpoint['babar_energy_pattern'].values():
                signs = summary_module.sign_pattern_by_range(pattern)
                self.assertEqual(''.join(signs), pattern['sign_sequence'])
                for item, sequence in zip(pattern['energy_ranges'], signs):
                    self.assertEqual((sequence.count('+'), sequence.count('-')),
                                     (item['positive'], item['negative']))

    def test_summary_is_deterministic_and_complete(self):
        first = summary_module.summary(*self.ledgers)
        self.assertEqual(first, summary_module.summary(*self.ledgers))
        for entry in self.ledgers[2]['entries']:
            self.assertIn(f"| {entry['name']} |", first)
        for entry in self.ledgers[1]['entries']:
            self.assertIn(f"| {entry['excluded_blocks'][0]} |", first)
        self.assertIn('225.726699', first)


if __name__ == '__main__':
    unittest.main()
