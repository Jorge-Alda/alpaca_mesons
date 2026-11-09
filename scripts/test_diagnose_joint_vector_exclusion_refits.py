"""Focused checks of eta-prime row removal and diagnostic reporting."""
import unittest
import numpy as np

from scripts import diagnose_joint_vector_exclusion_refits as diagnostic
from scripts import joint_vector_likelihood_options as options
from scripts import refit_joint_vector_sensitivity as refits


class ExclusionDiagnosticTests(unittest.TestCase):
    def test_eta_prime_mask_removes_only_active_radiative_rows(self):
        model = options.VariedLikelihood(1, 3, 4, bin_block_size=4)
        indices = diagnostic.exclude_eta_prime_radiative_rows(model)
        self.assertEqual(indices, [4, 11])
        self.assertEqual(np.flatnonzero(~model.row_mask).tolist(), indices)
        self.assertEqual([model.observations[index]['channel'] for index in indices],
                         ['etap_omega_gamma', 'etap_gamma_gamma'])
        self.assertTrue(model.row_mask[21])  # eta-prime total-width observation stays

    def test_signed_row_changes(self):
        evaluation = {'full_residuals': [3., -2.]}
        baseline = {'standardized_residuals': [{'pull': 1.}, {'pull': -3.}]}
        labels = [{'label': 'first'}, {'label': 'second'}]
        changes = diagnostic.row_changes(evaluation, baseline, labels)
        self.assertEqual([row['pull_change'] for row in changes], [2., 1.])

    def test_mass_exclusion_has_optimizer_status_and_high_order_rows(self):
        record = diagnostic.mass_exclusion_record()
        self.assertEqual(record['start_endpoint'], refits.BEST)
        self.assertIn('status', record['refit'])
        self.assertEqual([row['index'] for row in record['high_order']['removed_row_pulls']],
                         list(range(12, 19)))
        self.assertTrue(record['high_order']['valid'])

    def test_small_eta_prime_refit_records_status_and_removed_rows(self):
        context = refits.accepted_context()
        models = [options.VariedLikelihood(1, 3, 4, bin_block_size=4) for _ in range(2)]
        for model in models:
            diagnostic.exclude_eta_prime_radiative_rows(model)
        item = {'name': 'without_eta_prime_radiative_rows', 'excluded_blocks': ()}
        record = refits.run_exclusion(item, context, models=models,
                                      settings=dict(refits.OPTIMIZER_SETTINGS, max_nfev=2))
        self.assertEqual(record['refit']['state'], 'finished')
        self.assertIn('status', record['refit'])
        self.assertIn('success', record['refit'])
        self.assertEqual([row['index'] for row in record['high_order']['removed_row_pulls']], [4, 11])
        self.assertEqual(record['high_order']['blocks']['odd']['count'], 10)
        self.assertEqual(record['high_order']['blocks']['shared_widths']['count'], 7)


if __name__ == '__main__':
    unittest.main()
