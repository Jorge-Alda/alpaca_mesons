"""Tests for common-order comparison without a high-order calculation.

Run: .venv/bin/python -m unittest scripts.test_compare_joint_vector_endpoints
"""
import copy
import hashlib
import sys
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_joint_vector_endpoints as comparison
import fit_babar_kkpi_global_vector_sector as joint_fit


class FakeCommonOrderModel:
    """Use the real physical-domain check and a cheap common-order objective."""
    def __init__(self):
        self.model = joint_fit.GlobalLikelihood(1, 3, 4, bin_block_size=4)
        self.phase = SimpleNamespace(orders=(4, 128, 128))
        self.bin_block_size = 4
        self.odd_data = self.model.odd_data
        self.branch_sign = self.model.branch_sign
        self.evaluated = []

    def unpack(self, values):
        return self.model.unpack(values)

    def residual(self, values):
        self.evaluated.append(values.copy())
        return np.r_[values, np.zeros(43)]


class EndpointComparisonTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ledgers, cls.hashes = comparison.read_ledgers(comparison.DEFAULT_INPUTS)

    def test_common_objective_ranking_distances_and_acceptance(self):
        model = FakeCommonOrderModel()
        result = comparison.compare(self.ledgers, self.hashes, model)
        candidates = {candidate['id']: candidate for candidate in result['candidates']}
        self.assertEqual(len(candidates), 10)
        self.assertEqual(len(model.evaluated), 10)
        for candidate in candidates.values():
            vector = np.asarray(candidate['endpoint'])
            self.assertAlmostEqual(candidate['common_chi2'], float(vector @ vector))
            self.assertEqual(candidate['accepted'],
                             candidate['optimizer_success'] and candidate['endpoint_valid']
                             and candidate['quadrature_certified'])
        self.assertEqual(sum(candidate['accepted'] for candidate in candidates.values()), 4)
        self.assertFalse(candidates['continuation/stored_negative_c4']['accepted'])
        self.assertTrue(candidates['continuation_4128128/stored_negative_c4']['accepted'])
        self.assertEqual(result['ranking'], [candidate['id'] for candidate in sorted(
            candidates.values(), key=lambda candidate: (candidate['common_chi2'], candidate['id']))])
        self.assertEqual(len(result['pairwise_comparisons']), 45)
        pair = result['pairwise_comparisons'][0]
        first, second = candidates[pair['first']], candidates[pair['second']]
        expected = np.linalg.norm((np.asarray(first['endpoint']) - np.asarray(second['endpoint'])) /
                                  np.asarray(result['optimizer_coordinate_scales']))
        self.assertEqual(pair['scaled_distance'], expected)
        self.assertEqual(pair['within_both_tolerances'],
                         pair['scaled_distance'] <= .01 and pair['absolute_chi2_difference'] <= .01)
        self.assertEqual(result['comparison_tolerances'], {'scaled_distance': .01, 'absolute_chi2_difference': .01})

    def test_acceptance_requires_each_audit_and_optimizer_success(self):
        changed = copy.deepcopy(self.ledgers)
        changed['basin_probes']['attempts'][0]['quadrature_refinement']['passed'] = False
        changed['basin_probes']['attempts'][1]['success'] = False
        changed['basin_probes']['attempts'][1]['status'] = 0
        changed['continuation_audit']['attempts'][0]['quadrature_refinement']['refinements']['angle']['passed'] = False
        result = comparison.compare(changed, self.hashes, FakeCommonOrderModel())
        candidates = {candidate['id']: candidate for candidate in result['candidates']}
        self.assertFalse(candidates['basin_probes/even_shared_widths']['accepted'])
        self.assertFalse(candidates['basin_probes/odd']['accepted'])
        self.assertFalse(candidates['continuation_4128128/stored_negative_c4']['accepted'])
        self.assertTrue(candidates['basin_probes/babar']['accepted'])

    def test_rejects_inconsistent_scales_domain_and_provenance(self):
        for modification in ('scale', 'domain', 'lineage', 'names'):
            changed = copy.deepcopy(self.ledgers)
            if modification == 'scale':
                changed['basin_probes']['optimizer_coordinate_scales'][0] *= 2
            elif modification == 'domain':
                changed['basin_probes']['attempts'][0]['endpoint'][26] = 6
            elif modification == 'lineage':
                changed['basin_probes']['resume_source']['sha256'] = 'wrong'
            else:
                changed['basin_probes']['attempts'][0]['name'] = 'unknown'
            with self.subTest(modification=modification), self.assertRaises(ValueError):
                comparison.compare(changed, self.hashes, FakeCommonOrderModel())
        source = Path(comparison.DEFAULT_INPUTS['multistart'])
        self.assertEqual(self.hashes['multistart']['sha256'], hashlib.sha256(source.read_bytes()).hexdigest())


if __name__ == '__main__':
    unittest.main()
