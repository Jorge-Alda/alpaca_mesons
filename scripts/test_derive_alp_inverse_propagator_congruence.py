"""Tests for the congruence of the ALP-meson inverse propagator."""

import unittest

from derive_alp_inverse_propagator_congruence import derive


class InversePropagatorCongruenceTest(unittest.TestCase):
    def test_all_checks_for_several_field_changes(self):
        for seed in (20261111, 7, 123):
            result = derive(seed=seed)
            for label, passed in result["checks"].items():
                self.assertTrue(passed, f"seed {seed}: {label}")
            for entry in result["per_pole"]:
                for label, passed in entry.items():
                    if isinstance(passed, bool):
                        self.assertTrue(passed, f"seed {seed}, pole {entry['pole']}: {label}")


if __name__ == "__main__":
    unittest.main()
