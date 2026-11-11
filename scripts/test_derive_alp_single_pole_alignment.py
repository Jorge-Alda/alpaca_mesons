"""Tests for single-pole alignment of the ALP-meson inverse propagator."""

import unittest

from derive_alp_single_pole_alignment import (
    decoupled_meson_example,
    regular_pole_checks,
    singular_meson_block_determinant,
)


class SinglePoleAlignmentTest(unittest.TestCase):
    def test_regular_poles(self):
        for entry in regular_pole_checks():
            for label, passed in entry.items():
                if isinstance(passed, bool):
                    self.assertTrue(passed, f"pole {entry['pole']}: {label}")

    def test_singular_meson_block_needs_decoupled_direction(self):
        _, passed = singular_meson_block_determinant()
        self.assertTrue(passed)

    def test_pole_without_alp_component(self):
        result = decoupled_meson_example(degenerate=False)
        self.assertTrue(result["meson block singular"])
        self.assertEqual(result["null space dimension"], 1)
        self.assertTrue(result["meson rows solvable with free w_3"])
        self.assertTrue(result["no shift aligns this pole"])

    def test_degenerate_pole(self):
        result = decoupled_meson_example(degenerate=True)
        self.assertTrue(result["meson block singular"])
        self.assertEqual(result["null space dimension"], 2)
        self.assertTrue(result["shift aligns for every w_3"])
        self.assertEqual(result["residue rank"], 2)
        self.assertTrue(result["residue is sum over orthonormal pole vectors"])


if __name__ == "__main__":
    unittest.main()
