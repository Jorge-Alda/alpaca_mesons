"""Tests for the alignment of several ALP-meson poles."""

import unittest

import sympy as sp

from derive_alp_multi_pole_alignment import (
    floating_point_pencil_check,
    general_conversion_checks,
    indefinite_kinetic_counterexample,
    shift_cannot_diagonalize,
    shift_type_two_pole_counterexample,
)


def boolean_entries(result):
    return {label: value for label, value in result.items() if isinstance(value, bool)}


class MultiPoleAlignmentTest(unittest.TestCase):
    def test_shift_type_counterexample_for_every_pair(self):
        for first in range(4):
            for second in range(first + 1, 4):
                result = shift_type_two_pole_counterexample(first, second)
                for label, passed in boolean_entries(result).items():
                    self.assertTrue(passed, f"poles {first}, {second}: {label}")

    def test_shift_cannot_diagonalize(self):
        for label, passed in boolean_entries(shift_cannot_diagonalize()).items():
            self.assertTrue(passed, label)

    def test_general_conversion(self):
        for label, passed in boolean_entries(general_conversion_checks()).items():
            self.assertTrue(passed, label)

    def test_indefinite_kinetic_matrix_has_double_pole(self):
        result = indefinite_kinetic_counterexample()
        self.assertEqual(result["v^(k)T K v^(k)"], 0)
        self.assertEqual(result["coefficient of 1/s^2 in D^{-1}"], sp.Matrix([[0, 0], [0, 1]]))

    def test_floating_point_pencil(self):
        result = floating_point_pencil_check()
        self.assertLess(result["max |T^T K T - 1|"], 1e-12)
        self.assertLess(result["max |T^T M^2 T - diag|"], 1e-12)


if __name__ == "__main__":
    unittest.main()
