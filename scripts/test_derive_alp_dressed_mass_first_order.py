"""Unit tests for the first-order dressed quark mass source.

Run: .venv/bin/python -m unittest scripts.test_derive_alp_dressed_mass_first_order
"""

import unittest

from scripts.derive_alp_dressed_mass_first_order import (
    check_constant_redefinition_invariance,
    check_general_and_alda_form,
    check_hermitian_split,
    check_kernel_structure,
    check_numeric_rational,
    check_representative_consistency,
)


class DressedMassFirstOrderTests(unittest.TestCase):
    def test_general_and_alda_form(self):
        result = check_general_and_alda_form()
        for key in ("general_zero_residual", "general_one_residual",
                    "alda_zero_residual", "alda_one_residual"):
            self.assertEqual(result[key], 0, key)

    def test_hermitian_split(self):
        result = check_hermitian_split()
        self.assertEqual(result["scalar_residual"], 0)
        self.assertEqual(result["pseudoscalar_residual"], 0)

    def test_constant_redefinition_invariance(self):
        for key, value in check_constant_redefinition_invariance().items():
            self.assertEqual(value, 0, key)

    def test_numeric_rational(self):
        for seed in (20261111, 7, 12345):
            self.assertTrue(check_numeric_rational(seed=seed)["residual_is_zero"])

    def test_kernel_structure(self):
        result = check_kernel_structure()
        self.assertTrue(result["diagonal_scalar_one_is_zero"])
        self.assertEqual(result["commutator_kernel_nondegenerate"], 3)
        self.assertEqual(result["commutator_kernel_isospin_degenerate"], 5)
        self.assertEqual(result["commutator_kernel_generic_hermitian"], 3)
        self.assertEqual(result["anticommutator_kernel_nondegenerate"], 0)
        self.assertEqual(result["anticommutator_kernel_isospin_degenerate"], 0)
        self.assertEqual(result["anticommutator_kernel_generic_hermitian"], 0)

    def test_representative_consistency(self):
        result = check_representative_consistency()
        self.assertEqual(result["representative_residual"], 0)
        self.assertEqual(result["explicit_shift_residual"], 0)
        self.assertTrue(result["independent_of_axial_hat"])


if __name__ == "__main__":
    unittest.main()
