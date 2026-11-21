import unittest
import sympy as sp
from check_finite_pole_vector_normalization import (
    ALL_GENERATORS, NEUTRAL_GENERATORS, POLE_NORM, charged_source_residual, check_direct_quadratic_form,
    check_field_expansion, check_normalization, check_quadratic_form, normalization_with_free_norm, trace,
    u3_generators)


class TestFinitePoleVectorNormalization(unittest.TestCase):
    def test_generator_normalization(self):
        generators = u3_generators()
        for i in ALL_GENERATORS:
            for j in ALL_GENERATORS:
                self.assertEqual(sp.simplify(trace(generators[i] * generators[j])), 2 if i == j else 0)

    def test_quadratic_form_matches_left_side(self):
        for generator_indices in (ALL_GENERATORS, NEUTRAL_GENERATORS):
            self.assertEqual(check_quadratic_form(generator_indices), 0)

    def test_field_expansion(self):
        for generator_indices in (ALL_GENERATORS, NEUTRAL_GENERATORS):
            self.assertEqual(check_field_expansion(generator_indices), sp.zeros(3))

    def test_left_side_fixes_pole_norm(self):
        for generator_indices in (ALL_GENERATORS, NEUTRAL_GENERATORS):
            _, scaled_left_side, pole_norm_formula = normalization_with_free_norm(generator_indices)
            self.assertEqual(sp.simplify(scaled_left_side - pole_norm_formula), 0)
            self.assertNotIn(POLE_NORM, scaled_left_side.free_symbols)

    def test_normalization_and_zero_axial_source(self):
        for generator_indices in (ALL_GENERATORS, NEUTRAL_GENERATORS):
            self.assertEqual(check_normalization(generator_indices), 0)
            self.assertEqual(check_normalization(generator_indices, axial_is_zero=True), 0)
            self.assertEqual(check_direct_quadratic_form(generator_indices), 0)

    def test_charged_source_caveat(self):
        residual, expected = charged_source_residual()
        self.assertNotEqual(residual, 0)
        self.assertEqual(sp.simplify(residual - expected), 0)


if __name__ == '__main__':
    unittest.main()
