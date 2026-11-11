"""Tests for finite CCWZ source covariance."""

import unittest

from derive_alp_ccwz_redefinitions import axial_map_expansion, derive


class CcwzRedefinitionsTest(unittest.TestCase):
    def test_covariance_and_invariant_traces(self):
        *_, trace_vielbein, trace_mass = derive()
        self.assertEqual(trace_vielbein, 0)
        self.assertEqual(trace_mass, 0)

    def test_axial_map_first_nonlinear_term(self):
        quadratic, cubic = axial_map_expansion()
        self.assertEqual(quadratic, 0)
        self.assertNotEqual(cubic, 0)


if __name__ == "__main__":
    unittest.main()
