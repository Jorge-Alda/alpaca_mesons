"""Unit tests for derive_murayama_radial_frame_invariance.py.

Run: PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_murayama_radial_frame_invariance
"""
import unittest

import sympy as sp

from derive_murayama_radial_frame_invariance import (
    axial_exponential, contraction_checks, cyclic_trace, general_original_frame,
    series_product, series_sum, source_invariance, transformed_frame, transport)


class TransportTest(unittest.TestCase):
    def test_cyclic_trace(self):
        first, second, third = sp.symbols('A B C', commutative=False)
        self.assertEqual(cyclic_trace(first * second * third - third * first * second), 0)
        self.assertNotEqual(cyclic_trace(first * second * third - second * first * third), 0)

    def test_transport_rule_at_fixed_sigma(self):
        """W'_R = g_R W_R order by order in y and x, with g independent of sigma."""
        radial_order, alp_order = 2, 3
        generator = sp.Symbol('kappa', commutative=False)
        original = general_original_frame(radial_order, alp_order, with_connections=True)
        _, new_right, _, _ = transformed_frame(*original, generator, radial_order, alp_order)
        old_transport = transport(original[1], radial_order, alp_order)
        new_transport = transport(new_right, radial_order, alp_order)
        expected = series_product(axial_exponential(generator, +1, alp_order), old_transport,
                                  radial_order=radial_order, alp_order=alp_order)
        residual = series_sum(new_transport, {k: -v for k, v in expected.items()})
        self.assertTrue(all(value == 0 for value in residual.values()))


class DressedSourceTest(unittest.TestCase):
    def test_general_sources_and_connections(self):
        residuals = source_invariance(radial_order=2, alp_order=3, with_connections=True)
        self.assertEqual(len(residuals), 36)
        self.assertTrue(all(value == 0 for value in residuals.values()))

    def test_connection_free_sixth_order(self):
        residuals = source_invariance(radial_order=2, alp_order=6, with_connections=False)
        self.assertTrue(all(value == 0 for value in residuals.values()))


class ContractionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dressed_differences = contraction_checks()

    def test_dressed_contractions_and_pole_data(self):
        for name, value in self.dressed_differences.items():
            values = value.values() if isinstance(value, dict) else [value]
            for entry in values:
                self.assertEqual(entry, 0, name)


if __name__ == '__main__':
    unittest.main()
