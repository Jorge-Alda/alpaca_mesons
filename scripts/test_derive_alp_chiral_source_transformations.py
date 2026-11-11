"""Tests for the finite chiral source calculation.

Run: .venv/bin/python -m unittest discover -s scripts -p test_derive_alp_chiral_source_transformations.py
"""
import unittest

import sympy as sp

from derive_alp_chiral_source_transformations import (
    check_series, derive_anomaly, derive_source_laws, x,
)


class ChiralSourceTests(unittest.TestCase):
    def test_noncommuting_series(self):
        for name, residual in check_series().items():
            with self.subTest(name=name):
                self.assertEqual(residual, sp.zeros(2))

    def test_product_rule(self):
        self.assertEqual(derive_source_laws()['product rule residual'], 0)

    def test_singlet_jacobian_and_topological_source(self):
        result = derive_anomaly()
        flavor_number = sp.Symbol('N_f', integer=True, positive=True)
        theta = sp.Function('theta')(x)
        self.assertEqual(result['Jacobian shift of Theta'], -flavor_number * theta)
        self.assertEqual(result['shift of tr log u'], sp.I * flavor_number * theta)
        self.assertEqual(result['X_TC residual'], 0)
        self.assertEqual(result['normalization residual'], 0)


if __name__ == '__main__':
    unittest.main()
