#!/usr/bin/env python3
"""Derive the left-handed weak spurion transformation.

Invocation: .venv/bin/python scripts/weak_spurion_chiral_transformation.py
Tests: .venv/bin/python scripts/weak_spurion_chiral_transformation.py --test
Matrices are dimensionless and unitary where inverted; no expansion is used.
"""

import sys
import unittest

import sympy as sp


def derive():
    """Apply the chiral transformation to each dressed weak matrix."""
    dimension = 3
    u = sp.MatrixSymbol("u", dimension, dimension)
    w = sp.MatrixSymbol("w", dimension, dimension)
    w_dagger = sp.MatrixSymbol("w_dagger", dimension, dimension)
    g_L = sp.MatrixSymbol("g_L", dimension, dimension)
    h = sp.MatrixSymbol("h", dimension, dimension)

    # u' = h u g_L^{-1}, w' = g_L w g_L^{-1}; unitarity
    # identifies each inverse of a group element or u with its adjoint.
    transformed_u = h * u * g_L**-1
    transformed_w = g_L * w * g_L**-1
    transformed_w_dagger = g_L * w_dagger * g_L**-1
    transformed_plus = transformed_u * (transformed_w + transformed_w_dagger) * transformed_u**-1
    transformed_minus = transformed_u * (transformed_w - transformed_w_dagger) * transformed_u**-1
    plus = u * (w + w_dagger) * u**-1
    minus = u * (w - w_dagger) * u**-1
    return transformed_u, transformed_w, plus, minus, transformed_plus, transformed_minus, h


class TransformationTest(unittest.TestCase):
    def test_both_spurions_transform_as_matter(self):
        _, _, plus, minus, transformed_plus, transformed_minus, h = derive()
        self.assertEqual(sp.expand(transformed_plus - h * plus * h**-1), sp.ZeroMatrix(3, 3))
        self.assertEqual(sp.expand(transformed_minus - h * minus * h**-1), sp.ZeroMatrix(3, 3))


if __name__ == "__main__":
    if "--test" in sys.argv:
        sys.argv.remove("--test")
        unittest.main()
    else:
        transformed_u, transformed_w, plus, minus, transformed_plus, transformed_minus, h = derive()
        print("Assumptions: dimensionless unitary chiral matrices; no truncation.")
        for label, expression in (
            ("Transformed u", transformed_u),
            ("Transformed w", transformed_w),
            ("w plus", plus),
            ("w minus", minus),
            ("Transformed w plus", transformed_plus),
            ("Transformed w minus", transformed_minus),
            ("Expanded transformed w plus", sp.expand(transformed_plus)),
            ("Expanded transformed w minus", sp.expand(transformed_minus)),
        ):
            print(f"{label}:\n{sp.latex(expression)}\n")
