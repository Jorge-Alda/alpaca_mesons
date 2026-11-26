"""Exact checks of the compact trace notation and of its expansion into words.

Run:  PYTHONPATH=scripts .venv/bin/python -m unittest test_group_radial_loop_trace_permutations
"""
import itertools
import unittest

import sympy as sp

from group_radial_loop_trace_permutations import (SLOTS, Commutator, PermutationSum, SlotTrace, Trace,
                                                  canonical_word, matrix_words, polynomial_difference,
                                                  trace_polynomial)

X, Y, Z, W, SOURCE = sp.symbols('x y z w S_10', commutative=False)


class NotationChecks(unittest.TestCase):
    def test_trace_cyclicity_and_unit_trace(self):
        self.assertEqual(Trace(X * Y * Z), Trace(Z * X * Y))
        self.assertEqual(Trace(sp.Rational(1, 2) + 0 * X), sp.Rational(3, 2))
        self.assertEqual(Trace(2 * X * Y), 2 * Trace(Y * X))

    def test_permutation_sum_words(self):
        words = matrix_words(PermutationSum(X, Y, Z))
        expected = {tuple(m.name for m in ordering): 1 for ordering in itertools.permutations((X, Y, Z))}
        self.assertEqual(words, expected)
        self.assertEqual(matrix_words(PermutationSum(X, Y)), matrix_words(X * Y + Y * X))

    def test_permutation_sum_with_repeated_and_scalar_entries(self):
        self.assertEqual(matrix_words(PermutationSum(X, X, Y)), {('x', 'x', 'y'): 2, ('x', 'y', 'x'): 2, ('y', 'x', 'x'): 2})
        self.assertEqual(PermutationSum(X, Y, 1), 3 * PermutationSum(X, Y))

    def test_commutator(self):
        self.assertEqual(Commutator(X, X), 0)
        self.assertEqual(Commutator(Y, X), -Commutator(X, Y))
        self.assertEqual(matrix_words(Commutator(X, Y) * Z), {('x', 'y', 'z'): 1, ('y', 'x', 'z'): -1})

    def test_trace_polynomial_uses_cyclic_words(self):
        polynomial = trace_polynomial(Trace(SOURCE * PermutationSum(X, Y)) * Trace(Commutator(X, Y) * Z))
        expected = {}
        for first in (('S_10', 'x', 'y'), ('S_10', 'y', 'x')):
            for second, sign in ((('x', 'y', 'z'), 1), (('y', 'x', 'z'), -1)):
                key = tuple(sorted((canonical_word(first), canonical_word(second))))
                expected[key] = expected.get(key, 0) + sign
        self.assertEqual(polynomial_difference(polynomial, expected), {})

    def test_slot_trace(self):
        word = Commutator(SLOTS[0], X) * SLOTS[1]
        explicit = Trace(Commutator(Y, X) * Z) + Trace(Commutator(Z, X) * Y)
        self.assertEqual(polynomial_difference(trace_polynomial(SlotTrace(word, Y, Z)), trace_polynomial(explicit)), {})


if __name__ == '__main__':
    unittest.main()
