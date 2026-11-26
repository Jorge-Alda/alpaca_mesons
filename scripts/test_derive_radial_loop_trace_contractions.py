"""Exact checks of the U(3) completeness identities, of their compact application, and of the
vertex rules against the Lagrangian expansion.

The completeness test is the only place where explicit 3 x 3 matrices appear; nothing with
components is printed.

Run:  PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_radial_loop_trace_contractions
"""
import itertools
import unittest

import sympy as sp

import derive_radial_loop_trace_contractions as common
import expand_radial_loop_lagrangian as lagrangian
from group_radial_loop_trace_permutations import (Commutator, PermutationSum, Trace, polynomial_difference,
                                                  trace_polynomial)


def unitary_generators():
    """Nine Hermitian U(3) generators with <lambda_i lambda_j> = 2 delta_ij, lambda_0 = sqrt(2/3) 1."""
    i = sp.I
    generators = [sp.sqrt(sp.Rational(2, 3)) * sp.eye(3),
                  sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]]), sp.Matrix([[0, -i, 0], [i, 0, 0], [0, 0, 0]]),
                  sp.Matrix([[1, 0, 0], [0, -1, 0], [0, 0, 0]]), sp.Matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]]),
                  sp.Matrix([[0, 0, -i], [0, 0, 0], [i, 0, 0]]), sp.Matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]]),
                  sp.Matrix([[0, 0, 0], [0, 0, -i], [0, i, 0]]),
                  sp.diag(1, 1, -2) / sp.sqrt(3)]
    return generators


def general_matrix(name):
    return sp.Matrix(3, 3, sp.symbols(f'{name}_1:10'))


class CompletenessIdentities(unittest.TestCase):
    def setUp(self):
        self.generators = unitary_generators()
        self.first, self.second = general_matrix('a'), general_matrix('b')

    def test_normalization(self):
        for (k, first), (l, second) in itertools.product(enumerate(self.generators), repeat=2):
            self.assertEqual(sp.simplify((first * second).trace()), 2 if k == l else 0)

    def test_single_trace_identity(self):
        total = sum(((self.first * g * self.second * g).trace() for g in self.generators), sp.S.Zero)
        self.assertEqual(sp.expand(total - 2 * self.first.trace() * self.second.trace()), 0)

    def test_double_trace_identity(self):
        total = sum(((self.first * g).trace() * (g * self.second).trace() for g in self.generators), sp.S.Zero)
        self.assertEqual(sp.expand(total - 2 * (self.first * self.second).trace()), 0)

    def test_singlet_is_needed(self):
        octet = sum(((self.first * g * self.second * g).trace() for g in self.generators[1:]), sp.S.Zero)
        self.assertNotEqual(sp.expand(octet - 2 * self.first.trace() * self.second.trace()), 0)

    def test_no_floating_point_numbers(self):
        for generator in self.generators:
            self.assertFalse(generator.atoms(sp.Float))


class CompactContraction(unittest.TestCase):
    """contract_generator on compact traces equals the word rewrite <X L><L Y> -> 2<XY>."""

    def check(self, first_trace, second_trace):
        generator = GENERATOR
        compact = common.contract_generator(first_trace * second_trace, generator)
        words = trace_polynomial(first_trace * second_trace)
        expected = lagrangian.apply_completeness(words, generator.name)
        self.assertEqual(polynomial_difference(expected, trace_polynomial(compact)), {})
        self.assertFalse(compact.has(generator))

    def test_permutation_sums(self):
        self.check(Trace(SOURCE * PermutationSum(GENERATOR, X, Y)), Trace(OTHER * PermutationSum(GENERATOR, Z)))

    def test_commutators(self):
        self.check(Trace(Commutator(X, Y) * Commutator(Z, GENERATOR)),
                   Trace(Commutator(GENERATOR, W) * PermutationSum(SOURCE, X)))

    def test_nested_commutators(self):
        self.check(Trace(Commutator(Commutator(X, GENERATOR), Y) * Z), Trace(OTHER * PermutationSum(GENERATOR, W, X)))

    def test_nested_permutation_sum(self):
        self.check(Trace(PermutationSum(X, PermutationSum(GENERATOR, Y)) * SOURCE), Trace(GENERATOR * Z))


GENERATOR, SOURCE, OTHER, X, Y, Z, W = sp.symbols('G S_00 P_11 x y z w', commutative=False)


class VertexRules(unittest.TestCase):
    """Vertex rules of the compact module against the direct expansion of the Lagrangian."""

    def compare(self, compact_legs, word_legs, extra=None):
        compact = common.vertex(compact_legs, common.ma**2, extra)
        expected = lagrangian.vertex_polynomial(tuple(word_legs), common.ma**2, extra)
        self.assertEqual(polynomial_difference(expected, trace_polynomial(compact)), {})

    def general_legs(self, count, angular=True):
        matrices = sp.symbols(f'm_1:{count + 1}', commutative=False)
        components = sp.symbols(f'c_1:{count + 1}', real=True) if angular else (0,) * count
        labels = [f'q{k}' for k in range(count)]
        dots = {frozenset(pair): sp.Symbol(f'{pair[0]}.{pair[1]}') for pair in itertools.combinations(labels, 2)}
        compact = [(m, c, label) for m, c, label in zip(matrices, components, labels)]
        words = [(m.name, c, 0, label) for m, c, label in zip(matrices, components, labels)]
        return compact, words, dots

    def test_four_general_pole_legs(self):
        self.compare(*self.general_legs(4))

    def test_six_meson_legs(self):
        self.compare(*self.general_legs(6, angular=False))

    def test_radial_legs(self):
        compact, words, dots = self.general_legs(3)
        for radial_count in (1, 2):
            self.compare([None] * radial_count + compact[:2], [lagrangian.RADIAL] * radial_count + words[:2], dots)
            self.compare([None] * radial_count + compact, [lagrangian.RADIAL] * radial_count + words, dots)


if __name__ == '__main__':
    unittest.main()
