#!/usr/bin/env python3
"""Compact flavour-trace notation for the loop diagram products, and its expansion into words.

Objects (all exact SymPy expressions, matrices are noncommuting symbols):
    Trace(X)                 flavour trace <X>; cyclic rotations are identified and <1> = 3
    PermutationSum(X1..Xm)   sum over the m! orderings of the listed factors, X1...Xm included;
                             for m = 2 this is the anticommutator {X1, X2}
    Commutator(X, Y)         [X, Y] = XY - YX
    SlotTrace(word, M1..Mk)  sum over the k! orderings of M1..Mk placed in the k positions of
                             `word` held by the symbols in SLOTS (printed as a bullet)
No matrix indices or matrix components are used.  `trace_polynomial` expands any expression
built from these objects back into a sum of products of traces of plain words (tuples of
symbol names), with every word rotated to its cyclic representative; it is the reference
expansion used by the exact tests.

Usage (library):
    from group_radial_loop_trace_permutations import Trace, PermutationSum, trace_polynomial
    X, Y, Z = sympy.symbols('X Y Z', commutative=False)
    trace_polynomial(Trace(X*PermutationSum(Y, Z)))   # {(('X','Y','Z'),): 1, (('X','Z','Y'),): 1}
"""
import itertools
import math
from collections import defaultdict

import sympy as sp

# Number of quark flavours: the trace of the unit matrix.
FLAVOURS = 3
SLOTS = sp.symbols('_slot_1:7', commutative=False)


def sort_key(expression):
    """Deterministic ordering: plain symbols before composite matrix expressions."""
    return (0 if expression.is_Symbol else 1, str(expression))


def word_factors(expression):
    """Split a product into its commutative coefficient and its list of matrix factors."""
    commutative, noncommutative = expression.args_cnc()
    factors = []
    for factor in noncommutative:
        if factor.is_Pow and factor.exp.is_Integer and factor.exp > 0:
            factors.extend([factor.base] * int(factor.exp))
        else:
            factors.append(factor)
    return sp.Mul(*commutative), factors


def cyclic_representative(factors):
    """Rotation of a list of matrix factors that is smallest in the deterministic ordering."""
    if not factors:
        return []
    rotations = [factors[start:] + factors[:start] for start in range(len(factors))]
    return min(rotations, key=lambda rotation: tuple(sort_key(factor) for factor in rotation))


class Trace(sp.Function):
    """Flavour trace <X>: linear, cyclic, and <c 1> = 3 c."""
    is_commutative = True

    @classmethod
    def eval(cls, argument):
        if argument.is_Add:
            return sp.Add(*(cls(term) for term in argument.args))
        if argument == 0:
            return sp.S.Zero
        coefficient, factors = word_factors(argument)
        if not factors:
            return FLAVOURS * argument
        representative = cyclic_representative(factors)
        if coefficient != 1 or representative != factors:
            return coefficient * cls(sp.Mul(*representative))
        return None

    def _latex(self, printer):
        return r'\left\langle ' + printer._print(self.args[0]) + r'\right\rangle'


class PermutationSum(sp.Function):
    """Sum over all orderings of the listed matrix factors (labelled positions, m! terms)."""

    @classmethod
    def eval(cls, *arguments):
        coefficient = sp.S.One
        matrices = []
        scalars = 0
        for argument in arguments:
            if argument == 0:
                return sp.S.Zero
            factor, factors = word_factors(argument)
            coefficient *= factor
            if factors:
                matrices.append(sp.Mul(*factors))
            else:
                scalars += 1
        # A scalar entry contributes the same word in each of its m+1 positions.
        coefficient *= sp.Integer(math.factorial(len(matrices) + scalars) // math.factorial(len(matrices)))
        if len(matrices) <= 1:
            return coefficient * (matrices[0] if matrices else sp.S.One)
        ordered = sorted(matrices, key=sort_key)
        if coefficient != 1 or tuple(ordered) != tuple(arguments):
            return coefficient * cls(*ordered)
        return None

    def _latex(self, printer, exp=None):
        entries = [printer._print(argument) for argument in self.args]
        if len(entries) == 2:
            text = r'\left\{' + entries[0] + ',' + entries[1] + r'\right\}'
        else:
            text = r'\left(\sum_{\rm perm}' + r'\,'.join(entries) + r'\right)'
        return text if exp is None else text + '^{' + exp + '}'


class Commutator(sp.Function):
    """[X, Y] = XY - YX, with its arguments put in the deterministic order."""

    @classmethod
    def eval(cls, first, second):
        first_coefficient, first_factors = word_factors(first)
        second_coefficient, second_factors = word_factors(second)
        if not first_factors or not second_factors:
            return sp.S.Zero
        first_matrix, second_matrix = sp.Mul(*first_factors), sp.Mul(*second_factors)
        if first_matrix == second_matrix:
            return sp.S.Zero
        coefficient = first_coefficient * second_coefficient
        if sort_key(second_matrix) < sort_key(first_matrix):
            return -coefficient * cls(second_matrix, first_matrix)
        if coefficient != 1:
            return coefficient * cls(first_matrix, second_matrix)
        return None

    def _latex(self, printer, exp=None):
        text = r'\left[' + printer._print(self.args[0]) + ',' + printer._print(self.args[1]) + r'\right]'
        return text if exp is None else text + '^{' + exp + '}'


class SlotTrace(sp.Function):
    """Sum over the orderings of the movable matrices placed in the slots of a trace word."""
    is_commutative = True

    @classmethod
    def eval(cls, word, *movable):
        coefficient, factors = word_factors(word)
        ordered = sorted(movable, key=sort_key)
        if coefficient == 0:
            return sp.S.Zero
        if coefficient != 1 or tuple(ordered) != tuple(movable):
            return coefficient * cls(sp.Mul(*factors), *ordered)
        return None

    def _latex(self, printer):
        movable = r',\,'.join(printer._print(matrix) for matrix in self.args[1:])
        return (r'\sum_{\rm perm}^{(' + movable + r')}\left\langle ' + printer._print(self.args[0])
                + r'\right\rangle')


# ---------------------------------------------------------------- expansion into plain words
def matrix_words(expression):
    """{word: coefficient} of a matrix expression; a word is a tuple of symbol names."""
    if expression.is_commutative:
        return {(): expression}
    if expression.is_Symbol:
        return {(str(expression),): sp.S.One}
    if expression.is_Add:
        total = defaultdict(lambda: sp.S.Zero)
        for term in expression.args:
            for word, coefficient in matrix_words(term).items():
                total[word] += coefficient
        return dict(total)
    if expression.is_Mul:
        product = {(): sp.S.One}
        for factor in expression.args:
            product = multiply_words(product, matrix_words(factor))
        return product
    if expression.is_Pow and expression.exp.is_Integer and expression.exp > 0:
        product = {(): sp.S.One}
        for _ in range(int(expression.exp)):
            product = multiply_words(product, matrix_words(expression.base))
        return product
    if isinstance(expression, PermutationSum):
        total = defaultdict(lambda: sp.S.Zero)
        for ordering in itertools.permutations(expression.args):
            for word, coefficient in matrix_words(sp.Mul(*ordering)).items():
                total[word] += coefficient
        return dict(total)
    if isinstance(expression, Commutator):
        first, second = expression.args
        return add_words(matrix_words(first * second), matrix_words(second * first), -1)
    raise TypeError(expression)


def multiply_words(first, second):
    """Product of two {word: coefficient} dictionaries (word concatenation)."""
    product = defaultdict(lambda: sp.S.Zero)
    for first_word, first_coefficient in first.items():
        for second_word, second_coefficient in second.items():
            product[first_word + second_word] += first_coefficient * second_coefficient
    return dict(product)


def add_words(first, second, weight=1):
    """first + weight * second for {word: coefficient} dictionaries."""
    total = defaultdict(lambda: sp.S.Zero, first)
    for word, coefficient in second.items():
        total[word] += weight * coefficient
    return dict(total)


def canonical_word(word):
    """Smallest cyclic rotation of a tuple of symbol names."""
    if not word:
        return word
    return min(word[start:] + word[:start] for start in range(len(word)))


def single_trace_polynomial(words):
    """{(canonical word,): coefficient} of the trace of a {word: coefficient} dictionary."""
    total = defaultdict(lambda: sp.S.Zero)
    for word, coefficient in words.items():
        if word:
            total[(canonical_word(word),)] += coefficient
        else:
            total[()] += FLAVOURS * coefficient
    return dict(total)


def multiply_polynomials(first, second):
    """Product of trace polynomials, keys being sorted tuples of canonical words."""
    product = defaultdict(lambda: sp.S.Zero)
    for first_key, first_coefficient in first.items():
        for second_key, second_coefficient in second.items():
            product[tuple(sorted(first_key + second_key))] += first_coefficient * second_coefficient
    return dict(product)


def trace_polynomial(expression):
    """Expand an expression in Trace and SlotTrace into {sorted tuple of words: coefficient}."""
    if isinstance(expression, Trace):
        return single_trace_polynomial(matrix_words(expression.args[0]))
    if isinstance(expression, SlotTrace):
        word, movable = expression.args[0], expression.args[1:]
        slots = [slot for slot in SLOTS if expression.args[0].has(slot)]
        total = {}
        for ordering in itertools.permutations(movable):
            filled = word.xreplace(dict(zip(slots, ordering)))
            total = add_words(total, single_trace_polynomial(matrix_words(filled)))
        return total
    if expression.is_Add:
        total = {}
        for term in expression.args:
            total = add_words(total, trace_polynomial(term))
        return total
    if expression.is_Mul:
        product = {(): sp.S.One}
        for factor in expression.args:
            product = multiply_polynomials(product, trace_polynomial(factor))
        return product
    if expression.is_Pow and expression.exp.is_Integer and expression.exp > 0 and expression.base.has(Trace, SlotTrace):
        product = {(): sp.S.One}
        for _ in range(int(expression.exp)):
            product = multiply_polynomials(product, trace_polynomial(expression.base))
        return product
    if expression.has(Trace, SlotTrace):
        raise TypeError(expression)
    return {(): expression}


def polynomial_difference(first, second):
    """Nonzero coefficients of first - second after exact cancellation."""
    difference = {}
    for key in set(first) | set(second):
        value = sp.cancel(sp.expand(first.get(key, 0) - second.get(key, 0)))
        if value != 0:
            difference[key] = value
    return difference


# ---------------------------------------------------------------- printing
def trace_monomials(expression):
    """{product of traces: coefficient} of an expression linear in products of traces."""
    groups = defaultdict(lambda: sp.S.Zero)
    for term in sp.Add.make_args(sp.expand(expression)):
        traces = sp.Mul(*(factor for factor in sp.Mul.make_args(term) if factor.has(Trace, SlotTrace)))
        groups[traces] += term / traces
    return {traces: coefficient for traces, coefficient in groups.items() if coefficient != 0}


def compact_latex(expression, symbol_names):
    """LaTeX of a sum of trace products, one line per product, coefficients factorized."""
    names = dict(symbol_names)
    names.update({slot: r'\bullet' for slot in SLOTS})
    lines = []
    for traces, coefficient in sorted(trace_monomials(expression).items(), key=lambda item: str(item[0])):
        coefficient_text = sp.latex(sp.factor(coefficient), symbol_names=names)
        if traces == 1:
            lines.append(coefficient_text if coefficient_text.startswith('-') else '+ ' + coefficient_text)
        else:
            lines.append(r'+ \left(' + coefficient_text + r'\right)' + sp.latex(traces, symbol_names=names))
    return '\n'.join(lines) if lines else '0'
