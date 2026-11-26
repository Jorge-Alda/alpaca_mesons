"""Exact flavour-trace algebra for the size estimate of the one-loop sigma -> pi+ pi- amplitude.

This module is a library; it prints nothing.  It is used by the scripts
derive_sigma_pipi_loop_*.py.

Conventions (Alda:2026nrd, DiLuzio:2026poh; memory/alp/shared/
radial-decay-polar-cartesian-general-sources.md and
murayama-radial-up-only-estimate-result.md):

* Natural units, GeV.  epsilon = f_pi / f_a.  N_f = 3, <lambda_i lambda_j> = 2 delta_ij.
* Polar Lagrangian in the connection-free frame
      L = f_pi^2/4 <u_mu u^mu> + B_0 f_pi^2 sum_{r,s} (sigma/f_a)^r (a/f_a)^s
            [<S_rs cos(Phi/f_pi)> + <P_rs sin(Phi/f_pi)>] - m_1^2/12 <Phi>^2
          + 1/2 (d sigma)^2 + 1/2 (1 + sigma/f_a)^2 (d a)^2 - V_phi(sigma),
  V_phi = m_phi^2/2 (f_a + sigma)^2 + lambda/4 (f_a + sigma)^4.
* A vertex V(x_1, ..., x_n) is the n-th derivative of L with respect to the
  fields of the legs; derivatives become d_mu -> -i k_mu with all momenta
  incoming.  A leg is a pair (lambda_x, xi_x) (meson matrix, ALP component)
  or a radial leg sigma.
* A flavour trace is stored as a word of matrix letters (strings), compared
  modulo cyclic rotation only.  The trace of the empty word is N_f = 3.

Letters: 'S{r}{s}', 'P{r}{s}' are the dressed sources S_hat_{r,s}, P_hat_{r,s};
'l+' and 'l-' the charged-pion matrices lambda_{pi+}, lambda_{pi-}; 'Minv' is
M^{-1}_{m_a^2}(C) with C = B_0 P_hat_{0,1}; letters starting with 'G' are
internal U(3) generators that are summed by completeness.

The module never evaluates a loop integral.  It handles vertices, products of
vertices, completeness sums and bookkeeping of powers of epsilon only.
"""
import itertools
import math
from dataclasses import dataclass, field

import sympy as sp

# ---------------------------------------------------------------------------
# Symbols
# ---------------------------------------------------------------------------
EPSILON = sp.Symbol('epsilon', positive=True)
B0 = sp.Symbol('B_0', positive=True)
F_PI = sp.Symbol('f_pi', positive=True)
QUARTIC = sp.Symbol('lambda', positive=True)          # quartic coupling of V_phi
NORMALIZATION = sp.Symbol('N', positive=True)         # N(m_a^2)
M_SIGMA = sp.Symbol('m_sigma', positive=True)
M_A = sp.Symbol('m_a', positive=True)
M_PI = sp.Symbol('m_pi', positive=True)
# m_P: largest pole mass of the nine meson-like states; loop scale of the bare
# U(3) generator sum of a meson line (decision D3).  Always symbolic.
M_MESON = sp.Symbol('m_P', positive=True)
MOMENTUM_SQUARED = sp.Symbol('p^2', positive=True)    # external p^2 of a two-point function
NUMBER_OF_FLAVOURS = 3

LATEX_NAMES = {EPSILON: r'\epsilon', F_PI: r'f_\pi', NORMALIZATION: r'N(m_a^2)',
               M_SIGMA: r'm_\sigma', M_A: r'm_a', M_PI: r'm_{\pi^\pm}', M_MESON: r'm_P', QUARTIC: r'\lambda',
               MOMENTUM_SQUARED: r'p^2'}

CHARGES = {'l+': 1, 'l-': -1, 'lY+': 1, 'lY-': -1}


def is_generator(letter):
    return letter.startswith('G')


def is_source(letter):
    return letter[0] in 'SP' and letter[1:].isdigit()


def letter_charge(letter):
    """Electric charge of a letter (ad_Q eigenvalue); None for generators."""
    if is_generator(letter):
        return None
    return CHARGES.get(letter, 0)


def source_letter(kind, radial, alp):
    return f'{kind}{radial}{alp}'


# ---------------------------------------------------------------------------
# Words and sums of products of traces
# ---------------------------------------------------------------------------

def canonical_word(word):
    """Smallest cyclic rotation of a word (cyclicity of the trace only)."""
    word = tuple(word)
    if not word:
        return word
    return min(word[k:] + word[:k] for k in range(len(word)))


class TraceSum:
    """Linear combination of products of flavour traces.

    terms: {key: coefficient}; key is a sorted tuple of canonical words; the
    empty key is the number 1.  Traces of the empty word are replaced by N_f."""

    def __init__(self, terms=None):
        self.terms = {}
        for key, value in (terms or {}).items():
            self.add_term(key, value)

    @staticmethod
    def scalar(value):
        return TraceSum({(): value})

    @staticmethod
    def single(word, value=1):
        return TraceSum({(tuple(word),): value})

    def add_term(self, words, value):
        factor = 1
        canonical = []
        for word in words:
            if len(word) == 0:
                factor *= NUMBER_OF_FLAVOURS
            else:
                canonical.append(canonical_word(word))
        key = tuple(sorted(canonical))
        self.terms[key] = self.terms.get(key, 0) + factor * value

    def __add__(self, other):
        result = TraceSum(dict(self.terms))
        for key, value in other.terms.items():
            result.terms[key] = result.terms.get(key, 0) + value
        return result

    def __sub__(self, other):
        return self + other.scale(-1)

    def scale(self, factor):
        return TraceSum({key: factor * value for key, value in self.terms.items()})

    def __mul__(self, other):
        result = TraceSum()
        for key_one, value_one in self.terms.items():
            for key_two, value_two in other.terms.items():
                result.add_term(key_one + key_two, value_one * value_two)
        return result

    def map_coefficients(self, function):
        return TraceSum({key: function(value) for key, value in self.terms.items()})

    def simplified(self):
        """Expanded coefficients, zero terms removed."""
        result = {}
        for key, value in self.terms.items():
            value = sp.expand(value)
            if value != 0:
                result[key] = value
        out = TraceSum()
        out.terms = result
        return out

    def is_zero(self):
        return not self.simplified().terms

    def filtered(self, predicate):
        """Keep the terms whose key satisfies predicate(key)."""
        out = TraceSum()
        out.terms = {key: value for key, value in self.terms.items() if predicate(key)}
        return out

    def letters(self):
        return {letter for key in self.terms for word in key for letter in word}

    def __eq__(self, other):
        return (self - other).is_zero()

    def __repr__(self):
        return f'TraceSum({self.simplified().terms})'


def word_charge(word):
    charges = [letter_charge(letter) for letter in word]
    if any(charge is None for charge in charges):
        return None
    return sum(charges)


VANISHING_LETTERS = ('P00', 'P10')          # assumptions 1 and 2 (matrix identities)
VANISHING_SINGLE_TRACES = (('S01',), ('S11',))  # <S01> = 0, <S11> = 0


def drop_vanishing_letters(trace_sum):
    """P_hat_00 = 0 (no tadpole) and P_hat_10 = 0 (no radial mixing)."""
    return trace_sum.filtered(lambda key: not any(letter in VANISHING_LETTERS
                                                  for word in key for letter in word))


def apply_assumptions(trace_sum):
    """All tadpole and radial-mixing assumptions, including the trace conditions."""
    trace_sum = drop_vanishing_letters(trace_sum)
    return trace_sum.filtered(lambda key: not any(word in VANISHING_SINGLE_TRACES for word in key))


def apply_charge_rule(trace_sum):
    """A trace of ad_Q eigen-matrices with nonzero total charge vanishes."""
    return trace_sum.filtered(lambda key: all(word_charge(word) in (None, 0) for word in key))


def completeness(trace_sum, generator):
    """Sum over the nine U(3) generators of a letter that appears exactly twice:
    <A G B G> -> 2 <A><B>,  <A G><G B> -> 2 <A B>."""
    result = TraceSum()
    for key, value in trace_sum.terms.items():
        positions = [(index, position) for index, word in enumerate(key)
                     for position, letter in enumerate(word) if letter == generator]
        if len(positions) != 2:
            raise ValueError(f'generator {generator} appears {len(positions)} times in {key}')
        (first_word, first_position), (second_word, second_position) = positions
        others = [word for index, word in enumerate(key) if index not in (first_word, second_word)]
        if first_word == second_word:
            word = key[first_word]
            rotated = word[first_position:] + word[:first_position]
            second = second_position - first_position
            inner, outer = rotated[1:second], rotated[second + 1:]
            result.add_term(tuple(others) + (inner, outer), 2 * value)
        else:
            one, two = key[first_word], key[second_word]
            rest_one = one[first_position + 1:] + one[:first_position]
            rest_two = two[second_position + 1:] + two[:second_position]
            result.add_term(tuple(others) + (rest_one + rest_two,), 2 * value)
    return result


# ---------------------------------------------------------------------------
# Compact matrix expressions (trees)
#   ('L', name)            a letter
#   ('U',)                 the unit matrix
#   ('P', c1, ..., cn)     the product c1 ... cn
#   ('S', c1, ..., cn)     sum over all n! orderings of the product
#   ('C', a, b)            commutator [a, b]
#   ('A', a, b)            anticommutator {a, b}
#   ('+', ((k1, t1), ...)) linear combination with rational coefficients
# ---------------------------------------------------------------------------
UNIT = ('U',)
ZERO_TREE = ('+', ())


def letter(name):
    return ('L', name)


def tree_letters(tree):
    """Letters of a tree with multiplicity (first term for linear combinations)."""
    kind = tree[0]
    if kind == 'L':
        return [tree[1]]
    if kind == 'U':
        return []
    if kind == '+':
        return tree_letters(tree[1][0][1]) if tree[1] else []
    return [name for child in tree[1:] for name in tree_letters(child)]


def tree_count(tree, name):
    return tree_letters(tree).count(name)


def expand_tree(tree):
    """{word: coefficient} of a tree."""
    kind = tree[0]
    if kind == 'L':
        return {(tree[1],): 1}
    if kind == 'U':
        return {(): 1}
    if kind == '+':
        result = {}
        for coefficient, child in tree[1]:
            for word, value in expand_tree(child).items():
                result[word] = result.get(word, 0) + coefficient * value
        return result
    if kind == 'C' or kind == 'A':
        sign = -1 if kind == 'C' else 1
        first, second = expand_tree(tree[1]), expand_tree(tree[2])
        result = {}
        for word_one, value_one in first.items():
            for word_two, value_two in second.items():
                result[word_one + word_two] = result.get(word_one + word_two, 0) + value_one * value_two
                result[word_two + word_one] = result.get(word_two + word_one, 0) + sign * value_one * value_two
        return result
    orders = itertools.permutations(tree[1:]) if kind == 'S' else [tree[1:]]
    result = {}
    for order in orders:
        for word, value in _product_expand(order).items():
            result[word] = result.get(word, 0) + value
    return result


def _product_expand(children):
    partial = {(): 1}
    for child in children:
        new = {}
        for word_one, value_one in partial.items():
            for word_two, value_two in expand_tree(child).items():
                word = word_one + word_two
                new[word] = new.get(word, 0) + value_one * value_two
        partial = new
    return partial


def trace_of_tree(tree):
    """TraceSum of the single trace <tree>."""
    result = TraceSum()
    for word, value in expand_tree(tree).items():
        result.add_term((word,), value)
    return result


def cut(tree, rest, generator):
    """Matrix D with <tree rest> = <generator D>; tree is linear in generator."""
    kind = tree[0]
    if kind == 'L':
        if tree[1] != generator:
            raise ValueError('letter is not the generator')
        return rest
    if kind == '+':
        return ('+', tuple((coefficient, cut(child, rest, generator)) for coefficient, child in tree[1]))
    children = list(tree[1:])
    index = next(i for i, child in enumerate(children) if tree_count(child, generator))
    target = children[index]
    if kind == 'P':
        return cut(target, ('P', *children[index + 1:], rest, *children[:index]), generator)
    if kind == 'S':
        return cut(target, ('S', *(children[:index] + children[index + 1:]), rest), generator)
    if kind == 'C':
        if index == 0:
            return cut(target, ('C', children[1], rest), generator)
        return cut(target, ('C', rest, children[0]), generator)
    if kind == 'A':
        other = children[1 - index]
        return cut(target, ('A', other, rest), generator)
    raise ValueError(tree)


def split(tree, generator):
    """[(k, alpha, beta)] with tree = sum k alpha generator beta."""
    kind = tree[0]
    if kind == 'L':
        return [(1, UNIT, UNIT)]
    if kind == '+':
        return [(coefficient * k, alpha, beta) for coefficient, child in tree[1]
                for k, alpha, beta in split(child, generator)]
    children = list(tree[1:])
    index = next(i for i, child in enumerate(children) if tree_count(child, generator))
    pieces = split(children[index], generator)
    if kind == 'P':
        return [(k, ('P', *children[:index], alpha), ('P', beta, *children[index + 1:]))
                for k, alpha, beta in pieces]
    if kind == 'S':
        others = children[:index] + children[index + 1:]
        result = []
        for mask in range(2 ** len(others)):
            before = [child for bit, child in enumerate(others) if mask >> bit & 1]
            after = [child for bit, child in enumerate(others) if not mask >> bit & 1]
            left = ('S', *before) if before else UNIT
            right = ('S', *after) if after else UNIT
            result += [(k, ('P', left, alpha), ('P', beta, right)) for k, alpha, beta in pieces]
        return result
    sign = -1 if kind == 'C' else 1
    if index == 0:
        other = children[1]
        return ([(k, alpha, ('P', beta, other)) for k, alpha, beta in pieces]
                + [(sign * k, ('P', other, alpha), beta) for k, alpha, beta in pieces])
    other = children[0]
    return ([(k, ('P', other, alpha), beta) for k, alpha, beta in pieces]
            + [(sign * k, alpha, ('P', beta, other)) for k, alpha, beta in pieces])


def rename_first(tree, name, replacement):
    """Rename one occurrence of the letter `name` to `replacement` in every term of a tree."""
    kind = tree[0]
    if kind == 'L':
        return ('L', replacement) if tree[1] == name else tree
    if kind == 'U':
        return tree
    if kind == '+':
        return ('+', tuple((coefficient, rename_first(child, name, replacement)) for coefficient, child in tree[1]))
    children = list(tree[1:])
    index = next(i for i, child in enumerate(children) if tree_count(child, name))
    children[index] = rename_first(children[index], name, replacement)
    return (kind, *children)


def tree_key(tree):
    return repr(tree)


def simplify_tree(tree):
    """(coefficient, tree) with products flattened, units removed, short
    permutation sums written as anticommutators and children of symmetric
    operations sorted."""
    kind = tree[0]
    if kind in ('L', 'U'):
        return 1, tree
    if kind == '+':
        combined = {}
        for coefficient, child in tree[1]:
            factor, child = simplify_tree(child)
            if factor == 0:
                continue
            if child[0] == '+':
                for inner_coefficient, inner in child[1]:
                    combined[inner] = combined.get(inner, 0) + coefficient * factor * inner_coefficient
            else:
                combined[child] = combined.get(child, 0) + coefficient * factor
        items = sorted(((value, child) for child, value in combined.items() if value != 0),
                       key=lambda item: tree_key(item[1]))
        if not items:
            return 0, ZERO_TREE
        if len(items) == 1:
            return items[0]
        lead = items[0][0]
        return lead, ('+', tuple((sp.nsimplify(value / lead), child) for value, child in items))
    factor = 1
    children = []
    for child in tree[1:]:
        child_factor, child = simplify_tree(child)
        if child_factor == 0:
            return 0, ZERO_TREE
        factor *= child_factor
        children.append(child)
    if kind == 'P':
        flat = []
        for child in children:
            if child[0] == 'P':
                flat.extend(child[1:])
            elif child[0] != 'U':
                flat.append(child)
        if not flat:
            return factor, UNIT
        if len(flat) == 1:
            return factor, flat[0]
        return factor, ('P', *flat)
    if kind == 'S':
        units = sum(1 for child in children if child[0] == 'U')
        rest = sorted((child for child in children if child[0] != 'U'), key=tree_key)
        factor *= math.factorial(len(children)) // math.factorial(len(rest))
        if not rest:
            return factor, UNIT
        if len(rest) == 1:
            return factor, rest[0]
        if len(rest) == 2:
            return factor, ('A', rest[0], rest[1])
        return factor, ('S', *rest)
    first, second = children
    if kind == 'A':
        if first == second:
            return 2 * factor, ('P', first, first)
        if first[0] == 'U':
            return 2 * factor, second
        if second[0] == 'U':
            return 2 * factor, first
        if tree_key(second) < tree_key(first):
            first, second = second, first
        return factor, ('A', first, second)
    if kind == 'C':
        if first[0] == 'U' or second[0] == 'U' or first == second:
            return 0, ZERO_TREE
        if tree_key(second) < tree_key(first):
            return -factor, ('C', second, first)
        return factor, ('C', first, second)
    raise ValueError(tree)


def canonical_trace_tree(tree):
    """(coefficient, tree) for the trace <tree>: simplified, top product rotated
    to its smallest cyclic form, <[a,b]> = 0, <{a,b}> = 2<ab>, <1> = N_f."""
    factor, tree = simplify_tree(tree)
    if factor == 0:
        return 0, ZERO_TREE
    if tree[0] == 'U':
        return factor * NUMBER_OF_FLAVOURS, UNIT
    if tree[0] == 'C':
        return 0, ZERO_TREE
    if tree[0] == 'A':
        inner_factor, inner = canonical_trace_tree(('P', tree[1], tree[2]))
        return 2 * factor * inner_factor, inner
    if tree[0] == 'P':
        children = list(tree[1:])
        rotations = [tuple(children[k:] + children[:k]) for k in range(len(children))]
        best = min(rotations, key=lambda rotation: tuple(tree_key(child) for child in rotation))
        return factor, ('P', *best)
    return factor, tree


class CompactSum:
    """Linear combination of products of traces of trees.

    terms: {key: coefficient}; key is a sorted tuple of canonical trace trees."""

    def __init__(self, terms=None):
        self.terms = {}
        for key, value in (terms or {}).items():
            self.add_term(key, value)

    @staticmethod
    def scalar(value):
        return CompactSum({(): value})

    def add_term(self, traces, value):
        factor = 1
        canonical = []
        for tree in traces:
            trace_factor, tree = canonical_trace_tree(tree)
            if trace_factor == 0:
                return
            if tree[0] == 'U':
                factor *= trace_factor
                continue
            factor *= trace_factor
            canonical.append(tree)
        key = tuple(sorted(canonical, key=tree_key))
        self.terms[key] = self.terms.get(key, 0) + factor * value

    def __add__(self, other):
        result = CompactSum()
        result.terms = dict(self.terms)
        for key, value in other.terms.items():
            result.terms[key] = result.terms.get(key, 0) + value
        return result

    def scale(self, factor):
        result = CompactSum()
        result.terms = {key: factor * value for key, value in self.terms.items()}
        return result

    def __mul__(self, other):
        result = CompactSum()
        for key_one, value_one in self.terms.items():
            for key_two, value_two in other.terms.items():
                result.add_term(key_one + key_two, value_one * value_two)
        return result

    def map_coefficients(self, function):
        result = CompactSum()
        result.terms = {key: function(value) for key, value in self.terms.items()}
        return result

    def simplified(self):
        result = CompactSum()
        for key, value in self.terms.items():
            value = sp.factor(sp.expand(value))
            if value != 0:
                result.terms[key] = value
        return result

    def filtered(self, predicate):
        result = CompactSum()
        result.terms = {key: value for key, value in self.terms.items() if predicate(key)}
        return result

    def to_trace_sum(self):
        """Expansion into products of traces of words (independent of the trees)."""
        result = TraceSum()
        for key, value in self.terms.items():
            product = TraceSum.scalar(value)
            for tree in key:
                product = product * trace_of_tree(tree)
            result = result + product
        return result


def compact_completeness(compact_sum, generator):
    """Completeness on trees, applied before any expansion:
    <A G><G B> -> 2 <A B> and <A G B G> -> 2 <A><B>."""
    result = CompactSum()
    for key, value in compact_sum.terms.items():
        counts = [tree_count(tree, generator) for tree in key]
        holders = [index for index, count in enumerate(counts) if count]
        others = [tree for index, tree in enumerate(key) if index not in holders]
        if len(holders) == 2:
            first = cut(key[holders[0]], UNIT, generator)
            second = cut(key[holders[1]], UNIT, generator)
            result.add_term(tuple(others) + (('P', first, second),), 2 * value)
        elif len(holders) == 1 and counts[holders[0]] == 2:
            partner = generator + "'"
            tree = rename_first(key[holders[0]], generator, partner)
            remainder = cut(tree, UNIT, generator)
            for k, alpha, beta in split(remainder, partner):
                result.add_term(tuple(others) + (alpha, beta), 2 * k * value)
        else:
            raise ValueError(f'generator {generator} appears {counts} times')
    return result


def compact_assumptions(compact_sum):
    """Assumptions 1 and 2 and the charge rule on compact traces."""
    def keep(key):
        for tree in key:
            names = tree_letters(tree)
            if any(name in VANISHING_LETTERS for name in names):
                return False
            if tree == ('L', 'S01') or tree == ('L', 'S11'):
                return False
            charges = [letter_charge(name) for name in names]
            if None not in charges and sum(charges) != 0:
                return False
        return True
    return compact_sum.filtered(keep)


# ---------------------------------------------------------------------------
# LaTeX output (compact matrix notation only)
# ---------------------------------------------------------------------------

def letter_latex(name):
    if name in ('l+', 'l-'):
        return r'\lambda_{\pi^' + name[1] + '}'
    if name == 'Minv':
        return r'\mathbb M^{-1}_{m_a^2}(\mathcal C)'
    if is_source(name):
        return r'\widehat ' + name[0] + '_{' + name[1] + ',' + name[2:] + '}'
    if name == 'Pu':
        return r'P_u'
    if name == 'Pd':
        return r'P_d'
    return name


def tree_latex(tree):
    kind = tree[0]
    if kind == 'L':
        return letter_latex(tree[1])
    if kind == 'U':
        return r'\mathbf 1'
    if kind == 'C':
        return '[' + tree_latex(tree[1]) + ',' + tree_latex(tree[2]) + ']'
    if kind == 'A':
        return r'\{' + tree_latex(tree[1]) + ',' + tree_latex(tree[2]) + r'\}'
    if kind == 'S':
        return (r'\sum_{\rm perm}\big(' + r'\,'.join(tree_latex(child) for child in tree[1:]) + r'\big)')
    if kind == 'P':
        parts = []
        for child in tree[1:]:
            text = tree_latex(child)
            parts.append('(' + text + ')' if child[0] == '+' else text)
        return r'\,'.join(parts)
    if kind == '+':
        text = ''
        for coefficient, child in tree[1]:
            coefficient_text = sp.latex(coefficient)
            if coefficient == 1:
                term = tree_latex(child)
                text += ('+' if text else '') + term
            elif coefficient == -1:
                text += '-' + tree_latex(child)
            else:
                sign = '' if (not text or coefficient_text.startswith('-')) else '+'
                text += sign + coefficient_text + r'\,' + tree_latex(child)
        return text
    raise ValueError(tree)


def coefficient_latex(value):
    return sp.latex(value, symbol_names=LATEX_NAMES)


def compact_latex(compact_sum, separator='\n'):
    """One line per term: coefficient times the product of compact traces."""
    lines = []
    for key, value in sorted(compact_sum.terms.items(), key=lambda item: tuple(map(tree_key, item[0]))):
        traces = ''.join(r'\langle ' + tree_latex(tree) + r'\rangle' for tree in key)
        lines.append(r'\left(' + coefficient_latex(value) + r'\right)' + traces)
    return separator.join(lines) if lines else '0'


def assert_no_float(expression):
    """Exactness guard: no floating-point number in an expression."""
    if isinstance(expression, (TraceSum, CompactSum)):
        for value in expression.terms.values():
            assert_no_float(value)
        return
    assert not sp.sympify(expression).atoms(sp.Float), expression


# ---------------------------------------------------------------------------
# Legs, momenta and vertex rules
# ---------------------------------------------------------------------------

@dataclass
class Leg:
    """A vertex leg.  radial legs carry no meson matrix and no ALP component.
    meson: tuple of (coefficient, letter); alp: coefficient of the ALP field;
    momentum: {basis label: coefficient}, incoming."""
    radial: bool = False
    meson: tuple = ()
    alp: object = 0
    momentum: dict = field(default_factory=dict)


def radial_leg(momentum):
    return Leg(radial=True, momentum=dict(momentum))


def meson_leg(name, momentum, coefficient=1, alp=0):
    return Leg(meson=((coefficient, name),), alp=alp, momentum=dict(momentum))


def alp_like_leg(momentum):
    """The ALP-like pole vector v_a = (2 eps M^{-1}(C)/sqrt(N), 1/sqrt(N))."""
    return Leg(meson=((2 * EPSILON / sp.sqrt(NORMALIZATION), 'Minv'),),
               alp=1 / sp.sqrt(NORMALIZATION), momentum=dict(momentum))


def bare_angular_leg(momentum):
    return Leg(alp=sp.Integer(1), momentum=dict(momentum))


def negate(momentum):
    return {key: -value for key, value in momentum.items()}


def momentum_sum(*momenta):
    total = {}
    for momentum in momenta:
        for key, value in momentum.items():
            total[key] = total.get(key, 0) + value
    return {key: value for key, value in total.items() if value != 0}


def make_dot(basis_dot):
    """Bilinear dot product of momenta from the products of basis labels."""
    def dot(first, second):
        return sp.expand(sum(value_one * value_two * basis_dot(key_one, key_two)
                             for key_one, value_one in first.items()
                             for key_two, value_two in second.items()))
    return dot


def symbolic_dot(first, second):
    """Generic dot products d(i,j) of basis labels, used for the vertex checks."""
    def basis(one, two):
        one, two = sorted((one, two))
        return sp.Symbol(f'd_{one}_{two}', real=True)
    return make_dot(basis)(first, second)


def potential_vertex(legs):
    """Source part: B_0 f_pi^{2-n} r! s! eps^{r+s} (-1)^{m(m-1)/2}/m!
    prod xi  sum_{perms} <X_{rs} lambda...>, X = S (m even), P (m odd)."""
    n = len(legs)
    radial = sum(1 for leg in legs if leg.radial)
    pseudoscalar = [leg for leg in legs if not leg.radial]
    result = TraceSum()
    options = []
    for leg in pseudoscalar:
        roles = []
        if leg.alp != 0:
            roles.append(('alp', leg.alp, None))
        for coefficient, name in leg.meson:
            roles.append(('meson', coefficient, name))
        options.append(roles)
    for choice in itertools.product(*options):
        alp_count = sum(1 for role in choice if role[0] == 'alp')
        names = [role[2] for role in choice if role[0] == 'meson']
        m = len(names)
        kind = 'S' if m % 2 == 0 else 'P'
        coefficient = (B0 * F_PI ** (2 - n) * math.factorial(radial) * math.factorial(alp_count)
                       * EPSILON ** (radial + alp_count) * sp.Integer(-1) ** (m // 2) / math.factorial(m))
        for role in choice:
            coefficient *= role[1]
        source = source_letter(kind, radial, alp_count)
        for order in itertools.permutations(names):
            result.add_term(((source,) + order,), coefficient)
    return result


def meson_kinetic_vertex(legs, dot):
    """Four-meson part of f_pi^2/4 <u_mu u^mu>:
    1/(48 f_pi^2) sum_{perms} (-k_2.k_4) <[l_1,l_2][l_3,l_4]>."""
    if len(legs) != 4 or any(leg.radial for leg in legs):
        return TraceSum()
    result = TraceSum()
    for order in itertools.permutations(range(4)):
        momentum_factor = -dot(legs[order[1]].momentum, legs[order[3]].momentum) / (48 * F_PI ** 2)
        for choice in itertools.product(*(legs[i].meson for i in order)):
            coefficient = momentum_factor
            for value, _ in choice:
                coefficient *= value
            a, b, c, d = (name for _, name in choice)
            for word, sign in (((a, b, c, d), 1), ((a, b, d, c), -1), ((b, a, c, d), -1), ((b, a, d, c), 1)):
                result.add_term((word,), sign * coefficient)
    return result


def alp_kinetic_vertex(legs, dot):
    """1/2 (1 + sigma/f_a)^2 (d a)^2: -xi_x xi_y k_x.k_y C(sigma legs),
    C(sigma) = 2 eps/f_pi, C(sigma, sigma) = 2 eps^2/f_pi^2."""
    radial = [leg for leg in legs if leg.radial]
    pseudoscalar = [leg for leg in legs if not leg.radial]
    if len(pseudoscalar) != 2 or len(radial) not in (1, 2):
        return TraceSum()
    first, second = pseudoscalar
    if first.alp == 0 or second.alp == 0:
        return TraceSum()
    factor = {1: 2 * EPSILON / F_PI, 2: 2 * EPSILON ** 2 / F_PI ** 2}[len(radial)]
    return TraceSum.scalar(-first.alp * second.alp * dot(first.momentum, second.momentum) * factor)


def radial_lambda_times_fa():
    """lambda f_a = eps [m_sigma^2 - eps^2 B_0(<S10> - 2<S20>)]/(2 f_pi), from
    m_sigma^2 = 2 lambda f_a^2 + B_0 eps^2 (<S10> - 2<S20>)."""
    return (TraceSum.scalar(EPSILON * M_SIGMA ** 2 / (2 * F_PI))
            + TraceSum.single(('S10',), -EPSILON ** 3 * B0 / (2 * F_PI))
            + TraceSum.single(('S20',), EPSILON ** 3 * B0 / F_PI))


def radial_self_vertex(legs, eliminate_lambda=True):
    """V_phi part: sigma^3 -> -6 lambda f_a, sigma^4 -> -6 lambda."""
    if not legs or not all(leg.radial for leg in legs) or len(legs) not in (3, 4):
        return TraceSum()
    if eliminate_lambda:
        lambda_fa = radial_lambda_times_fa()
    else:
        lambda_fa = TraceSum.scalar(QUARTIC * F_PI / EPSILON)
    if len(legs) == 3:
        return lambda_fa.scale(-6)
    return lambda_fa.scale(-6 * EPSILON / F_PI)


def vertex(legs, dot, eliminate_lambda=True):
    """Full vertex rule V(legs) (word level)."""
    return (potential_vertex(legs) + meson_kinetic_vertex(legs, dot)
            + alp_kinetic_vertex(legs, dot) + radial_self_vertex(legs, eliminate_lambda))


# Compact versions of the same rules, built from trees.

def compact_vertex(legs, dot, eliminate_lambda=True):
    n = len(legs)
    radial = sum(1 for leg in legs if leg.radial)
    pseudoscalar = [leg for leg in legs if not leg.radial]
    result = CompactSum()
    options = []
    for leg in pseudoscalar:
        roles = []
        if leg.alp != 0:
            roles.append(('alp', leg.alp, None))
        for coefficient, name in leg.meson:
            roles.append(('meson', coefficient, name))
        options.append(roles)
    for choice in itertools.product(*options):
        alp_count = sum(1 for role in choice if role[0] == 'alp')
        names = [role[2] for role in choice if role[0] == 'meson']
        m = len(names)
        kind = 'S' if m % 2 == 0 else 'P'
        coefficient = (B0 * F_PI ** (2 - n) * math.factorial(radial) * math.factorial(alp_count)
                       * EPSILON ** (radial + alp_count) * sp.Integer(-1) ** (m // 2) / math.factorial(m))
        for role in choice:
            coefficient *= role[1]
        source = letter(source_letter(kind, radial, alp_count))
        tree = ('P', source, ('S', *(letter(name) for name in names))) if names else source
        result.add_term((tree,), coefficient)
    if n == 4 and radial == 0:
        for order in itertools.permutations(range(4)):
            momentum_factor = -dot(legs[order[1]].momentum, legs[order[3]].momentum) / (48 * F_PI ** 2)
            for choice in itertools.product(*(legs[i].meson for i in order)):
                coefficient = momentum_factor
                for value, _ in choice:
                    coefficient *= value
                a, b, c, d = (letter(name) for _, name in choice)
                result.add_term((('P', ('C', a, b), ('C', c, d)),), coefficient)
    kinetic = alp_kinetic_vertex(legs, dot)
    for value in kinetic.terms.values():
        result.add_term((), value)
    for key, value in radial_self_vertex(legs, eliminate_lambda).terms.items():
        result.add_term(tuple(('L', word[0]) for word in key), value)
    return result


# ---------------------------------------------------------------------------
# Exact zero test with the normalization N(m_a^2) written through q
# ---------------------------------------------------------------------------
# N(m_a^2) = 1 + 2 eps^2 q with q = <(M^{-1}_{m_a^2}(C))^2> (exact identity,
# murayama-sigma-pipi-loop-large-epsilon-scaling.md, step 2).
TRACE_MINV_SQUARED = sp.Symbol('q', positive=True)
MINV_SQUARED_WORD = ('Minv', 'Minv')


def merge_minv_squared(trace_sum):
    """Replace every trace <M^{-1} M^{-1}> in a TraceSum by the symbol q, so that
    terms that differ only by this trace are combined."""
    result = TraceSum()
    for key, value in trace_sum.terms.items():
        count = sum(1 for word in key if word == MINV_SQUARED_WORD)
        rest = tuple(word for word in key if word != MINV_SQUARED_WORD)
        result.add_term(rest, value * TRACE_MINV_SQUARED ** count)
    return result.simplified()


def coefficient_is_exactly_zero(value, normalization=NORMALIZATION):
    """True if a coefficient that depends on eps and on N^{-1/2} vanishes after
    the exact substitution N = 1 + 2 eps^2 q.  The coefficient is split as
    E(eps, 1/N) + N^{-1/2} O(eps, 1/N); E and O become rational functions of
    eps and q, and the two parts cannot cancel each other (one is rational in
    eps, the other is not), so both must vanish."""
    value = sp.expand(value)
    if value == 0:
        return True
    root = sp.Dummy('root_inverse_normalization', positive=True)   # N^{-1/2}
    value = sp.expand(value.subs(normalization, root ** -2))
    even, odd = 0, 0
    for term in sp.Add.make_args(value):
        exponent = term.as_powers_dict().get(root, 0)
        if exponent % 2 == 0:
            even += term
        else:
            odd += term / root
    exact = 1 + 2 * EPSILON ** 2 * TRACE_MINV_SQUARED
    return all(sp.cancel(sp.together(sp.sympify(part).subs(root, 1 / sp.sqrt(exact)))) == 0
               for part in (even, odd))


def is_exactly_zero(trace_sum):
    """True if a TraceSum vanishes after the exact substitution
    N = 1 + 2 eps^2 q, every other product of traces being independent."""
    return all(coefficient_is_exactly_zero(value) for value in merge_minv_squared(trace_sum).terms.values())


def source_count(key):
    return sum(1 for word in key for name in word if is_source(name))
