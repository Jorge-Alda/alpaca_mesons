#!/usr/bin/env python3
"""Independent expansion of the Lagrangian into flavour-trace words, for the exact checks; GeV units.

This module does not use the vertex rules of derive_radial_loop_trace_contractions.  It expands

    L = (f_pi^2/4)<chi_+> + (f_pi^2/4)<d U d U^dag> + (1/2)(1 + sigma/f_a)^2 (d a)^2 - lambda f_a sigma^3,
    chi_+ = u^dag chi u^dag + u chi^dag u,  U = u u,  u = exp(i Phi/(2 f_pi)),
    chi = 2 B_0 sum_{r,s} (S_rs + i P_rs) (sigma/f_a)^r (a/f_a)^s,

as truncated power series in noncommuting symbols, with one marker per leg; the vertex is the
coefficient multilinear in all markers (no factor i, incoming momenta, d -> -i k).  A leg is
(meson matrix name or None, ALP component, radial component, momentum label).  Every product of
two momenta containing a loop momentum is the loop scale squared.  Internal lines are pairs of
legs: ALP line = ALP pole leg at both ends; radial line = sigma at both ends; meson line =
1_10 - v_a v_a^T in the ten-field space: a placeholder pair L_k ... L_k (removed by
<A L B L> -> 2<A><B> and <A L><L B> -> 2<AB>), plus the bare ALP at both ends, minus the
ALP pole leg at both ends.  Results are trace polynomials {sorted tuple of words: coefficient},
words being cyclically rotated tuples of symbol names, <1> = 3.

Library module, used by the tests.  Example:
    PYTHONPATH=scripts .venv/bin/python -c "import expand_radial_loop_lagrangian as e; \
print(e.diagram_polynomial((('A',), ('+', '-', '0')), ('S', 'A')))"
"""
import itertools
import math
from collections import defaultdict
from functools import lru_cache

import sympy as sp

from derive_radial_loop_trace_contractions import (B0, P, S, external_dots, fa, fpi, la, l0, lam,
                                                   lm, lp, ma, ms, s, xi_0, xi_a)
from group_radial_loop_trace_permutations import (FLAVOURS, canonical_word, multiply_polynomials,
                                                  polynomial_difference, trace_polynomial)

ALP_POLE = (str(la), xi_a, 0, 'a')
EXTERNAL = {'A': ALP_POLE, '0': (str(l0), xi_0, 0, '0'), '+': (str(lp), 0, 0, '+'), '-': (str(lm), 0, 0, '-')}
RADIAL = (None, 0, 1, 'sigma')
INTERNAL_ALP = (str(la), xi_a, 0, 'loop')
BARE_ALP = (None, 1, 0, 'loop')
BRANCH_PAIRS = {'beta': ('+', '-'), 'gamma': ('A', '0')}


def multiply_series(first, second):
    """Product of {(mask, word): coefficient} series, keeping each leg marker at most once."""
    product = defaultdict(lambda: sp.S.Zero)
    for (first_mask, first_word), first_coefficient in first.items():
        for (second_mask, second_word), second_coefficient in second.items():
            if first_mask & second_mask:
                continue
            product[first_mask | second_mask, first_word + second_word] += first_coefficient * second_coefficient
    return dict(product)


def exponential(generator_series, order, factor):
    """sum_k (factor)^k X^k / k! up to `order`, X a series without constant term."""
    total = {(0, ()): sp.S.One}
    power = {(0, ()): sp.S.One}
    for k in range(1, order + 1):
        power = multiply_series(power, generator_series)
        for key, coefficient in power.items():
            total[key] = total.get(key, 0) + factor**k / math.factorial(k) * coefficient
    return total


def vertex_polynomial(legs, loop_scale, extra_dots=None):
    """Multilinear coefficient of L in the leg markers, as a trace polynomial."""
    count = len(legs)
    full = (1 << count) - 1

    def dot(first, second):
        labels = (legs[first][3], legs[second][3])
        if 'loop' in labels:
            return loop_scale
        pair = frozenset(labels)
        return (extra_dots or {}).get(pair, external_dots.get(pair))

    phi = {(1 << index, (leg[0],)): sp.S.One for index, leg in enumerate(legs) if leg[0] is not None}
    sigma = {(1 << index, ()): sp.sympify(leg[2]) for index, leg in enumerate(legs) if leg[2] != 0}
    alp = {(1 << index, ()): sp.sympify(leg[1]) for index, leg in enumerate(legs) if leg[1] != 0}
    u = exponential(phi, count, sp.I / (2 * fpi))
    u_dagger = exponential(phi, count, -sp.I / (2 * fpi))
    result = defaultdict(lambda: sp.S.Zero)

    def add_trace(series, weight):
        for (mask, word), coefficient in series.items():
            if mask == full:
                key = ((canonical_word(word),) if word else ())
                result[key] += weight * coefficient * (1 if word else FLAVOURS)

    # (f^2/4)<u^dag chi u^dag + u chi^dag u>
    sigma_powers = [{(0, ()): sp.S.One}]
    alp_powers = [{(0, ()): sp.S.One}]
    for _ in range(count):
        sigma_powers.append(multiply_series(sigma_powers[-1], sigma))
        alp_powers.append(multiply_series(alp_powers[-1], alp))
    for r in range(count + 1):
        for t in range(count + 1 - r):
            fields = multiply_series(sigma_powers[r], alp_powers[t])
            if not fields:
                continue
            scalar = {(0, (str(S[r, t]),)): sp.S.One}
            pseudoscalar = {(0, (str(P[r, t]),)): sp.S.One}
            weight = fpi**2 / 4 * 2 * B0 / fa**(r + t)
            for source, phase in ((scalar, 1), (pseudoscalar, sp.I)):
                with_fields = multiply_series(fields, source)
                add_trace(multiply_series(multiply_series(u_dagger, with_fields), u_dagger), weight * phase)
                add_trace(multiply_series(multiply_series(u, with_fields), u), weight * sp.conjugate(phase))
    # (f^2/4)<dU dU^dag>, d acting on the markers of a word gives -i times their momentum sum.
    big_u = multiply_series(u, u)
    big_u_dagger = multiply_series(u_dagger, u_dagger)
    for (first_mask, first_word), first_coefficient in big_u.items():
        for (second_mask, second_word), second_coefficient in big_u_dagger.items():
            if first_mask & second_mask or first_mask | second_mask != full or not first_mask or not second_mask:
                continue
            momenta = sum(dot(i, j) for i in range(count) if first_mask >> i & 1
                          for j in range(count) if second_mask >> j & 1)
            word = first_word + second_word
            key = ((canonical_word(word),) if word else ())
            result[key] += fpi**2 / 4 * (-momenta) * first_coefficient * second_coefficient * (1 if word else FLAVOURS)
    # (1/2)(1 + sigma/f_a)^2 (d a)^2 and -lambda f_a sigma^3
    alp_derivatives = defaultdict(lambda: sp.S.Zero)
    for (first_mask, _), first_coefficient in alp.items():
        for (second_mask, _), second_coefficient in alp.items():
            if first_mask != second_mask:
                first, second = first_mask.bit_length() - 1, second_mask.bit_length() - 1
                alp_derivatives[first_mask | second_mask, ()] += -dot(first, second) * first_coefficient * second_coefficient
    dressing = {(0, ()): sp.S.One}
    for key, coefficient in sigma.items():
        dressing[key] = 2 * coefficient / fa
    for key, coefficient in sigma_powers[2].items():
        dressing[key] = dressing.get(key, 0) + coefficient / fa**2
    for (mask, _), coefficient in multiply_series(dressing, dict(alp_derivatives)).items():
        if mask == full:
            result[()] += coefficient / 2
    for (mask, _), coefficient in sigma_powers[3].items():
        if mask == full:
            result[()] += -lam * fa * coefficient
    return {key: sp.expand(value) for key, value in result.items() if sp.expand(value) != 0}


@lru_cache(maxsize=None)
def cached_vertex(legs, loop_scale):
    return vertex_polynomial(legs, loop_scale)


def apply_completeness(polynomial, placeholder):
    """<A L B L> -> 2<A><B> and <A L><L B> -> 2<AB> for the placeholder L."""
    result = defaultdict(lambda: sp.S.Zero)
    for key, coefficient in polynomial.items():
        holding = [word for word in key if placeholder in word]
        others = list(key)
        for word in holding:
            others.remove(word)
        new_words, factor = [], 2 * coefficient
        if len(holding) == 2 and all(word.count(placeholder) == 1 for word in holding):
            joined = ()
            for word in holding:
                position = word.index(placeholder)
                joined += word[position + 1:] + word[:position]
            new_words.append(joined)
        elif len(holding) == 1 and holding[0].count(placeholder) == 2:
            word = holding[0]
            first = word.index(placeholder)
            second = word.index(placeholder, first + 1)
            new_words += [word[first + 1:second], word[second + 1:] + word[:first]]
        else:
            raise ValueError((placeholder, key))
        for word in new_words:
            if word:
                others.append(canonical_word(word))
            else:
                factor *= FLAVOURS
        result[tuple(sorted(others))] += factor
    return dict(result)


def scale_polynomial(polynomial, factor):
    return {key: factor * value for key, value in polynomial.items()}


def add_polynomials(first, second):
    total = defaultdict(lambda: sp.S.Zero, first)
    for key, value in second.items():
        total[key] += value
    return dict(total)


def branch_polynomial(name):
    """V(sigma; x, y)/(m_sigma^2 - s_{+-}) for the radial exchange `name`."""
    first, second = BRANCH_PAIRS[name]
    return scale_polynomial(cached_vertex((RADIAL, EXTERNAL[first], EXTERNAL[second]), 0), 1 / (ms**2 - s))


def diagram_polynomial(blocks, lines):
    """Diagram product (symmetry and loop factors excluded) as a trace polynomial."""
    scale = ma**2 if 'A' in lines else ms**2
    meson_lines = [index for index, line in enumerate(lines) if line == 'M']
    total = {}
    for choice in itertools.product(('placeholder', 'bare', 'pole'), repeat=len(meson_lines)):
        ends, sign, placeholders = [], 1, []
        for index, line in enumerate(lines):
            if line == 'A':
                ends.append(INTERNAL_ALP)
            elif line == 'S':
                ends.append((None, 0, 1, 'loop'))
            else:
                term = choice[meson_lines.index(index)]
                if term == 'placeholder':
                    placeholders.append(f'L_{index}')
                    ends.append((f'L_{index}', 0, 0, 'loop'))
                elif term == 'bare':
                    ends.append(BARE_ALP)
                else:
                    ends.append(INTERNAL_ALP)
                    sign = -sign
        product = {(): sp.S(sign)}
        for index, block in enumerate(blocks):
            fixed = tuple(RADIAL if name in BRANCH_PAIRS else EXTERNAL[name] for name in block)
            product = multiply_polynomials(product, cached_vertex(fixed + (ends[index - 1], ends[index]), scale))
            for name in block:
                if name in BRANCH_PAIRS:
                    product = multiply_polynomials(product, branch_polynomial(name))
        for placeholder in placeholders:
            product = apply_completeness(product, placeholder)
        total = add_polynomials(total, product)
    return total


def compact_difference(blocks, lines, compact_expression):
    """Nonzero coefficients of (Lagrangian expansion) - (compact expression expanded into words)."""
    return polynomial_difference(diagram_polynomial(blocks, lines), trace_polynomial(compact_expression))
