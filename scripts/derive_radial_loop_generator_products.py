#!/usr/bin/env python3
"""Loop diagram products of a -> pi+ pi- pi0 in compact flavour traces, completeness applied; GeV units.

For each retained diagram (list DIAGRAMS, classes 1-6 of the estimate record) the vertices of its
blocks are multiplied around the loop.  Internal lines: ALP line = ALP pole leg at both ends;
radial line = sigma at both ends; meson line = 1_10 - v_a v_a^T in the ten-field space, i.e. the
U(3) completeness term (applied before any expansion, so no generator remains), plus the bare
angular term, minus the ALP pole term.  beta and gamma are the radial tree exchanges into
(pi+, pi-) and (a, pi0), with factor V(sigma; x, y)/(m_sigma^2 - s_{+-}).  The loop scale is m_a
for diagrams with an ALP line and m_sigma otherwise.  Symmetry and loop factors are excluded.
No loop integral is evaluated.

Run:  PYTHONPATH=scripts .venv/bin/python scripts/derive_radial_loop_generator_products.py
"""
import itertools
from functools import lru_cache

import sympy as sp

import derive_radial_loop_trace_contractions as common

# (label, symmetry factor, blocks of external legs at each vertex, internal lines); line j joins
# block j and block j+1 cyclically.  Same list as the estimate record.
DIAGRAMS = (
    ('1.1', sp.Rational(1, 2), (('A', '+', '-', '0'),), ('A',)),
    ('1.2', sp.Rational(1, 2), (('beta', 'A', '0'),), ('A',)),
    ('1.3', sp.Rational(1, 2), (('gamma', '+', '-'),), ('A',)),
    ('2.1', 1, (('A', '+'), ('-', '0')), ('A', 'M')),
    ('2.2', 1, (('A', '-'), ('+', '0')), ('A', 'M')),
    ('2.3', 1, (('A', '0'), ('+', '-')), ('A', 'M')),
    ('2.4', 1, (('beta',), ('A', '0')), ('A', 'M')),
    ('2.5', 1, (('gamma',), ('+', '-')), ('A', 'M')),
    ('3.1', sp.Rational(1, 2), (('A', '0'), ('+', '-')), ('A', 'A')),
    ('3.2', sp.Rational(1, 2), (('beta',), ('A', '0')), ('A', 'A')),
    ('3.3', sp.Rational(1, 2), (('gamma',), ('+', '-')), ('A', 'A')),
    ('4.1', 1, (('A',), ('+', '-', '0')), ('S', 'M')),
    ('4.2', 1, (('+',), ('A', '-', '0')), ('S', 'M')),
    ('4.3', 1, (('-',), ('A', '+', '0')), ('S', 'M')),
    ('4.4', 1, (('0',), ('A', '+', '-')), ('S', 'M')),
    ('4.5', 1, (('A',), ('beta', '0')), ('S', 'M')),
    ('4.6', 1, (('0',), ('beta', 'A')), ('S', 'M')),
    ('4.7', 1, (('A', '+'), ('-',), ('0',)), ('M', 'S', 'M')),
    ('4.8', 1, (('A', '-'), ('+',), ('0',)), ('M', 'S', 'M')),
    ('4.9', 1, (('A', '0'), ('+',), ('-',)), ('M', 'S', 'M')),
    ('4.10', 1, (('+', '-'), ('A',), ('0',)), ('M', 'S', 'M')),
    ('4.11', 1, (('+', '0'), ('A',), ('-',)), ('M', 'S', 'M')),
    ('4.12', 1, (('-', '0'), ('A',), ('+',)), ('M', 'S', 'M')),
    ('4.13', 1, (('beta',), ('A',), ('0',)), ('M', 'S', 'M')),
    ('4.14', 1, (('gamma',), ('+',), ('-',)), ('M', 'S', 'M')),
    ('5.1', 1, (('A',), ('+', '-', '0')), ('S', 'A')),
    ('5.2', 1, (('0',), ('A', '+', '-')), ('S', 'A')),
    ('5.3', 1, (('A',), ('beta', '0')), ('S', 'A')),
    ('5.4', 1, (('0',), ('beta', 'A')), ('S', 'A')),
    ('5.5', 1, (('A', '+'), ('-',), ('0',)), ('M', 'S', 'A')),
    ('5.6', 1, (('A', '-'), ('+',), ('0',)), ('M', 'S', 'A')),
    ('5.7', 1, (('+', '-'), ('A',), ('0',)), ('A', 'S', 'M')),
    ('5.8', 1, (('+', '-'), ('A',), ('0',)), ('M', 'S', 'A')),
    ('5.9', 1, (('+', '0'), ('A',), ('-',)), ('A', 'S', 'M')),
    ('5.10', 1, (('-', '0'), ('A',), ('+',)), ('A', 'S', 'M')),
    ('5.11', 1, (('beta',), ('A',), ('0',)), ('A', 'S', 'M')),
    ('5.12', 1, (('beta',), ('A',), ('0',)), ('M', 'S', 'A')),
    ('6.1', sp.Rational(1, 2), (('beta',), ('A', '0')), ('S', 'S')),
    ('6.2', 1, (('A', '0'), ('+',), ('-',)), ('S', 'M', 'S')),
    ('6.3', 1, (('beta',), ('A',), ('0',)), ('S', 'M', 'S')),
    ('6.4', 1, (('gamma',), ('+',), ('-',)), ('S', 'M', 'S')),
    ('6.5', 1, (('A',), ('+',), ('-',), ('0',)), ('S', 'M', 'S', 'M')),
    ('6.6', 1, (('+',), ('-',), ('A',), ('0',)), ('M', 'S', 'M', 'S')),
    ('6.7', sp.Rational(1, 2), (('beta',), ('gamma',)), ('S', 'S')),
)
DIAGRAM_DATA = {label: (symmetry, blocks, lines) for label, symmetry, blocks, lines in DIAGRAMS}

EXTERNAL_LEGS = {'A': common.a, '0': common.zero, '+': common.plus, '-': common.minus}
BRANCH_PAIRS = {'beta': (common.plus, common.minus), 'gamma': (common.a, common.zero)}
MESON_LINE_TERMS = ('completeness', 'bare angular', 'ALP pole')
GENERATORS = sp.symbols('G_1:5', commutative=False)   # placeholders, removed by completeness


def radial_branch(name):
    """V(sigma; x, y)/(m_sigma^2 - s_{+-}) for the radial tree exchange `name`."""
    return common.vertex([None, *BRANCH_PAIRS[name]]) / (common.ms**2 - common.s)


def loop_scale(lines):
    """m_a^2 for diagrams with an ALP line, m_sigma^2 otherwise."""
    return common.ma**2 if 'A' in lines else common.ms**2


@lru_cache(maxsize=None)
def diagram_terms(label):
    """[(meson-line terms, sign, factors)] of a diagram; the product of the factors times the sign,
    summed over the list, is the diagram product without symmetry and loop factors."""
    _, blocks, lines = DIAGRAM_DATA[label]
    meson_lines = [index for index, line in enumerate(lines) if line == 'M']
    scale = loop_scale(lines)
    branches = [radial_branch(name) for block in blocks for name in block if name in BRANCH_PAIRS]
    result = []
    for choice in itertools.product(MESON_LINE_TERMS, repeat=len(meson_lines)):
        ends, sign, generators = [], 1, []
        for index, line in enumerate(lines):
            if line == 'A':
                ends.append(common.internal)
            elif line == 'S':
                ends.append(None)
            else:
                term = choice[meson_lines.index(index)]
                if term == 'completeness':
                    generators.append(GENERATORS[index])
                    ends.append((GENERATORS[index], sp.S.Zero, 'loop'))
                elif term == 'bare angular':
                    ends.append(common.loop_bare)
                else:
                    ends.append(common.internal)
                    sign = -sign
        factors = []
        for index, block in enumerate(blocks):
            fixed = [None if name in BRANCH_PAIRS else EXTERNAL_LEGS[name] for name in block]
            factors.append(common.vertex(fixed + [ends[index - 1], ends[index]], scale))
        if any(factor == 0 for factor in factors):
            continue
        for generator in generators:
            joined = [factor for factor in factors if factor.has(generator)]
            factors = [factor for factor in factors if not factor.has(generator)]
            factors.append(common.contract_generator(sp.Mul(*joined), generator))
        if any(factor == 0 for factor in factors):
            continue
        result.append((choice, sign, tuple(factors) + tuple(branches)))
    return tuple(result)


def diagram_total(label):
    """Unexpanded diagram product: sum over meson-line terms of sign times product of factors."""
    return sp.Add(*(sign * sp.Mul(*factors) for _, sign, factors in diagram_terms(label)))


def describe_choice(label, choice):
    """Text naming the meson-line term used on each meson line of the diagram."""
    _, _, lines = DIAGRAM_DATA[label]
    meson_lines = [index + 1 for index, line in enumerate(lines) if line == 'M']
    if not meson_lines:
        return 'no meson line'
    parts = []
    for number, term in zip(meson_lines, choice):
        text = {'completeness': 'U(3) completeness term', 'bare angular': 'bare angular term',
                'ALP pole': 'ALP pole term (subtracted)'}[term]
        parts.append(f'line {number}: {text}')
    return '; '.join(parts)


def print_diagrams(labels):
    """Print the compact products of the listed diagrams; factors shared by all meson-line terms
    are printed once."""
    print(common.NOTATION + '\n')
    for label in labels:
        symmetry, blocks, lines = DIAGRAM_DATA[label]
        terms = diagram_terms(label)
        shared = [factor for factor in terms[0][2] if all(factor in factors for _, _, factors in terms)]
        print(f'Diagram {label}: blocks {blocks}, lines {lines}, symmetry factor {symmetry} (not included)')
        print('  product of the shared factors and of the sum over meson-line terms below')
        for factor in shared:
            print('  shared factor:\n' + common.to_latex(factor))
        for choice, sign, factors in terms:
            print(f'  [{describe_choice(label, choice)}] overall sign {sign:+d}, product of the factors:')
            remaining = [factor for factor in factors if factor not in shared]
            for factor in remaining or [sp.S.One]:
                print('  factor:\n' + common.to_latex(factor))
        print()


if __name__ == '__main__':
    print_diagrams([label for label, *_ in DIAGRAMS])
