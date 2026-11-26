#!/usr/bin/env python3
"""Vertex rules in compact flavour traces and the U(3) completeness contraction; GeV units.

Conventions (Alda:2026nrd, model of DiLuzio:2026poh, polar fields, connection-free frame):
    L = (f_pi^2/4)<d U d U^dag> + (f_pi^2/4)<chi U^dag + U chi^dag> + (1/2)(1 + sigma/f_a)^2 (d a)^2
        - lambda f_a sigma^3,   U = exp(i Phi/f_pi),  <lambda_i lambda_j> = 2 delta_ij,
    chi = 2 B_0 sum_{r,s} (S_rs + i P_rs) (sigma/f_a)^r (a/f_a)^s  (S_rs, P_rs: hatted sources).
A leg is (meson matrix, ALP component, momentum label); a radial leg is None.  Vertices are
multilinear coefficients of L (no factor i), momenta incoming.  Every product of two momenta
that contains the loop momentum is set to the loop scale squared (sizing prescription).
Potential part: B_0 f_pi^(2-n) r! s! (f_pi/f_a)^(r+s) (-1)^floor(m/2)/m! prod xi <X_rs PermutationSum>.
Four-meson kinetic part: -(1/(24 f_pi^2)) sum over the three pairings (k_1-k_2).(k_3-k_4)
<[L_1,L_2][L_3,L_4]>, from L_4 = <[Phi,dPhi][Phi,dPhi]>/(48 f_pi^2).
Six-meson kinetic part: from L_6 = <[Phi,[Phi,dPhi]][Phi,[Phi,dPhi]]>/(1440 f_pi^4).
Completeness, applied before any expansion: sum_i <A lambda_i B lambda_i> = 2<A><B> and
sum_i <A lambda_i><lambda_i B> = 2<AB>.  The meson line is 1_10 - v_a v_a^T in the ten-field
space (nine meson fields and the bare ALP): completeness term + bare angular term - ALP pole.

Library module; no output.  Example:
    PYTHONPATH=scripts .venv/bin/python -c "import derive_radial_loop_trace_contractions as c; \
print(c.to_latex(c.vertex([None, c.plus, c.minus])))"
"""
import itertools

import sympy as sp

from group_radial_loop_trace_permutations import (SLOTS, Commutator, PermutationSum, SlotTrace, Trace,
                                                  compact_latex, sort_key, word_factors)

B0, fpi, fa, ma, ms, s = sp.symbols('B_0 f_pi f_a m_a m_sigma s_pm', positive=True)
lam = sp.Symbol('lambda', positive=True)
xi_a, xi_0 = sp.symbols('xi_a xi_pi0', real=True)
la, l0, lp, lm = sp.symbols('lambda_a lambda_pi0 lambda_pi+ lambda_pi-', commutative=False)
S = {(r, t): sp.Symbol(f'S_{r}{t}', commutative=False) for r in range(4) for t in range(7)}
P = {(r, t): sp.Symbol(f'P_{r}{t}', commutative=False) for r in range(4) for t in range(7)}
tr = Trace

MOMENTUM_NAMES = {'a': 'a', '0': r'\pi^0', '+': r'\pi^+', '-': r'\pi^-'}
external_dots = {frozenset(pair): sp.Symbol(f'k_{pair[0]}.k_{pair[1]}', real=True)
                 for pair in itertools.combinations(('a', '0', '+', '-'), 2)}
k_a0 = external_dots[frozenset(('a', '0'))]
k_pm = external_dots[frozenset(('+', '-'))]

# Legs: (meson matrix, ALP component, momentum label).
a = (la, xi_a, 'a')
zero = (l0, xi_0, '0')
plus = (lp, sp.S.Zero, '+')
minus = (lm, sp.S.Zero, '-')
internal = (la, xi_a, 'loop')           # ALP pole leg on an internal line
loop_bare = (sp.S.Zero, sp.S.One, 'loop')  # bare angular field on a meson line
loop_alp = internal

LATEX_NAMES = {B0: 'B_0', fpi: r'f_\pi', fa: 'f_a', ma: 'm_a', ms: r'm_\sigma', s: r's_{+-}', lam: r'\lambda',
               xi_a: r'\xi_a', xi_0: r'\xi_{\pi^0}', la: r'\lambda_a', l0: r'\lambda_{\pi^0}',
               lp: r'\lambda_{\pi^+}', lm: r'\lambda_{\pi^-}'}
LATEX_NAMES.update({symbol: r'\widehat S_{%d%d}' % key for key, symbol in S.items()})
LATEX_NAMES.update({symbol: r'\widehat P_{%d%d}' % key for key, symbol in P.items()})
LATEX_NAMES.update({symbol: r'k_{%s}\cdot k_{%s}' % tuple(MOMENTUM_NAMES[label] for label in sorted(pair, key='a0+-'.index))
                    for pair, symbol in external_dots.items()})

NOTATION = (r'Notation: \langle X\rangle is the flavour trace; (\sum_{\rm perm}X_1\cdots X_m) is the sum over the m! '
            r'orderings of the listed factors and \{X,Y\}=XY+YX is its m=2 case; [X,Y]=XY-YX; '
            r'\sum_{\rm perm}^{(M_1,\ldots,M_k)}\langle\cdots\rangle places M_1..M_k in all k! orders at the positions '
            r'marked \bullet; \widehat S_{rs}, \widehat P_{rs} are the dressed sources; momenta incoming.')


def to_latex(expression, names=None):
    """Compact LaTeX of a sum of trace products with the record notation."""
    return compact_latex(expression, {**LATEX_NAMES, **(names or {})})


def momentum_product(first, second, loop_scale=ma**2, extra=None):
    """k_first . k_second: the loop scale squared if one leg is internal, else the external symbol."""
    if 'loop' in (first[2], second[2]):
        return loop_scale
    pair = frozenset((first[2], second[2]))
    if extra and pair in extra:
        return extra[pair]
    return external_dots[pair]


def potential_vertex(legs):
    """Multilinear coefficient of (f_pi^2/4)<chi U^dag + U chi^dag> for the given legs."""
    radial = sum(leg is None for leg in legs)
    pseudoscalars = [leg for leg in legs if leg is not None]
    total = sp.S.Zero
    for roles in itertools.product((0, 1), repeat=len(pseudoscalars)):
        angular = [leg for leg, role in zip(pseudoscalars, roles) if role]
        mesons = [leg for leg, role in zip(pseudoscalars, roles) if not role]
        if any(leg[1] == 0 for leg in angular) or any(leg[0] == 0 for leg in mesons):
            continue
        m, t = len(mesons), len(angular)
        source = (S if m % 2 == 0 else P)[radial, t]
        coefficient = (B0 * fpi**(2 - len(legs)) * sp.factorial(radial) * sp.factorial(t)
                       * (fpi / fa)**(radial + t) * (-1)**(m // 2) / sp.factorial(m))
        coefficient *= sp.prod([leg[1] for leg in angular])
        total += coefficient * Trace(source * PermutationSum(*(leg[0] for leg in mesons)))
    return total


def four_meson_kinetic(legs, dot):
    """-(1/(24 f_pi^2)) sum over pairings (k_1-k_2).(k_3-k_4) <[L_1,L_2][L_3,L_4]>."""
    total = sp.S.Zero
    for first, second, third, fourth in ((0, 1, 2, 3), (0, 2, 1, 3), (0, 3, 1, 2)):
        momenta = (dot(legs[first], legs[third]) - dot(legs[first], legs[fourth])
                   - dot(legs[second], legs[third]) + dot(legs[second], legs[fourth]))
        total += -momenta / (24 * fpi**2) * Trace(Commutator(legs[first][0], legs[second][0])
                                                  * Commutator(legs[third][0], legs[fourth][0]))
    return total


def six_meson_kinetic(legs, dot):
    """Sum over the derivative pair {c, f} of -(k_c.k_f)/(720 f_pi^4) times the slot sum of
    <[.,[.,L_c]][.,[.,L_f]]> over the other four meson matrices."""
    total = sp.S.Zero
    for pair in itertools.combinations(range(6), 2):
        # The two commutator factors may be exchanged inside the trace; order them by their matrices.
        first, second = sorted(pair, key=lambda index: sort_key(legs[index][0]))
        others = [legs[index][0] for index in range(6) if index not in (first, second)]
        word = (Commutator(SLOTS[0], Commutator(SLOTS[1], legs[first][0]))
                * Commutator(SLOTS[2], Commutator(SLOTS[3], legs[second][0])))
        total += -dot(legs[first], legs[second]) / (720 * fpi**4) * SlotTrace(word, *others)
    return total


def vertex(legs, loop_scale=ma**2, extra_dots=None):
    """Vertex V(legs): potential, meson kinetic, ALP kinetic and sigma^3 parts."""
    def dot(first, second):
        return momentum_product(first, second, loop_scale, extra_dots)
    total = potential_vertex(legs)
    radial = sum(leg is None for leg in legs)
    pseudoscalars = [leg for leg in legs if leg is not None]
    if len(pseudoscalars) == 2 and radial in (1, 2) and all(leg[1] != 0 for leg in pseudoscalars):
        total -= 2 * dot(*pseudoscalars) * pseudoscalars[0][1] * pseudoscalars[1][1] / fa**radial
    mesonic = radial == 0 and all(leg[0] != 0 for leg in legs)
    if mesonic and len(legs) == 4:
        total += four_meson_kinetic(legs, dot)
    if mesonic and len(legs) == 6:
        total += six_meson_kinetic(legs, dot)
    if radial == 3 and len(legs) == 3:
        total -= 6 * lam * fa
    return total


# ---------------------------------------------------------------- completeness contraction
def line_end_cofactor(rest, matrix, generator):
    """Y with <rest matrix> = <Y generator>, for `matrix` linear in `generator`."""
    if matrix == generator:
        return rest
    coefficient, factors = word_factors(matrix)
    if coefficient != 1 or len(factors) > 1:
        position = next(index for index, factor in enumerate(factors) if factor.has(generator))
        rotated = sp.Mul(*factors[position + 1:]) * rest * sp.Mul(*factors[:position])
        return coefficient * line_end_cofactor(rotated, factors[position], generator)
    if matrix.is_Add:
        return sp.Add(*(line_end_cofactor(rest, term, generator) for term in matrix.args))
    if isinstance(matrix, PermutationSum):
        if generator in matrix.args:
            others = list(matrix.args)
            others.remove(generator)
            # <W sum_perm(M, G)> = <sum_perm(M, W) G>
            return PermutationSum(*others, rest)
        expanded = sp.Add(*(sp.Mul(*ordering) for ordering in itertools.permutations(matrix.args)))
        return line_end_cofactor(rest, expanded, generator)
    if isinstance(matrix, Commutator):
        first, second = matrix.args
        # <W [A,B]> = <[W,A] B> = -<[W,B] A>
        if second.has(generator):
            return line_end_cofactor(Commutator(rest, first), second, generator)
        return -line_end_cofactor(Commutator(rest, second), first, generator)
    raise ValueError((matrix, generator))


def trace_cofactor(trace, generator):
    """Y with trace = <Y generator>."""
    coefficient, factors = word_factors(trace.args[0])
    position = next(index for index, factor in enumerate(factors) if factor.has(generator))
    rest = sp.Mul(*factors[position + 1:]) * sp.Mul(*factors[:position])
    return coefficient * line_end_cofactor(rest, factors[position], generator)


def trace_factors(term):
    """Commutative coefficient and list of trace factors (repeated according to powers)."""
    coefficient, traces = sp.S.One, []
    for factor in sp.Mul.make_args(term):
        if isinstance(factor, Trace):
            traces.append(factor)
        elif factor.is_Pow and isinstance(factor.base, Trace):
            traces.extend([factor.base] * int(factor.exp))
        else:
            coefficient *= factor
    return coefficient, traces


def contract_generator(expression, generator):
    """Sum over the nine U(3) generators of the placeholder `generator`, which appears once in
    each of two traces of every term: <Y_1 G><Y_2 G> -> 2<Y_1 Y_2>."""
    result = sp.S.Zero
    for term in sp.Add.make_args(sp.expand(expression)):
        coefficient, traces = trace_factors(term)
        selected = [trace for trace in traces if trace.has(generator)]
        if len(selected) != 2 or coefficient.has(generator):
            raise ValueError((generator, term))
        for trace in selected:
            traces.remove(trace)
        merged = Trace(trace_cofactor(selected[0], generator) * trace_cofactor(selected[1], generator))
        result += 2 * coefficient * sp.Mul(*traces) * merged
    return result
