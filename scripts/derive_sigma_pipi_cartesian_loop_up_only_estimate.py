"""Cartesian frame: up-only sources - leading large-epsilon term of the size estimate M_1^est/M_0 of the one-loop sigma -> pi+ pi- amplitude (task 4.3, up-only part).

Run:
    PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_cartesian_loop_up_only_estimate.py
    PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_cartesian_loop_up_only_estimate.py --markdown
    PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_cartesian_loop_up_only_estimate

OpenSpec change estimate-sigma-pipi-one-loop-murayama-cartesian.  Units GeV for
B_0, m_u, f_pi, m_1, m_P, m_pi; eps = f_pi/f_a and u = 1/eps are dimensionless.
Conventions of sigma_pipi_cartesian_loop_algebra.py.

Up-only sources (decision D-A4; model assumption): M_hat = m_u P_u phi^*/f_a, i.e.
S^C_00 = S^C_10 = m_u P_u, P^C_01 = -m_u P_u, every other S^C, P^C zero.  Then
A = kappa_sigma = B_0 m_u (boundary r = kappa_sigma/A = 1 of the regions of the
large-eps routing choice) and C = -B_0 m_u P_u.

Relations of the matrices (no matrix components): P_u^2 = P_u, P_d^2 = P_d,
P_u P_d = P_d P_u = 0, <P_u> = <P_d> = 1, <1> = 3; pi+ = u dbar:
lambda_{pi+} lambda_{pi-} = 2 P_u, lambda_{pi-} lambda_{pi+} = 2 P_d,
P_u lambda_{pi+} = lambda_{pi+} P_d = lambda_{pi+}, lambda_{pi+} P_u = P_d lambda_{pi+} = 0,
the same with u <-> d for lambda_{pi-}, lambda_{pi+-}^2 = 0.  A trace that these
relations do not reduce raises an error.

Method (nothing is truncated in the estimate; series only read leading powers):
 1. M^{-1}_{m_a^2}(C) is solved exactly on the span of P_u and 1, which
    L(X) = B_0 {S^C_00, X} + (m_1^2/3)<X> 1 maps into itself:
    M^{-1}_s(C) = a_P(s) P_u + a_1(s) 1 with a_P (coefficient of P_u) and a_1
    (coefficient of the unit matrix) exact rational functions of s (code names only).
    N(s) = 1 + 2 eps^2 <(a_P P_u + a_1 1)^2> = 1 + 2 eps^2 (a_P^2 + 2 a_P a_1 + 3 a_1^2).
 2. m_a^2 is the root of the exact up-only D_a(s) = s - eps^2 A + 2 eps^2 <C M_s^{-1}(C)>,
    as a series in u^2 obtained by iteration through relative order u^(2K+2);
    m_sigma^2 = 2 lambda f_pi^2/eps^2 + eps^2 B_0 m_u (exact).
 3. Every kept graph is taken with the large-eps routing choice of
    derive_sigma_pipi_cartesian_loop_large_epsilon_graphs.py (D12 at large eps),
    in the generic-source basis, then the up-only substitution is made and the
    traces are reduced.  max(m_P, ...) -> the largest of m_a, m_sigma (m_P fixed);
    max(m_a, m_sigma) -> m_a, because m_a^2 - m_sigma^2 -> 2 B_0 m_u > 0.
 4. Contributions to M_1^est/M_0: 1PI graph X_Gamma/M_0, s self-energy dZ/2,
    pi self-energy dZ (dZ_pi+ = dZ_pi-, both legs), M_0 = -B_0 m_u eps/f_pi.
 5. The leading term is read with K and K + 2 and must agree (exactness check).
 6. The same is done with the routing choices of the two neighbouring regions
    r < 1 and r > 1 (in the latter also with max(m_a, m_sigma) -> m_sigma, the
    branch of that region), to see whether the limit agrees from both sides.

No loop integral, scalar loop function, loop coefficient or loop amplitude is
computed.  The result is an order of magnitude; graphs with the same epsilon
scaling can cancel in the exact loop.
"""
import functools
import sys

import sympy as sp

import sigma_pipi_cartesian_loop_algebra as cartesian
import derive_sigma_pipi_cartesian_loop_traces as traces
import derive_sigma_pipi_cartesian_loop_large_epsilon_graphs as large
import derive_sigma_pipi_cartesian_loop_large_epsilon_vertices as vertices_script
import derive_sigma_pipi_cartesian_loop_tree as cartesian_tree
import derive_sigma_pipi_loop_tree as generic_tree_tools        # generic word rewriting only
from sigma_pipi_loop_algebra import (B0, EPSILON, F_PI, M_A, M_MESON, M_SIGMA, NORMALIZATION, QUARTIC,
                                     TraceSum, canonical_word)
from sigma_pipi_cartesian_loop_algebra import MINV_LETTER, SINGLET_MASS_SQUARED, UP_QUARK_MASS

U = vertices_script.INVERSE_EPSILON
SQUARED = vertices_script.SQUARED_MOMENTUM
# coefficients of P_u and of the unit matrix in M_s^{-1}(C) = (coefficient of P_u) P_u + (coefficient of 1) 1
PROJECTOR_COEFFICIENT, IDENTITY_COEFFICIENT = sp.symbols('projector_coefficient identity_coefficient')
ROOT_ORDER = 3                       # K: m_a^2 exact through u^(2K) (relative order u^(2K+2))

ALP_COEFFICIENT_UP_ONLY = B0 * UP_QUARK_MASS     # A = kappa_sigma = B_0 m_u
TREE_UP_ONLY = -B0 * UP_QUARK_MASS * EPSILON / F_PI

UP_ONLY_RULES = {**generic_tree_tools.UP_ONLY_RULES,
                 ('Pu', 'l+'): [(1, ('l+',))], ('l+', 'Pd'): [(1, ('l+',))], ('l+', 'Pu'): [], ('Pd', 'l+'): [],
                 ('Pd', 'l-'): [(1, ('l-',))], ('l-', 'Pu'): [(1, ('l-',))], ('l-', 'Pd'): [], ('Pu', 'l-'): [],
                 ('l+', 'l+'): [], ('l-', 'l-'): []}


class UnreducedTrace(Exception):
    """A trace that the up-only relations do not reduce to a number."""


# ---------------------------------------------------------------------------
# 1. Pole data, exact
# ---------------------------------------------------------------------------

@functools.lru_cache(maxsize=None)
def inverse_components():
    """(coefficient of P_u, coefficient of 1) of M_s^{-1}(C), exact functions of s, from M_s(X) = L(X) - s X = C = -B_0 m_u P_u,
    L(P_u) = 2 B_0 m_u P_u + (m_1^2/3) 1, L(1) = 2 B_0 m_u P_u + m_1^2 1."""
    m_u, m1 = UP_QUARK_MASS, SINGLET_MASS_SQUARED
    s = SQUARED
    equations = [2 * B0 * m_u * (PROJECTOR_COEFFICIENT + IDENTITY_COEFFICIENT) - s * PROJECTOR_COEFFICIENT + B0 * m_u,     # P_u component
                 m1 / 3 * PROJECTOR_COEFFICIENT + m1 * IDENTITY_COEFFICIENT - s * IDENTITY_COEFFICIENT]                    # 1 component
    solution = sp.solve(equations, [PROJECTOR_COEFFICIENT, IDENTITY_COEFFICIENT], dict=True)[0]
    return sp.factor(solution[PROJECTOR_COEFFICIENT]), sp.factor(solution[IDENTITY_COEFFICIENT])


def trace_c_inverse():
    """<C M_s^{-1}(C)> = -B_0 m_u (coefficient of P_u + coefficient of 1) (exact)."""
    projector_part, identity_part = inverse_components()
    return sp.factor(-B0 * UP_QUARK_MASS * (projector_part + identity_part))


def trace_inverse_squared():
    """<(M_s^{-1}(C))^2> = a_P^2 + 2 a_P a_1 + 3 a_1^2, a_P, a_1 the coefficients of P_u and 1 (exact)."""
    projector_part, identity_part = inverse_components()
    return sp.factor(projector_part ** 2 + 2 * projector_part * identity_part + 3 * identity_part ** 2)


def alp_pole_function():
    """u^2 D_a(s) = u^2 s - A + 2 <C M_s^{-1}(C)> (exact up-only)."""
    return U ** 2 * SQUARED - ALP_COEFFICIENT_UP_ONLY + 2 * trace_c_inverse()


@functools.lru_cache(maxsize=None)
def alp_mass_squared(order=ROOT_ORDER):
    """m_a^2 = A/u^2 + sum_{k<=order} b_k u^(2k): fixed-point iteration of s = (A - 2 <C M_s^{-1}(C)>)/u^2,
    each step exact through one more power of u^2."""
    root = ALP_COEFFICIENT_UP_ONLY / U ** 2
    for _ in range(order + 2):
        update = (ALP_COEFFICIENT_UP_ONLY - 2 * trace_c_inverse().subs(SQUARED, root)) / U ** 2
        root = sp.expand(sp.series(update, U, 0, 2 * order + 1).removeO())
    return root


def root_residual(order=ROOT_ORDER):
    """Leading power of u of u^2 D_a(m_a^2) at the truncated root (it must be >= 2 order + 4)."""
    value = alp_pole_function().subs(SQUARED, alp_mass_squared(order))
    series = sp.expand(sp.series(value, U, 0, 2 * order + 8).removeO())
    return min(sp.Poly(series, U).monoms())[0] if series != 0 else sp.oo


def radial_mass_squared():
    return 2 * QUARTIC * F_PI ** 2 * U ** 2 + ALP_COEFFICIENT_UP_ONLY / U ** 2


# ---------------------------------------------------------------------------
# 2. Up-only value of a TraceSum
# ---------------------------------------------------------------------------

def up_only_trace(word):
    """<word> for a word in P_u, lambda_pi+, lambda_pi- (the empty word is <1> = 3)."""
    if not word:
        return sp.Integer(3)
    total = sp.Integer(0)
    for reduced, factor in generic_tree_tools._rewrite({tuple(word): 1}, UP_ONLY_RULES).items():
        rotated = generic_tree_tools._rewrite({canonical_word(reduced): 1}, UP_ONLY_RULES) if reduced else {(): 1}
        for final, final_factor in rotated.items():
            if final == ():
                total += factor * final_factor * 3
            elif final in generic_tree_tools.UP_ONLY_SINGLE_TRACES:
                total += factor * final_factor * generic_tree_tools.UP_ONLY_SINGLE_TRACES[final]
            else:
                raise UnreducedTrace(final)
    return sp.expand(total)


def word_options(word):
    """[(letters, factor)]: source letters -> +-m_u P_u or zero, M^{-1}_{m_a^2}(C) -> (coefficient of P_u) P_u + (coefficient of 1) 1."""
    options = [((), sp.Integer(1))]
    for name in word:
        if name == MINV_LETTER:
            options = [(letters + ('Pu',), factor * PROJECTOR_COEFFICIENT) for letters, factor in options] + \
                      [(letters, factor * IDENTITY_COEFFICIENT) for letters, factor in options]
        elif cartesian.is_cartesian_source(name):
            replacement = cartesian.up_only_letter(name)
            if replacement is None:
                return []
            options = [(letters + (replacement[1],), factor * replacement[0]) for letters, factor in options]
        else:
            options = [(letters + (name,), factor) for letters, factor in options]
    return options


def up_only_scalar(trace_sum):
    """Up-only value of a TraceSum as a scalar in the coefficients of P_u and 1 of M^{-1}_{m_a^2}(C) and the other symbols."""
    total = sp.Integer(0)
    for key, coefficient in trace_sum.terms.items():
        product = sp.sympify(coefficient)
        for word in key:
            product *= sum((factor * up_only_trace(letters) for letters, factor in word_options(word)),
                           sp.Integer(0))
            if product == 0:
                break
        total += product
    return sp.expand(total)


# ---------------------------------------------------------------------------
# 3. Large-eps series
# ---------------------------------------------------------------------------

def heavier_mass(variant):
    return M_SIGMA if variant == 'r > 1' else M_A


def replace_maximum(expression, variant):
    """max(m_P, ...) -> the largest of m_a, m_sigma (m_P fixed); max(m_a, m_sigma) -> m_a (or m_sigma for the
    branch of the region r > 1)."""
    def choose(*arguments):
        rest = set(arguments) - {M_MESON}
        if rest == {M_A}:
            return M_A
        if rest == {M_SIGMA}:
            return M_SIGMA
        if rest == {M_A, M_SIGMA}:
            return heavier_mass(variant)
        raise ValueError(arguments)
    return sp.sympify(expression).replace(lambda node: isinstance(node, sp.Max), lambda node: choose(*node.args))


def laurent_of(expression, exact_relative_order):
    """Laurent series in u of an expression in u (sympy series), exact for the powers below
    lowest power + exact_relative_order."""
    series = sp.expand(sp.series(expression, U, 0, exact_relative_order + 8).removeO())
    terms = {}
    for term in sp.Add.make_args(series):
        coefficient, power = term.as_coeff_exponent(U)
        terms[int(power)] = terms.get(int(power), 0) + coefficient
    terms = {power: value for power, value in terms.items() if sp.expand(value) != 0}
    return large.Laurent(terms, min(terms) + exact_relative_order)


@functools.lru_cache(maxsize=None)
def pole_block(name, order):
    """Series of m_a^2, of the coefficients of P_u and 1 of M_s^{-1}(C), and of N at s = m_a^2.  m_a^2 is exact through u^(2 order), i.e. to relative
    order u^(2 order + 2); a power-law function f(s) of s then has the relative error u^(2 order + 4) at most,
    since f(s + ds) - f(s) ~ f ds/s and ds/s ~ u^(2 order + 4).  Each block is cut there."""
    root = alp_mass_squared(order)
    relative = 2 * order + 4
    if name == 'm_a^2':
        return laurent_of(root, relative)          # lowest power u^-2: exact below u^(2 order + 2)
    projector_part, identity_part = inverse_components()
    functions = {'projector coefficient': projector_part, 'identity coefficient': identity_part, 'N': 1 + 2 * trace_inverse_squared() / U ** 2}
    return laurent_of(functions[name].subs(SQUARED, root), relative)


def block_series(base, exponent, order):
    if base == EPSILON:
        return large.Laurent({-exponent: 1})
    if base == M_A:
        assert exponent % 2 == 0, exponent
        return pole_block('m_a^2', order).power(exponent // 2)
    if base == M_SIGMA:
        assert exponent % 2 == 0, exponent
        return large.Laurent({2: 2 * QUARTIC * F_PI ** 2, -2: ALP_COEFFICIENT_UP_ONLY}).power(exponent // 2)
    names = {PROJECTOR_COEFFICIENT: 'projector coefficient', IDENTITY_COEFFICIENT: 'identity coefficient', NORMALIZATION: 'N'}
    if base in names:
        return pole_block(names[base], order).power(exponent)
    return None


def series_in_u(expression, order, variant):
    """Laurent series in u of an exact up-only scalar (monomials in eps, m_a, m_sigma, the coefficients of P_u and 1, N times
    fixed symbols), with max(...) replaced as in replace_maximum.  Returns the Laurent series; its order is the
    first power of u that is not exact."""
    expression = sp.expand(replace_maximum(expression, variant))
    total = large.Laurent({}, large.EXACT)
    for term in sp.Add.make_args(expression):
        constant = sp.Integer(1)
        result = large.Laurent.constant(1)
        for factor in sp.Mul.make_args(term):
            base, exponent = factor.as_base_exp()
            series = block_series(base, int(exponent), order) if exponent.is_integer else None
            if series is None:
                constant *= factor
            else:
                result = result * series
        total = total + result.scale(constant)
    return total


def leading(series):
    """(power of eps, coefficient) of the lowest power of u of a Laurent series, None if it is exactly zero;
    raises if every power read below the exactness order vanishes but the series is not exact."""
    terms = {power: sp.factor(value) for power, value in series.terms.items() if sp.simplify(value) != 0}
    if not terms:
        if series.order < large.EXACT:
            raise ValueError('leading term beyond the exact range of the series')
        return None
    lowest = min(terms)
    return -lowest, terms[lowest]


# ---------------------------------------------------------------------------
# 4. Graphs, routing choices and contributions
# ---------------------------------------------------------------------------

VARIANTS = ('r < 1', 'r = 1', 'r > 1')


@functools.lru_cache(maxsize=None)
def routing_choices(variant):
    """Large-eps routing choices: at r = 1, and at the sample point closest to r = 1 of each neighbouring region
    (the choice is constant in each region)."""
    if variant == 'r = 1':
        return large.choices_at(sp.Integer(1))
    merged, _ = large.regions()
    for interval in merged:
        if variant == 'r < 1' and interval['upper'] == 1:
            return large.choices_at(interval['points'][-1])
        if variant == 'r > 1' and interval['lower'] == 1:
            return large.choices_at(interval['points'][0])
    raise KeyError(variant)


def contribution_factor(group):
    """Factor that turns X_Gamma (1PI) or dZ into its share of M_1^est/M_0."""
    if group == 'one_particle_irreducible':
        return 1 / TREE_UP_ONLY
    if group == 'sigma':
        return sp.Rational(1, 2)
    return sp.Integer(1)          # (dZ_pi+ + dZ_pi-)/2 = dZ_pi


@functools.lru_cache(maxsize=None)
def graph_scalars(variant):
    """{label: up-only exact share of M_1^est/M_0 (scalar, before the pole data are inserted)}."""
    choices = routing_choices(variant)
    result = {}
    for entry in traces.entries():
        value = large.selected_value(entry, choices)
        result[entry['label']] = sp.expand(up_only_scalar(value) * contribution_factor(entry['group']))
    return result


@functools.lru_cache(maxsize=None)
def graph_series(variant, order=ROOT_ORDER):
    """{label: Laurent series in u of the share of M_1^est/M_0} (exactly zero shares give an empty exact series)."""
    return {label: series_in_u(value, order, variant) for label, value in graph_scalars(variant).items()}


@functools.lru_cache(maxsize=None)
def graph_leading(variant, order=ROOT_ORDER):
    """{label: (power of eps, coefficient) or None} of each share of M_1^est/M_0."""
    return {label: leading(series) for label, series in graph_series(variant, order).items()}


@functools.lru_cache(maxsize=None)
def total_leading(variant, order=ROOT_ORDER):
    """(power of eps, coefficient, graphs at that power) of M_1^est/M_0, read on the sum of the series of all
    graphs (the sum is exact below the smallest order of its parts)."""
    total = large.Laurent({}, large.EXACT)
    for series in graph_series(variant, order).values():
        total = total + series
    power, coefficient = leading(total)
    labels = sorted(label for label, value in graph_leading(variant, order).items()
                    if value is not None and value[0] == power)
    return power, sp.factor(coefficient), labels


def exactness_check(variant):
    """The leading term of M_1^est/M_0 is the same with the pole root through K and K + 2."""
    return total_leading(variant, ROOT_ORDER)[:2] == total_leading(variant, ROOT_ORDER + 2)[:2]


def sum_of_graphs_check(variant):
    power, coefficient, labels = total_leading(variant)
    leading_terms = graph_leading(variant)
    return sp.simplify(sum(leading_terms[label][1] for label in labels) - coefficient) == 0


def mass_difference():
    """(power of eps, coefficient) of m_a^2 - m_sigma^2 (up-only)."""
    difference = sp.expand(alp_mass_squared() - radial_mass_squared())
    terms = {}
    for term in sp.Add.make_args(difference):
        coefficient, power = term.as_coeff_exponent(U)
        terms[int(power)] = terms.get(int(power), 0) + coefficient
    terms = {power: value for power, value in terms.items() if sp.expand(value) != 0}
    return -min(terms), sp.factor(terms[min(terms)])


def zero_graphs(variant):
    return sorted(label for label, value in graph_scalars(variant).items() if value == 0)


def generic_up_only_value(numerator):
    """Up-only value of a generic leading numerator divided by the tree trace <S^C_10 {lambda_pi+, lambda_pi-}>."""
    value = up_only_scalar(numerator.map_coefficients(
        lambda c: sp.sympify(c).subs({large.KAPPA: ALP_COEFFICIENT_UP_ONLY, large.A: ALP_COEFFICIENT_UP_ONLY})))
    return sp.factor(value / up_only_scalar(large.TREE_TRACE))


@functools.lru_cache(maxsize=None)
def generic_comparison(variant):
    """Cross-check with the generic large-eps results of the same region (task 4.3):
    'total': (generic leading power, up-only value of the generic leading coefficient);
    'graphs': {label: up-only value of the generic leading coefficient = up-only leading coefficient} for the
    graphs whose up-only leading power equals their generic leading power."""
    power, numerator, _ = large.total_leading(variant)
    results = large.graph_results(variant)
    agreement = {}
    for label, value in graph_leading(variant).items():
        generic_power, generic_numerator = results[label]['contribution']
        if value is None or value[0] != generic_power:
            continue
        agreement[label] = sp.simplify(generic_up_only_value(generic_numerator) - value[1]) == 0
    return {'total': (power, generic_up_only_value(numerator)), 'graphs': agreement}


# ---------------------------------------------------------------------------
# 5. Output
# ---------------------------------------------------------------------------

def latex(value):
    return cartesian.coefficient_latex(value)


def markdown():
    power, coefficient, labels = total_leading('r = 1')
    lines = ['| graph | share of $M_1^{\\rm est}/M_0$: leading power | leading coefficient |', '|---|---|---|']
    for label, value in graph_leading('r = 1').items():
        if value is None:
            continue
        lines.append(f"| ${traces.label_latex(label)}$ | $\\epsilon^{{{value[0]}}}$ | ${latex(value[1])}$ |")
    lines.append('')
    lines.append('Exactly zero: ' + ', '.join(f'${traces.label_latex(label)}$' for label in zero_graphs('r = 1'))
                 + '.')
    lines.append('')
    lines.append(f"$$\n\\frac{{M_1^{{\\rm est}}}}{{M_0}}\\to {latex(coefficient)}\\,\\epsilon^{{{power}}},"
                 f"\\qquad R\\to {latex(sp.Abs(coefficient))}\\,\\epsilon^{{{power}}}\n$$")
    lines.append('')
    lines.append('Graphs at the leading power: ' + ', '.join(f'${traces.label_latex(label)}$' for label in labels)
                 + '.')
    lines.append('')
    power_difference, coefficient_difference = mass_difference()
    lines.append(f'$m_a^2-m_\\sigma^2\\to{latex(coefficient_difference)}\\,\\epsilon^{{{power_difference}}}$, '
                 f'so $\\max(m_a,m_\\sigma)=m_a$.')
    lines.append('')
    lines.append('| routing choices and branch | leading power | leading coefficient | graphs | generic leading '
                 'power of the region | its coefficient for up-only sources | graphs where the up-only and generic '
                 'leading terms coincide |')
    lines.append('|---|---|---|---|---|---|---|')
    for variant in VARIANTS:
        p, c, graph_labels = total_leading(variant)
        comparison = generic_comparison(variant)
        generic_power, generic_value = comparison['total']
        agreeing = [label for label, ok in comparison['graphs'].items() if ok]
        assert all(comparison['graphs'].values())
        lines.append(f'| ${variant}$ | $\\epsilon^{{{p}}}$ | ${latex(c)}$ | '
                     + ', '.join(f'${traces.label_latex(label)}$' for label in graph_labels)
                     + f' | $\\epsilon^{{{generic_power}}}$ | ${latex(generic_value)}$ | '
                     + ', '.join(f'${traces.label_latex(label)}$' for label in agreeing) + ' |')
    return '\n'.join(lines)


def main():
    if '--markdown' in sys.argv:
        print(markdown())
        return
    projector_part, identity_part = inverse_components()
    print('M_s^{-1}(C): coefficient of P_u =', projector_part, ';  coefficient of 1 =', identity_part)
    print('m_a^2 (series) =', alp_mass_squared())
    print('Residual of D_a at the root: lowest power of u =', root_residual(), '(needs >=', 2 * ROOT_ORDER + 4, ')')
    print('m_a^2 - m_sigma^2: eps power and coefficient', mass_difference())
    print('Tree (up-only, from the Cartesian tree):', cartesian_tree.up_only_tree()[0])
    for variant in VARIANTS:
        print(f'\n== routing choices of {variant}')
        for label, value in graph_leading(variant).items():
            print(f'  {label}: {value}')
        power, coefficient, labels = total_leading(variant)
        print(f'  M_1^est/M_0 -> {coefficient} eps^{power}; graphs {labels}')
        print('  exact through the series:', exactness_check(variant), '; sum of graphs:', sum_of_graphs_check(variant))
        comparison = generic_comparison(variant)
        print('  generic leading term of the region (power, up-only value):', comparison['total'])
        print('  graphs with the same up-only and generic power, coefficients agree:', comparison['graphs'])


if __name__ == '__main__':
    main()
