"""Cartesian frame: size estimate M_1^est of the one-loop sigma -> pi+ pi- amplitude and the ratio R (task 4.1).

Run:
    PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_cartesian_loop_estimate.py
    PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_cartesian_loop_estimate.py --markdown
    PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_cartesian_loop_estimate

OpenSpec change estimate-sigma-pipi-one-loop-murayama-cartesian.  Units GeV;
conventions of sigma_pipi_cartesian_loop_algebra.py (Cartesian fields
phi = f_a + s + i p, eps = f_pi/f_a, M = V for a contact graph).

    M_1^est = sum_Gamma X_Gamma^max + M_0 (dZ_s^max + dZ_pi+^max + dZ_pi-^max)/2,
    R = |M_1^est| / |M_0|,   M_0 = -(B_0 eps/(2 f_pi)) <S^C_10 {lambda_pi+, lambda_pi-}>.

X_Gamma^max and dZ^max are the termwise routing maxima (D12) with the per-term
loop scales (D3) of derive_sigma_pipi_cartesian_loop_traces.py (stage B);
they are taken from estimate_inputs() and entries() of that script.  Where
the tree trace <S^C_10 {lambda_pi+, lambda_pi-}> vanishes, R is undefined
(no division) and only M_1^est is reported.

No loop integral, scalar loop function, loop coefficient or loop amplitude is
computed: every loop enters through S_Gamma (m_term^2)^{2-r}/(16 pi^2) per term.
The result is an order of magnitude per graph class and can miss
cancellations between graphs with the same epsilon scaling.
"""
import functools
import itertools
import sys

import sympy as sp

import sigma_pipi_cartesian_loop_algebra as cartesian
import derive_sigma_pipi_cartesian_loop_traces as traces
from sigma_pipi_loop_algebra import (B0, F_PI, M_A, M_MESON, M_PI, M_SIGMA, TraceSum, assert_no_float)
from sigma_pipi_cartesian_loop_algebra import SINGLET_MASS_SQUARED

TREE_TRACE = TraceSum.single(('SC10', 'l+', 'l-')) + TraceSum.single(('SC10', 'l-', 'l+'))   # <S^C_10 {l+, l-}>
SCALE_FACTOR = sp.Symbol('c', positive=True)                    # common rescaling of every mass
# Mass dimension of the quantities that are rescaled (sources are rescaled letter by letter).
DIMENSIONFUL = {M_SIGMA: 1, M_A: 1, M_PI: 1, M_MESON: 1, B0: 1, F_PI: 1, SINGLET_MASS_SQUARED: 2}


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------

@functools.lru_cache(maxsize=None)
def _estimate():
    inputs = traces.estimate_inputs()
    pi_minus = TraceSum()
    for value in traces.pi_minus_residue_corrections().values():
        pi_minus = pi_minus + value
    parts = {
        'one particle irreducible': dict(inputs['one particle irreducible']),
        'delta Z sigma': inputs['delta Z sigma'].simplified(),
        'delta Z pi+': inputs['delta Z pi'].simplified(),
        'delta Z pi-': pi_minus.simplified(),
        'tree': inputs['tree'].simplified(),
    }
    external = (parts['tree'] * (parts['delta Z sigma'] + parts['delta Z pi+'] + parts['delta Z pi-'])).scale(
        sp.Rational(1, 2)).simplified()
    one_particle_irreducible = TraceSum()
    for value in parts['one particle irreducible'].values():
        one_particle_irreducible = one_particle_irreducible + value
    parts['one particle irreducible sum'] = one_particle_irreducible.simplified()
    parts['external legs'] = external
    parts['total'] = (one_particle_irreducible + external).simplified()
    return parts


def estimate():
    """{'one particle irreducible': {label: X^max}, 'one particle irreducible sum', 'delta Z sigma',
    'delta Z pi+', 'delta Z pi-', 'tree': M_0, 'external legs': M_0 (dZ_s + dZ_pi+ + dZ_pi-)/2,
    'total': M_1^est}, all TraceSums of words (cached)."""
    return _estimate()


def sum_of_parts_check():
    """M_1^est - sum_Gamma X_Gamma^max - M_0 (dZ_s + dZ_pi+ + dZ_pi-)/2 = 0, computed term by term."""
    parts = estimate()
    difference = parts['total']
    for value in parts['one particle irreducible'].values():
        difference = difference - value
    for name in ('delta Z sigma', 'delta Z pi+', 'delta Z pi-'):
        difference = difference - (parts['tree'] * parts[name]).scale(sp.Rational(1, 2))
    return difference.is_zero()


def pion_residues_equal():
    return (estimate()['delta Z pi+'] - estimate()['delta Z pi-']).is_zero()


# ---------------------------------------------------------------------------
# Ratio R
# ---------------------------------------------------------------------------

def trace_symbol(word):
    """Commuting symbol standing for the trace of a canonical word (generic-source basis)."""
    return sp.Symbol('<' + ' '.join(word) + '>')


def as_expression(trace_sum):
    """A TraceSum as a sympy expression in the trace symbols (products of traces are products of symbols)."""
    total = sp.Integer(0)
    for key, value in trace_sum.terms.items():
        total += value * sp.Mul(*(trace_symbol(word) for word in key))
    return total


def ratio(total=None, tree=None):
    """R = |M_1^est|/|M_0| as a sympy expression, or None if the tree vanishes (R undefined, no division)."""
    parts = estimate()
    total = parts['total'] if total is None else total
    tree = parts['tree'] if tree is None else tree
    tree_value = sp.expand(as_expression(tree.simplified()))
    if tree_value == 0:
        return None
    return sp.Abs(as_expression(total)) / sp.Abs(tree_value)


def set_letter_to_zero(trace_sum, name):
    """The TraceSum with the matrix `name` set to zero (every product containing it removed)."""
    return trace_sum.filtered(lambda key: not any(name in word for word in key)).simplified()


def vanishing_tree_case():
    """(tree, R, M_1^est) for S^C_10 = 0: the tree vanishes, R is None, M_1^est is still nonzero."""
    parts = estimate()
    tree = set_letter_to_zero(parts['tree'], 'SC10')
    total = set_letter_to_zero(parts['total'], 'SC10')
    return tree, ratio(total, tree), total


# ---------------------------------------------------------------------------
# Dimensions
# ---------------------------------------------------------------------------

def factor_out_scale(expression, scale=SCALE_FACTOR):
    """max(c x, c y, ...) -> c max(x, y, ...) for c > 0 (checked separately by maximum_is_homogeneous)."""
    return sp.sympify(expression).replace(
        lambda node: isinstance(node, sp.Max),
        lambda node: scale * sp.Max(*(sp.cancel(argument / scale) for argument in node.args)))


def maximum_is_homogeneous():
    """max(c x_1, ..., c x_n) = c max(x_1, ..., x_n) for c > 0, checked for n = 2, 3 in every ordering
    of the arguments, written with positive gaps (x_(1) = y, x_(2) = y + d_1, x_(3) = y + d_1 + d_2)."""
    c = SCALE_FACTOR
    smallest = sp.Symbol('y', positive=True)
    gaps = sp.symbols('d1 d2', positive=True)
    for count in (2, 3):
        ordered = [smallest + sum(gaps[:k], sp.Integer(0)) for k in range(count)]
        for arguments in itertools.permutations(ordered):
            if sp.Max(*(c * argument for argument in arguments)) - c * sp.Max(*arguments) != 0:
                return False
            if sp.Max(*arguments) != ordered[-1]:
                return False
    return True


def rescaled(trace_sum):
    """Every mass, B_0, f_pi and every source letter multiplied by c (m_1^2 by c^2);
    eps, N(m_a^2), M^{-1}_{m_a^2}(C) and lambda_pi+- are dimensionless."""
    c = SCALE_FACTOR
    substitution = {symbol: c ** power * symbol for symbol, power in DIMENSIONFUL.items()}
    result = TraceSum()
    for key, value in trace_sum.terms.items():
        value = factor_out_scale(sp.sympify(value).subs(substitution, simultaneous=True))
        result.add_term(key, value * c ** cartesian.source_count(key))
    return result


def piecewise_conditions(trace_sum):
    conditions = set()
    for value in trace_sum.terms.values():
        for piecewise in sp.sympify(value).atoms(sp.Piecewise):
            conditions |= {condition for _, condition in piecewise.args if condition is not sp.true}
    return conditions


def conditions_are_homogeneous(trace_sum):
    """Each Piecewise condition lhs >= rhs has lhs - rhs -> c^k (lhs - rhs) under the rescaling."""
    c = SCALE_FACTOR
    substitution = {symbol: c ** power * symbol for symbol, power in DIMENSIONFUL.items()}
    for condition in piecewise_conditions(trace_sum):
        difference = condition.lhs - condition.rhs
        scaled = factor_out_scale(difference.subs(substitution, simultaneous=True))
        quotient = sp.simplify(scaled / difference)
        if quotient.free_symbols - {c} or not (quotient.is_Pow and quotient.base == c or quotient == c):
            return False
    return True


def branch_values(trace_sum):
    """The TraceSum in every combination of Piecewise branches (one branch per distinct condition).
    The conditions are homogeneous, so selecting a branch commutes with the rescaling."""
    conditions = sorted(piecewise_conditions(trace_sum), key=sp.default_sort_key)
    results = []
    for truth in itertools.product((True, False), repeat=len(conditions)):
        assignment = dict(zip(conditions, truth))

        def choose(*arguments):
            for expression, condition in arguments:
                if condition is sp.true or assignment.get(condition, False):
                    return expression
            raise ValueError(arguments)
        results.append(trace_sum.map_coefficients(
            lambda value: sp.sympify(value).replace(lambda node: isinstance(node, sp.Piecewise),
                                                    lambda node: choose(*node.args))))
    return results


def has_dimension(trace_sum, dimension):
    """X(c x) = c^dimension X(x) in every branch."""
    for value in branch_values(trace_sum):
        if not (rescaled(value) - value.scale(SCALE_FACTOR ** dimension)).is_zero():
            return False
    return True


def dimension_checks():
    """{name: True/False}: X_Gamma^max ~ c, dZ^max ~ c^0, M_0 ~ c, M_1^est ~ c, conditions homogeneous."""
    parts = estimate()
    checks = {'max(c x, ...) = c max(x, ...)': maximum_is_homogeneous()}
    for label, value in parts['one particle irreducible'].items():
        checks[f'{label}: dimension 1'] = has_dimension(value, 1)
        checks[f'{label}: homogeneous conditions'] = conditions_are_homogeneous(value)
    for name in ('delta Z sigma', 'delta Z pi+', 'delta Z pi-'):
        checks[f'{name}: dimension 0'] = has_dimension(parts[name], 0)
    checks['M_0: dimension 1'] = has_dimension(parts['tree'], 1)
    checks['M_1^est: dimension 1'] = has_dimension(parts['total'], 1)
    checks['M_1^est: homogeneous conditions'] = conditions_are_homogeneous(parts['total'])
    return checks


def no_float():
    parts = estimate()
    for name in ('one particle irreducible sum', 'delta Z sigma', 'delta Z pi+', 'delta Z pi-', 'tree', 'total'):
        assert_no_float(parts[name])
    return True


# ---------------------------------------------------------------------------
# Summary: letters, number of terms, S^C_00
# ---------------------------------------------------------------------------

def source_letters(trace_sum):
    return sorted((name for name in trace_sum.letters() if cartesian.is_cartesian_source(name)),
                  key=cartesian.source_indices)


def letters_by_graph():
    """{label: source letters} of the 1PI graphs and of the kept self-energies."""
    result = {}
    for entry in traces.entries():
        result[entry['label']] = source_letters(entry['maximum']['value'])
    return result


def products_with(trace_sum, name):
    return [key for key in trace_sum.simplified().terms if any(name in word for word in key)]


def letter_part(trace_sum, name):
    """The products of traces of a TraceSum that contain the letter `name`."""
    return trace_sum.filtered(lambda key: any(name in word for word in key)).simplified()


def source_by_scale(name='SC00'):
    """Where a source letter enters the kept graphs, loop scale by loop scale.

    Uses the displayed form of the termwise maximum (base routing per scale plus the
    routing-dependent corrections).  Returns (table, vanish, absent):
      table  {scale: {label: number of products with the letter}} in the base routing;
      vanish True if, at every scale, the sum over the graphs of the products with the letter is zero;
      absent True if the letter does not enter any routing-dependent correction."""
    table = {}
    sums = {}
    absent_from_corrections = True
    for entry in traces.entries():
        display = entry['display']
        for scale, value in display['base by scale'].items():
            part = letter_part(value.to_trace_sum(), name)
            if part.terms:
                table.setdefault(scale, {})[entry['label']] = len(part.terms)
                sums[scale] = sums.get(scale, TraceSum()) + part
        for _, _, difference in display['corrections']:
            if letter_part(difference.to_trace_sum(), name).terms:
                absent_from_corrections = False
    vanish = all(value.simplified().is_zero() for value in sums.values())
    return table, vanish, absent_from_corrections


def summary():
    """Counts and letters of M_1^est and the status of S^C_00."""
    parts = estimate()
    graphs_with_s00 = [label for label, letters in letters_by_graph().items() if 'SC00' in letters]
    s00_one_particle_irreducible = products_with(parts['one particle irreducible sum'], 'SC00')
    s00_external = products_with(parts['external legs'], 'SC00')
    return {
        'products in M_1^est': len(parts['total'].terms),
        'products in the 1PI sum': len(parts['one particle irreducible sum'].terms),
        'products in the external-leg term': len(parts['external legs'].terms),
        'products in dZ_s^max': len(parts['delta Z sigma'].terms),
        'products in dZ_pi^max': len(parts['delta Z pi+'].terms),
        'kept 1PI graphs': len(parts['one particle irreducible']),
        'kept s self-energies': sum(1 for entry in traces.entries() if entry['group'] == 'sigma'),
        'kept pi self-energies (each charged leg)': sum(1 for entry in traces.entries() if entry['group'] == 'pi'),
        'source letters in M_1^est': source_letters(parts['total']),
        'graphs with S^C_00': graphs_with_s00,
        'products with S^C_00 in M_1^est': len(products_with(parts['total'], 'SC00')),
        'products with S^C_00 in the 1PI sum': len(s00_one_particle_irreducible),
        'products with S^C_00 in the external-leg term': len(s00_external),
        'Piecewise products in M_1^est': sum(1 for value in parts['total'].terms.values()
                                             if sp.sympify(value).has(sp.Piecewise)),
    }


def letter_text(name):
    return '$' + cartesian.letter_latex(name) + '$'


def markdown_summary():
    """Script output for the estimate record (the test checks that the record contains it)."""
    data = summary()
    lines = ['| quantity | value |', '|---|---|']
    for name in ('kept 1PI graphs', 'kept s self-energies', 'kept pi self-energies (each charged leg)',
                 'products in the 1PI sum', 'products in dZ_s^max', 'products in dZ_pi^max',
                 'products in the external-leg term', 'products in M_1^est', 'Piecewise products in M_1^est',
                 'products with S^C_00 in the 1PI sum', 'products with S^C_00 in the external-leg term',
                 'products with S^C_00 in M_1^est'):
        lines.append(f'| {name} | {data[name]} |')
    lines.append('| source letters in $M_1^{\\rm est}$ | ' + ', '.join(letter_text(name) for name in
                                                                      data['source letters in M_1^est']) + ' |')
    lines.append('| graphs with $\\widehat S^C_{0,0}$ | ' + ', '.join('$' + traces.label_latex(label) + '$'
                                                                    for label in data['graphs with S^C_00']) + ' |')
    table, vanish, absent = source_by_scale('SC00')
    for scale in sorted(table, key=sp.default_sort_key):
        lines.append(f'| products with $\\widehat S^C_{{0,0}}$ at $m_{{\\rm term}}={traces.scale_latex(scale)}$ | '
                     + ', '.join(f'${traces.label_latex(label)}$: {count}' for label, count in sorted(table[scale].items()))
                     + ' |')
    lines.append(f'| $\\widehat S^C_{{0,0}}$ part cancels at each $m_{{\\rm term}}$ | {vanish} |')
    lines.append(f'| $\\widehat S^C_{{0,0}}$ absent from the routing-dependent branches | {absent} |')
    return '\n'.join(lines)


def main():
    if '--markdown' in sys.argv:
        print(markdown_summary())
        return
    print('M_1^est = sum of its parts:', sum_of_parts_check())
    print('dZ_pi+^max = dZ_pi-^max:', pion_residues_equal())
    print('No float:', no_float())
    for name, value in dimension_checks().items():
        if not value:
            print('  FAILED:', name)
    print('All dimension checks:', all(dimension_checks().values()))
    tree, value, total = vanishing_tree_case()
    print('S^C_10 = 0: tree zero:', tree.is_zero(), '; R:', value, '; M_1^est products:', len(total.terms))
    for name, value in summary().items():
        print(f'{name}: {value}')
    print('M_0 =', cartesian.trace_sum_latex(estimate()['tree']))


if __name__ == '__main__':
    main()
