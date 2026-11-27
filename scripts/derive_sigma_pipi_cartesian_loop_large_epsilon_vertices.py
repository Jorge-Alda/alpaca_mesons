"""Large epsilon at fixed Cartesian sources: pole data, meson-line residue pieces, leg factors and vertex powers (task 2.3).

Run:
    PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_cartesian_loop_large_epsilon_vertices.py
    PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_cartesian_loop_large_epsilon_vertices

OpenSpec change estimate-sigma-pipi-one-loop-murayama-cartesian.  Units GeV
(s, m_a^2, m_sigma^2, m_1^2 in GeV^2; B_0, f_pi and the sources in GeV).
Conventions of sigma_pipi_cartesian_loop_algebra.py.

Regime (D10, restated in the Cartesian invariants; assumptions):
 * held fixed as eps = f_pi/f_a grows: f_pi, B_0, m_1^2, m_P, every source
   letter S^C_{r,q}, P^C_{r,q}, the matrices lambda_{pi+-}, and the quartic
   coupling lambda of V_phi; f_a = f_pi/eps varies;
 * A = B_0(<S^C_10> - 2<S^C_02>) > 0 and kappa_sigma = B_0(<S^C_10> - 2<S^C_20>) > 0;
 * eps^2 A lies above the spectrum of L(X) = B_0 {S^C_00, X} + (m_1^2/3) <X> 1, so
   that M_s^{-1} = (L - s)^{-1} = -sum_k L^k / s^{k+1} at s = m_a^2;
 * generic sources: every independent product of traces is nonzero, so the
   powers are upper bounds unless a trace combination vanishes.

Every quantity is exact (derive_sigma_pipi_cartesian_loop_tree.py):
  D_a(s) = s - eps^2 A + 2 eps^2 <C M_s^{-1}(C)>,  N(s) = 1 + 2 eps^2 <(M_s^{-1}(C))^2>,
  m_sigma^2 = 2 lambda f_pi^2/eps^2 + eps^2 kappa_sigma,
  v_a = (2 eps M^{-1}_{m_a^2}(C)/sqrt(N), 1/sqrt(N)),  C = B_0 P^C_01.
Series in 1/eps are used only to read the leading power and coefficient;
c_k = <C L^k(C)> are kept up to k = 3 and the omitted terms are of higher
order (c_k enters with eps^{-2k} relative to c_0).  In the vertices lambda is
first eliminated through m_sigma^2 (D-A3); the leading powers are read after
m_sigma^2 = 2 lambda f_pi^2/eps^2 + eps^2 kappa_sigma and the exact rewriting
<S^C_10> = kappa_sigma/B_0 + 2 <S^C_20>, <S^C_02> = (kappa_sigma - A)/(2 B_0) + <S^C_20>,
with lambda fixed.  No loop is evaluated.
"""
import itertools

import sympy as sp

import sigma_pipi_cartesian_loop_algebra as cartesian
import derive_sigma_pipi_cartesian_loop_tree as cartesian_tree
from sigma_pipi_loop_algebra import B0, EPSILON, F_PI, M_SIGMA, QUARTIC, TraceSum, symbolic_dot
from sigma_pipi_cartesian_loop_algebra import (ALP_MASS_COEFFICIENT, KAPPA_SIGMA, SINGLET_MASS_SQUARED,
                                               UP_QUARK_MASS)

INVERSE_EPSILON = sp.Symbol('u', positive=True)              # u = 1/eps
SQUARED_MOMENTUM = sp.Symbol('s', positive=True)
LOOP_OPERATOR = sp.Symbol('L', commutative=False)
ORDER = 3                                                    # c_0 ... c_3
TRACE_C_POWER = sp.symbols(f'c_0:{ORDER + 1}', real=True)     # c_k = <C L^k(C)>
MATRIX_C_POWER = sp.symbols(f'LC_0:{ORDER + 1}')              # the matrices L^k(C) (linear combinations only)
MATRIX_C_POWER_RIGHT = sp.symbols(f"LC'_0:{ORDER + 1}")       # second factor of a tensor product
TRACE_C_SQUARED = sp.Symbol('<PC01 PC01>', positive=True)     # <P^C_01 P^C_01>
P01 = sp.Symbol('P^C_01')                                     # the matrix P^C_01 (linear combinations only)
P01_RIGHT = sp.Symbol("P^C_01'")
SERIES_ORDER = 6                                              # powers of u kept when reading leading terms


# ---------------------------------------------------------------------------
# 1. Regime: large-s expansion and self-adjointness of L
# ---------------------------------------------------------------------------

def resolvent_series_check(order=ORDER):
    """(L - s)(-(1/s)) sum_{k<=K} (L/s)^k = 1 - (L/s)^{K+1}."""
    s = SQUARED_MOMENTUM
    series = -sum((LOOP_OPERATOR / s) ** k for k in range(order + 1)) / s
    return sp.expand((LOOP_OPERATOR - s) * series - (1 - (LOOP_OPERATOR / s) ** (order + 1))) == 0


def loop_operator_self_adjoint():
    """<X L(Y)> = <L(X) Y>, so <L^j(C) L^k(C)> = c_{j+k}."""
    def trace_with_operator(left, argument):
        unit = TraceSum()
        unit.add_term(((left,), (argument,)), SINGLET_MASS_SQUARED / 3)
        return (TraceSum.single((left, 'SC00', argument), B0) + TraceSum.single((left, argument, 'SC00'), B0)
                + unit)
    return (trace_with_operator('X', 'Y') - trace_with_operator('Y', 'X')).is_zero()


def large_s_traces(s):
    """<C M_s^{-1}(C)> = -sum c_k/s^{k+1};  <(M_s^{-1}(C))^2> = sum (n+1) c_n/s^{n+2}."""
    c_minv = -sum(TRACE_C_POWER[k] / s ** (k + 1) for k in range(ORDER + 1))
    minv_squared = sum((n + 1) * TRACE_C_POWER[n] / s ** (n + 2) for n in range(ORDER + 1))
    return c_minv, minv_squared


# ---------------------------------------------------------------------------
# 2. Pole data
# ---------------------------------------------------------------------------

def alp_pole_root():
    """Root of D_a(s) = 0 as a series in u = 1/eps: m_a^2 = A/u^2 + b_0 + b_1 u^2 + b_2 u^4."""
    u = INVERSE_EPSILON
    s = SQUARED_MOMENTUM
    unknowns = sp.symbols('b_0:3')
    root = ALP_MASS_COEFFICIENT / u ** 2 + sum(unknowns[k] * u ** (2 * k) for k in range(3))
    c_minv, _ = large_s_traces(s)
    # u^2 D_a(s) = u^2 s - A + 2 <C M_s^{-1}(C)>
    equation = sp.series((u ** 2 * s - ALP_MASS_COEFFICIENT + 2 * c_minv).subs(s, root), u, 0, 8).removeO()
    equation = sp.expand(equation)
    solution = sp.solve([equation.coeff(u, 2 * k) for k in range(1, 4)], unknowns, dict=True)[0]
    return sp.expand(root.subs(solution))


def radial_mass_squared():
    return 2 * QUARTIC * F_PI ** 2 * INVERSE_EPSILON ** 2 + KAPPA_SIGMA / INVERSE_EPSILON ** 2


def series_in_u(expression, order=SERIES_ORDER):
    return sp.expand(sp.series(expression, INVERSE_EPSILON, 0, order).removeO())


def leading_term(expression):
    """(power of eps, coefficient) of the leading term of a Laurent polynomial in u = 1/eps."""
    expression = sp.expand(expression)
    terms = {}
    for term in sp.Add.make_args(expression):
        coefficient, exponent = term.as_coeff_exponent(INVERSE_EPSILON)
        terms[exponent] = terms.get(exponent, 0) + coefficient
    terms = {exponent: sp.factor(value) for exponent, value in terms.items() if sp.simplify(value) != 0}
    lowest = min(terms)
    return -lowest, terms[lowest]


def pole_quantities():
    """Series in u of every pole quantity, before reading the leading term."""
    u = INVERSE_EPSILON
    root = alp_pole_root()
    c_minv, minv_squared = large_s_traces(SQUARED_MOMENTUM)
    normalization = series_in_u((1 + 2 * minv_squared / u ** 2).subs(SQUARED_MOMENTUM, root))
    # M^{-1}_{m_a^2}(C) = -sum_k L^k(C)/m_a^{2(k+1)} as a linear combination of the matrices L^k(C)
    minv = series_in_u(sum(-MATRIX_C_POWER[k] / root ** (k + 1) for k in range(ORDER + 1)))
    minv_right = minv.subs(dict(zip(MATRIX_C_POWER, MATRIX_C_POWER_RIGHT)))
    inverse_normalization = series_in_u(1 / normalization)
    root_inverse_normalization = series_in_u(1 / sp.sqrt(normalization))
    derivative = sp.diff(SQUARED_MOMENTUM - (ALP_MASS_COEFFICIENT - 2 * c_minv) / u ** 2, SQUARED_MOMENTUM)
    return {
        'm_a^2': root,
        'm_sigma^2': radial_mass_squared(),
        'N': normalization,
        "N from D_a'": series_in_u(derivative.subs(SQUARED_MOMENTUM, root)),
        'M^{-1}(C)': minv,
        '1/N': inverse_normalization,
        'N - 1': series_in_u(normalization - 1),
        'xi_a = 1/sqrt(N)': root_inverse_normalization,
        'meson matrix of v_a': series_in_u(2 / u * minv * root_inverse_normalization),
        '(i) 1_9': sp.Integer(1),
        '(ii) 1 - 1/N': series_in_u(1 - inverse_normalization),
        '(iii) -2 eps M^{-1}(C)/N': series_in_u(-2 / u * minv * inverse_normalization),
        '(iv) -4 eps^2 M^{-1}(C) x M^{-1}(C)/N': series_in_u(-4 / u ** 2 * minv * minv_right * inverse_normalization),
    }


def in_sources(value):
    """Write C = B_0 P^C_01 and c_0 = B_0^2 <P^C_01 P^C_01>."""
    return sp.factor(sp.sympify(value).subs({MATRIX_C_POWER[0]: B0 * P01, MATRIX_C_POWER_RIGHT[0]: B0 * P01_RIGHT,
                                             TRACE_C_POWER[0]: B0 ** 2 * TRACE_C_SQUARED}))


def pole_data_table():
    """{name: (leading power of eps, leading coefficient in the sources)} for generic sources."""
    quantities = pole_quantities()
    table = {}
    for name, value in quantities.items():
        if name == "N from D_a'":
            continue
        power, coefficient = leading_term(value)
        table[name] = (power, in_sources(coefficient))
    difference = quantities['m_a^2'] - quantities['m_sigma^2']
    table['m_a^2 - m_sigma^2 (A != kappa_sigma)'] = leading_term(difference)
    power, coefficient = leading_term(difference.subs(KAPPA_SIGMA, ALP_MASS_COEFFICIENT))
    table['m_a^2 - m_sigma^2 (A = kappa_sigma)'] = (power, in_sources(coefficient))
    return table


# ---------------------------------------------------------------------------
# 3. Up-only pole data (D-A4: S^C_00 = S^C_10 = m_u P_u, P^C_01 = -m_u P_u)
# ---------------------------------------------------------------------------

def operator_terms(terms):
    """L(X) = B_0 {S^C_00, X} + (m_1^2/3) <X> 1 on a matrix written as [(coefficient, word, trace words)]."""
    result = []
    for coefficient, word, traces in terms:
        result.append((B0 * coefficient, ('SC00',) + word, traces))
        result.append((B0 * coefficient, word + ('SC00',), traces))
        result.append((SINGLET_MASS_SQUARED / 3 * coefficient, (), traces + (word,)))
    return result


def trace_of_product(first, second):
    """<first second> as a TraceSum."""
    result = TraceSum()
    for (one, word_one, traces_one), (two, word_two, traces_two) in itertools.product(first, second):
        result.add_term(traces_one + traces_two + (word_one + word_two,), one * two)
    return result


def c_power_trace_sums():
    """c_0, c_1, c_2 as TraceSums in the Cartesian letters (C = B_0 P^C_01)."""
    c_matrix = [(B0, ('PC01',), ())]
    once = operator_terms(c_matrix)
    return [trace_of_product(c_matrix, c_matrix), trace_of_product(c_matrix, once), trace_of_product(once, once)]


def up_only_values():
    """A, kappa_sigma, c_0, c_1, c_2 in the up-only limit (projector relations, <1> = 3, <P_u> = 1)."""
    def scalar(trace_sum):
        reduced = cartesian.up_only_substitution(trace_sum)
        assert all(key == () for key in reduced.terms), reduced
        return sp.expand(reduced.terms.get((), 0))
    values = {'A': scalar(cartesian_tree.alp_mass_coefficient()), 'kappa_sigma': scalar(cartesian_tree.kappa_sigma())}
    for k, trace_sum in enumerate(c_power_trace_sums()):
        values[f'c_{k}'] = scalar(trace_sum)
    return values


def up_only_mass_difference():
    """Series of m_a^2 - m_sigma^2 in the up-only limit and its first nonvanishing term."""
    values = up_only_values()
    substitution = {ALP_MASS_COEFFICIENT: values['A'], KAPPA_SIGMA: values['kappa_sigma'],
                    TRACE_C_POWER[0]: values['c_0'], TRACE_C_POWER[1]: values['c_1'],
                    TRACE_C_POWER[2]: values['c_2']}
    difference = sp.expand((alp_pole_root() - radial_mass_squared()).subs(substitution))
    return difference, leading_term(difference)


def up_only_exact_alp_pole_function():
    """Exact D_a(s) in the up-only limit.  The span of P_u and 1 is invariant under L:
    L(P_u) = 2 B_0 m_u P_u + (m_1^2/3) 1, L(1) = 2 B_0 m_u P_u + m_1^2 1 (S^C_00 = m_u P_u, <P_u> = 1,
    <1> = 3).  M_s(X) = C = -B_0 m_u P_u is solved with X = alpha P_u + beta 1, and
    <C X> = -B_0 m_u (alpha + beta)."""
    alpha, beta = sp.symbols('alpha beta')
    s = SQUARED_MOMENTUM
    m_u, m1 = UP_QUARK_MASS, SINGLET_MASS_SQUARED
    equations = [2 * B0 * m_u * (alpha + beta) - s * alpha + B0 * m_u,
                 m1 / 3 * alpha + m1 * beta - s * beta]
    solution = sp.solve(equations, [alpha, beta], dict=True)[0]
    trace_c_x = -B0 * m_u * (solution[alpha] + solution[beta])
    return sp.together(s - B0 * m_u / INVERSE_EPSILON ** 2 + 2 * trace_c_x / INVERSE_EPSILON ** 2)


def up_only_root_check():
    """The up-only series root (from the c_k) is a zero of the exact D_a(s) through relative order u^6."""
    values = up_only_values()
    substitution = {ALP_MASS_COEFFICIENT: values['A'], TRACE_C_POWER[0]: values['c_0'],
                    TRACE_C_POWER[1]: values['c_1'], TRACE_C_POWER[2]: values['c_2']}
    root = alp_pole_root().subs(substitution)
    value = (INVERSE_EPSILON ** 2 * up_only_exact_alp_pole_function()).subs(SQUARED_MOMENTUM, root)
    return sp.expand(sp.series(value, INVERSE_EPSILON, 0, 8).removeO()) == 0


def up_only_pole_table():
    """Leading terms of the pole data in the up-only limit (P^C_01 = -m_u P_u)."""
    values = up_only_values()
    substitution = {ALP_MASS_COEFFICIENT: values['A'], KAPPA_SIGMA: values['kappa_sigma'],
                    TRACE_C_POWER[0]: values['c_0'], TRACE_C_POWER[1]: values['c_1'],
                    TRACE_C_POWER[2]: values['c_2']}
    table = {}
    for name, value in pole_quantities().items():
        if name in ("N from D_a'",):
            continue
        power, coefficient = leading_term(sp.sympify(value).subs(substitution))
        coefficient = sp.sympify(coefficient).subs({MATRIX_C_POWER[0]: -B0 * UP_QUARK_MASS * sp.Symbol('P_u'),
                                                    MATRIX_C_POWER_RIGHT[0]: -B0 * UP_QUARK_MASS * sp.Symbol("P_u'")})
        table[name] = (power, sp.factor(coefficient))
    return table


# ---------------------------------------------------------------------------
# 4. Vertices by field-end pattern
# ---------------------------------------------------------------------------
FIELD_ENDS = ('radial', 'ALP field', 'meson')


def field_end_leg(kind, index):
    momentum = {f'q{index}': 1}
    if kind == 'radial':
        return cartesian.radial_leg(momentum)
    if kind == 'ALP field':
        return cartesian.alp_field_leg(momentum)
    return cartesian.meson_leg(f'G{index}', momentum)


def field_end_patterns():
    return [pattern for n in (3, 4, 5) for pattern in itertools.combinations_with_replacement(FIELD_ENDS, n)]


def large_epsilon_rewriting(trace_sum):
    """m_sigma^2 = 2 lambda f_pi^2/eps^2 + eps^2 kappa_sigma, <S^C_10> = kappa_sigma/B_0 + 2<S^C_20>,
    <S^C_02> = (kappa_sigma - A)/(2 B_0) + <S^C_20> (exact; lambda fixed)."""
    replacements = {('SC10',): TraceSum.scalar(KAPPA_SIGMA / B0) + TraceSum.single(('SC20',), 2),
                    ('SC02',): TraceSum.scalar((KAPPA_SIGMA - ALP_MASS_COEFFICIENT) / (2 * B0))
                    + TraceSum.single(('SC20',))}
    mass = 2 * QUARTIC * F_PI ** 2 / EPSILON ** 2 + EPSILON ** 2 * KAPPA_SIGMA
    result = TraceSum()
    for key, value in trace_sum.terms.items():
        product = TraceSum.scalar(sp.sympify(value).subs(M_SIGMA, sp.sqrt(mass)))
        for word in key:
            product = product * replacements.get(word, TraceSum.single(word))
        result = result + product
    return result.simplified()


def up_only_rewriting(trace_sum):
    """Up-only limit with m_sigma^2 = 2 lambda f_pi^2/eps^2 + eps^2 B_0 m_u."""
    mass = 2 * QUARTIC * F_PI ** 2 / EPSILON ** 2 + EPSILON ** 2 * B0 * UP_QUARK_MASS
    reduced = cartesian.up_only_substitution(trace_sum)
    return reduced.map_coefficients(lambda value: sp.sympify(value).subs(M_SIGMA, sp.sqrt(mass))).simplified()


def epsilon_exponents(trace_sum, loop_momenta=False):
    """Set of powers of eps over all products of traces (k.k -> eps^2 k.k for loop momenta)."""
    exponents = set()
    for value in trace_sum.terms.values():
        value = sp.sympify(value)
        if loop_momenta:
            value = value.subs({symbol: EPSILON ** 2 * symbol for symbol in value.free_symbols
                                if symbol.name.startswith('d_')}, simultaneous=True)
        for term in sp.Add.make_args(sp.expand(value)):
            exponents.add(term.as_coeff_exponent(EPSILON)[1])
    return exponents


def pattern_vertex(pattern):
    legs = [field_end_leg(kind, index + 1) for index, kind in enumerate(pattern)]
    value = cartesian.vertex(legs, symbolic_dot, eliminate_lambda=True)
    return cartesian.apply_assumptions(value).simplified()


def vertex_power_table():
    """{pattern: {'generic': set, 'generic, loop momenta': set, 'up-only': set, 'up-only, loop momenta': set,
    'rule': set}}; 'rule' is eps^{r+q} for the source part (r radial ends, q ALP-field ends)."""
    table = {}
    for pattern in field_end_patterns():
        value = pattern_vertex(pattern)
        generic_value = large_epsilon_rewriting(value)
        up_only_value = up_only_rewriting(value)
        source = cartesian.apply_assumptions(cartesian.source_vertex(
            [field_end_leg(kind, index + 1) for index, kind in enumerate(pattern)])).simplified()
        ends = pattern.count('radial') + pattern.count('ALP field')
        table[pattern] = {'generic': epsilon_exponents(generic_value),
                          'generic, loop momenta': epsilon_exponents(generic_value, loop_momenta=True),
                          'up-only': epsilon_exponents(up_only_value),
                          'up-only, loop momenta': epsilon_exponents(up_only_value, loop_momenta=True),
                          'rule': {ends} if source.terms else set()}
    return table


def exceptions_to_rule(table):
    """Patterns whose powers differ from eps^{r+q} (generic, fixed momenta)."""
    return {pattern: row for pattern, row in table.items() if row['generic'] != row['rule']}


def lambda_check():
    """V(s,s,s) and V(s;p,p) with lambda fixed: -6 lambda f_pi/eps + 6 B_0 eps^3 <S^C_30>/f_pi and
    -2 lambda f_pi/eps + 2 B_0 eps^3 <S^C_12>/f_pi (generic)."""
    sss = large_epsilon_rewriting(pattern_vertex(('radial',) * 3))
    spp = large_epsilon_rewriting(pattern_vertex(('radial', 'ALP field', 'ALP field')))
    expected_sss = (TraceSum.scalar(-6 * QUARTIC * F_PI / EPSILON) + TraceSum.single(('SC30',), 6 * B0 * EPSILON ** 3 / F_PI))
    expected_spp = (TraceSum.scalar(-2 * QUARTIC * F_PI / EPSILON) + TraceSum.single(('SC12',), 2 * B0 * EPSILON ** 3 / F_PI))
    return (sss - expected_sss).is_zero() and (spp - expected_spp).is_zero()


# ---------------------------------------------------------------------------
# 5. Leg factors with placement, and tree
# ---------------------------------------------------------------------------
# Each ALP-field end at a source vertex carries eps (vertex rule eps^{r+q}); at a V_phi vertex the
# power of the vertex (eps^-1 cubic, eps^0 quartic) does not depend on the number of ALP-field ends.
LEG_FACTORS = (('xi_a = 1/sqrt(N)', 'xi_a = 1/sqrt(N)', 1),
               ('meson matrix of v_a', 'meson matrix of v_a', 0),
               ('(i) 1_9', '(i) 1_9', 0),
               ('(ii) 1 - 1/N', '(ii) 1 - 1/N', 2),
               ('(iii) mixed', '(iii) -2 eps M^{-1}(C)/N', 1),
               ('(iv) meson block', '(iv) -4 eps^2 M^{-1}(C) x M^{-1}(C)/N', 0))


def leg_factor_table(pole_table=None):
    """{factor: (power, {placement: net power})}; placement = number of ALP-field ends at source
    vertices (the other ALP-field ends sit at V_phi vertices)."""
    pole_table = pole_table or pole_data_table()
    table = {}
    for name, key, alp_ends in LEG_FACTORS:
        power = pole_table[key][0]
        placements = {f'{at_source} of {alp_ends} ALP-field ends at source vertices': power + at_source
                      for at_source in range(alp_ends, -1, -1)} if alp_ends else {'': power}
        table[name] = (power, placements)
    return table


def tree_leading():
    """(powers of eps of M_0 generic, M_0 up-only)."""
    tree = cartesian_tree.tree_amplitude()
    up_only, _ = cartesian_tree.up_only_tree()
    return epsilon_exponents(tree), sp.factor(up_only)


def pattern_text(pattern):
    names = {'radial': 's', 'ALP field': 'p', 'meson': 'M'}
    return ''.join(names[kind] for kind in pattern)


def main():
    print('Regime: resolvent series identity:', resolvent_series_check(),
          '; L self-adjoint:', loop_operator_self_adjoint())
    quantities = pole_quantities()
    print("N from <(M^{-1})^2> equals N from D_a'(m_a^2) (series):",
          sp.expand(quantities['N'] - quantities["N from D_a'"]) == 0)
    print('m_a^2 series:', sp.latex(quantities['m_a^2']))
    print('\nPole data and residue pieces (generic): quantity: (power of eps, leading coefficient)')
    pole_table = pole_data_table()
    for name, (power, coefficient) in pole_table.items():
        print(f'  {name}: eps^{power}, {sp.latex(coefficient)}')
    print('\nUp-only values:', up_only_values())
    difference, (power, coefficient) = up_only_mass_difference()
    print(f'Up-only m_a^2 - m_sigma^2 = {sp.latex(difference)}; first nonvanishing: eps^{power} {sp.latex(coefficient)}')
    print('Up-only series root is a zero of the exact up-only D_a(s):', up_only_root_check())
    print('Up-only pole data:')
    for name, (power, coefficient) in up_only_pole_table().items():
        print(f'  {name}: eps^{power}, {sp.latex(coefficient)}')
    print('\nLeg factors (net power by placement):')
    for name, (power, placements) in leg_factor_table(pole_table).items():
        print(f'  {name}: eps^{power}; {placements}')
    table = vertex_power_table()
    print(f'\nVertices by field-end pattern ({len(table)} patterns): generic | generic loop | up-only | up-only loop | rule')
    for pattern, row in table.items():
        print(f"  {pattern_text(pattern)}: {sorted(row['generic'])} | {sorted(row['generic, loop momenta'])} | "
              f"{sorted(row['up-only'])} | {sorted(row['up-only, loop momenta'])} | {sorted(row['rule'])}")
    print('Exceptions to eps^(r+q):', [pattern_text(p) for p in exceptions_to_rule(table)])
    print('lambda handling (V(sss), V(spp)):', lambda_check())
    powers, up_only = tree_leading()
    print(f'\nTree: powers {sorted(powers)}; up-only M_0 = {sp.latex(up_only)}')


if __name__ == '__main__':
    main()
