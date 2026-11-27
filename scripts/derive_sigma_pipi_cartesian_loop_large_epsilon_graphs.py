"""Cartesian frame: large epsilon at fixed sources - leading powers of every kept graph, of M_1^est/M_0 and of R (task 4.3, generic sources).

Run:
    PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_cartesian_loop_large_epsilon_graphs.py
    PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_cartesian_loop_large_epsilon_graphs.py --markdown
    PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_cartesian_loop_large_epsilon_graphs

OpenSpec change estimate-sigma-pipi-one-loop-murayama-cartesian.  Units GeV.
Conventions of sigma_pipi_cartesian_loop_algebra.py; u = 1/eps.

Regime (D10 restated in the Cartesian invariants; record
cartesian-sigma-pipi-loop-large-epsilon-inputs.md; assumptions): eps grows at
fixed S^C_hat, P^C_hat, lambda, f_pi, B_0, m_1^2, m_P, m_pi; A = B_0(<S^C_10> -
2<S^C_02>) > 0, kappa_sigma = B_0(<S^C_10> - 2<S^C_20>) > 0;
m_sigma^2 = 2 lambda f_pi^2/eps^2 + eps^2 kappa_sigma (exact),
m_a^2 = eps^2 A + 2 c_0/A + ... (root of D_a, series), N(m_a^2) = 1 + 2 eps^2
<(M^{-1}_{m_a^2}(C))^2>, M^{-1}_{m_a^2}(C) = -sum_k L^k(C)/m_a^{2(k+1)} with
L(X) = B_0 {S^C_00, X} + (m_1^2/3)<X> 1 and c_k = <C L^k(C)>, C = B_0 P^C_01.
eps is large enough that m_a, m_sigma > m_P, so max(m_P, m_a) = m_a,
max(m_P, m_sigma) = m_sigma.

Method.  Nothing is truncated in the estimate: the series in u are used only
to read the leading power and its coefficient.
 1. Routing at large eps (D12 restated): for every product of traces t of a
    graph the routing is chosen by the largest leading power of eps of C_t^rho,
    then the largest |leading coefficient| (exact ratio kappa_sigma = r A),
    then the next coefficients; exact ties follow D12 (sign of the sum over
    the routings, else the first routing, marked).
 2. Regions in r = kappa_sigma/A: bounded by r = 1 (where max(m_a, m_sigma)
    and the Piecewise conditions of stage B change) and by the positive roots
    of |lead_rho(r)| = |lead_rho'(r)|; each open interval is checked at three
    points, the boundary points separately.
 3. The selected X_Gamma (and dZ) are rewritten exactly: M^{-1}_{m_a^2}(C) by its
    series in the matrices L^k(C) (words of P^C_01, S^C_00 and <P^C_01>), then
    <S^C_10> = kappa_sigma/B_0 + 2<S^C_20>, <S^C_02> = (kappa_sigma - A)/(2 B_0) + <S^C_20>
    (stage-A convention), then each scalar coefficient is expanded in u and
    c_k is written back as traces.  Distinct canonical words are independent
    (generic sources, as in stage B; no Cayley-Hamilton relation is used).
 4. The leading power is read on the sum of all pieces of a quantity; the
    script checks that it lies inside the range where the truncated series
    is exact (otherwise it raises).

M_1^est/M_0 = sum_Gamma X_Gamma/M_0 + (dZ_s + dZ_pi+ + dZ_pi-)/2, with
M_0 = -(B_0 eps/(2 f_pi)) <S^C_10 {lambda_pi+, lambda_pi-}> (the tree trace).
Every quantity Q/M_0 is written as eps^p (numerator)/<S^C_10 {lambda_pi+, lambda_pi-}>;
no symbol is introduced for the tree trace or the numerator in the output.

No loop integral, scalar loop function, loop coefficient or loop amplitude is
computed.  The powers are labels of an order-of-magnitude estimate; graphs
with the same epsilon scaling can cancel in the exact loop.
"""
import functools
import itertools
import sys
from dataclasses import dataclass

import sympy as sp

import sigma_pipi_cartesian_loop_algebra as cartesian
import derive_sigma_pipi_cartesian_loop_graphs as graphs
import derive_sigma_pipi_cartesian_loop_traces as traces
import derive_sigma_pipi_cartesian_loop_large_epsilon_vertices as vertices_script
from sigma_pipi_loop_algebra import (B0, EPSILON, F_PI, M_A, M_MESON, M_PI, M_SIGMA, NORMALIZATION, QUARTIC,
                                     TraceSum)
from sigma_pipi_cartesian_loop_algebra import ALP_MASS_COEFFICIENT, KAPPA_SIGMA, MINV_LETTER, SINGLET_MASS_SQUARED

U = vertices_script.INVERSE_EPSILON            # u = 1/eps
A = ALP_MASS_COEFFICIENT
KAPPA = KAPPA_SIGMA
RATIO = sp.Symbol('r', positive=True)          # r = kappa_sigma / A
TRACE_C = sp.symbols('c_0:4', real=True)       # c_k = <C L^k(C)>
TRACE_C_SQUARED = sp.Symbol('w', positive=True)    # <P^C_01 P^C_01> > 0 (Hermitian, nonzero), used for signs
EXACT = 10 ** 6                                 # order of an exact (untruncated) series
INVERSE_DEPTH = 12                              # relative powers of u kept in the inverse of an infinite series
MINV_TERMS = 3                                  # L^0(C), L^1(C), L^2(C)
DEPTHS = (2, 4)                                 # relative depths in u tried when reading a leading term
TREE_TRACE = TraceSum.single(('SC10', 'l+', 'l-')) + TraceSum.single(('SC10', 'l-', 'l+'))
TREE_TRACE_LATEX = '\\langle\\widehat S^C_{1,0}\\{\\lambda_{\\pi^+},\\lambda_{\\pi^-}\\}\\rangle'


# ---------------------------------------------------------------------------
# 1. Truncated Laurent series in u
# ---------------------------------------------------------------------------

class Laurent:
    """sum_n c_n u^n, exact for n < order."""

    def __init__(self, terms, order=EXACT):
        self.order = order
        self.terms = {}
        for power, value in terms.items():
            if power < order:
                value = sp.expand(value)
                if value != 0:
                    self.terms[power] = value

    @staticmethod
    def constant(value):
        return Laurent({0: value})

    def lowest(self):
        return min(self.terms) if self.terms else self.order

    def __add__(self, other):
        terms = dict(self.terms)
        for power, value in other.terms.items():
            terms[power] = terms.get(power, 0) + value
        return Laurent(terms, min(self.order, other.order))

    def __mul__(self, other):
        if not self.terms or not other.terms:
            return Laurent({}, min(self.order + other.lowest(), other.order + self.lowest()))
        order = min(self.order + other.lowest(), other.order + self.lowest())
        terms = {}
        for (p1, v1), (p2, v2) in itertools.product(self.terms.items(), other.terms.items()):
            if p1 + p2 < order:
                terms[p1 + p2] = terms.get(p1 + p2, 0) + v1 * v2
        return Laurent(terms, order)

    def scale(self, factor):
        return Laurent({p: factor * v for p, v in self.terms.items()}, self.order)

    def inverse(self):
        """1/series, exact for relative powers below min(order - lowest, INVERSE_DEPTH): the inverse of
        an exact series with two or more terms is infinite and is cut at INVERSE_DEPTH (the cut is recorded
        in the order, so a leading term read beyond it raises in leading_of)."""
        low = self.lowest()
        lead = self.terms[low]
        depth = min(self.order - low, INVERSE_DEPTH)
        rest = Laurent({p - low: v / lead for p, v in self.terms.items() if p != low}, depth)
        result = Laurent.constant(1)
        power = Laurent.constant(1)
        step = min(rest.terms) if rest.terms else depth
        for _ in range(depth // max(step, 1) + 1):
            power = power * rest.scale(-1)
            result = result + power
        return Laurent({p - low: v / lead for p, v in result.terms.items()}, depth - low)

    def power(self, exponent):
        if exponent == 0:
            return Laurent.constant(1)
        base = self if exponent > 0 else self.inverse()
        result = Laurent.constant(1)
        for _ in range(abs(exponent)):
            result = result * base
        return result

    def substitute(self, substitution):
        return Laurent({p: sp.sympify(v).subs(substitution) for p, v in self.terms.items()}, self.order)


# ---------------------------------------------------------------------------
# 2. Building blocks
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Region:
    """A region of r = kappa_sigma/A: kappa is KAPPA (open region) or A (boundary r = 1);
    heavier is M_A or M_SIGMA (the larger of m_a and m_sigma at large eps)."""
    name: str
    kappa: object
    heavier: object


REGION_BELOW = Region('r < 1', KAPPA, M_A)
REGION_ABOVE = Region('r > 1', KAPPA, M_SIGMA)
REGION_BOUNDARY = Region('r = 1', A, M_A)


@functools.lru_cache(maxsize=None)
def alp_mass_squared_series():
    """m_a^2 = A/u^2 + b_0 + b_1 u^2 + b_2 u^4 (root of D_a, exact through u^4), in c_0, c_1, c_2."""
    root = sp.expand(vertices_script.alp_pole_root().subs(dict(zip(vertices_script.TRACE_C_POWER, TRACE_C))))
    terms = {}
    for term in sp.Add.make_args(root):
        coefficient, power = term.as_coeff_exponent(U)
        terms[int(power)] = terms.get(int(power), 0) + coefficient
    return Laurent(terms, 6)


def radial_mass_squared_series(region):
    return Laurent({2: 2 * QUARTIC * F_PI ** 2, -2: region.kappa})


@functools.lru_cache(maxsize=None)
def normalization_series():
    """N(m_a^2) = 1 + 2 eps^2 <(M^{-1}_{m_a^2}(C))^2> = 1 + (2/u^2) sum_n (n+1) c_n / m_a^{2(n+2)}."""
    inverse_mass = alp_mass_squared_series().inverse()
    total = Laurent.constant(1)
    for n in range(4):
        total = total + inverse_mass.power(n + 2).scale(2 * (n + 1) * TRACE_C[n]) * Laurent({-2: 1})
    return total


def maximum_value(node, region):
    """max(...) at large eps: m_P is fixed, so it is never the largest; max(m_a, m_sigma) -> region.heavier."""
    arguments = set(node.args) - {M_MESON}
    if arguments == {M_A}:
        return M_A
    if arguments == {M_SIGMA}:
        return M_SIGMA
    if arguments == {M_A, M_SIGMA}:
        return region.heavier
    raise ValueError(node)


def block_series(base, exponent, region):
    """Series of base**exponent for the eps-dependent symbols."""
    if base == EPSILON:
        return Laurent({-exponent: 1})
    if base in (M_A, M_SIGMA):
        assert exponent % 2 == 0, (base, exponent)
        squared = alp_mass_squared_series() if base == M_A else radial_mass_squared_series(region)
        return squared.power(exponent // 2)
    if base == NORMALIZATION:
        return normalization_series().power(exponent)
    return None


@functools.lru_cache(maxsize=None)
def _monomial_series(factors, region):
    result = Laurent.constant(1)
    for base, exponent in factors:
        result = result * block_series(base, exponent, region)
    return result


def scalar_series(expression, region):
    """Laurent series in u of a scalar coefficient (masses, N, eps, max)."""
    expression = sp.sympify(expression).replace(lambda node: isinstance(node, sp.Max),
                                                lambda node: maximum_value(node, region))
    if region.kappa != KAPPA:
        expression = expression.subs(KAPPA, region.kappa)
    total = Laurent({}, EXACT)
    for term in sp.Add.make_args(sp.expand(expression)):
        constant = sp.Integer(1)
        factors = []
        for factor in sp.Mul.make_args(term):
            base, exponent = factor.as_base_exp()
            if base in (EPSILON, M_A, M_SIGMA, NORMALIZATION):
                factors.append((base, int(exponent)))
            else:
                constant *= factor
        total = total + _monomial_series(tuple(sorted(factors, key=str)), region).scale(constant)
    return total


# Matrices L^k(C) as lists of (coefficient, word, trace words).
def apply_loop_operator(matrix):
    """L(X) = B_0 {S^C_00, X} + (m_1^2/3) <X> 1."""
    result = []
    for coefficient, word, trace_words in matrix:
        result.append((B0 * coefficient, ('SC00',) + word, trace_words))
        result.append((B0 * coefficient, word + ('SC00',), trace_words))
        result.append((SINGLET_MASS_SQUARED / 3 * coefficient, (), trace_words + (word,)))
    return result


@functools.lru_cache(maxsize=None)
def loop_operator_powers():
    """[L^0(C), L^1(C), ...] with C = B_0 P^C_01."""
    powers = [[(B0, ('PC01',), ())]]
    for _ in range(MINV_TERMS):
        powers.append(apply_loop_operator(powers[-1]))
    return powers


@functools.lru_cache(maxsize=None)
def trace_c_sums():
    """c_k = <C L^k(C)> as TraceSums, k = 0..3."""
    powers = loop_operator_powers()
    result = []
    for k in range(4):
        total = TraceSum()
        for coefficient, word, trace_words in powers[k]:
            total.add_term(trace_words + (('PC01',) + word,), B0 * coefficient)
        result.append(total.simplified())
    return result


def expand_inverse(trace_sum, depth):
    """M^{-1}_{m_a^2}(C) = -sum_k L^k(C) m_a^{-2(k+1)} in every word; combinations with sum_i 2 k_i > depth
    are dropped (they start at relative order u^{depth+2}).  Returns (TraceSum, {key: number of letters})."""
    powers = loop_operator_powers()
    result = TraceSum()
    for key, value in trace_sum.terms.items():
        positions = [(i, j) for i, word in enumerate(key) for j, name in enumerate(word) if name == MINV_LETTER]
        for orders in itertools.product(range(MINV_TERMS), repeat=len(positions)):
            if 2 * sum(orders) > depth:
                continue
            for choice in itertools.product(*(powers[k] for k in orders)):
                coefficient = value
                extra_traces = ()
                substitution = {}
                for position, k, (factor, word, trace_words) in zip(positions, orders, choice):
                    coefficient *= -factor * M_A ** (-2 * (k + 1))
                    extra_traces += trace_words
                    substitution[position] = word
                words = []
                for i, word in enumerate(key):
                    new_word = ()
                    for j, name in enumerate(word):
                        new_word += substitution.get((i, j), (name,))
                    words.append(new_word)
                result.add_term(tuple(words) + extra_traces, coefficient)
    return result


def rewrite_single_traces(trace_sum):
    """<S^C_10> = kappa_sigma/B_0 + 2<S^C_20>, <S^C_02> = (kappa_sigma - A)/(2 B_0) + <S^C_20> (exact)."""
    replacements = {('SC10',): TraceSum.scalar(KAPPA / B0) + TraceSum.single(('SC20',), 2),
                    ('SC02',): TraceSum.scalar((KAPPA - A) / (2 * B0)) + TraceSum.single(('SC20',))}
    result = TraceSum()
    for key, value in trace_sum.terms.items():
        product = TraceSum.scalar(value)
        for word in key:
            product = product * replacements.get(word, TraceSum.single(word))
        result = result + product
    return result


def traces_from_c(expression, key):
    """Write the polynomial in c_k of a series coefficient back as traces, multiplied by the product `key`."""
    result = TraceSum()
    polynomial = sp.Poly(sp.expand(expression), *TRACE_C)
    c_sums = trace_c_sums()
    for monomial, coefficient in polynomial.terms():
        product = TraceSum({key: coefficient})
        for k, power in enumerate(monomial):
            for _ in range(power):
                product = product * c_sums[k]
        result = result + product
    return result


def term_series(trace_sum, region, depth):
    """([(key, Laurent series of its coefficient)], exactness bound) after the exact rewriting
    (M^{-1}_{m_a^2}(C) series, <S^C_10>, <S^C_02> in A and kappa_sigma).  Powers of u below the bound are exact:
    the bound is the smallest series order and, for every original product with n_M letters
    M^{-1}_{m_a^2}(C), its leading power + 2 n_M + depth + 2 (first dropped combination)."""
    expanded = expand_inverse(trace_sum, depth)
    expanded = rewrite_single_traces(expanded) if region.kappa == KAPPA else _rewrite_boundary(expanded)
    bound = EXACT
    for key, value in trace_sum.terms.items():
        count = sum(word.count(MINV_LETTER) for word in key)
        if count:
            bound = min(bound, scalar_series(value, region).lowest() + 2 * count + depth + 2)
    items = []
    for key, value in expanded.simplified().terms.items():
        series = scalar_series(value, region)
        bound = min(bound, series.order)
        if series.terms:
            items.append((key, series))
    return items, bound


def assemble(items, power):
    """TraceSum of the coefficient of u^power (c_k written back as traces)."""
    total = TraceSum()
    for key, series in items:
        if power in series.terms:
            total = total + traces_from_c(series.terms[power], key)
    return total.simplified()


def lowest_nonzero(items, bound):
    """(power of u, TraceSum) of the lowest nonzero power below the bound, or None."""
    for power in sorted({p for _, series in items for p in series.terms if p < bound}):
        value = assemble(items, power)
        if value.terms:
            return power, value
    return None


def _rewrite_boundary(trace_sum):
    """At r = 1 (kappa_sigma = A) the same rewriting with kappa_sigma replaced by A."""
    return rewrite_single_traces(trace_sum).map_coefficients(lambda value: sp.sympify(value).subs(KAPPA, A))


def leading_of(trace_sum, region):
    """(power of eps, leading TraceSum, relative depth used).  Tries the relative depths of DEPTHS and raises if
    the leading term is not inside the exact range."""
    for depth in DEPTHS:
        items, bound = term_series(trace_sum, region, depth)
        found = lowest_nonzero(items, bound)
        if found is not None:
            return -found[0], found[1], depth
    raise ValueError('leading term beyond the exact range of the series')


# ---------------------------------------------------------------------------
# 3. Routing at large eps and regions of r
# ---------------------------------------------------------------------------

def sign_of(expression):
    """Sign of an expression in positive symbols (A, B_0, f_pi, lambda, m_pi, m_P, m_1^2, w, m_u); r given."""
    expression = sp.factor(sp.expand(expression))
    if expression == 0:
        return 0
    if expression.is_positive:
        return 1
    if expression.is_negative:
        return -1
    raise ValueError(f'sign not determined: {expression}')


SIGN_SUBSTITUTION = {TRACE_C[0]: B0 ** 2 * TRACE_C_SQUARED}


def magnitude_key(series, ratio_value):
    """The series of |C| at kappa_sigma = ratio_value A, as a list of (power of eps, coefficient) from the
    leading term on, for the lexicographic comparison."""
    series = series.substitute({KAPPA: ratio_value * A, **SIGN_SUBSTITUTION})
    if not series.terms:
        return None
    sign = sign_of(series.terms[series.lowest()])
    return [(-power, sign * series.terms[power]) for power in sorted(series.terms)], series.order


def compare_magnitudes(first, second):
    """+1 if |first| > |second| at large eps, -1 if smaller, 0 if equal through the exact range."""
    if first is None or second is None:
        return (first is not None) - (second is not None)
    (terms_one, order_one), (terms_two, order_two) = first, second
    limit = -min(order_one, order_two)
    one, two = dict(terms_one), dict(terms_two)
    powers = sorted(set(one) | set(two), reverse=True)
    for power in powers:
        if power <= limit:
            break
        difference = sp.expand(one.get(power, 0) - two.get(power, 0))
        if difference != 0:
            return sign_of(difference)
    return 0


@functools.lru_cache(maxsize=None)
def routing_dependent_products():
    """{label: {key: [C_t^rho for every routing]}} for the products whose coefficient depends on the routing."""
    result = {}
    for entry in traces.entries():
        keys = set()
        for value in entry['words']:
            keys |= set(value.terms)
        dependent = {}
        for key in keys:
            coefficients = [sp.expand(value.terms.get(key, 0)) for value in entry['words']]
            if any(sp.expand(c - coefficients[0]) != 0 for c in coefficients):
                dependent[key] = coefficients
        if dependent:
            result[entry['label']] = dependent
    return result


@functools.lru_cache(maxsize=None)
def _coefficient_series(coefficients, region):
    return tuple(scalar_series(value, region) for value in coefficients)


def coefficient_series(coefficients, region):
    return list(_coefficient_series(tuple(coefficients), region))


def chosen_routing(coefficients, region, ratio_value):
    """Large-eps routing choice for one product at r = ratio_value: (index, marked)."""
    series = coefficient_series(coefficients, region)
    keys = [magnitude_key(value, ratio_value) for value in series]
    best = 0
    for index in range(1, len(keys)):
        if compare_magnitudes(keys[index], keys[best]) > 0:
            best = index
    tied = [index for index in range(len(keys)) if compare_magnitudes(keys[index], keys[best]) == 0]
    if len(tied) == 1:
        return best, False
    substitution = {KAPPA: ratio_value * A}
    exact = [sp.expand(coefficients[index].subs(substitution)) for index in tied]
    if not all(sp.expand(value - exact[0]) == 0 or sp.expand(value + exact[0]) == 0 for value in exact):
        raise ValueError('tie through the exact range without an exact tie')
    if all(sp.expand(value - exact[0]) == 0 for value in exact):
        return tied[0], False
    # D12: equal magnitude, opposite signs -> sign of the sum over all routings; zero sum -> first, marked.
    total = sum(series[1:], series[0]).substitute({**substitution, **SIGN_SUBSTITUTION})
    if not total.terms:
        return tied[0], True
    total_sign = sign_of(total.terms[total.lowest()])
    for index in tied:
        value = series[index].substitute({**substitution, **SIGN_SUBSTITUTION})
        if sign_of(value.terms[value.lowest()]) == total_sign:
            return index, False
    return tied[0], True


def region_for(ratio_value):
    if ratio_value < 1:
        return REGION_BELOW
    if ratio_value > 1:
        return REGION_ABOVE
    return REGION_BOUNDARY


def choices_at(ratio_value):
    """{(label, key): (routing index, marked)} for every routing-dependent product at r = ratio_value."""
    region = region_for(ratio_value)
    result = {}
    for label, products in routing_dependent_products().items():
        for key, coefficients in products.items():
            result[(label, key)] = chosen_routing(tuple(coefficients), region, ratio_value)
    return result


def comparison_roots():
    """Positive roots in r of |lead_rho(r)| = |lead_rho'(r)| for every pair of routings with the same
    leading power (open regions, kappa_sigma = r A), and of lead_rho(r) = 0."""
    roots = set()
    for region in (REGION_BELOW, REGION_ABOVE):
        for products in routing_dependent_products().values():
            for coefficients in products.values():
                leads = []
                for value in coefficient_series(coefficients, region):
                    low = value.lowest()
                    leads.append((low, sp.expand(value.terms[low].subs(KAPPA, RATIO * A)) if value.terms else 0))
                for (p1, l1), (p2, l2) in itertools.combinations(leads, 2):
                    if p1 != p2:
                        continue
                    for equation in (l1 - l2, l1 + l2):
                        if sp.expand(equation) == 0:
                            continue
                        roots |= {root for root in sp.solve(equation, RATIO) if root.is_positive}
    return sorted(roots)


def sample_points(lower, upper):
    if upper == sp.oo:
        return [lower + sp.Rational(1, 2), lower + 1, 4 * lower]
    return [lower + (upper - lower) * sp.Rational(k, 4) for k in (1, 2, 3)]


@functools.lru_cache(maxsize=None)
def regions():
    """(merged regions, choices at the boundary points).  The open intervals bounded by r = 1 and the
    comparison roots are checked at three points each; adjacent intervals with the same choices (and the
    same larger of m_a, m_sigma) are merged when the choice at the root between them is also the same.
    A merged region is {'lower', 'upper', 'points', 'choices', 'constant', 'interior roots'}."""
    boundaries = sorted(set(comparison_roots()) | {sp.Integer(1)})
    edges = [sp.Integer(0)] + boundaries + [sp.oo]
    intervals = []
    for lower, upper in zip(edges[:-1], edges[1:]):
        points = sample_points(lower, upper)
        choices = [choices_at(point) for point in points]
        constant = all(choice == choices[0] for choice in choices) and len({region_for(p) for p in points}) == 1
        intervals.append({'lower': lower, 'upper': upper, 'points': points, 'choices': choices[0],
                          'constant': constant, 'interior roots': []})
    boundary_choices = {point: choices_at(point) for point in boundaries}
    merged = [intervals[0]]
    for interval in intervals[1:]:
        previous = merged[-1]
        root = interval['lower']
        if (interval['choices'] == previous['choices'] and boundary_choices[root] == previous['choices']
                and region_for(interval['points'][0]) == region_for(previous['points'][0])):
            previous.update(upper=interval['upper'], points=previous['points'] + interval['points'],
                            constant=previous['constant'] and interval['constant'],
                            **{'interior roots': previous['interior roots'] + [root]})
        else:
            merged.append(dict(interval))
    return merged, boundary_choices


def stage_b_branch_agrees(ratio_value):
    """The large-eps choice equals the stage-B termwise maximum (Piecewise in the exact masses) with its
    conditions evaluated at large eps: m_a^2 >= m_sigma^2 - m_pi^2 holds for r <= 1 and fails for r > 1."""
    region = region_for(ratio_value)
    holds = ratio_value <= 1
    for (label, key), (index, _) in choices_at(ratio_value).items():
        entry = traces.entry_by_label(label)
        maximum = sp.sympify(entry['maximum']['value'].terms[key])
        chosen = maximum.replace(lambda node: isinstance(node, sp.Piecewise),
                                 lambda node: node.args[0].expr if holds else node.args[1].expr)
        selected = entry['words'][index].terms.get(key, 0)
        if sp.expand(chosen - selected) != 0:
            return False
    return True


def large_epsilon_conditions(region):
    """Leading terms of m_a^2 - m_sigma^2 and of m_a^2 - (m_sigma^2 - m_pi^2) (the stage-B condition)."""
    difference = alp_mass_squared_series() + radial_mass_squared_series(region).scale(-1)
    shifted = difference + Laurent.constant(M_PI ** 2)
    result = []
    for series in (difference, shifted):
        low = series.lowest()
        result.append((-low, sp.factor(series.terms[low])))
    return result


# ---------------------------------------------------------------------------
# 4. Selected values and leading terms
# ---------------------------------------------------------------------------

def selected_value(entry, choices):
    """X_Gamma (or dZ) with the large-eps routing choice of every product."""
    result = TraceSum()
    keys = set()
    for value in entry['words']:
        keys |= set(value.terms)
    for key in keys:
        index = choices.get((entry['label'], key), (0, False))[0]
        result.add_term(key, entry['words'][index].terms.get(key, 0))
    return result.simplified()


def contribution(group, power, value):
    """Share of M_1^est/M_0 = (numerator)/<S^C_10 {lambda_pi+, lambda_pi-}> of a leading term eps^power value:
    1PI: X/M_0 = X (-2 f_pi/(B_0 eps))/T_0; s leg: dZ/2 = (dZ T_0/2)/T_0; pi legs: (dZ_pi+ + dZ_pi-)/2 = dZ."""
    if group == 'one_particle_irreducible':
        return power - 1, value.scale(-2 * F_PI / B0).simplified()
    if group == 'sigma':
        return power, (value * TREE_TRACE).scale(sp.Rational(1, 2)).simplified()
    return power, (value * TREE_TRACE).simplified()


def region_data(region_name):
    """(Region, routing choices) for 'r < 1', 'r > 1' or 'r = 1'."""
    merged, boundary_choices = regions()
    if region_name == 'r = 1':
        return REGION_BOUNDARY, boundary_choices[sp.Integer(1)]
    for interval in merged:
        if region_for(interval['points'][0]).name == region_name:
            return region_for(interval['points'][0]), interval['choices']
    raise KeyError(region_name)


def region_names():
    return [region_for(interval['points'][0]).name for interval in regions()[0]]


@functools.lru_cache(maxsize=None)
def graph_results(region_name):
    """{label: {'own': (power, TraceSum), 'contribution': (power, numerator TraceSum), 'group'}}.
    'own' is the leading term of X_Gamma^max (1PI) or dZ^max; 'contribution' the numerator of its share of
    M_1^est/M_0 over the tree trace <S^C_10 {lambda_pi+, lambda_pi-}>."""
    region, choices = region_data(region_name)
    result = {}
    for entry in traces.entries():
        power, value, _ = leading_of(selected_value(entry, choices), region)
        result[entry['label']] = {'own': (power, value), 'contribution': contribution(entry['group'], power, value),
                                  'group': entry['group']}
    return result


@functools.lru_cache(maxsize=None)
def total_leading(region_name):
    """(power of eps, numerator TraceSum, labels of the graphs whose share has that power) of the leading term of
    M_1^est/M_0 = numerator/<S^C_10 {lambda_pi+, lambda_pi-}>, read on the sum of all graphs."""
    region, choices = region_data(region_name)
    for depth in DEPTHS:
        items, bound = [], EXACT
        for entry in traces.entries():
            entry_items, entry_bound = term_series(selected_value(entry, choices), region, depth)
            bound = min(bound, entry_bound + (1 if entry['group'] == 'one_particle_irreducible' else 0))
            for key, series in entry_items:
                if entry['group'] == 'one_particle_irreducible':
                    items.append((key, series * Laurent({1: -2 * F_PI / B0})))
                else:
                    factor = sp.Rational(1, 2) if entry['group'] == 'sigma' else 1
                    for tree_key, tree_value in TREE_TRACE.terms.items():
                        items.append((key + tree_key, series.scale(factor * tree_value)))
        found = lowest_nonzero(items, bound)
        if found is not None:
            power = -found[0]
            labels = sorted(label for label, data in graph_results(region_name).items()
                            if data['contribution'][0] == power)
            return power, found[1], labels
    raise ValueError('leading term of M_1^est/M_0 beyond the exact range')


def graphs_sum_check(region_name):
    """The leading numerator of M_1^est/M_0 equals the sum of the leading numerators of the graphs at that power."""
    power, numerator, labels = total_leading(region_name)
    total = TraceSum()
    for label in labels:
        total = total + graph_results(region_name)[label]['contribution'][1]
    return (total - numerator).is_zero()


# ---------------------------------------------------------------------------
# 5. D4 terms and groups of graphs with the same leading power
# ---------------------------------------------------------------------------

def d4_powers():
    """{label: largest power of eps over the routings and regions of (D4 term)/M_0 (for the pi leg: of its share
    dZ of M_1/M_0)}.  The D4 terms depend on eps only through m_sigma^2."""
    result = {}
    for group in graphs.GROUPS:
        for row in graphs.table(group):
            if row['D4 status'] != 'excluded (D4)':
                continue
            powers = set()
            for region in (REGION_BELOW, REGION_ABOVE, REGION_BOUNDARY):
                for value in row['D4 term']:
                    power, leading, _ = leading_of(value, region)
                    powers.add(contribution(group, power, leading)[0])
            result[row['label']] = max(powers)
    return result


def coefficient_signs(trace_sum, ratio_values):
    """{key: set of signs of the coefficient at the sample points r}."""
    result = {}
    for key, value in trace_sum.terms.items():
        result[key] = {sign_of(sp.sympify(value).subs({KAPPA: point * A, **SIGN_SUBSTITUTION}))
                       for point in ratio_values}
    return result


def proportionality(first, second):
    """c if first = c * second with c a scalar, else None."""
    if set(first.terms) != set(second.terms):
        return None
    ratios = {sp.factor(sp.sympify(first.terms[key]) / sp.sympify(second.terms[key])) for key in first.terms}
    return ratios.pop() if len(ratios) == 1 else None


def cancellation_status(numerators, ratio_values):
    """Statement for a group of graphs with the same leading power:
    'single graph'; 'same sign' / 'could cancel' if all numerators are proportional to one trace polynomial
    with ratios of fixed sign (all positive / some negative); otherwise 'depends on sources'.
    Also the products of traces shared by two graphs with coefficients of opposite sign."""
    labels = list(numerators)
    opposite = []
    for one, two in itertools.combinations(labels, 2):
        signs_one = coefficient_signs(numerators[one], ratio_values)
        signs_two = coefficient_signs(numerators[two], ratio_values)
        for key in set(signs_one) & set(signs_two):
            if any(a * b < 0 for a in signs_one[key] for b in signs_two[key]):
                opposite.append((one, two, key))
    if len(labels) == 1:
        return 'single graph', opposite
    reference = numerators[labels[0]]
    ratios = [proportionality(numerators[label], reference) for label in labels[1:]]
    if all(ratio is not None for ratio in ratios):
        signs = set()
        for ratio in ratios:
            signs |= {sign_of(sp.sympify(ratio).subs({KAPPA: point * A, **SIGN_SUBSTITUTION}))
                      for point in ratio_values}
        if signs == {1}:
            return 'same sign', opposite
        if -1 in signs and len(signs) == 1 or signs == {1, -1}:
            return 'could cancel', opposite
    return 'depends on sources', opposite


def region_sample_points(region_name):
    intervals, _ = regions()
    if region_name == 'r = 1':
        return [sp.Integer(1)]
    for interval in intervals:
        if region_for(interval['points'][0]).name == region_name:
            return interval['points']
    raise KeyError(region_name)


@functools.lru_cache(maxsize=None)
def same_power_groups(region_name):
    """{power: (labels, status, opposite-sign shared products)} for the contributions to M_1^est/M_0."""
    results = graph_results(region_name)
    powers = sorted({data['contribution'][0] for data in results.values()}, reverse=True)
    groups = {}
    for power in powers:
        numerators = {label: data['contribution'][1] for label, data in results.items()
                      if data['contribution'][0] == power}
        status, opposite = cancellation_status(numerators, region_sample_points(region_name))
        groups[power] = (sorted(numerators), status, opposite)
    return groups


def exactly_opposite_pairs(region_name):
    """Pairs of graphs with the same leading power whose leading numerators add up to zero exactly."""
    results = graph_results(region_name)
    pairs = []
    for (one, first), (two, second) in itertools.combinations(results.items(), 2):
        if first['contribution'][0] != second['contribution'][0]:
            continue
        if (first['contribution'][1] + second['contribution'][1]).simplified().is_zero():
            pairs.append((one, two))
    return pairs


def common_roots_in_ratio(trace_sum):
    """Positive r at which every coefficient of the TraceSum vanishes (kappa_sigma = r A)."""
    polynomials = [sp.numer(sp.factor(sp.sympify(value).subs(KAPPA, RATIO * A))) for value in trace_sum.terms.values()]
    common = polynomials[0]
    for polynomial in polynomials[1:]:
        common = sp.gcd(common, polynomial)
    return sorted(root for root in sp.solve(common, RATIO) if root.is_positive) if common.has(RATIO) else []


# ---------------------------------------------------------------------------
# 6. Output
# ---------------------------------------------------------------------------

def power_text(power):
    return f'\\epsilon^{{{power}}}'


def numerator_latex(trace_sum):
    terms = TraceSum()
    terms.terms = {key: sp.factor(value) for key, value in trace_sum.terms.items()}
    return cartesian.trace_sum_latex(terms)


def label_text(label):
    return '$' + traces.label_latex(label) + '$'


def markdown_graph_table():
    """Leading power of each contribution to M_1^est/M_0 per region, and of X^max (1PI) or dZ^max."""
    names = region_names() + ['r = 1']
    lines = ['| graph | ' + ' | '.join(f'${name}$: own / share of $M_1^{{\\rm est}}/M_0$' for name in names) + ' |',
             '|---|' + '---|' * len(names)]
    for entry in traces.entries():
        cells = []
        for name in names:
            data = graph_results(name)[entry['label']]
            cells.append(f"${power_text(data['own'][0])}$ / ${power_text(data['contribution'][0])}$")
        lines.append(f"| {label_text(entry['label'])} | " + ' | '.join(cells) + ' |')
    return '\n'.join(lines)


def graph_coefficient_text(label, region_name):
    """Leading term of X^max/M_0 (1PI graph, numerator over the tree trace) or of dZ^max (self-energy)."""
    data = graph_results(region_name)[label]
    if data['group'] == 'one_particle_irreducible':
        power, numerator = data['contribution']
        return f"$X^{{\\max}}/M_0\\to{power_text(power)}\\,\\big[" + numerator_latex(numerator) + '\\big]/' + TREE_TRACE_LATEX + '$'
    power, value = data['own']
    return f"$\\delta Z^{{\\max}}\\to{power_text(power)}\\,\\big[" + numerator_latex(value) + '\\big]$'


def coefficient_differences(label):
    """Regions (other than r < 1) in which the leading term of the graph differs from the one of r < 1."""
    reference = graph_results('r < 1')[label]
    names = [name for name in region_names() + ['r = 1'] if name != 'r < 1']
    different = []
    for name in names:
        data = graph_results(name)[label]
        if (data['own'][0] != reference['own'][0]
                or not (data['own'][1] - reference['own'][1]).simplified().is_zero()):
            different.append(name)
    return different


def markdown_graph_coefficients():
    """Leading coefficient of every kept graph in r < 1, and in the other regions where it differs."""
    blocks = []
    for entry in traces.entries():
        label = entry['label']
        text = f"- {label_text(label)} ($r<1$): " + graph_coefficient_text(label, 'r < 1')
        for name in coefficient_differences(label):
            text += f"\n  - ${name}$: " + graph_coefficient_text(label, name)
        blocks.append(text)
    return '\n'.join(blocks)


def markdown_d4():
    """Largest power of each excluded D4 term relative to M_0 and its distance to the leading power of R."""
    leading_powers = {total_leading(name)[0] for name in region_names() + ['r = 1']}
    assert len(leading_powers) == 1
    power_of_ratio = leading_powers.pop()
    lines = ['| D4 term | largest power relative to $M_0$ | powers below the leading power of $R$ |', '|---|---|---|']
    for label, power in d4_powers().items():
        lines.append(f'| {label_text(label)} | ${power_text(power)}$ | {power_of_ratio - power} |')
    return '\n'.join(lines)


def markdown_leading(region_name):
    """Leading term of M_1^est/M_0 in a region, written out explicitly (R is its absolute value)."""
    power, numerator, labels = total_leading(region_name)
    groups = same_power_groups(region_name)
    text = (f"**${region_name}$.** Graphs at the leading power: "
            + ', '.join(label_text(label) for label in labels) + '. Status of this group: '
            + groups[power][1] + '. $R=|M_1^{\\rm est}|/|M_0|$ tends to the absolute value of\n$$\n'
            + '\\frac{M_1^{\\rm est}}{M_0}\\to\\frac{' + power_text(power) + '}{' + TREE_TRACE_LATEX + '}\\Big['
            + numerator_latex(numerator) + '\\Big]\n$$')
    return text


def markdown_contributions(region_name):
    """Leading numerators of the graphs at the leading power of M_1^est/M_0."""
    power, _, labels = total_leading(region_name)
    blocks = []
    for label in labels:
        numerator = graph_results(region_name)[label]['contribution'][1]
        blocks.append(f"{label_text(label)}: $" + numerator_latex(numerator) + '$')
    return '\n\n'.join(blocks)


def markdown_groups(region_name):
    lines = ['| power | graphs | status | shared products with opposite signs | leading terms exactly opposite |',
             '|---|---|---|---|---|']
    exact = exactly_opposite_pairs(region_name)
    for power, (labels, status, opposite) in same_power_groups(region_name).items():
        pairs = sorted({(one, two) for one, two, _ in opposite})
        exact_here = [(one, two) for one, two in exact if one in labels]
        lines.append(f"| ${power_text(power)}$ | " + ', '.join(label_text(label) for label in labels)
                     + f' | {status} | ' + (', '.join(f'{label_text(a)}–{label_text(b)}' for a, b in pairs) or '–')
                     + ' | ' + (', '.join(f'{label_text(a)}–{label_text(b)}' for a, b in exact_here) or '–') + ' |')
    return '\n'.join(lines)


def main():
    if '--markdown' in sys.argv:
        print(markdown_graph_table() + '\n')
        print(markdown_graph_coefficients() + '\n')
        for name in region_names() + ['r = 1']:
            print(markdown_leading(name) + '\n')
            print(markdown_contributions(name) + '\n')
            print(markdown_groups(name) + '\n')
        print(markdown_d4())
        return
    print('Large-eps conditions (power, coefficient) of m_a^2 - m_sigma^2 and m_a^2 - m_sigma^2 + m_pi^2:')
    for region in (REGION_BELOW, REGION_BOUNDARY):
        print('   ', region.name, large_epsilon_conditions(region))
    print('Routing-dependent products:', {label: len(products) for label, products in routing_dependent_products().items()})
    print('Comparison roots in r:', comparison_roots())
    intervals, boundary_choices = regions()
    for interval in intervals:
        counts = {}
        for (label, _), (index, marked) in interval['choices'].items():
            counts.setdefault(label, set()).add((index + 1, marked))
        print(f"Region ({interval['lower']}, {interval['upper']}): constant at {interval['points']}: "
              f"{interval['constant']}; routings {counts}; stage-B branch agrees: "
              f"{all(stage_b_branch_agrees(point) for point in interval['points'])}")
    for point, choices in boundary_choices.items():
        counts = {}
        for (label, _), (index, marked) in choices.items():
            counts.setdefault(label, set()).add((index + 1, marked))
        print(f'Boundary r = {point}: routings {counts}; stage-B branch agrees: {stage_b_branch_agrees(point)}')
    for name in region_names() + ['r = 1']:
        print(f'\n== {name}')
        for label, data in graph_results(name).items():
            print(f"  {label}: own eps^{data['own'][0]}, share eps^{data['contribution'][0]}")
        power, numerator, labels = total_leading(name)
        print(f'  M_1/M_0 leading: eps^{power}; graphs {labels}; sum check {graphs_sum_check(name)}')
        print('  numerator:', numerator_latex(numerator))
        print('  common roots in r of the leading coefficient:', common_roots_in_ratio(numerator))
        for group_power, (labels, status, opposite) in same_power_groups(name).items():
            print(f'  group eps^{group_power}: {labels}: {status}; opposite shared: '
                  f'{sorted({(a, b) for a, b, _ in opposite})}')
    print('\nD4 terms, largest power relative to M_0:', d4_powers())


if __name__ == '__main__':
    main()
