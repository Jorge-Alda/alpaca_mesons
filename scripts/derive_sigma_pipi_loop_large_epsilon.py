"""Behaviour at large epsilon = f_pi/f_a of the ALP-like pole data and of every quantity of the one-loop sigma -> pi+ pi- estimate (task 2.3).

Run:
    PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_loop_large_epsilon.py
    PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_loop_large_epsilon.py --markdown
    PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_loop_large_epsilon
The first line prints sections 1-4 (a few seconds); the second prints the
markdown blocks of sections 1-10, labelled by the memory record and section
they belong to (about two minutes).

Units GeV (s, m_a^2, m_1^2 in GeV^2; B_0 and the sources in GeV); conventions
of sigma_pipi_loop_algebra.py (polar fields, connection-free frame,
<lambda_i lambda_j> = 2 delta_ij, C = B_0 P_hat_01, Theta_hat_1 = 0).

Contents (all exact symbolic, operator level: matrices are noncommuting
letters inside traces; no components):
 1. Quadratic part of the Lagrangian from the vertex rules (two legs) plus the
    two-point kinetic and anomalous terms read off the Lagrangian:
      D_PP(X, Y) = -(1/2) <X M_s(Y)>,  D_aP(X) = eps <C X>,
      D_aa = s + 2 B_0 eps^2 <S_02>    (ALP mass entry -2 B_0 eps^2 <S_02>),
    with M_s(X) = B_0 {S_00, X} + (m_1^2/3) <X> 1 - s X.
 2. Elimination of the meson field at stationarity, Phi = 2 eps a M_s^{-1}(C),
    and the ALP inverse propagator after the mesons are integrated out
      D_a(s) = s + 2 B_0 eps^2 <S_02> + 2 eps^2 <C M_s^{-1}(C)>,
    with N(s) = D_a'(s) = 1 + 2 eps^2 <(M_s^{-1}(C))^2> equal to the norm of
    the pole vector (2 eps M_s^{-1}(C), 1).
 3. Exact identities of the meson-line residue (1_9 (+) 1) - v_a v_a^T, valid
    at any eps: v_a^T v_a = 1, so the residue is a projector; the meson-block
    pole is idempotent and its trace over the nine mesons is one.
 4. Large-eps behaviour (fixed sources and fixed lambda) of the pole data:
    ALP mass entry eps^2 A with A = -2 B_0 <S_02> > 0 fixed.  Large-s
    expansion of M_s^{-1}, root m_a^2 = eps^2 A + O(1), N -> 1; check with the
    up-only sources.
 5. Every quantity of the estimate at large eps.  Assumptions:
    * held fixed (independent of eps): f_pi, B_0, m_1^2, every source letter
      S_hat_{r,s}, P_hat_{r,s} (generic: all independent products of traces are
      nonzero), lambda_{pi+-}, and the quartic coupling lambda of V_phi;
      f_a = f_pi/eps varies.  m_pi is not an independent input (it comes
      from S_00 only) and is eps-independent; it is kept as a fixed symbol in
      the external kinematics p1.p2 = (m_sigma^2 - 2 m_pi^2)/2,
      p_sigma.p_i = -m_sigma^2/2; no relation between m_pi and S_00 is used.
    * A = -2 B_0 <S_02> > 0; m_a^2 is the root of D_a(s) = s - eps^2 F(s),
      m_a^2 = eps^2 A + 2 c_0/A + ..., c_k = <C L^k(C)>, C = B_0 P_01.
    * m_sigma^2 = 2 lambda f_pi^2/eps^2 + eps^2 kappa_sigma with
      kappa_sigma = B_0 (<S_10> - 2 <S_20>) > 0 (assumption) and lambda fixed.
      In V(sigma^3), V(sigma^4) the coupling lambda enters, i.e.
      m_sigma^2 - eps^2 kappa_sigma = 2 lambda f_pi^2/eps^2 exactly; this is
      implemented by the exact rewriting <S_10> = kappa_sigma/B_0 + 2 <S_20>
      (and <S_02> = -A/(2 B_0)) in every trace product before the series, so
      that the lambda-eliminated vertices of the graph scripts cancel exactly.
      A single trace <S_20> left in a leading coefficient (meson-only sigma
      self-energy) comes from a genuine <S_10> of a completeness sum written
      in these variables.
    * Loop scales per term as in derive_sigma_pipi_loop_graphs.py (D3):
      m_sigma, m_a, m_P and max(...) of them; m_P (largest pole mass of the
      nine meson-like states) is fixed (eps^0), so at large eps it drops out
      of every max containing m_a or m_sigma (resolve_maxima).
      max(m_sigma, m_a) is m_a for A > kappa_sigma and m_sigma for
      A < kappa_sigma; for A = kappa_sigma, m_a^2 - m_sigma^2 -> 2 c_0/A > 0,
      i.e. the first branch.
    * Routing (D12) at large eps: for every product of traces t of the
      canonical generic-source basis (before any rewriting), the coefficient
      C_t^rho of each routing is expanded at large eps; the routing with the
      largest leading power is chosen, then the one with the largest |leading
      coefficient| (exact ratios, with kappa_sigma = r A), then, for equal
      leading terms, the next coefficients times the sign of the leading one.
      Ties of equal magnitude and opposite sign take the sign of the sum of
      the coefficients over all routings; if this sum vanishes, the first
      routing in the canonical order is used and the product is marked (a
      conventional choice at a degenerate mass relation).  The choice is constant in
      regions of r = kappa_sigma/A bounded by the roots of L_i = +-L_j, L_i
      and by r = 1 (break points, i.e. degenerate mass relations, are not
      resolved); it is checked at three points per region.
    * M^{-1}_{m_a^2}(C) = -sum_k s^{-(k+1)} L^k(C), s = m_a^2 (series), and every
      other series, truncated by total order in 1/eps^2 (not by k); the
      truncation order is tracked exactly (TruncatedSeries.precision).
      N = 1 + 2 eps^2 <(M^{-1})^2>.  All pieces of a quantity are summed before
      the leading power is taken, so 1 (ALP-field entry) and -1/N (ALP-field
      component of -v_a v_a^T) combine into 1 - 1/N = O(eps^-2), and explicit
      M^{-1} letters and direct P_01 letters become the same traces.  If every
      coefficient cancels up to the highest order tried, InsufficientOrder is
      raised and no power is reported.
    * Leading power: over the independent products of traces (cyclic words,
      after the assumptions and the charge rule of the graph scripts), in the
      fixed parameters A, kappa_sigma, lambda, B_0, f_pi, m_1^2, m_pi.
 6-10. Outputs: pole and input data, vertices by field-end pattern (lambda
    fixed), leg factors, tree, every kept graph per region (T_Gamma,
    T_Gamma/M_0, X_Gamma^max/M_0, dZ^max), the total R = M_1^est/M_0, the
    contribution of the four pieces of each meson line (bare generator sum
    1_9, ALP-field entry 1 - 1/N, mixed entries -2 eps M^{-1}/N, meson-block
    pole -4 eps^2 M^{-1} (x) M^{-1}/N); each for generic sources and for the
    "r <= 1" sources S_hat_{r>=2,s} = P_hat_{r>=2,s} = 0 (letter substitution
    only); the D4 terms (scale m_P) with the restated criterion.
No loop integral, loop function or loop amplitude is evaluated.
"""
import itertools
import math
import sys
import time

import sympy as sp

import sigma_pipi_loop_algebra as algebra
from sigma_pipi_loop_algebra import (B0, EPSILON, M_A, TraceSum, bare_angular_leg, completeness,
                                     meson_leg, potential_vertex)

SQUARED_MOMENTUM = sp.Symbol('s', real=True)
SINGLET_MASS_SQUARED = sp.Symbol('m_1^2', positive=True)
TRACE_MINV_SQUARED = algebra.TRACE_MINV_SQUARED                  # q = <(M^{-1}(C))^2>
TRACE_S02 = sp.Symbol('<S_02>', real=True)                       # <S_hat_{0,2}>
TRACE_C_MINV = sp.Symbol('<C Minv>', real=True)                  # <C M^{-1}(C)>
ALP_MASS_COEFFICIENT = sp.Symbol('A', positive=True)             # ALP mass entry eps^2 A
UP_QUARK_MASS = sp.Symbol('m_u', positive=True)


# ---------------------------------------------------------------------------
# 1. Quadratic part of the Lagrangian
# ---------------------------------------------------------------------------

def unit_word_product(first, second):
    """TraceSum of <first><second> for two words."""
    result = TraceSum()
    result.add_term((tuple(first), tuple(second)), 1)
    return result


def meson_meson_two_point(first, second, first_coefficient=1, second_coefficient=1):
    """Second derivative of L with respect to the meson fields along the
    matrices `first` and `second` at momenta p, -p (p^2 = s): potential part
    from the vertex rule, kinetic part (s/2)<first second> from
    f_pi^2/4 <u_mu u^mu>, anomalous part -(m_1^2/6)<first><second>."""
    legs = [meson_leg(first, {}, first_coefficient), meson_leg(second, {}, second_coefficient)]
    product = first_coefficient * second_coefficient
    kinetic = TraceSum.single((first, second), SQUARED_MOMENTUM * product / 2)
    anomalous = unit_word_product((first,), (second,)).scale(-SINGLET_MASS_SQUARED * product / 6)
    return (potential_vertex(legs) + kinetic + anomalous).simplified()


def meson_alp_two_point(first, first_coefficient=1):
    """Mixed second derivative (meson matrix `first`, ALP field): vertex rule only;
    there is no kinetic mixing in the connection-free frame."""
    return potential_vertex([meson_leg(first, {}, first_coefficient), bare_angular_leg({})]).simplified()


def alp_alp_two_point():
    """Second derivative with respect to the ALP field: s from 1/2 (d a)^2 plus
    the vertex-rule term 2 B_0 eps^2 <S_02> (minus the ALP mass entry)."""
    return (TraceSum.scalar(SQUARED_MOMENTUM)
            + potential_vertex([bare_angular_leg({}), bare_angular_leg({})])).simplified()


def trace_with_mass_operator(left, argument, argument_coefficient=1):
    """<left M_s(argument)> with M_s(X) = B_0 {S_00, X} + (m_1^2/3) <X> 1 - s X;
    `left` is a word, `argument` a letter."""
    left = tuple(left)
    c = argument_coefficient
    return (TraceSum.single(left + ('S00', argument), B0 * c)
            + TraceSum.single(left + (argument, 'S00'), B0 * c)
            + unit_word_product(left, (argument,)).scale(SINGLET_MASS_SQUARED * c / 3)
            + TraceSum.single(left + (argument,), -SQUARED_MOMENTUM * c)).simplified()


def quadratic_form_checks():
    """Operator-level identities of the quadratic Lagrangian."""
    checks = {}
    checks['D_PP(X,Y) = -(1/2)<X M_s(Y)>'] = (
        meson_meson_two_point('X', 'Y') + trace_with_mass_operator(('X',), 'Y').scale(sp.Rational(1, 2))).is_zero()
    checks['D_aP(X) = eps <C X>, C = B_0 P_01'] = (
        meson_alp_two_point('X') - TraceSum.single(('P01', 'X'), B0 * EPSILON)).is_zero()
    checks['D_aa = s + 2 B_0 eps^2 <S_02>'] = (
        alp_alp_two_point() - TraceSum.scalar(SQUARED_MOMENTUM)
        - TraceSum.single(('S02',), 2 * B0 * EPSILON ** 2)).is_zero()
    checks['M_s self-adjoint: <X M_s(Y)> = <M_s(X) Y>'] = (
        trace_with_mass_operator(('X',), 'Y') - trace_with_mass_operator(('Y',), 'X')).is_zero()
    expansion = trace_with_mass_operator(('X',), 'Y')
    checks['M_s and C contain no epsilon'] = all(not sp.sympify(value).has(EPSILON)
                                                  for value in expansion.terms.values())
    return checks


# ---------------------------------------------------------------------------
# 2. Elimination of the mesons: D_a(s) and N(s)
# ---------------------------------------------------------------------------
# 'Minv' denotes M_s^{-1}(C), defined by M_s(Minv) = C = B_0 P_01, i.e.
# <X M_s(Minv)> = B_0 <X P_01> for every X (M_s invertible, assumption 5).

def defining_relation(left):
    """<left M_s(Minv)> - B_0 <left P_01>, which vanishes by definition of Minv."""
    return trace_with_mass_operator(left, 'Minv') - TraceSum.single(tuple(left) + ('P01',), B0)


def stationarity_check():
    """dL_2/dPhi along a probe matrix X at Phi = a Z, Z = 2 eps Minv:
    D_PP(X, Z) + D_aP(X) = -eps (<X M_s(Minv)> - B_0 <X P_01>) = 0."""
    variation = meson_meson_two_point('X', 'Minv', second_coefficient=2 * EPSILON) + meson_alp_two_point('X')
    return (variation + defining_relation(('X',)).scale(EPSILON)).is_zero()


def alp_inverse_propagator():
    """D_a(s) = D_aa + 2 D_aP(Z) + D_PP(Z, Z) at Z = 2 eps Minv, with
    <Minv M_s(Minv)> replaced by B_0 <Minv P_01> through the defining relation.
    Returns (TraceSum, sympy expression in TRACE_S02 and TRACE_C_MINV)."""
    coefficient = 2 * EPSILON
    meson_part = meson_meson_two_point('Minv', 'Minv', coefficient, coefficient)
    # D_PP(Z,Z) = -(1/2) <Z M_s(Z)> = -2 eps^2 <Minv M_s(Minv)>
    assert (meson_part + trace_with_mass_operator(('Minv',), 'Minv').scale(2 * EPSILON ** 2)).is_zero()
    meson_part = meson_part + defining_relation(('Minv',)).scale(2 * EPSILON ** 2)
    total = (alp_alp_two_point() + meson_alp_two_point('Minv', coefficient).scale(2) + meson_part).simplified()
    expression = 0
    for key, value in total.terms.items():
        if key == ():
            expression += value
        elif key == (('S02',),):
            expression += value * TRACE_S02
        elif key == (algebra.canonical_word(('P01', 'Minv')),):
            expression += value * TRACE_C_MINV / B0
        else:
            raise ValueError(key)
    return total, sp.expand(expression)


def pole_vector_norm():
    """Kinetic norm of the pole vector (Z, 1): (1/2) <Z Z> + 1 with Z = 2 eps Minv,
    as a function of q = <Minv Minv>."""
    meson = TraceSum.single(('Minv', 'Minv'), (2 * EPSILON) ** 2 / 2)
    return sp.expand(algebra.merge_minv_squared(meson + TraceSum.scalar(1)).terms[()])


# ---------------------------------------------------------------------------
# 3. Exact identities of the meson-line residue (any eps)
# ---------------------------------------------------------------------------

NORMALIZATION_EXACT = 1 + 2 * EPSILON ** 2 * TRACE_MINV_SQUARED


def residue_identities():
    """Exact identities of the meson-line residue (1_9 (+) 1) - v_a v_a^T at any
    eps, with N = 1 + 2 eps^2 q and q = <(M^{-1}_{m_a^2}(C))^2>."""
    n = NORMALIZATION_EXACT
    q = TRACE_MINV_SQUARED
    checks = {}
    # v_a^T v_a = (1/2) <(2 eps Minv/sqrt N)^2> + 1/N = 1  => the residue is a projector of rank 9
    vector_norm = algebra.merge_minv_squared(
        TraceSum.single(('Minv', 'Minv'), (2 * EPSILON) ** 2 / (2 * algebra.NORMALIZATION))
        + TraceSum.scalar(1 / algebra.NORMALIZATION)).terms[()]
    checks['v_a^T v_a = 1'] = sp.simplify(vector_norm.subs(algebra.NORMALIZATION, n) - 1) == 0
    checks['N = norm of (2 eps Minv, 1)'] = sp.expand(pole_vector_norm() - n) == 0
    checks['4 eps^2/N = (2/q)(1 - 1/N)'] = sp.simplify(4 * EPSILON ** 2 / n - (2 / q) * (1 - 1 / n)) == 0
    # O_W(X) = W <W X>/2 with W = Minv; P_W = (2/q) O_W.  <Y P_W(P_W(X))> = <Y P_W(X)>:
    once = TraceSum({(('Minv', 'Y'), ('Minv', 'X')): 2 / q / 2})
    twice = TraceSum({(('Minv', 'Y'), ('Minv', 'Minv'), ('Minv', 'X')): (2 / q) ** 2 / 4})
    checks['P_W idempotent'] = (algebra.merge_minv_squared(twice) - algebra.merge_minv_squared(once)).is_zero()
    # trace over the meson space: sum_i (1/2) <lambda_i P_W(lambda_i)> = 1 (completeness)
    trace_value = completeness(TraceSum({(('G1', 'Minv'), ('Minv', 'G1')): (2 / q) / 4}), 'G1')
    checks['trace of P_W over the nine mesons = 1'] = algebra.merge_minv_squared(trace_value).terms == {(): 1}
    # trace of the generator sum 1_9 over the meson space: sum_i (1/2)<lambda_i lambda_i> = 9
    unit_trace = completeness(TraceSum.single(('G1', 'G1'), sp.Rational(1, 2)), 'G1')
    checks['trace of 1_9 = 9'] = unit_trace.simplified().terms == {(): 9}
    return checks


# ---------------------------------------------------------------------------
# 4. Large eps with fixed sources: ALP mass entry eps^2 A
# ---------------------------------------------------------------------------
LOOP_OPERATOR = sp.Symbol('L', commutative=False)     # L(X) = B_0 {S_00, X} + (m_1^2/3) <X> 1
TRACE_C_POWER = sp.symbols('c_0:8', real=True)        # c_k = <C L^k(C)>
INVERSE_EPSILON_SQUARED = sp.Symbol('t', positive=True)   # t = 1/eps^2


def resolvent_series_check(order):
    """(L - s)(-(1/s)) sum_{k<=K} (L/s)^k = 1 - (L/s)^{K+1}, with M_s = L - s."""
    s = SQUARED_MOMENTUM
    series = -sum((LOOP_OPERATOR / s) ** k for k in range(order + 1)) / s
    difference = sp.expand((LOOP_OPERATOR - s) * series - (1 - (LOOP_OPERATOR / s) ** (order + 1)))
    return difference == 0


def pole_function_with_alp_mass_entry():
    """D_a(s) = s - eps^2 F(s), F(s) = A - 2 <C M_s^{-1}(C)>, A = -2 B_0 <S_02>."""
    _, pole_function = alp_inverse_propagator()
    function_f = ALP_MASS_COEFFICIENT - 2 * TRACE_C_MINV
    rewritten = SQUARED_MOMENTUM - EPSILON ** 2 * function_f
    difference = sp.expand(pole_function.subs(TRACE_S02, -ALP_MASS_COEFFICIENT / (2 * B0)) - rewritten)
    return difference == 0


def large_s_traces(order):
    """<C M_s^{-1}(C)> = -sum_k c_k/s^{k+1} and <(M_s^{-1}(C))^2> = sum_n (n+1) c_n/s^{n+2}
    (self-adjoint L: <L^j(C) L^k(C)> = c_{j+k}); truncated at c_order."""
    s = SQUARED_MOMENTUM
    c_minv = -sum(TRACE_C_POWER[k] / s ** (k + 1) for k in range(order + 1))
    minv_squared = sum((n + 1) * TRACE_C_POWER[n] / s ** (n + 2) for n in range(order + 1))
    return c_minv, minv_squared


def alp_pole_series(order=2):
    """Root of D_a(s) = 0 as a series in t = 1/eps^2, and N, xi_a, M^{-1}(C) at the root."""
    t = INVERSE_EPSILON_SQUARED
    s = SQUARED_MOMENTUM
    a = ALP_MASS_COEFFICIENT
    c_minv, minv_squared = large_s_traces(order)
    function_f = a - 2 * c_minv
    unknowns = sp.symbols(f'b_0:{order}')
    root = a / t + sum(unknowns[k] * t ** k for k in range(order))
    # t D_a(root) = t root - F(root)
    equation = sp.series(t * root - function_f.subs(s, root), t, 0, order + 1).removeO()
    solution = sp.solve([sp.expand(equation).coeff(t, k) for k in range(1, order + 1)], unknowns, dict=True)[0]
    root = sp.expand(root.subs(solution))
    normalization_from_trace = 1 + (2 * minv_squared / t).subs(s, root)
    normalization_from_derivative = sp.diff(s - function_f / t, s).subs(s, root)
    n_trace = sp.series(normalization_from_trace, t, 0, order).removeO()
    n_derivative = sp.series(normalization_from_derivative, t, 0, order).removeO()
    minv_leading = sp.series((-1 / s).subs(s, root), t, 0, 2).removeO()   # coefficient of C in M^{-1}(C)
    return {'root': root, 'solution': solution,
            'N from <(M^{-1})^2>': sp.expand(n_trace), 'N from D_a\'': sp.expand(n_derivative),
            'coefficient of C in M^{-1}(C) (leading)': sp.expand(minv_leading),
            'xi_a': sp.expand(sp.series(1 / sp.sqrt(normalization_from_trace), t, 0, order).removeO())}


def up_only_reduce(word):
    """Reduce a trace of a word in the letters 'Pu' (up projector) and 'S00' with
    the projector relations Pu Pu = Pu, S00 Pu = Pu S00 = m_u Pu, <Pu> = 1.
    Words without Pu are returned unreduced (as a TraceSum)."""
    word = tuple(word)
    if 'Pu' not in word:
        return TraceSum.single(word) if word else TraceSum.scalar(algebra.NUMBER_OF_FLAVOURS)
    return TraceSum.scalar(UP_QUARK_MASS ** word.count('S00'))


def up_only_check():
    """Up-only sources S_rs + i P_rs = m_u (-i)^s P_u/s!: S_02 = -m_u P_u/2,
    P_01 = -m_u P_u, so A = B_0 m_u, C = -B_0 m_u P_u; c_0 and c_1 from the
    projector relations; resulting root and N."""
    trace_s02 = -UP_QUARK_MASS / 2            # <S_02> = -m_u <P_u>/2
    a_value = sp.expand(-2 * B0 * trace_s02)
    c_factor = -B0 * UP_QUARK_MASS             # C = c_factor P_u
    # c_0 = <C C>; c_1 = <C L(C)> = c_factor^2 (B_0 <Pu {S00, Pu}> + (m_1^2/3) <Pu><Pu>)
    c0 = c_factor ** 2 * up_only_reduce(('Pu', 'Pu')).terms[()]
    c1 = c_factor ** 2 * (B0 * (up_only_reduce(('Pu', 'S00', 'Pu')).terms[()]
                                + up_only_reduce(('Pu', 'Pu', 'S00')).terms[()])
                          + SINGLET_MASS_SQUARED / 3 * up_only_reduce(('Pu',)).terms[()] ** 2)
    generic = alp_pole_series(order=2)
    substitution = {ALP_MASS_COEFFICIENT: a_value, TRACE_C_POWER[0]: sp.expand(c0), TRACE_C_POWER[1]: sp.expand(c1)}
    return {'A': a_value, 'c_0': sp.expand(c0), 'c_1': sp.expand(c1),
            'root': sp.expand(generic['root'].subs(substitution)),
            'N': sp.expand(generic['N from <(M^{-1})^2>'].subs(substitution)),
            # coefficient of P_u in the meson matrix 2 eps M^{-1}(C)/sqrt(N) of v_a, leading term (eps = t^{-1/2})
            'coefficient of P_u in the meson matrix of v_a (leading)': sp.expand(
                2 / sp.sqrt(INVERSE_EPSILON_SQUARED)
                * generic['coefficient of C in M^{-1}(C) (leading)'].subs(substitution) * c_factor)}


# ---------------------------------------------------------------------------
# 5. Every quantity of the estimate at large eps: truncated series in u = 1/eps
# ---------------------------------------------------------------------------
INVERSE_EPSILON = sp.Symbol('u', positive=True)            # u = 1/eps, only used in printed series
KAPPA_SIGMA = sp.Symbol('kappa_sigma', positive=True)      # B_0 (<S_10> - 2 <S_20>), assumed > 0
QUARTIC = algebra.QUARTIC
F_PI = algebra.F_PI
M_PI = algebra.M_PI
M_SIGMA = algebra.M_SIGMA
NORMALIZATION = algebra.NORMALIZATION
BRANCHES = ('A > kappa_sigma', 'A < kappa_sigma')
SOURCE_SETS = ('generic', 'r <= 1')
ORDERS = (2, 3, 4)                                         # relative orders in u^2 tried in turn
LATEX_NAMES_LARGE_EPSILON = {**algebra.LATEX_NAMES, KAPPA_SIGMA: r'\kappa_\sigma', ALP_MASS_COEFFICIENT: 'A',
                 SINGLET_MASS_SQUARED: 'm_1^2', B0: 'B_0', UP_QUARK_MASS: 'm_u'}


M_MESON = algebra.M_MESON


def resolve_maxima(expression, branch):
    """Large-eps form of every max(...) of loop scales: m_P is
    fixed (eps^0) while m_a and m_sigma grow as eps, so m_P drops out of any
    max that contains m_a or m_sigma (exact for eps large enough);
    max(m_a, m_sigma) is m_a in the branch A > kappa_sigma and m_sigma in the
    branch A < kappa_sigma."""
    larger = {BRANCHES[0]: M_A, BRANCHES[1]: M_SIGMA}.get(branch)

    def replacement(*arguments):
        growing = [argument for argument in arguments if argument != M_MESON]
        if not growing:
            return M_MESON
        if len(set(growing)) == 1:
            return growing[0]
        assert set(growing) == {M_A, M_SIGMA} and larger is not None, 'branch of max(m_sigma, m_a) needed'
        return larger
    return expression.replace(sp.Max, replacement)


class InsufficientOrder(Exception):
    """Every coefficient of a series vanishes up to its truncation order."""


def merge_keys(first, second):
    return tuple(sorted(first + second))


def word_key(word):
    """(key, factor) of the trace of one word: the empty word gives N_f = 3."""
    word = tuple(word)
    if not word:
        return (), algebra.NUMBER_OF_FLAVOURS
    return (algebra.canonical_word(word),), 1


class TruncatedSeries:
    """Laurent series in u = 1/eps whose coefficients are linear combinations of
    products of traces (keys as in TraceSum), known exactly up to and including
    u^precision (precision = inf for an exact expression).
    terms: {key: {exponent of u: coefficient independent of eps}}."""

    def __init__(self, precision=math.inf):
        self.terms = {}
        self.precision = precision

    @staticmethod
    def monomial(value=1, exponent=0, key=(), precision=math.inf):
        result = TruncatedSeries(precision)
        result.add(key, exponent, value)
        return result

    @staticmethod
    def from_trace_sum(trace_sum, exponent=0):
        result = TruncatedSeries()
        for key, value in trace_sum.terms.items():
            result.add(key, exponent, value)
        return result

    def add(self, key, exponent, value):
        if exponent > self.precision:
            return
        inner = self.terms.setdefault(key, {})
        inner[exponent] = inner.get(exponent, 0) + value

    def accumulate(self, other):
        """In-place sum."""
        self.precision = min(self.precision, other.precision)
        for key, inner in other.terms.items():
            for exponent, value in inner.items():
                self.add(key, exponent, value)
        self.terms = {key: {e: v for e, v in inner.items() if e <= self.precision}
                      for key, inner in self.terms.items()}
        return self

    def copy(self):
        result = TruncatedSeries(self.precision)
        result.terms = {key: dict(inner) for key, inner in self.terms.items()}
        return result

    def lowest(self):
        exponents = [exponent for inner in self.terms.values() for exponent in inner]
        return min(exponents) if exponents else self.precision + 1

    def truncated(self, precision):
        result = TruncatedSeries(min(precision, self.precision))
        for key, inner in self.terms.items():
            for exponent, value in inner.items():
                result.add(key, exponent, value)
        return result

    def __add__(self, other):
        return self.copy().accumulate(other)

    def __sub__(self, other):
        return self.copy().accumulate(other.scale(-1))

    def scale(self, value):
        """Multiply by an eps-independent factor."""
        result = TruncatedSeries(self.precision)
        result.terms = {key: {e: value * v for e, v in inner.items()} for key, inner in self.terms.items()}
        return result

    def shift(self, exponent):
        """Multiply by u^exponent."""
        result = TruncatedSeries(self.precision + exponent)
        result.terms = {key: {e + exponent: v for e, v in inner.items()} for key, inner in self.terms.items()}
        return result

    def with_key(self, key, factor=1):
        """Multiply by the product of traces `key` (and a number)."""
        result = TruncatedSeries(self.precision)
        for own, inner in self.terms.items():
            merged = merge_keys(own, key)
            for exponent, value in inner.items():
                result.add(merged, exponent, factor * value)
        return result

    def __mul__(self, other):
        precision = min(self.precision + other.lowest(), other.precision + self.lowest())
        result = TruncatedSeries(precision)
        for key_one, inner_one in self.terms.items():
            for key_two, inner_two in other.terms.items():
                key = merge_keys(key_one, key_two)
                for exponent_one, value_one in inner_one.items():
                    for exponent_two, value_two in inner_two.items():
                        if exponent_one + exponent_two <= precision:
                            result.add(key, exponent_one + exponent_two, value_one * value_two)
        return result.expanded()

    def power(self, exponent):
        result = TruncatedSeries.monomial(1)
        for _ in range(exponent):
            result = result * self
        return result

    def expanded(self):
        result = TruncatedSeries(self.precision)
        for key, inner in self.terms.items():
            for exponent, value in inner.items():
                value = sp.expand(value)
                if value != 0:
                    result.terms.setdefault(key, {})[exponent] = value
        return result

    def filtered(self, predicate):
        result = TruncatedSeries(self.precision)
        result.terms = {key: dict(inner) for key, inner in self.terms.items() if predicate(key)}
        return result

    def coefficient(self, exponent):
        """TraceSum multiplying u^exponent."""
        result = TraceSum()
        for key, inner in self.terms.items():
            if exponent in inner:
                result.add_term(key, inner[exponent])
        return result.simplified()

    def leading(self):
        """(power of eps, TraceSum of the leading coefficient); raises
        InsufficientOrder if every coefficient vanishes up to the precision."""
        exponents = sorted({e for inner in self.expanded().terms.values() for e in inner})
        if not exponents:
            raise InsufficientOrder
        return -exponents[0], self.coefficient(exponents[0])

    def as_expression(self):
        """sympy expression in u with one symbol per trace (for printing and tests)."""
        total = 0
        for key, inner in self.expanded().terms.items():
            traces = sp.Mul(*(trace_symbol(word) for word in key))
            for exponent, value in inner.items():
                total += value * traces * INVERSE_EPSILON ** exponent
        return total


def trace_symbol(word):
    return sp.Symbol(r'\langle ' + r'\,'.join(algebra.letter_latex(name) for name in word) + r'\rangle')


def one_plus_power(correction, exponent):
    """(1 + correction)^exponent as a binomial series; correction has lowest
    power of u >= 1 and a finite precision."""
    correction = correction.expanded()
    assert correction.lowest() >= 1 and correction.precision < math.inf
    result = TruncatedSeries.monomial(1, precision=correction.precision)
    term = TruncatedSeries.monomial(1)
    order = 0
    while True:
        order += 1
        term = term * correction
        if term.lowest() > correction.precision:
            break
        result.accumulate(term.scale(sp.binomial(exponent, order)))
    return result.expanded()


def inverse(series, leading_value, leading_exponent):
    """1/series, for series = leading_value u^leading_exponent (1 + O(u))."""
    delta = (series - TruncatedSeries.monomial(leading_value, leading_exponent)).scale(1 / leading_value)
    delta = delta.shift(-leading_exponent)
    return one_plus_power(delta, -1).scale(1 / leading_value).shift(-leading_exponent)


def apply_loop_operator(matrix):
    """L(X) = B_0 {S_00, X} + (m_1^2/3) <X> 1 for X = {(word, key): coefficient}
    (X = sum coefficient * word * product of traces `key`)."""
    result = {}

    def add(word, key, value):
        result[(word, key)] = result.get((word, key), 0) + value
    for (word, key), value in matrix.items():
        add(('S00',) + word, key, B0 * value)
        add(word + ('S00',), key, B0 * value)
        trace, factor = word_key(word)
        add((), merge_keys(key, trace), factor * SINGLET_MASS_SQUARED * value / 3)
    return result


def substitute_mass_traces(trace_sum):
    """Exact rewriting in the fixed parameters A and kappa_sigma:
    <S_02> = -A/(2 B_0) and <S_10> = kappa_sigma/B_0 + 2 <S_20>."""
    replacements = {('S02',): TraceSum.scalar(-ALP_MASS_COEFFICIENT / (2 * B0)),
                    ('S10',): TraceSum.scalar(KAPPA_SIGMA / B0) + TraceSum.single(('S20',), 2)}
    result = TraceSum()
    for key, value in trace_sum.terms.items():
        product = TraceSum.scalar(value)
        for word in key:
            product = product * replacements.get(word, TraceSum.single(word))
        result = result + product
    return result.simplified()


def radial_index_at_most_one(key):
    return all(not (algebra.is_source(name) and int(name[1]) >= 2) for word in key for name in word)


def restrict_sources(trace_sum, sources):
    """'generic' or 'r <= 1' (S_hat_{r,s} = P_hat_{r,s} = 0 for r >= 2, by letter only)."""
    if sources == 'generic':
        return trace_sum
    return trace_sum.filtered(radial_index_at_most_one)


class FixedSourceSeries:
    """Exact series of the pole data at large eps (fixed sources, fixed lambda), each known to relative
    order u^(2 order) beyond its leading power:
      m_a^2 from the fixed-point form s = eps^2 (A + 2 sum_k c_k / s^{k+1}),
      M^{-1}_{m_a^2}(C) = -sum_k s^{-(k+1)} L^k(C)   (truncated by total order in u),
      N = 1 + 2 eps^2 <(M^{-1}(C))^2>,
      m_sigma^2 = 2 lambda f_pi^2 u^2 + kappa_sigma u^{-2}  (exact)."""

    def __init__(self, order):
        self.order = order
        powers = [{(('P01',), ()): B0}]                       # C = B_0 P_01
        for _ in range(order + 1):
            powers.append(apply_loop_operator(powers[-1]))
        self.operator_powers = powers
        self.trace_powers = []
        for matrix in powers:
            trace = TraceSum()
            for (word, key), value in matrix.items():
                trace.add_term(key + (('P01',) + word,), B0 * value)
            self.trace_powers.append(trace.simplified())
        self.alp_mass_squared = self._root()
        self.inverse_alp_mass_squared = inverse(self.alp_mass_squared, ALP_MASS_COEFFICIENT, -2)
        self.inverse_matrix = self._inverse_matrix()
        self.word_cache, self.key_cache, self.factor_cache = {}, {}, {}
        self.normalization = (TruncatedSeries.monomial(1)
                              + self.key_series((('Minv', 'Minv'),)).scale(2).shift(-2)).expanded()
        radial = TruncatedSeries()
        radial.add((), -2, KAPPA_SIGMA)
        radial.add((), 2, 2 * QUARTIC * F_PI ** 2)
        self.radial_mass_squared = radial

    def _root(self):
        a = ALP_MASS_COEFFICIENT
        root = TruncatedSeries.monomial(a, -2, precision=-2)
        for _ in range(self.order):
            inverse_root = inverse(root, a, -2)
            total = TruncatedSeries.monomial(a)
            inverse_power = inverse_root
            for k in range(self.order + 1):
                total.accumulate((inverse_power * TruncatedSeries.from_trace_sum(self.trace_powers[k])).scale(2))
                inverse_power = inverse_power * inverse_root
            root = total.shift(-2).expanded()
        return root

    def _inverse_matrix(self):
        precision = 2 + 2 * self.order
        matrix = {}
        inverse_power = self.inverse_alp_mass_squared
        for k in range(self.order + 1):
            for (word, key), value in self.operator_powers[k].items():
                piece = inverse_power.with_key(key, -value).truncated(precision)
                if word in matrix:
                    matrix[word].accumulate(piece)
                else:
                    matrix[word] = piece
            inverse_power = inverse_power * self.inverse_alp_mass_squared
        return {word: series.expanded() for word, series in matrix.items()}

    def word_series(self, word):
        """Series of the trace of one word, every letter 'Minv' replaced by its series."""
        if word in self.word_cache:
            return self.word_cache[word]
        if 'Minv' not in word:
            key, factor = word_key(word)
            result = TruncatedSeries.monomial(factor, key=key)
        else:
            partial = {(): TruncatedSeries.monomial(1)}
            target = 2 * word.count('Minv') + 2 * self.order
            for name in word:
                if name != 'Minv':
                    partial = {prefix + (name,): value for prefix, value in partial.items()}
                    continue
                new = {}
                for prefix, coefficient in partial.items():
                    for piece_word, piece in self.inverse_matrix.items():
                        if coefficient.lowest() + piece.lowest() > target:
                            continue
                        product = coefficient * piece
                        if not product.terms:
                            continue
                        if prefix + piece_word in new:
                            new[prefix + piece_word].accumulate(product)
                        else:
                            new[prefix + piece_word] = product
                partial = new
            result = TruncatedSeries(target)
            for prefix, coefficient in partial.items():
                key, factor = word_key(prefix)
                result.accumulate(coefficient.with_key(key, factor).truncated(target))
            result = result.expanded()
        self.word_cache[word] = result
        return result

    def key_series(self, key):
        if key not in self.key_cache:
            result = TruncatedSeries.monomial(1)
            for word in key:
                result = result * self.word_series(word)
            self.key_cache[key] = result
        return self.key_cache[key]

    def factor_series(self, exponents):
        """u^{-a} N^b (m_a^2)^{c/2} (m_sigma^2)^{d/2} for exponents (a, b, c, d)."""
        if exponents in self.factor_cache:
            return self.factor_cache[exponents]
        a, b, c, d = exponents
        assert c % 2 == 0 and d % 2 == 0, exponents
        result = TruncatedSeries.monomial(1, -a)
        if b != 0:
            result = result * one_plus_power(self.normalization - TruncatedSeries.monomial(1), b)
        if c > 0:
            result = result * self.alp_mass_squared.power(c // 2)
        elif c < 0:
            result = result * self.inverse_alp_mass_squared.power(-c // 2)
        if d > 0:
            result = result * self.radial_mass_squared.power(d // 2)
        elif d < 0:
            inverse_radial = inverse(self.radial_mass_squared.truncated(-2 + 2 * self.order), KAPPA_SIGMA, -2)
            result = result * inverse_radial.power(-d // 2)
        self.factor_cache[exponents] = result
        return result

    def coefficient_series(self, coefficient):
        groups = {}
        for term in sp.Add.make_args(sp.expand(coefficient)):
            powers = term.as_powers_dict()
            exponents = tuple(powers.get(symbol, 0) for symbol in (EPSILON, NORMALIZATION, M_A, M_SIGMA))
            rest = (term * EPSILON ** -exponents[0] * NORMALIZATION ** -exponents[1]
                    * M_A ** -exponents[2] * M_SIGMA ** -exponents[3])
            assert not rest.has(EPSILON, NORMALIZATION, M_A, M_SIGMA, sp.Max), term
            groups[exponents] = groups.get(exponents, 0) + rest
        total = TruncatedSeries()
        for exponents, rest in groups.items():
            total.accumulate(self.factor_series(exponents).scale(rest))
        return total

    def series(self, trace_sum, branch=None):
        """Series of a TraceSum of the graph scripts (words with 'Minv', coefficients
        in eps, N, m_a, m_sigma, m_P and max(...) of them); branch selects
        max(m_sigma, m_a) (resolve_maxima)."""
        trace_sum = substitute_mass_traces(trace_sum)
        total = TruncatedSeries()
        for key, coefficient in trace_sum.terms.items():
            coefficient = resolve_maxima(sp.sympify(coefficient), branch)
            total.accumulate(self.coefficient_series(coefficient) * self.key_series(key))
        return total.expanded()


def series_identities(order):
    """Exact checks of the large-eps series up to their precision:
    pole condition m_a^2 - eps^2 (A - 2 <C M^{-1}(C)>) = 0, defining relation
    <Y M_{m_a^2}(M^{-1}(C))> - B_0 <Y P_01> = 0 for a probe Y, and agreement of
    the root with alp_pole_series(order=2) after c_k -> traces."""
    ctx = FixedSourceSeries(order)
    s = ctx.alp_mass_squared
    pole = (s - (TruncatedSeries.monomial(ALP_MASS_COEFFICIENT)
                 + ctx.word_series(("P01", "Minv")).scale(-2 * B0)).shift(-2)).expanded()
    relation = (ctx.word_series(('Y', 'S00', 'Minv')).scale(B0) + ctx.word_series(('Y', 'Minv', 'S00')).scale(B0)
                + (ctx.word_series(('Y',)) * ctx.word_series(('Minv',))).scale(SINGLET_MASS_SQUARED / 3)
                - s * ctx.word_series(('Y', 'Minv'))
                - TruncatedSeries.monomial(B0, key=(('P01', 'Y'),))).expanded()
    reference = alp_pole_series(order=2)['root'].subs(INVERSE_EPSILON_SQUARED, INVERSE_EPSILON ** 2)
    traces = {TRACE_C_POWER[k]: as_trace_expression(ctx.trace_powers[k]) for k in range(2)}
    root = s.truncated(2).as_expression()
    return {'pole condition': not pole.terms, 'defining relation of M^{-1}(C)': not relation.terms,
            'root agrees with alp_pole_series(order=2)': sp.expand(root - reference.subs(traces)) == 0,
            'precisions': (s.precision, relation.precision)}


CONTEXTS = {}


def context(order):
    if order not in CONTEXTS:
        CONTEXTS[order] = FixedSourceSeries(order)
    return CONTEXTS[order]


def leading_of(build, nonzero=True):
    """Leading power and coefficient of the series build(context), trying the
    orders in ORDERS in turn.  Raises InsufficientOrder if every coefficient
    cancels up to the highest order tried."""
    if not nonzero:
        return None
    for order in ORDERS:
        series = build(context(order))
        try:
            power, coefficient = series.leading()
        except InsufficientOrder:
            continue
        return {'power': power, 'coefficient': coefficient, 'order': order}
    raise InsufficientOrder('cancellation persists up to the highest order tried')


def leading_trace_sum(trace_sum, branch=None, sources='generic', multiplier=None):
    """Leading behaviour at large eps of a TraceSum of the graph scripts;
    None if it vanishes exactly for the chosen sources."""
    restricted = restrict_sources(trace_sum, sources)
    if sources != 'generic':
        restricted = restrict_sources(substitute_mass_traces(restricted), sources)

    def build(ctx):
        series = ctx.series(restricted, branch)
        if sources != 'generic':
            series = series.filtered(radial_index_at_most_one)
        if multiplier is not None:
            series = series * multiplier(ctx)
        return series
    return leading_of(build, not restricted.is_zero())


# ---------------------------------------------------------------------------
# 6. Item 1: pole and input data
# ---------------------------------------------------------------------------

def probe(letter_name):
    return TraceSum.single((letter_name, 'Minv'))


def pole_data_series(ctx):
    """{name: (series, kind)}; kind 'scalar', 'matrix' (probe Y) or 'tensor' (probes Y1, Y2)."""
    one = TruncatedSeries.monomial(1)
    n = ctx.normalization
    inverse_n = one_plus_power(n - one, -1)
    root_inverse_n = one_plus_power(n - one, sp.Rational(-1, 2))
    matrix = ctx.word_series(('Y', 'Minv'))
    tensor = ctx.word_series(('Y1', 'Minv')) * ctx.word_series(('Y2', 'Minv'))
    return {
        'm_a^2': (ctx.alp_mass_squared, 'scalar'),
        'm_sigma^2': (ctx.radial_mass_squared, 'scalar'),
        'm_a^2 - m_sigma^2': ((ctx.alp_mass_squared - ctx.radial_mass_squared).truncated(
            ctx.alp_mass_squared.precision), 'scalar'),
        'M^{-1}_{m_a^2}(C)': (matrix, 'matrix'),
        'N(m_a^2)': (n, 'scalar'),
        'N - 1': (n - one, 'scalar'),
        '1/N': (inverse_n, 'scalar'),
        'xi_a = 1/sqrt(N)': (root_inverse_n, 'scalar'),
        'v_a, meson matrix 2 eps M^{-1}/sqrt(N)': ((matrix * root_inverse_n).scale(2).shift(-1), 'matrix'),
        'v_a, ALP-field component 1/sqrt(N)': (root_inverse_n, 'scalar'),
        'residue (i): bare generator sum 1_9': (one, 'scalar'),
        'residue (ii): ALP-field entry 1 - 1/N': (one - inverse_n, 'scalar'),
        'residue (iii): mixed entry -2 eps M^{-1}/N': ((matrix * inverse_n).scale(-2).shift(-1), 'matrix'),
        'residue (iv): meson-block pole -4 eps^2 M^{-1} (x) M^{-1}/N':
            ((tensor * inverse_n).scale(-4).shift(-2), 'tensor'),
    }


def pole_data_table():
    """Leading power and coefficient of every pole and input quantity (item 1)."""
    rows = {}
    names = list(pole_data_series(context(ORDERS[0])))
    for name in names:
        result = leading_of(lambda ctx, name=name: pole_data_series(ctx)[name][0])
        series, kind = pole_data_series(context(result['order']))[name]
        rows[name] = {**result, 'kind': kind, 'series': series}
    # degenerate case A = kappa_sigma
    difference = rows['m_a^2 - m_sigma^2']['series']
    degenerate = TruncatedSeries(difference.precision)
    for key, inner in difference.terms.items():
        for exponent, value in inner.items():
            degenerate.add(key, exponent, value.subs(KAPPA_SIGMA, ALP_MASS_COEFFICIENT))
    power, coefficient = degenerate.leading()
    rows['m_a^2 - m_sigma^2 (A = kappa_sigma)'] = {'power': power, 'coefficient': coefficient,
                                                   'kind': 'scalar', 'series': degenerate.expanded()}
    rows['m_a^2 - m_sigma^2 (up-only)'] = up_only_mass_difference(difference)
    return rows


def up_only_value(trace_sum):
    """Up-only sources: P_01 = -m_u P_u, S_10 = m_u P_u (kappa_sigma = B_0 m_u),
    with the projector relations of up_only_reduce."""
    total = 0
    for key, value in trace_sum.terms.items():
        product = value
        for word in key:
            factor = (-UP_QUARK_MASS) ** word.count('P01')
            product *= factor * up_only_reduce(tuple('Pu' if name == 'P01' else name for name in word)).terms[()]
        total += product
    return sp.expand(total.subs({ALP_MASS_COEFFICIENT: B0 * UP_QUARK_MASS, KAPPA_SIGMA: B0 * UP_QUARK_MASS}))


def up_only_mass_difference(difference):
    exponents = sorted({e for inner in difference.terms.values() for e in inner})
    values = {exponent: up_only_value(difference.coefficient(exponent)) for exponent in exponents}
    nonzero = [exponent for exponent in exponents if values[exponent] != 0]
    series = sum(values[exponent] * INVERSE_EPSILON ** exponent for exponent in exponents)
    return {'power': -nonzero[0], 'coefficient': TraceSum.scalar(values[nonzero[0]]), 'kind': 'scalar',
            'up-only series': sp.expand(series), 'precision': difference.precision}


# ---------------------------------------------------------------------------
# 7. Item 2: vertices by field-end pattern, lambda fixed
# ---------------------------------------------------------------------------

def epsilon_laurent_exponents(trace_sum, loop_momenta=False):
    """Exponents of eps in a vertex whose coefficients are Laurent polynomials in
    eps; with loop_momenta every product of momenta k_i.k_j counts as eps^2."""
    exponents = set()
    for value in trace_sum.simplified().terms.values():
        for term in sp.Add.make_args(sp.expand(value)):
            powers = term.as_powers_dict()
            exponent = powers.get(EPSILON, 0)
            if loop_momenta:
                exponent += 2 * sum(power for symbol, power in powers.items()
                                    if isinstance(symbol, sp.Symbol) and symbol.name.startswith('d_'))
            exponents.add(exponent)
    return sorted(exponents)


def vertex_power_table():
    """{field ends: {'generic', 'r <= 1'} x {'fixed momenta', 'loop momenta'}: exponents}
    for the 46 patterns n = 3, 4, 5, with lambda kept as a fixed coupling."""
    import derive_sigma_pipi_loop_vertices as vertices
    table = {}
    for n in range(3, vertices.MAX_LEGS + 1):
        for kinds in itertools.combinations_with_replacement(vertices.FIELD_ENDS, n):
            legs = [vertices.field_end_leg(kind, index + 1) for index, kind in enumerate(kinds)]
            value = algebra.apply_assumptions(algebra.vertex(legs, algebra.symbolic_dot, eliminate_lambda=False))
            row = {}
            for sources in SOURCE_SETS:
                restricted = restrict_sources(value, sources)
                for momenta in ('fixed momenta', 'loop momenta'):
                    row[(sources, momenta)] = epsilon_laurent_exponents(restricted, momenta == 'loop momenta')
            table[kinds] = row
    return table


def radial_cubic_vertex_forms():
    """V(sigma^3) with lambda f_a = lambda f_pi/eps, and the same vertex built from the
    lambda-eliminated rule with m_sigma^2 = 2 lambda f_pi^2/eps^2 + eps^2 kappa_sigma
    substituted (they must agree exactly)."""
    legs = [algebra.radial_leg({f'q{i}': 1}) for i in range(1, 4)]
    with_lambda = algebra.vertex(legs, algebra.symbolic_dot, eliminate_lambda=False)
    eliminated = substitute_mass_traces(algebra.vertex(legs, algebra.symbolic_dot, eliminate_lambda=True))
    eliminated = eliminated.map_coefficients(lambda value: sp.expand(
        sp.sympify(value).subs(M_SIGMA, sp.sqrt(2 * QUARTIC * F_PI ** 2 / EPSILON ** 2
                                                 + EPSILON ** 2 * KAPPA_SIGMA))))
    return with_lambda.simplified(), eliminated.simplified()


# ---------------------------------------------------------------------------
# 8. Item 3: leg factors
# ---------------------------------------------------------------------------
LEG_FACTORS = (
    ('ALP-field component of an ALP-like end, xi_a', 'xi_a = 1/sqrt(N)', 1),
    ('meson component of v_a', 'v_a, meson matrix 2 eps M^{-1}/sqrt(N)', 0),
    ('meson-line residue (i): generator-sum entry 1_9', 'residue (i): bare generator sum 1_9', 0),
    ('meson-line residue (ii): ALP-field entry 1 - 1/N', 'residue (ii): ALP-field entry 1 - 1/N', 2),
    ('meson-line residue (iii): mixed entry', 'residue (iii): mixed entry -2 eps M^{-1}/N', 1),
    ('meson-line residue (iv): meson-block pole', 'residue (iv): meson-block pole -4 eps^2 M^{-1} (x) M^{-1}/N', 0),
)


def leg_factor_table(pole_rows=None):
    """Leading power of each leg factor at large eps and the net power for each
    placement of its ALP-field ends: +1 for an ALP-field end at a source vertex,
    0 at a kinetic sigma-a-a or sigma-sigma-a-a vertex (whose k_x.k_y is counted
    with the vertex)."""
    pole_rows = pole_data_table() if pole_rows is None else pole_rows
    rows = []
    for name, key, alp_field_ends in LEG_FACTORS:
        power = pole_rows[key]['power']
        net = {f'{alp_field_ends - kinetic} at source, {kinetic} at kinetic vertices': power + alp_field_ends - kinetic
               for kinetic in range(alp_field_ends + 1)}
        rows.append({'name': name, 'alp_field_ends': alp_field_ends, 'power': power,
                     'coefficient': pole_rows[key]['coefficient'], 'kind': pole_rows[key]['kind'], 'net': net})
    return rows


# ---------------------------------------------------------------------------
# 9. Items 4-6: tree, graphs, total
# ---------------------------------------------------------------------------

def tree_leading():
    import derive_sigma_pipi_loop_tree as tree_script
    tree = tree_script.tree_amplitude()
    return tree, leading_trace_sum(tree)


def graph_setups():
    import derive_sigma_pipi_loop_graphs as graphs
    return {'one_particle_irreducible': graphs.three_point_setup(), 'sigma': graphs.two_point_setup('sigma'),
            'pi+': graphs.two_point_setup('pi+'), 'pi-': graphs.two_point_setup('pi-')}


def graph_entries():
    """Kept graphs with the size-rule value of every routing (per-term scales of
    D3; 1PI: S (m_term^2)^{2-r} T/(16 pi^2), self-energies: dZ) and the unsized
    T of every routing; from derive_sigma_pipi_loop_graphs.py."""
    import derive_sigma_pipi_loop_graphs as graphs
    data = graphs.summary()
    setups = graph_setups()
    entries = []
    number = {'one_particle_irreducible': 0, 'sigma': 0, 'pi+': 0, 'pi-': 0}
    prefix = {'one_particle_irreducible': '', 'sigma': 'sigma', 'pi+': 'pi', 'pi-': 'pi-'}
    for group in ('one_particle_irreducible', 'sigma', 'pi+', 'pi-'):
        legs, basis_dot = setups[group]
        for row in graphs.kept(data[group]):
            number[group] += 1
            unsized = [graphs.graph_value(row, legs, basis_dot, start=start, orientation=orientation,
                                          pieces=row['pieces'])
                       for start, orientation in graphs.all_routings(row)]
            entries.append({**row, 'group': group, 'name': f'{prefix[group]}{number[group]}',
                            'r': len(row['blocks']), 'unsized routing values': unsized})
    return entries, data


# ---------------------------------------------------------------------------
# Routing maximum at large eps (D12): largest |C_t^rho| at eps -> oo
# ---------------------------------------------------------------------------
RATIO = sp.Symbol('r', positive=True)          # r = kappa_sigma / A
DECISION_ORDER = ORDERS[0]


def branch_of(ratio):
    """Branch of max(m_sigma, m_a) for r = kappa_sigma/A (r = 1 belongs to the first branch)."""
    return BRANCHES[0] if ratio <= 1 else BRANCHES[1]


EXPANSION_CACHE = {}


def symbolic_expansion(ctx, coefficient, branch):
    """([(power of eps, expression in A, kappa_sigma, ...)] in decreasing power,
    lowest known power) of a scalar coefficient at large eps (cached)."""
    key = (ctx.order, coefficient, branch)
    if key not in EXPANSION_CACHE:
        series = ctx.coefficient_series(resolve_maxima(sp.sympify(coefficient), branch)).expanded()
        terms = []
        for exponent in sorted({e for inner in series.terms.values() for e in inner}):
            value = sp.expand(as_trace_expression(series.coefficient(exponent)))
            if value != 0:
                terms.append((-exponent, value))
        EXPANSION_CACHE[key] = (terms, -series.precision)
    return EXPANSION_CACHE[key]


def expansion(ctx, coefficient, branch, ratio):
    """symbolic_expansion with kappa_sigma = ratio * A."""
    coefficient = sp.expand(coefficient)
    if coefficient == 0:
        return [], -math.inf
    terms, lowest = symbolic_expansion(ctx, coefficient, branch)
    substituted = []
    for power, value in terms:
        value = sp.expand(value.subs(KAPPA_SIGMA, ratio * ALP_MASS_COEFFICIENT))
        if value != 0:
            substituted.append((power, value))
    return substituted, lowest


def sign_of(expression):
    expression = sp.factor(expression)
    if expression.is_positive:
        return 1
    if expression.is_negative:
        return -1
    return None


def largest_routing_at_large_epsilon(coefficients, ctx, branch, ratio):
    """(index, mark) of the routing with the largest |C_t^rho| at eps -> oo:
    largest leading power of eps, then largest |leading coefficient| (exact
    ratios), then, for equal leading terms, the next coefficients times the sign
    of the leading one.  Ties of equal magnitude and opposite sign take the sign
    of the sum of the coefficients over all routings; if this sum vanishes, the
    first routing is used with the mark 'sum over routings zero'.  mark is None
    or a string."""
    expansions = [expansion(ctx, c, branch, ratio) for c in coefficients]
    candidates = [i for i, (terms, _) in enumerate(expansions) if terms]
    if not candidates:
        return 0, None
    top = max(expansions[i][0][0][0] for i in candidates)
    candidates = [i for i in candidates if expansions[i][0][0][0] == top]
    reference = expansions[candidates[0]][0][0][1]
    ratios = {i: sp.cancel(expansions[i][0][0][1] / reference) for i in candidates}
    if any(sp.sympify(value).free_symbols for value in ratios.values()):
        raise ValueError(f'leading coefficients not comparable: {ratios}')
    largest = max(abs(value) for value in ratios.values())
    group = [i for i in candidates if abs(ratios[i]) == largest]
    if len({ratios[i] for i in group}) > 1:
        routing_sum, _ = expansion(ctx, sum(coefficients), branch, ratio)
        if not routing_sum:
            return group[0], 'sum over routings zero'
        sum_over_reference = sp.cancel(routing_sum[0][1] / reference)
        if sum_over_reference.free_symbols:
            return group[0], 'sign of the sum over routings undecided'
        return next(i for i in group if ratios[i] * sum_over_reference > 0), None
    lead_sign = sign_of(reference * ratios[group[0]])
    known = max(expansions[i][1] for i in group)
    powers = sorted({power for i in group for power, _ in expansions[i][0][1:] if power >= known}, reverse=True)
    remaining = list(group)
    for power in powers:
        values = {i: dict(expansions[i][0]).get(power, 0) for i in remaining}
        if len({sp.expand(value) for value in values.values()}) == 1:
            continue
        if lead_sign is None:
            return remaining[0], 'undecided (sign of the leading coefficient)'
        best = []
        for i in remaining:
            signs = [sign_of(lead_sign * (values[i] - values[j])) for j in remaining if j != i
                     and sp.expand(values[i] - values[j]) != 0]
            if None in signs:
                return remaining[0], 'undecided (subleading coefficients)'
            if all(sign > 0 for sign in signs):
                best.append(i)
        tied = [i for i in remaining if sp.expand(values[i] - values[best[0]]) == 0] if best else remaining
        remaining = tied
        if len(remaining) == 1:
            return remaining[0], None
    return remaining[0], None


def routing_dependent_keys(per_routing):
    keys = set()
    for value in per_routing:
        keys |= set(value.terms)
    result = []
    for key in sorted(keys, key=repr):
        coefficients = [sp.expand(value.terms.get(key, 0)) for value in per_routing]
        if any(sp.expand(c - coefficients[0]) != 0 for c in coefficients):
            result.append((key, coefficients))
    return result


DECISION_CACHE = {}


def routing_choices(per_routing, ratio, cache_key=None):
    """({t: routing index}, {t: mark}) for the routing-dependent products of a
    list of per-routing values, in the region of r = kappa_sigma/A containing
    `ratio` (an exact number)."""
    if cache_key is not None and (cache_key, ratio) in DECISION_CACHE:
        return DECISION_CACHE[(cache_key, ratio)]
    ctx = context(DECISION_ORDER)
    choices, marks = {}, {}
    for key, coefficients in routing_dependent_keys(per_routing):
        index, mark = largest_routing_at_large_epsilon(coefficients, ctx, branch_of(ratio), ratio)
        choices[key] = index
        if mark:
            marks[key] = mark
    if cache_key is not None:
        DECISION_CACHE[(cache_key, ratio)] = (choices, marks)
    return choices, marks


def decision_breakpoints(per_routing):
    """Positive values of r = kappa_sigma/A at which the comparison of the
    leading coefficients can change: roots of L_i = +-L_j and of L_i."""
    ctx = context(DECISION_ORDER)
    points = set()
    for key, coefficients in routing_dependent_keys(per_routing):
        for branch in BRANCHES:
            leads = []
            expansions = [expansion(ctx, c, branch, RATIO)[0] for c in coefficients]
            top = max((terms[0][0] for terms in expansions if terms), default=None)
            if top is None:
                continue
            leads = [terms[0][1] for terms in expansions if terms and terms[0][0] == top]
            equations = list(leads) + [one - sign * two for i, one in enumerate(leads)
                                       for two in leads[i + 1:] for sign in (1, -1)]
            for equation in equations:
                equation = sp.numer(sp.together(sp.expand(equation)))
                if equation.has(RATIO):
                    points |= {root for root in sp.solve(equation, RATIO) if root.is_real and root > 0}
    return points


def chosen_value(per_routing, choices):
    """sum_t C_t^{rho*_t} t, with rho*_t from `choices` (routing 0 for the
    routing-independent products)."""
    keys = set()
    for value in per_routing:
        keys |= set(value.terms)
    result = TraceSum()
    for key in keys:
        result.add_term(key, per_routing[choices.get(key, 0)].terms.get(key, 0))
    return result.simplified()


def regions(entries, extra_routing_lists=()):
    """Regions of r = kappa_sigma/A in which every routing choice at large eps is
    constant: open intervals between the break points (the break points
    themselves, degenerate mass relations, are not resolved).  r = 1 (the
    branch point of max(m_sigma, m_a)) is always a break point.  Each region
    carries a test value, the branch and the choices per graph name; the
    choices are checked at two more points of each interval."""
    points = {sp.Integer(1)}
    lists = [(entry['name'], entry['routing_values']) for entry in entries] + list(extra_routing_lists)
    for _, per_routing in lists:
        points |= decision_breakpoints(per_routing)
    points = sorted(points)
    edges = [sp.Integer(0)] + points + [sp.oo]
    result = []
    for lower, upper in zip(edges[:-1], edges[1:]):
        if upper == sp.oo:
            tests = [lower + 1, lower + 2, lower + 3]
        else:
            tests = [lower + (upper - lower) * fraction for fraction in (sp.Rational(1, 2), sp.Rational(1, 4),
                                                                        sp.Rational(3, 4))]
        tests = [sp.nsimplify(test) for test in tests]
        choices = {}
        for name, per_routing in lists:
            reference = routing_choices(per_routing, tests[0], name)
            for test in tests[1:]:
                assert routing_choices(per_routing, test, name)[0] == reference[0], (name, lower, upper)
            choices[name] = reference
        region = {'lower': lower, 'upper': upper, 'test': tests[0], 'branch': branch_of(tests[0]),
                  'choices': choices}
        last = result[-1] if result else None
        if (last is not None and last['branch'] == region['branch']
                and all(last['choices'][name][0] == choices[name][0] for name in choices)):
            last['upper'] = upper
        else:
            result.append(region)
    return result


def region_text(region):
    lower, upper = region['lower'], region['upper']
    parts = []
    if lower != 0:
        parts.append(f'{sp.latex(lower)} A<\\kappa_\\sigma' if lower != 1 else 'A<\\kappa_\\sigma')
    if upper != sp.oo:
        parts.append(f'\\kappa_\\sigma<{sp.latex(upper)} A' if upper != 1 else '\\kappa_\\sigma<A')
    return ',\\ '.join(parts) if parts else '\\text{all } \\kappa_\\sigma/A'


LEADING_CACHE = {}


def cached_leading(name, kind, value, branch, sources, choices):
    key = (name, kind, branch, sources, tuple(sorted((repr(k), v) for k, v in choices.items())))
    if key not in LEADING_CACHE:
        LEADING_CACHE[key] = leading_trace_sum(value, branch, sources)
    return LEADING_CACHE[key]


def graph_leading_table(entries, tree_lead, region):
    """Per kept graph and source set in one region: leading power and coefficient
    of T_Gamma (unsized, same routing choices), T_Gamma/M_0, the contribution
    X_Gamma^max/M_0 (1PI) and dZ^max (self-energies)."""
    rows = []
    branch = region['branch']
    for entry in entries:
        choices = region['choices'][entry['name']][0]
        value = chosen_value(entry['routing_values'], choices)
        unsized = chosen_value(entry['unsized routing values'], choices)
        for sources in SOURCE_SETS:
            row = {'name': entry['name'], 'group': entry['group'], 'label': entry['label'],
                   'sources': sources, 'r': entry['r'], 'marks': region['choices'][entry['name']][1]}
            row['T'] = cached_leading(entry['name'], 'T', unsized, branch, sources, choices)
            if entry['group'] == 'one_particle_irreducible':
                contribution = cached_leading(entry['name'], 'X', value, branch, sources, choices)
                row['contribution'] = contribution
                row['T/M0'] = None if row['T'] is None else row['T']['power'] - tree_lead['power']
                row['contribution/M0'] = (None if contribution is None
                                          else contribution['power'] - tree_lead['power'])
            else:
                row['dZ'] = cached_leading(entry['name'], 'dZ', value, branch, sources, choices)
            rows.append(row)
    return rows


def total_estimate(entries, tree, region, sources):
    """Leading power and coefficient of M_1^est and of R = M_1^est/M_0 in one
    region, with the graphs that attain the leading power and the leading
    coefficient of each of them ('by graph').  M_1^est = sum of the 1PI X^max
    + M_0 (dZ_sigma^max + dZ_pi+^max + dZ_pi-^max)/2, summed before the leading
    power is taken."""
    branch = region['branch']
    values = {entry['name']: (entry['group'], chosen_value(entry['routing_values'],
                                                            region['choices'][entry['name']][0]))
              for entry in entries}

    def pieces(ctx):
        result = {}
        tree_series = ctx.series(restrict_sources(tree, sources))
        for name, (group, value) in values.items():
            value = restrict_sources(value, sources)
            if group == 'one_particle_irreducible':
                series = ctx.series(value, branch)
            else:
                series = (tree_series * ctx.series(value, branch)).scale(sp.Rational(1, 2))
            if sources != 'generic':
                series = series.filtered(radial_index_at_most_one)
            result[name] = series
        return result

    for order in ORDERS:
        parts = pieces(context(order))
        total = TruncatedSeries()
        for series in parts.values():
            total.accumulate(series)
        try:
            power, coefficient = total.leading()
        except InsufficientOrder:
            continue
        by_graph = {}
        for name, series in parts.items():
            if series.precision < -power:
                raise InsufficientOrder(name)
            piece = series.coefficient(-power)
            if not piece.is_zero():
                by_graph[name] = piece
        tree_power, tree_coefficient = context(order).series(restrict_sources(tree, sources)).leading()
        ratio = sp.expand(anticommutator_expression(coefficient) / anticommutator_expression(tree_coefficient))
        return {'power': power, 'coefficient': coefficient, 'R power': power - tree_power,
                'R coefficient': ratio, 'attaining': list(by_graph), 'by graph': by_graph, 'order': order}
    raise InsufficientOrder('total estimate cancels up to the highest order tried')


def as_trace_expression(trace_sum):
    total = 0
    for key, value in trace_sum.terms.items():
        total += value * sp.Mul(*(trace_symbol(word) for word in key))
    return total


# ---------------------------------------------------------------------------
# 10. Item 7: survival of the four pieces of each meson line
# ---------------------------------------------------------------------------
RESIDUE_PIECES = ('(i) generator sum', '(ii) ALP-field entry', '(iii) mixed', '(iv) meson-block pole')


def four_piece_line(line_type, generator):
    """Meson line split into the four entries of (1_9 + 1) - v_a v_a^T:
    (i) 1_9 (scale m_P); (ii) (1 - 1/N) on the ALP field at both ends; (iii)
    -2 eps/N with M^{-1}(C) at one end and the ALP field at the other (both
    orientations); (iv) -4 eps^2/N with M^{-1}(C) at both ends; (ii)-(iv) have
    the scale m_a (D3).  Other lines as in line_pieces.  The fifth entry of
    each piece is its loop scale."""
    import derive_sigma_pipi_loop_graphs as graphs
    if line_type != 'M':
        return graphs.line_pieces(line_type, generator)
    n = NORMALIZATION

    def inverse_leg(momentum):
        return meson_leg('Minv', momentum)
    return [(1, lambda momentum: meson_leg(generator, momentum), generator, graphs.GENERATOR_SUM, M_MESON),
            (1 - 1 / n, bare_angular_leg, None, RESIDUE_PIECES[1], M_A),
            (-2 * EPSILON / n, (inverse_leg, bare_angular_leg), None, RESIDUE_PIECES[2], M_A),
            (-2 * EPSILON / n, (bare_angular_leg, inverse_leg), None, RESIDUE_PIECES[2], M_A),
            (-4 * EPSILON ** 2 / n, inverse_leg, None, RESIDUE_PIECES[3], M_A)]


def piece_name(piece):
    import derive_sigma_pipi_loop_graphs as graphs
    return RESIDUE_PIECES[0] if piece[3] == graphs.GENERATOR_SUM else piece[3]


def exactly_nonzero(per_routing):
    """Nonzero in some routing after the exact substitution N = 1 + 2 eps^2 <(M^{-1})^2>."""
    return not all(algebra.is_exactly_zero(value) for value in per_routing)


SURVIVAL_RATIO = sp.Rational(1, 2)      # region kappa_sigma < A used for the powers of the survival table


def maximum_power(per_routing, name):
    """Leading large-eps power of the termwise routing maximum (region kappa_sigma < A)."""
    choices, _ = routing_choices(per_routing, SURVIVAL_RATIO, name)
    return leading_trace_sum(chosen_value(per_routing, choices), branch_of(SURVIVAL_RATIO))['power']


def survival_table(entries, with_powers=True):
    """For every kept graph with a meson line and every meson line in it: which of
    the four pieces gives a nonzero contribution (exact, in some routing; D4 term
    removed as in the estimate; the measure is the size-rule value: X for 1PI
    graphs and dZ for self-energies) and the leading large-eps power of its
    termwise routing maximum (generic sources, region kappa_sigma < A)."""
    import derive_sigma_pipi_loop_graphs as graphs
    setups = graph_setups()
    rows = []
    for entry in entries:
        if entry['group'] == 'pi-' or 'M' not in entry['lines']:
            continue
        legs, basis_dot = setups[entry['group']]

        def piece_values(choice_filter, piece_function=four_piece_line):
            return graphs.routing_values(entry, legs, basis_dot, entry['group'], pieces=entry['pieces'],
                                         piece_function=piece_function, choice_filter=choice_filter)
        total = None
        lines = []
        for index, line in enumerate(entry['lines']):
            if line != 'M':
                continue
            pieces = {}
            for name in RESIDUE_PIECES:
                values = piece_values(lambda choice, i=index, n=name: piece_name(choice[i]) == n)
                nonzero = exactly_nonzero(values)
                result = {'nonzero': nonzero, 'values': values}
                if with_powers and nonzero:
                    result['power'] = maximum_power(values, (entry['name'], index, name))
                pieces[name] = result
                if not lines:
                    total = values if total is None else [one + two for one, two in zip(total, values)]
            survivors = [name for name in RESIDUE_PIECES if pieces[name]['nonzero']]
            separate = {}
            if with_powers:
                for name in (graphs.ALP_FIELD, graphs.ALP_LIKE_POLE):
                    values = piece_values(lambda choice, i=index, n=name: choice[i][3] == n, graphs.line_pieces)
                    separate[name] = (maximum_power(values, (entry['name'], index, 'separate', name))
                                      if exactly_nonzero(values) else None)
            lines.append({'line': index, 'pieces': pieces, 'survivors': survivors,
                          'only ALP-field entry': survivors == [RESIDUE_PIECES[1]],
                          'separate pieces': separate})
        meson_lines = [line['line'] for line in lines]
        combinations = {}
        for combination in itertools.product(RESIDUE_PIECES, repeat=len(meson_lines)):
            if len(meson_lines) == 1:
                values = lines[0]['pieces'][combination[0]]['values']
            else:
                def matches(choice, combination=combination):
                    return all(piece_name(choice[i]) == name for i, name in zip(meson_lines, combination))
                values = piece_values(matches)
            if exactly_nonzero(values):
                combinations[combination] = (maximum_power(values, (entry['name'], 'combination', combination))
                                             if with_powers else None)
        exact = all((one - two).is_zero() for one, two in zip(total, entry['routing_values']))
        rows.append({'name': entry['name'], 'label': entry['label'], 'lines': lines,
                     'combinations': combinations, 'decomposition exact': exact,
                     'D4 term removed': entry['pieces'] == 'rest'})
    return rows


def excluded_terms(entries, tree_lead):
    """Large-eps behaviour of the D4 terms (all meson lines on the generator
    sum, scale m_P) removed from the estimate: restated D4 criterion, power of
    the measure (X for 1PI graphs, dZ for self-energies; termwise routing
    maximum) and of its ratio to M_0."""
    import derive_sigma_pipi_loop_graphs as graphs
    result = {}
    for entry in entries:
        if entry['excluded'] is None or entry['group'] == 'pi-':
            continue
        values = entry['excluded_routings']
        choices, _ = routing_choices(values, SURVIVAL_RATIO, (entry['name'], 'D4'))
        measure = leading_trace_sum(chosen_value(values, choices), branch_of(SURVIVAL_RATIO))
        relative = measure['power'] - tree_lead['power'] if entry['group'] == 'one_particle_irreducible' \
            else measure['power']
        letters = set()
        for value in values:
            letters |= value.letters()
        result[entry['name']] = {'measure': measure, 'relative': relative,
                                 'criterion': all(graphs.satisfies_restated_d4(value, entry['group'])
                                                  for value in values),
                                 'letters': sorted(letters),
                                 'scales': sorted({str(atom) for value in values for c in value.terms.values()
                                                   for atom in sp.sympify(c).free_symbols
                                                   if atom in (M_A, M_SIGMA, M_MESON)})}
    return result


# ---------------------------------------------------------------------------
# 11. Markdown output
# ---------------------------------------------------------------------------

def coefficient_text(value):
    return sp.latex(sp.factor(value), symbol_names=LATEX_NAMES_LARGE_EPSILON)


def anticommutator_terms(trace_sum):
    """[(coefficient, key, paired)] where two products that differ only by
    lambda_{pi+} lambda_{pi-} -> lambda_{pi-} lambda_{pi+} in one trace and have
    equal coefficients are merged (paired = True) into <W {lambda_{pi+}, lambda_{pi-}}>;
    for paired terms the word in `key` is W followed by the marker letter '{+-}'."""
    remaining = dict(trace_sum.simplified().terms)
    result = []
    for key in sorted(remaining, key=repr):
        if key not in remaining:
            continue
        value = remaining.pop(key)
        paired = False
        for index, word in enumerate(key):
            if word.count('l+') != 1 or word.count('l-') != 1:
                continue
            position = word.index('l+')
            rotated = word[position:] + word[:position]
            if rotated[1] != 'l-':
                continue
            rest = rotated[2:]
            others = key[:index] + key[index + 1:]
            partner = tuple(sorted(others + (algebra.canonical_word(('l-', 'l+') + rest),)))
            if partner in remaining and sp.expand(remaining[partner] - value) == 0:
                remaining.pop(partner)
                result.append((value, others + (rest + ('{+-}',),), True))
                paired = True
                break
        if not paired:
            result.append((value, key, False))
    return result


def word_latex(word):
    names = []
    for name in word:
        if name == '{+-}':
            names.append(r'\{\lambda_{\pi^+},\lambda_{\pi^-}\}')
        else:
            names.append(algebra.letter_latex(name))
    return r'\,'.join(names)


def trace_sum_latex(trace_sum, kind='scalar', limit=12):
    """LaTeX of a leading coefficient.  Matrices are shown by removing the probe
    letters Y (matrix) and Y1, Y2 (tensor product); pairs of traces with
    lambda_{pi+} lambda_{pi-} and lambda_{pi-} lambda_{pi+} and equal coefficients
    are written with the anticommutator."""
    parts = []
    items = anticommutator_terms(trace_sum)
    for value, key, _ in items[:limit]:
        factors = []
        matrices = []
        for word in key:
            probes = [name for name in word if name in ('Y', 'Y1', 'Y2')]
            if probes:
                position = word.index(probes[0])
                rest = word[position + 1:] + word[:position]
                matrices.append((probes[0], word_latex(rest) or r'\mathbf 1'))
            else:
                factors.append(r'\langle ' + word_latex(word) + r'\rangle')
        matrices.sort()
        matrix_text = r'\otimes '.join(text for _, text in matrices)
        parts.append(r'\left(' + coefficient_text(value) + r'\right)' + ''.join(factors) + matrix_text)
    text = '+'.join(parts) if parts else '0'
    if len(items) > limit:
        text += rf'+\dots\ ({len(items)}\ \text{{terms}})'
    return '$' + text + '$'


def anticommutator_expression(trace_sum):
    """sympy expression with one symbol per (possibly anticommutator) trace."""
    total = 0
    for value, key, _ in anticommutator_terms(trace_sum):
        total += value * sp.Mul(*(sp.Symbol(r'\langle ' + word_latex(word) + r'\rangle') for word in key))
    return total


def up_only_ratio(total):
    """Leading coefficient of R for the up-only sources S_rs + i P_rs = m_u (-i)^s P_u/s!
    (r <= 1), lambda_{pi+} lambda_{pi-} = 2 P_u, lambda_{pi-} lambda_{pi+} = 2 P_d,
    A = kappa_sigma = B_0 m_u."""
    import derive_sigma_pipi_loop_vertices as vertices
    import derive_sigma_pipi_loop_tree as tree_script
    tree = tree_script.tree_amplitude()

    def value(trace_sum):
        reduced = tree_script.up_only_trace_value(vertices.up_only_substitution(trace_sum))
        return sp.expand(sp.sympify(reduced).subs({ALP_MASS_COEFFICIENT: B0 * UP_QUARK_MASS,
                                                   KAPPA_SIGMA: B0 * UP_QUARK_MASS,
                                                   sp.Symbol('m_u', positive=True): UP_QUARK_MASS}))
    tree_value = value(tree.map_coefficients(lambda c: sp.sympify(c).subs(EPSILON, 1)))
    return sp.factor(value(total['coefficient']) / tree_value)


def power_text(result):
    if result is None:
        return '0 (exact)'
    return f'$\\epsilon^{{{result["power"]}}}$'


RECORDS = {
    'scaling': 'murayama-sigma-pipi-loop-large-epsilon-scaling.md',
    'inputs': 'murayama-sigma-pipi-loop-large-epsilon-inputs.md',
    'vertices': 'murayama-sigma-pipi-loop-large-epsilon-vertices-tree.md',
    'graphs': 'murayama-sigma-pipi-loop-large-epsilon-graphs.md',
    'pieces': 'murayama-sigma-pipi-loop-meson-line-pieces.md',
}


def record_marker(record, section):
    """Comment line naming the memory record and section of the following block."""
    return f'\n<!-- {RECORDS[record]}: {section} -->\n'


def pole_data_markdown():
    """Sections 1-4 (exact checks and pole data at large eps) as markdown."""
    lines = [record_marker('scaling', 'Steps 1-4'), '### Steps 1-2. Quadratic Lagrangian and elimination of the mesons\n']
    for name, value in quadratic_form_checks().items():
        lines.append(f'- {name}: {value}')
    lines.append(f'- stationarity solved by $\\Phi=2\\epsilon a\\,\\mathbb M_s^{{-1}}(\\mathcal C)$: {stationarity_check()}')
    _, pole_function = alp_inverse_propagator()
    trace_names = {**LATEX_NAMES_LARGE_EPSILON, TRACE_S02: r'\langle\widehat S_{0,2}\rangle',
                   TRACE_C_MINV: r'\langle\mathcal C\,\mathbb M_s^{-1}(\mathcal C)\rangle'}
    lines.append(f'- $D_a(s)={sp.latex(pole_function, symbol_names=trace_names)}$')
    lines.append(f'- norm of $(2\\epsilon\\,\\mathbb M_s^{{-1}}(\\mathcal C),1)$: '
                 f'${sp.latex(pole_vector_norm(), symbol_names=LATEX_NAMES_LARGE_EPSILON)}$ with '
                 r'$q=\langle(\mathbb M_s^{-1}(\mathcal C))^2\rangle$')
    lines.append('\n### Step 3. Exact identities of the meson-line residue (any $\\epsilon$)\n')
    for name, value in residue_identities().items():
        lines.append(f'- {name}: {value}')
    lines.append('\n### Step 4. Pole data at large $\\epsilon$ (fixed sources, ALP mass entry $\\epsilon^2A$)\n')
    lines.append(f'- $D_a(s)=s-\\epsilon^2F(s)$, $F(s)=A-2\\langle\\mathcal C\\,\\mathbb M_s^{{-1}}(\\mathcal C)\\rangle$: '
                 f'{pole_function_with_alp_mass_entry()}')
    lines.append(f'- resolvent series identities $K=0,\\dots,4$: {all(resolvent_series_check(k) for k in range(5))}')
    result = alp_pole_series(order=2)
    names = {**LATEX_NAMES_LARGE_EPSILON, INVERSE_EPSILON_SQUARED: 't'}
    lines.append(f'- $m_a^2={sp.latex(result["root"], symbol_names=names)}+O(t^2)$, $t=1/\\epsilon^2$, '
                 r'$c_k=\langle\mathcal C\,L^k(\mathcal C)\rangle$')
    normalization = result['N from <(M^{-1})^2>']
    agree = sp.expand(normalization - result["N from D_a'"]) == 0
    lines.append(f'- $N(m_a^2)={sp.latex(normalization, symbol_names=names)}+O(t^2)$ '
                 f'(equal to $D_a\'(m_a^2)$: {agree})')
    lines.append(f'- $\\mathbb M^{{-1}}_{{m_a^2}}(\\mathcal C)=\\big({sp.latex(result["coefficient of C in M^{-1}(C) (leading)"], symbol_names=names)}\\big)\\,\\mathcal C+\\dots$; '
                 f'$\\xi_a={sp.latex(result["xi_a"], symbol_names=names)}+O(t^2)$')
    up = up_only_check()
    lines.append(f'- up-only sources: $A={sp.latex(up["A"], symbol_names=names)}$, $c_0={sp.latex(up["c_0"], symbol_names=names)}$, '
                 f'$c_1={sp.latex(up["c_1"], symbol_names=names)}$, $m_a^2={sp.latex(up["root"], symbol_names=names)}+O(t^2)$, '
                 f'$N={sp.latex(up["N"], symbol_names=names)}+O(t^2)$')
    return '\n'.join(lines)


def markdown(results):
    """Markdown blocks of sections 1-10, each preceded by a comment that names the
    memory record and section it belongs to."""
    lines = [pole_data_markdown()]
    lines.append(record_marker('inputs', 'Result: pole data and inputs; Result: pieces of the meson-line residue'))
    lines.append('### Pole and input data at large $\\epsilon$\n')
    lines.append('| quantity | leading power | leading coefficient |\n|---|---|---|')
    for name, row in results['pole'].items():
        lines.append(f'| {name} | $\\epsilon^{{{row["power"]}}}$ | {trace_sum_latex(row["coefficient"], row["kind"])} |')
    up = results['pole']['m_a^2 - m_sigma^2 (up-only)']
    lines.append(f'\nUp-only $m_a^2-m_\\sigma^2$ with $u=1/\\epsilon$: ${sp.latex(up["up-only series"], symbol_names=LATEX_NAMES_LARGE_EPSILON)}+O(u^{{{up["precision"] + 1}}})$')
    lines.append(record_marker('inputs', 'Result: leg factors with their placement'))
    lines.append('### Leg factors at large $\\epsilon$\n')
    lines.append('| leg factor | ALP-field ends | leading power | leading coefficient | net power by placement |\n|---|---|---|---|---|')
    for row in results['legs']:
        net = '; '.join(f'{k}: $\\epsilon^{{{v}}}$' for k, v in row['net'].items())
        lines.append(f'| {row["name"]} | {row["alp_field_ends"]} | $\\epsilon^{{{row["power"]}}}$ | '
                     f'{trace_sum_latex(row["coefficient"], row["kind"])} | {net} |')
    lines.append(record_marker('vertices', 'Vertices by field-end pattern'))
    lines.append('### Vertices by field-end pattern ($\\lambda$ fixed)\n')
    lines.append('| field ends | generic, fixed momenta | generic, loop momenta | r<=1, fixed momenta | r<=1, loop momenta | leading (generic, loop momenta) |')
    lines.append('|---|---|---|---|---|---|')
    for kinds, row in results['vertices'].items():
        cells = [str(row[(s, m)]) if row[(s, m)] else 'vanishes'
                 for s in SOURCE_SETS for m in ('fixed momenta', 'loop momenta')]
        lead = max(row[('generic', 'loop momenta')]) if row[('generic', 'loop momenta')] else '-'
        lines.append(f'| {" ".join(kinds)} | ' + ' | '.join(cells) + f' | {lead} |')
    tree = results['tree']
    lines.append(record_marker('vertices', 'Tree'))
    lines.append(f'### Tree\n\n$M_0$: {power_text(tree)}, coefficient {trace_sum_latex(tree["coefficient"])}')
    lines.append(record_marker('graphs', 'Result: 1PI graphs; Result: self-energies'))
    lines.append('### Graphs at large $\\epsilon$\n')
    lines.append('Regions of $r=\\kappa_\\sigma/A$ in which every routing choice of D12 is constant: '
                 + '; '.join(f'region {k + 1}: ${region_text(region)}$' for k, region in enumerate(results['regions']))
                 + '.\n')
    lines.append('Powers of $\\epsilon$ in region 1; the last column states whether the other regions give the same '
                 'powers. "contribution" is $X_\\Gamma^{\\max}/M_0$ for 1PI graphs and $\\delta Z^{\\max}$ for '
                 'self-energies.\n')
    lines.append('| graph | lines | $T_\\Gamma$ | $T_\\Gamma/M_0$ | contribution | $T_\\Gamma$ ($r\\le1$) | '
                 'contribution ($r\\le1$) | regions agree |')
    lines.append('|---|---|---|---|---|---|---|---|')
    tables = [{(row['name'], row['sources']): row for row in rows} for rows in results['graphs']]

    def contribution(row):
        if row['group'] == 'one_particle_irreducible':
            return '-' if row['contribution/M0'] is None else f'$\\epsilon^{{{row["contribution/M0"]}}}$'
        return power_text(row['dZ'])

    def powers(row):
        return (power_text(row['T']), contribution(row))
    names = []
    for row in results['graphs'][0]:
        if row['name'] not in names and row['group'] != 'pi-':
            names.append(row['name'])
    for name in names:
        generic, restricted = tables[0][(name, 'generic')], tables[0][(name, 'r <= 1')]
        same = all(powers(table[(name, sources)]) == powers(tables[0][(name, sources)])
                   for table in tables for sources in SOURCE_SETS)
        ratio = '-' if generic.get('T/M0') is None else f'$\\epsilon^{{{generic["T/M0"]}}}$'
        lines.append(f'| {name} | {generic["label"].split()[-1]} | {power_text(generic["T"])} | {ratio} | '
                     f'{contribution(generic)} | {power_text(restricted["T"])} | {contribution(restricted)} | '
                     f'{"yes" if same else "no"} |')
    lines.append('\nLeading coefficients (contribution $X_\\Gamma^{\\max}$ for 1PI graphs, $\\delta Z^{\\max}$ for '
                 'self-energies; generic sources; per region where they differ):\n')
    for name in names:
        texts = []
        for k, table in enumerate(tables):
            row = table[(name, 'generic')]
            measure = row['contribution'] if row['group'] == 'one_particle_irreducible' else row['dZ']
            text = 'exactly zero' if measure is None else (f'$\\epsilon^{{{measure["power"]}}}$ times '
                                                          f'{trace_sum_latex(measure["coefficient"], limit=40)}')
            texts.append(text)
        if len(set(texts)) == 1:
            lines.append(f'- {name}: {texts[0]}')
        else:
            lines.append(f'- {name}: ' + '; '.join(f'region {k + 1}: {text}' for k, text in enumerate(texts)))
        marks = {mark for table in tables for mark in table[(name, 'generic')]['marks'].values()}
        if marks:
            lines.append(f'  - routing choice marked: {", ".join(sorted(marks))}')
    lines.append(record_marker('graphs', 'Result: graphs that set the leading power; Result: M1_est and R'))
    lines.append('### Total $R=M_1^{\\rm est}/M_0$\n')
    for sources in SOURCE_SETS:
        for indices in grouped_regions(results['total'], sources):
            total = results['total'][indices[0]][sources]
            region_names = ', '.join(f'{k + 1} (${region_text(results["regions"][k])}$)' for k in indices)
            lines.append(f'- sources {sources}, regions {region_names}: $R\\sim\\epsilon^{{{total["R power"]}}}$, '
                         f'attained by {", ".join(total["attaining"])}.')
            lines.append(f'  - $M_1^{{\\rm est}}/M_0\\to\\epsilon^{{{total["R power"]}}}\\big('
                         f'{sp.latex(total["R coefficient"], symbol_names=LATEX_NAMES_LARGE_EPSILON)}\\big)$')
            if sources == 'r <= 1':
                lines.append(f'  - up-only sources: $M_1^{{\\rm est}}/M_0\\to\\epsilon^{{{total["R power"]}}}\\big('
                             f'{sp.latex(up_only_ratio(total), symbol_names=LATEX_NAMES_LARGE_EPSILON)}\\big)$, from '
                             f'{", ".join(up_only_contributors(total))}')
    lines.append(record_marker('graphs', 'Result: D4 terms (restated criterion)'))
    lines.append('### D4 terms (all meson lines on the generator sum, scale $m_P$, removed from the estimate)\n')
    for name, result in results['D4 terms'].items():
        lines.append(f'- {name}: restated criterion satisfied: {result["criterion"]}; letters '
                     f'{", ".join(result["letters"])}; masses {", ".join(result["scales"])}; '
                     f'{power_text(result["measure"])} for the measure; relative to $M_0$: '
                     f'$\\epsilon^{{{result["relative"]}}}$; coefficient {trace_sum_latex(result["measure"]["coefficient"])}')
    if 'survival' in results:
        lines.append(record_marker('pieces', 'Result'))
        lines.append('### Contribution of the four pieces of each meson line (exact, generic sources)\n')
        lines.append('Measure: $X_\\Gamma$ (1PI) or $\\delta Z$ (self-energies), nonzero in some routing, D4 term '
                     'removed. Line $i$ joins vertex $i$ and vertex $i+1$ in the label. In parentheses the leading '
                     'power at large $\\epsilon$ of the termwise routing maximum (region $\\kappa_\\sigma<A$). The '
                     'last column gives the power of the ALP-field entry $1$ alone and of the ALP-like pole '
                     '$-v_av_a^{\\rm T}$ alone.\n')
        lines.append('| graph | line | ' + ' | '.join(RESIDUE_PIECES) + ' | only (i), (ii) | D4 term removed | '
                     'entry 1 alone / $-v_av_a^{\\rm T}$ alone |')
        lines.append('|---|---|' + '---|' * len(RESIDUE_PIECES) + '---|---|---|')
        for row in results['survival']:
            for line in row['lines']:
                cells = []
                for name in RESIDUE_PIECES:
                    piece = line['pieces'][name]
                    cells.append(f'yes ($\\epsilon^{{{piece["power"]}}}$)' if piece['nonzero'] and 'power' in piece
                                 else ('yes' if piece['nonzero'] else 'no'))
                description = 'yes' if set(line['survivors']) <= set(RESIDUE_PIECES[:2]) else 'no'
                separate = ' / '.join('0' if power is None else f'$\\epsilon^{{{power}}}$'
                                      for power in line['separate pieces'].values())
                lines.append(f'| {row["name"]} | {line["line"]} | ' + ' | '.join(cells)
                             + f' | {description} | {"yes" if row["D4 term removed"] else "no"} | {separate} |')
        lines.append('\nNonzero combinations of pieces, one per meson line in the order of the lines (leading power):\n')
        for row in results['survival']:
            combinations = ', '.join('(' + ', '.join(name.split()[0] for name in combination)
                                     + f') $\\epsilon^{{{power}}}$' for combination, power in row['combinations'].items())
            lines.append(f'- {row["name"]} (lines {", ".join(str(line["line"]) for line in row["lines"])}): {combinations}')
    return '\n'.join(lines)


def same_total(one, two):
    return one['R power'] == two['R power'] and sp.expand(one['R coefficient'] - two['R coefficient']) == 0


def grouped_regions(totals, sources):
    """Lists of region indices with the same leading power and coefficient of R."""
    groups = []
    for index, total in enumerate(totals):
        for group in groups:
            if same_total(totals[group[0]][sources], total[sources]):
                group.append(index)
                break
        else:
            groups.append([index])
    return groups


def up_only_contributors(total):
    """Graphs whose leading term is nonzero for the up-only sources."""
    import derive_sigma_pipi_loop_vertices as vertices
    import derive_sigma_pipi_loop_tree as tree_script
    names = []
    for name, coefficient in total['by graph'].items():
        reduced = tree_script.up_only_trace_value(vertices.up_only_substitution(coefficient))
        value = sp.expand(sp.sympify(reduced).subs({ALP_MASS_COEFFICIENT: B0 * UP_QUARK_MASS,
                                                    KAPPA_SIGMA: B0 * UP_QUARK_MASS}))
        if value != 0:
            names.append(name)
    return names


def all_results(with_survival=True):
    timings = {}
    start = time.time()
    results = {'pole': pole_data_table()}
    timings['pole data'] = time.time() - start
    start = time.time()
    results['vertices'] = vertex_power_table()
    results['legs'] = leg_factor_table(results['pole'])
    tree, tree_lead = tree_leading()
    results['tree'] = tree_lead
    timings['vertices, legs, tree'] = time.time() - start
    start = time.time()
    entries, _ = graph_entries()
    results['entries'] = entries
    timings['graph values'] = time.time() - start
    start = time.time()
    results['regions'] = regions(entries)
    timings['routing decisions'] = time.time() - start
    start = time.time()
    results['graphs'] = [graph_leading_table(entries, tree_lead, region) for region in results['regions']]
    timings['graph series'] = time.time() - start
    start = time.time()
    results['total'] = [{sources: total_estimate(entries, tree, region, sources) for sources in SOURCE_SETS}
                        for region in results['regions']]
    timings['total'] = time.time() - start
    results['D4 terms'] = excluded_terms(entries, tree_lead)
    if with_survival:
        start = time.time()
        results['survival'] = survival_table(entries)
        timings['survival'] = time.time() - start
    results['timings'] = timings
    return results


def main():
    if '--markdown' in sys.argv:
        results = all_results()
        print(markdown(results))
        print('\nRuntime (s):', {key: round(value, 1) for key, value in results['timings'].items()})
        return
    print('1. Quadratic Lagrangian from the vertex rules (operator level):')
    for name, value in quadratic_form_checks().items():
        print(f'   {name}: {value}')
    print('\n2. Elimination of the mesons:')
    print('   stationarity solved by Phi = 2 eps a M_s^{-1}(C):', stationarity_check())
    _, pole_function = alp_inverse_propagator()
    print('   D_a(s) =', pole_function, '   [<C Minv> = <C M_s^{-1}(C)>]')
    print('   N(s) = D_a\'(s) with d<C M_s^{-1}(C)>/ds = <(M_s^{-1}(C))^2>; norm of (2 eps Minv, 1) =',
          pole_vector_norm())
    print('\n3. Exact identities of the meson-line residue (1_9 (+) 1) - v_a v_a^T, any eps '
          '(N = 1 + 2 eps^2 q, q = <(M^{-1}(C))^2>):')
    for name, value in residue_identities().items():
        print(f'   {name}: {value}')
    print('\n4. Large-eps behaviour (fixed sources, ALP mass entry eps^2 A, A > 0):')
    print('   D_a(s) = s - eps^2 F(s), F(s) = A - 2<C M_s^{-1}(C)>:', pole_function_with_alp_mass_entry())
    print('   resolvent series identities K = 0..4:', all(resolvent_series_check(k) for k in range(5)))
    result = alp_pole_series(order=2)
    print('   root m_a^2 =', result['root'], '  (t = 1/eps^2, c_k = <C L^k(C)>)')
    print('   N from <(M^{-1})^2> =', result['N from <(M^{-1})^2>'], '; N from D_a\' =', result['N from D_a\''])
    print('   M^{-1}_{m_a^2}(C) = (', result['coefficient of C in M^{-1}(C) (leading)'], ') C + ...')
    print('   xi_a =', result['xi_a'])
    print('   up-only sources:', up_only_check())


if __name__ == '__main__':
    main()
