"""Cartesian frame: assumptions, quadratic Lagrangian, pole data and tree amplitude of sigma -> pi+ pi- (task 1.1).

Run:
    PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_cartesian_loop_tree.py
    PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_cartesian_loop_tree

OpenSpec change estimate-sigma-pipi-one-loop-murayama-cartesian.  Units GeV;
conventions of sigma_pipi_cartesian_loop_algebra.py (Cartesian fields
phi = f_a + s + i p, eps = f_pi/f_a, all momenta incoming, d_mu -> -i k_mu,
vertex = n-th derivative of the Lagrangian).

Everything is derived from the Cartesian Lagrangian by the direct
noncommutative expansion (construction 1 of the module):
 1. one-point functions -> no-tadpole conditions P^C_00 = 0, <S^C_01> = 0 and
    the stationarity condition m_phi^2 + lambda f_a^2 = eps^2 B_0 <S^C_10>;
 2. two-point functions -> no s-p mixing (<S^C_11> = 0), no s-meson mixing
    (P^C_10 = 0), and, after stationarity,
      D_ss = k^2 - m_sigma^2,  m_sigma^2 = 2 lambda f_a^2 + eps^2 kappa_sigma,
      D_pp = k^2 - eps^2 A,
      D_pPhi(X) = eps <C X>,  C = B_0 P^C_01,
      D_PhiPhi(X, Y) = -(1/2) <X M_s(Y)>,  M_s(X) = B_0 {S^C_00, X} + (m_1^2/3) <X> 1 - s X,
    with A = B_0(<S^C_10> - 2<S^C_02>), kappa_sigma = B_0(<S^C_10> - 2<S^C_20>);
 3. elimination of the mesons -> D_a(s) = s - eps^2 A + 2 eps^2 <C M_s^{-1}(C)>,
    N(s) = D_a'(s) = 1 + 2 eps^2 <(M_s^{-1}(C))^2> = kinetic norm of (2 eps M_s^{-1}(C), 1);
 4. CP (s even, p and Phi odd, Phi -> -Phi^T) -> S^C_{r,q} = (-1)^q (S^C_{r,q})^T,
    P^C_{r,q} = (-1)^{q+1} (P^C_{r,q})^T;
 5. charge rule (generic proof reused) and Q-neutrality of M_s;
 6. tree M_0 = V(s; pi+, pi-) = -(B_0 eps/(2 f_pi)) <S^C_10 {lambda_{pi+}, lambda_{pi-}}>,
    and its up-only value -B_0 m_u eps/f_pi = -B_0 m_u/f_a.
No loop is evaluated.
"""
import sympy as sp

import sigma_pipi_cartesian_loop_algebra as cartesian
import derive_sigma_pipi_loop_tree as generic_tree_tools   # generic charge-rule proof and projector reduction
from sigma_pipi_loop_algebra import B0, EPSILON, F_PI, M_SIGMA, QUARTIC, TraceSum, CompactSum, letter
from sigma_pipi_cartesian_loop_algebra import (SINGLET_MASS_SQUARED, RADIAL_MASS_PARAMETER, SQUARED_MOMENTUM,
                                               UP_QUARK_MASS, MINV_LETTER)

F_A = F_PI / EPSILON
TRACE_SYMBOLS = {('SC10',): sp.Symbol('<SC10>', real=True), ('SC20',): sp.Symbol('<SC20>', real=True),
                 ('SC02',): sp.Symbol('<SC02>', real=True), ('SC01',): sp.Symbol('<SC01>', real=True),
                 ('SC11',): sp.Symbol('<SC11>', real=True)}


def as_expression(trace_sum):
    """A TraceSum whose keys are single traces of the letters in TRACE_SYMBOLS, as a sympy expression."""
    total = 0
    for key, value in trace_sum.simplified().terms.items():
        product = value
        for word in key:
            product *= TRACE_SYMBOLS[word]
        total += product
    return sp.expand(total)


def substitute_scalar(trace_sum, symbol, replacement):
    """Replace a scalar symbol entering linearly in the coefficients by a TraceSum."""
    result = TraceSum()
    for key, value in trace_sum.terms.items():
        value = sp.expand(value)
        linear = value.coeff(symbol, 1)
        constant = value.coeff(symbol, 0)
        if sp.expand(value - constant - linear * symbol) != 0:
            raise ValueError('symbol does not enter linearly')
        result = result + TraceSum({key: constant}) + TraceSum({key: 1}) * replacement.scale(linear)
    return result.simplified()


def two_point_momentum(trace_sum):
    """Momenta k and -k on the two legs: k_1.k_2 = -k^2."""
    dot = sp.Symbol('d_q1_q2', real=True)
    return trace_sum.map_coefficients(lambda value: sp.sympify(value).subs(dot, -SQUARED_MOMENTUM)).simplified()


# ---------------------------------------------------------------------------
# 1. One-point functions
# ---------------------------------------------------------------------------

def one_point_functions():
    """{'radial': ..., 'pseudoscalar': ...} from the Lagrangian with m_phi^2 kept."""
    return {'radial': cartesian.direct_vertex(1, 1, include_quadratic=True).simplified(),
            'pseudoscalar': cartesian.direct_vertex(1, 0, include_quadratic=True).simplified()}


def stationarity_solution():
    """m_phi^2 from the vanishing of the s tadpole."""
    expression = as_expression(one_point_functions()['radial'])
    solution, = sp.solve(expression, RADIAL_MASS_PARAMETER)
    return sp.expand(solution)


def stationarity_trace_sum():
    """m_phi^2 = eps^2 B_0 <S^C_10> - lambda f_a^2 as a TraceSum."""
    return TraceSum.single(('SC10',), EPSILON ** 2 * B0) + TraceSum.scalar(-QUARTIC * F_A ** 2)


# ---------------------------------------------------------------------------
# 2. Two-point functions
# ---------------------------------------------------------------------------

def two_point_functions():
    """Second derivatives with momenta k, -k, after stationarity:
    'radial-radial', 'radial-pseudoscalar' (leg 2 = (x2, xi2)), 'pseudoscalar-pseudoscalar'."""
    stationary = stationarity_trace_sum()
    result = {}
    for name, radial in (('radial-radial', 2), ('radial-pseudoscalar', 1), ('pseudoscalar-pseudoscalar', 0)):
        value = two_point_momentum(cartesian.direct_vertex(2, radial, include_quadratic=True))
        result[name] = substitute_scalar(value, RADIAL_MASS_PARAMETER, stationary)
    return result


def radial_mass_squared():
    """m_sigma^2 = 2 lambda f_a^2 + eps^2 B_0 (<S^C_10> - 2 <S^C_20>), as a TraceSum."""
    return (TraceSum.scalar(2 * QUARTIC * F_A ** 2)
            + TraceSum.single(('SC10',), EPSILON ** 2 * B0) + TraceSum.single(('SC20',), -2 * EPSILON ** 2 * B0))


def alp_mass_coefficient():
    """A = B_0 (<S^C_10> - 2 <S^C_02>) (the p mass entry is eps^2 A)."""
    return TraceSum.single(('SC10',), B0) + TraceSum.single(('SC02',), -2 * B0)


def kappa_sigma():
    """kappa_sigma = B_0 (<S^C_10> - 2 <S^C_20>)."""
    return TraceSum.single(('SC10',), B0) + TraceSum.single(('SC20',), -2 * B0)


def split_pseudoscalar_two_point():
    """Split D(x1 xi1, x2 xi2) = xi1 xi2 D_pp + xi1 D_pPhi(x2) + xi2 D_pPhi(x1) + D_PhiPhi(x1, x2)."""
    value = two_point_functions()['pseudoscalar-pseudoscalar']
    xi1, xi2 = sp.symbols('xi1 xi2', real=True)
    parts = {'pp': TraceSum(), 'p-meson 2': TraceSum(), 'p-meson 1': TraceSum(), 'meson-meson': TraceSum()}
    for key, coefficient in value.terms.items():
        polynomial = sp.Poly(sp.expand(coefficient), xi1, xi2)
        for (power_one, power_two), part in polynomial.terms():
            name = {(1, 1): 'pp', (1, 0): 'p-meson 2', (0, 1): 'p-meson 1', (0, 0): 'meson-meson'}[(power_one,
                                                                                                    power_two)]
            parts[name] = parts[name] + TraceSum({key: part})
    return {name: part.simplified() for name, part in parts.items()}


def mass_operator_trace(left, argument, coefficient=1):
    """<left M_s(argument)>, M_s(X) = B_0 {S^C_00, X} + (m_1^2/3) <X> 1 - s X, s = k^2."""
    left = tuple(left)
    unit = TraceSum()
    unit.add_term((left, (argument,)), SINGLET_MASS_SQUARED * coefficient / 3)
    return (TraceSum.single(left + ('SC00', argument), B0 * coefficient)
            + TraceSum.single(left + (argument, 'SC00'), B0 * coefficient)
            + unit + TraceSum.single(left + (argument,), -SQUARED_MOMENTUM * coefficient)).simplified()


def quadratic_lagrangian_checks():
    """Identities of the quadratic Cartesian Lagrangian after stationarity."""
    functions = two_point_functions()
    parts = split_pseudoscalar_two_point()
    checks = {}
    checks['s tadpole -> m_phi^2 = eps^2 B_0 <S^C_10> - lambda f_a^2'] = sp.expand(
        stationarity_solution() - as_expression(stationarity_trace_sum())) == 0
    checks['D_ss = k^2 - m_sigma^2, m_sigma^2 = 2 lambda f_a^2 + eps^2 kappa_sigma'] = (
        functions['radial-radial'] - TraceSum.scalar(SQUARED_MOMENTUM) + radial_mass_squared()).is_zero()
    checks['D_s(x, xi) = eps B_0 <P^C_10 x> + eps^2 B_0 xi <S^C_11>'] = (
        functions['radial-pseudoscalar'] - TraceSum.single(('PC10', 'x2'), EPSILON * B0)
        - TraceSum.single(('SC11',), EPSILON ** 2 * B0 * sp.Symbol('xi2', real=True))).is_zero()
    checks['D_pp = k^2 - eps^2 A'] = (
        parts['pp'] - TraceSum.scalar(SQUARED_MOMENTUM) + alp_mass_coefficient().scale(EPSILON ** 2)).is_zero()
    checks['D_pPhi(X) = eps <C X>, C = B_0 P^C_01'] = (
        (parts['p-meson 2'] - TraceSum.single(('PC01', 'x2'), EPSILON * B0)).is_zero()
        and (parts['p-meson 1'] - TraceSum.single(('PC01', 'x1'), EPSILON * B0)).is_zero())
    checks['D_PhiPhi(X, Y) = -(1/2) <X M_s(Y)>'] = (
        parts['meson-meson'] + mass_operator_trace(('x1',), 'x2').scale(sp.Rational(1, 2))).is_zero()
    checks['M_s self-adjoint: <X M_s(Y)> = <M_s(X) Y>'] = (
        mass_operator_trace(('X',), 'Y') - mass_operator_trace(('Y',), 'X')).is_zero()
    return checks


# ---------------------------------------------------------------------------
# 3. Elimination of the mesons: pole equation and N
# ---------------------------------------------------------------------------
# 'MinvC' = M_s^{-1}(C), defined by M_s(MinvC) = C = B_0 P^C_01 (assumption 5).
DERIVATIVE_LETTER = 'dMinvC'     # d/ds M_s^{-1}(C)


def defining_relation(left):
    """<left M_s(MinvC)> - B_0 <left P^C_01> (zero by definition)."""
    return mass_operator_trace(left, MINV_LETTER) - TraceSum.single(tuple(left) + ('PC01',), B0)


def derivative_relation(left):
    """d/ds of the defining relation: <left M_s(dMinvC)> - <left MinvC> (zero)."""
    return mass_operator_trace(left, DERIVATIVE_LETTER) - TraceSum.single(tuple(left) + (MINV_LETTER,))


def evaluate_two_point(two_point, first, second):
    """Evaluate D(x1 xi1, x2 xi2) on two vectors (meson coefficient, meson letter, ALP-field component)."""
    xi1, xi2 = sp.symbols('xi1 xi2', real=True)
    result = TraceSum()
    for key, value in two_point.terms.items():
        factor = sp.sympify(value).subs({xi1: first[2], xi2: second[2]})
        words = []
        for word in key:
            new_word = []
            for name in word:
                if name == 'x1':
                    factor *= first[0]
                    new_word.append(first[1])
                elif name == 'x2':
                    factor *= second[0]
                    new_word.append(second[1])
                else:
                    new_word.append(name)
            words.append(tuple(new_word))
        result.add_term(tuple(words), factor)
    return result.simplified()


def meson_stationarity():
    """dL_2/dPhi along X at Phi = 2 eps p MinvC: D_PhiPhi(X, Z) + D_pPhi(X) = -eps (defining relation)."""
    value = two_point_functions()['pseudoscalar-pseudoscalar']
    variation = evaluate_two_point(value, (1, 'X', 0), (2 * EPSILON, MINV_LETTER, 1))
    return (variation + defining_relation(('X',)).scale(EPSILON)).is_zero()


def alp_inverse_propagator():
    """D_a(s) = D(v, v) at v = (2 eps MinvC, 1), with <MinvC M_s(MinvC)> replaced by B_0 <MinvC P^C_01>."""
    value = two_point_functions()['pseudoscalar-pseudoscalar']
    vector = (2 * EPSILON, MINV_LETTER, 1)
    total = evaluate_two_point(value, vector, vector) + defining_relation((MINV_LETTER,)).scale(2 * EPSILON ** 2)
    return total.simplified()


def pole_data_checks():
    checks = {}
    checks['meson stationarity: Phi = 2 eps p M_s^{-1}(C)'] = meson_stationarity()
    expected = (TraceSum.scalar(SQUARED_MOMENTUM) - alp_mass_coefficient().scale(EPSILON ** 2)
                + TraceSum.single(('PC01', MINV_LETTER), 2 * EPSILON ** 2 * B0))
    checks['D_a(s) = s - eps^2 A + 2 eps^2 <C M_s^{-1}(C)>'] = (alp_inverse_propagator() - expected).is_zero()
    # N(s) = D_a'(s) = 1 + 2 eps^2 B_0 <P^C_01 dMinvC>, and B_0 <P^C_01 dMinvC> = <MinvC MinvC>:
    identity = (TraceSum.single(('PC01', DERIVATIVE_LETTER), B0) - TraceSum.single((MINV_LETTER, MINV_LETTER)))
    decomposition = (defining_relation((DERIVATIVE_LETTER,)).scale(-1)
                     + mass_operator_trace((DERIVATIVE_LETTER,), MINV_LETTER)
                     - mass_operator_trace((MINV_LETTER,), DERIVATIVE_LETTER)
                     + derivative_relation((MINV_LETTER,)))
    checks["N(s) = D_a'(s) = 1 + 2 eps^2 <(M_s^{-1}(C))^2>"] = (identity - decomposition).is_zero()
    # kinetic norm of (Z, 1): the coefficient of k^2 in D(v, v) before the defining relation
    vector = (2 * EPSILON, MINV_LETTER, 1)
    raw = evaluate_two_point(two_point_functions()['pseudoscalar-pseudoscalar'], vector, vector)
    norm = raw.map_coefficients(lambda value: sp.diff(value, SQUARED_MOMENTUM)).simplified()
    checks['kinetic norm of (2 eps M_s^{-1}(C), 1) = N(s)'] = (
        norm - TraceSum.scalar(1) - TraceSum.single((MINV_LETTER, MINV_LETTER), 2 * EPSILON ** 2)).is_zero()
    return checks


# ---------------------------------------------------------------------------
# 4. CP
# ---------------------------------------------------------------------------

def transpose_letter(name):
    return name[:-2] if name.endswith('^T') else name + '^T'


def cp_image(trace_sum):
    """Vertex of the CP-transformed Lagrangian L(s, -p, -Phi^T): legs (x_j, xi_j) -> (-x_j^T, -xi_j),
    then <W> = <W^T> (reverse the word, transpose each letter).  Source letters end up transposed."""
    result = TraceSum()
    for key, value in trace_sum.terms.items():
        factor = sp.sympify(value).subs({symbol: -symbol for symbol in sp.sympify(value).free_symbols
                                         if symbol.name.startswith('xi')}, simultaneous=True)
        words = []
        for word in key:
            meson_count = sum(1 for name in word if name.startswith('x'))
            factor *= (-1) ** meson_count
            transformed = [transpose_letter(name) if name.startswith('x') else name for name in word]
            words.append(tuple(transpose_letter(name) for name in reversed(transformed)))
        result.add_term(tuple(words), factor)
    return result.simplified()


def cp_sign(name):
    """eta in X = eta X^T: (-1)^q for S^C_{r,q}, (-1)^{q+1} for P^C_{r,q}."""
    kind, _, q = cartesian.source_indices(name)
    return (-1) ** q if kind == 'S' else (-1) ** (q + 1)


def cp_condition_substitution(trace_sum, sign=cp_sign):
    """The vertex with every source letter X replaced by eta X^T (eta = sign(X))."""
    result = TraceSum()
    for key, value in trace_sum.terms.items():
        factor = value
        words = []
        for word in key:
            new_word = []
            for name in word:
                if cartesian.is_cartesian_source(name):
                    factor *= sign(name)
                    new_word.append(name + '^T')
                else:
                    new_word.append(name)
            words.append(tuple(new_word))
        result.add_term(tuple(words), factor)
    return result.simplified()


def cp_checks(max_legs=4, sign=cp_sign):
    """For every pattern with n <= max_legs: the vertex of L(s, -p, -Phi^T) (direct expansion,
    all terms of L) equals the vertex with X -> eta X^T.  Since the meson matrices span all
    Hermitian matrices, CP invariance of L is then equivalent to X = eta X^T for every source letter."""
    results = {}
    for n in range(1, max_legs + 1):
        for radial in range(n + 1):
            value = cartesian.direct_vertex(n, radial, include_quadratic=True)
            results[(n, radial)] = (cp_image(value) - cp_condition_substitution(value, sign)).is_zero()
    return results


# ---------------------------------------------------------------------------
# 5. Charge rule
# ---------------------------------------------------------------------------

def mass_operator_commutes_with_charge():
    """[Q, M_s(X)] - M_s([Q, X]) with [Q, S^C_00] = 0: matrix part (words) and trace part."""
    pole_mass = sp.Symbol('s')
    rules = {('Q', 'SC00'): [(1, ('SC00', 'Q'))]}

    def apply_mass_operator(words):
        out = {}
        for word, value in words.items():
            for new_word, factor in ((('SC00',) + word, B0), (word + ('SC00',), B0), (word, -pole_mass)):
                out[new_word] = out.get(new_word, 0) + factor * value
        return out

    def commutator_with_charge(words):
        out = {}
        for word, value in words.items():
            out[('Q',) + word] = out.get(('Q',) + word, 0) + value
            out[word + ('Q',)] = out.get(word + ('Q',), 0) - value
        return out

    probe = {('X',): 1}
    left = commutator_with_charge(apply_mass_operator(probe))
    right = apply_mass_operator(commutator_with_charge(probe))
    difference = {word: left.get(word, 0) - right.get(word, 0) for word in set(left) | set(right)}
    difference = generic_tree_tools._rewrite({w: v for w, v in difference.items() if v != 0}, rules)
    trace_part = (TraceSum.single(('Q', 'X')) - TraceSum.single(('X', 'Q'))).scale(SINGLET_MASS_SQUARED / 3)
    return difference, trace_part.simplified().terms


# ---------------------------------------------------------------------------
# 6. Tree amplitude
# ---------------------------------------------------------------------------

def external_legs():
    """s (radial pole, no meson or p component: P^C_10 = 0, <S^C_11> = 0), pi+ = (lambda_{pi+}, 0),
    pi- = (lambda_{pi-}, 0); incoming momenta ps, p1, p2."""
    return {'sigma': cartesian.radial_leg({'ps': 1}),
            'pi+': cartesian.meson_leg('l+', {'p1': 1}),
            'pi-': cartesian.meson_leg('l-', {'p2': 1})}


def tree_amplitude():
    """M_0 = V(s; pi+, pi-) from the closed vertex rule, assumptions and charge rule applied."""
    legs = external_legs()
    value = cartesian.vertex([legs['sigma'], legs['pi+'], legs['pi-']], cartesian.symbolic_dot)
    return cartesian.apply_charge_rule(cartesian.apply_assumptions(value)).simplified()


def rename_legs(trace_sum, mapping, xi_values):
    """Rename leg letters and set the ALP-field components of a direct-expansion vertex."""
    result = TraceSum()
    for key, value in trace_sum.terms.items():
        factor = sp.sympify(value).subs(xi_values)
        result.add_term(tuple(tuple(mapping.get(name, name) for name in word) for word in key), factor)
    return result.simplified()


def tree_amplitude_direct(up_only=False):
    """M_0 from construction 1 (direct expansion of the Lagrangian)."""
    value = cartesian.direct_vertex(3, 1, up_only=up_only)
    xi2, xi3 = sp.symbols('xi2 xi3', real=True)
    value = rename_legs(value, {'x2': 'l+', 'x3': 'l-'}, {xi2: 0, xi3: 0})
    return cartesian.apply_charge_rule(cartesian.apply_assumptions(value)).simplified()


def tree_amplitude_compact():
    """M_0 = -(B_0 eps/(2 f_pi)) <S^C_10 {lambda_{pi+}, lambda_{pi-}}>."""
    result = CompactSum()
    result.add_term((('P', letter('SC10'), ('A', letter('l+'), letter('l-'))),), -B0 * EPSILON / (2 * F_PI))
    return result


TREE_LATEX = (r'M_0=V(s;\pi^+,\pi^-)=-\frac{B_0\epsilon}{2f_\pi}'
              r'\big\langle\widehat S^C_{1,0}\{\lambda_{\pi^+},\lambda_{\pi^-}\}\big\rangle')


def up_only_tree():
    """(M_0 after the D-A4 substitution, M_0 from the direct up-only expansion), both reduced with
    lambda_{pi+} lambda_{pi-} = 2 P_u, lambda_{pi-} lambda_{pi+} = 2 P_d, P_u^2 = P_u, P_d^2 = P_d,
    P_u P_d = 0, <P_u> = <P_d> = 1 (no explicit matrix)."""
    substituted = cartesian.up_only_substitution(tree_amplitude())
    direct = tree_amplitude_direct(up_only=True)
    return (generic_tree_tools.up_only_trace_value(substituted),
            generic_tree_tools.up_only_trace_value(direct))


def up_only_quadratic_data():
    """A, kappa_sigma and C in the up-only limit (<P_u> = 1)."""
    def value(trace_sum):
        reduced = cartesian.up_only_substitution(trace_sum)
        return sp.expand(sum(coefficient for key, coefficient in reduced.terms.items() if key == ()))
    return {'A': value(alp_mass_coefficient()), 'kappa_sigma': value(kappa_sigma()),
            'C / P_u': -B0 * UP_QUARK_MASS}


def tree_has_no_momentum():
    return all(not any(symbol.name.startswith('d_') for symbol in sp.sympify(value).free_symbols)
               for value in tree_amplitude().terms.values())


# ---------------------------------------------------------------------------
# Assumptions in the Cartesian invariants, linked to the polar assumption numbers
# ---------------------------------------------------------------------------
ASSUMPTIONS = [
    (1, 'No tree tadpoles: P^C_00 = 0, <S^C_01> = 0; s tadpole removed by m_phi^2 + lambda f_a^2 = eps^2 B_0 <S^C_10>',
     'derived: one_point_functions, stationarity_solution'),
    (2, 'No quadratic s-p mixing: <S^C_11> = 0; no s-meson mixing: P^C_10 = 0', 'derived: two_point_functions'),
    (3, 'Anomalous singlet term -m_1^2/12 <Phi>^2 independent of s and p', 'assumption (Lagrangian)'),
    (4, 'Leading U(3) chiral order; Cartesian frame with no derivative coupling of s, p to the mesons',
     'assumption (Lagrangian); with 3 it is the Cartesian form of Theta_hat_1 = 0'),
    (5, 'Regular simple poles; M_{m_a^2} invertible', 'assumption'),
    (6, 'p-meson quadratic entry only from the source: C = B_0 P^C_01', 'derived: D_pPhi(X) = eps <C X>'),
    (7, 'All S^C, P^C commute with Q', 'assumption; charge rule derived (generic proof)'),
    (8, 'f_pi independent of s and p: no s-meson or p-meson kinetic coupling', 'assumption; checked: tree has no momentum'),
    (9, 'One-loop tadpoles: pi0, eta, eta\', p vanish by CP; s tadpole set to zero (<Re phi> = f_a at one loop)',
     'assumption'),
    (10, 'CP: S^C_rq = (-1)^q S^C_rq^T, P^C_rq = (-1)^(q+1) P^C_rq^T; flavour-diagonal sources and quark masses',
     'derived: cp_checks; flavour-diagonality is an assumption'),
]


def main():
    print('1. One-point functions (m_phi^2 kept):')
    for name, value in one_point_functions().items():
        print(f'   {name}: {cartesian.trace_sum_latex(value)}')
    print(f'   stationarity: m_phi^2 = {sp.latex(stationarity_solution())}')
    print('\n2. Two-point functions after stationarity (momenta k, -k; s_k = k^2):')
    for name, value in two_point_functions().items():
        print(f'   {name}: {cartesian.trace_sum_latex(value)}')
    for name, ok in quadratic_lagrangian_checks().items():
        print(f'   {name}: {ok}')
    print('\n3. Pole data:')
    for name, ok in pole_data_checks().items():
        print(f'   {name}: {ok}')
    print('\n4. CP: vertex of L(s,-p,-Phi^T) equals the vertex with X -> eta X^T (n <= 4):',
          all(cp_checks().values()))
    derivation, trace = generic_tree_tools.charge_rule_proof()
    difference, trace_part = mass_operator_commutes_with_charge()
    print('\n5. Charge rule: [Q, W] - (sum q) W =', derivation or 0, '; <QW> - <WQ> =', trace or 0,
          '; [Q, M_s(X)] - M_s([Q, X]) =', difference or 0, trace_part or 0)
    tree = tree_amplitude()
    print('\n6. Tree from the vertex rule:', cartesian.trace_sum_latex(tree))
    print('   equal to direct expansion:', (tree - tree_amplitude_direct()).is_zero(),
          '; equal to compact form:', (tree_amplitude_compact().to_trace_sum() - tree).is_zero(),
          '; no momentum dependence:', tree_has_no_momentum())
    print('   ' + TREE_LATEX)
    substituted, direct = up_only_tree()
    print(f'   up-only: M_0 = {sp.latex(substituted)} (substitution), {sp.latex(direct)} (direct expansion)')
    print('   up-only quadratic data:', up_only_quadratic_data())
    print('\nAssumptions (polar number, Cartesian statement, status):')
    for number, text, status in ASSUMPTIONS:
        print(f'   {number}. {text}  [{status}]')


if __name__ == '__main__':
    main()
