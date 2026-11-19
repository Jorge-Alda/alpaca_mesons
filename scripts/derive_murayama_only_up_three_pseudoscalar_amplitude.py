#!/usr/bin/env python3
"""a_phys -> P_i P_j P_k in the model of DiLuzio:2026poh: cubic-vertex selection rule, radial-mode
(sigma) exchange channel by channel, the quartic vertex in invariant form (original fields and
original source coefficients) and the vanishing contact for three distinct neutral daughters.

Conventions: those of Alda:2026nrd and of scripts/derive_alp_three_pseudoscalar_radial_exchange.py
(polar coordinates, connection-free representative):
  * M(a, sigma) = S + i P = sum_{r,s} (S_rs + i P_rs) (sigma/f_a)^r (a/f_a)^s, Hermitian S_rs, P_rs.
  * Phi = sum_k P_k lambda_k, <lambda_k lambda_l> = 2 delta_kl, U = exp(i Phi/f_pi), chi = 2 B0 M,
    L = (f_pi^2/4) <dU dU^dagger + chi U^dagger + chi^dagger U> - (m_1^2/12) <Phi>^2
        + (1/2)(d sigma)^2 + (1/2)(1 + sigma/f_a)^2 (d a)^2 - V_phi.
  * M_s(X) = B0 {S_00, X} + (m_1^2/3) <X> 1 - s X;  N(m^2) = 1 + (2 B0^2 f_pi^2/f_a^2) <M^{-1}_{m^2}(P_01)^2>.
  * A pole x with an ALP component has, in the original fields, the meson matrix
    (2 B0 f_pi/f_a) M^{-1}_{m_x^2}(P_01)/sqrt(N(m_x^2)) and the ALP component 1/sqrt(N(m_x^2)); a
    charged or K^0-type pole has its meson generator and ALP component 0.
  * Vertices omit the factor i; all momenta outgoing (k_a = -p for the parent); the tree amplitude
    is M = V4 + sum_channels V3 (M^2 - s K)^{-1} V3 (S = 1 + i T).
  * PQ-charge models: M = diag[m_q (1 + sigma/f_a)^{n_q} exp(-i n_q a/f_a)]; only-up: n = (1, 0, 0).
Units: GeV for masses, fields, B0, f_pi, f_a and source coefficients; GeV^2 for s_uv, m^2 and k.k.

Assumptions: tree level, LO U(3) chiral order, no expansion in f_pi/f_a; P_00 = 0, <S_01> = 0,
P_10 = 0, <S_11> = 0; anomalous singlet term independent of sigma and a; stationary radial
potential; stable simple poles, no widths.

Content:
  (1) sympy: parity of the PQ-charge sources; cubic trace words at sigma = 0; the symmetric chiral
      map exp(X) exp(Y) exp(X) has no even terms (basis independence of on-shell vertices);
      quartic Lagrangian and the fully symmetric quartic vertex in original fields (brute-force
      fourth derivative); its form on pole legs; diagonal neutral and charged-pair forms; the
      neutral equations of motion, the secular product and the three-term identity behind the
      vanishing contact; closed form of the only-up sigma vertex.
  (2) numerics: cubic vertices of the exact Lagrangian (original fields and basis without
      ALP-meson mixing); direct amplitudes with complex (charged) legs; sigma exchange per channel;
      neutral-sector cosine model; quartic vertex of the original fields against the basis without
      ALP-meson mixing, on shell and off shell.

Run:
    .venv/bin/python scripts/derive_murayama_only_up_three_pseudoscalar_amplitude.py
Tests:
    .venv/bin/python -m unittest scripts.test_derive_murayama_only_up_three_pseudoscalar_amplitude
"""
import itertools
import math
import sys
from pathlib import Path

import numpy as np
import sympy as sp

_SCRIPTS = str(Path(__file__).resolve().parent)
if _SCRIPTS not in sys.path:
    sys.path.insert(0, _SCRIPTS)

from derive_alp_three_pseudoscalar_radial_exchange import (  # noqa: E402
    LAMBDA, NUMBER_OF_FIELDS, PAIRINGS, B_0, f_a, f_pi, dalitz_invariants, decay_gram, derivative_product,
    direct_amplitude, field_point, impose_stationarity, inverse_mass_operator, matrix_exponential,
    pole_normalization, radial_mass_squared_trace_form, radial_vertex_trace_form, random_parameters, source_matrices,
    spectrum, vertex)

MAXIMUM_SOURCE_ORDER = 4
LEG_LABELS = ('a', 'i', 'j', 'k')

# ---------------------------------------------------------------------------
# (1) Symbolic derivations
# ---------------------------------------------------------------------------


def derive_charge_source_parity(maximum_radial=2, maximum_angular=4):
    """Taylor coefficients of m (1 + rho)^n exp(-i n x) for symbolic real n.

    Returns (table {(r, s): (S_rs/m, P_rs/m)}, list of parity violations).  A violation is a
    nonzero S_rs with odd s or a nonzero P_rs with even s.
    """
    rho, x, charge = sp.symbols('rho x n', real=True)
    function = (1 + rho) ** charge * sp.exp(-sp.I * charge * x)
    table, violations = {}, []
    for r in range(maximum_radial + 1):
        for s in range(maximum_angular + 1):
            value = sp.diff(function, rho, r, x, s).subs({rho: 0, x: 0}) / (sp.factorial(r) * sp.factorial(s))
            scalar, pseudoscalar = sp.simplify(sp.re(value)), sp.simplify(sp.im(value))
            table[(r, s)] = (scalar, pseudoscalar)
            if (s % 2 == 1 and scalar != 0) or (s % 2 == 0 and pseudoscalar != 0):
                violations.append((r, s))
    return table, violations


def _words(expression):
    """Terms of an expanded noncommutative polynomial as (commutative coefficient, word)."""
    result = []
    for term in sp.Add.make_args(sp.expand(expression)):
        commutative, noncommutative = term.args_cnc()
        word = []
        for factor in noncommutative:
            if isinstance(factor, sp.Pow):
                word.extend([factor.base] * int(factor.exp))
            else:
                word.append(factor)
        result.append((sp.Mul(*commutative), tuple(word)))
    return result


def cyclic_trace_classes(expression):
    """Coefficients of the trace of a noncommutative polynomial, words identified up to rotation."""
    classes = {}
    for coefficient, word in _words(expression):
        names = tuple(str(symbol) for symbol in word)
        canonical = min(names[k:] + names[:k] for k in range(len(names))) if names else ()
        classes[canonical] = sp.expand(classes.get(canonical, 0) + coefficient)
    return {word: value for word, value in classes.items() if value != 0}


def _part_of_order(expression, counter, order):
    """Terms of a noncommutative polynomial with counter^order (counter commutative)."""
    return sp.Add(*[coefficient.coeff(counter, order) * sp.Mul(*word) for coefficient, word in _words(expression)])


def derive_cubic_trace_parity():
    """Cubic part of the polar Lagrangian at sigma = 0, as trace words.

    Noncommutative Hermitian symbols: Phi, D (= d_mu Phi, contracted with itself), S_n, P_n
    (Taylor coefficients of the source in x = a/f_a).  t counts fields (x, Phi and D).
    Returns (words of <chi U^dagger + chi^dagger U> with general sources,
             the same with S_1 = S_3 = P_0 = P_2 = 0, words of <dU dU^dagger>).
    """
    meson, derivative = sp.symbols('Phi D', commutative=False)
    scalar = sp.symbols('S0:4', commutative=False)
    pseudoscalar = sp.symbols('P0:4', commutative=False)
    count, angle = sp.symbols('t x', real=True)

    def exponential(sign):
        return sum((sign * sp.I * count * meson) ** n / sp.factorial(n) for n in range(4))

    def potential(scalar_list, pseudoscalar_list):
        chi = sum((count * angle) ** n * (scalar_list[n] + sp.I * pseudoscalar_list[n]) for n in range(4))
        chi_dagger = sum((count * angle) ** n * (scalar_list[n] - sp.I * pseudoscalar_list[n]) for n in range(4))
        return _part_of_order(chi * exponential(-1) + chi_dagger * exponential(1), count, 3)

    general = cyclic_trace_classes(potential(scalar, pseudoscalar))
    parity = cyclic_trace_classes(potential([scalar[0], 0, scalar[2], 0], [0, pseudoscalar[1], 0, pseudoscalar[3]]))

    def derivative_of_exponential(sign):
        total = 0
        for n in range(1, 4):
            total += (sign * sp.I * count) ** n / sp.factorial(n) * sum(
                meson ** k * derivative * meson ** (n - 1 - k) for k in range(n))
        return total

    kinetic = cyclic_trace_classes(_part_of_order(sp.expand(derivative_of_exponential(1) * derivative_of_exponential(-1)),
                                                  count, 3))
    return general, parity, kinetic


def derive_symmetric_chiral_map(order=4):
    """log(e^{tX} e^{tY} e^{tX}) by order in t, for noncommutative X, Y.

    With X = i a B0 M^{-1}_{m_a^2}(P_01)/f_a and Y = i Phi'/f_pi this is i Phi/f_pi of the chiral
    shift: the fields change as Phi = Phi' + a (2 B0 f_pi/f_a) M^{-1}_{m_a^2}(P_01) + (cubic) + ...  Returns
    {n: coefficient of t^n}; the even orders vanish, so the map has no quadratic part.
    """
    first, second = sp.symbols('X Y', commutative=False)
    t = sp.Symbol('t')

    def exponential(x):
        return sum((t * x) ** n / sp.factorial(n) for n in range(order + 1))

    product =sp.expand(exponential(first) * exponential(second) * exponential(first))
    deviation = sum(_part_of_order(product, t, n) * t ** n for n in range(1, order + 1))
    logarithm, power = 0, 1
    for n in range(1, order + 1):
        power = sp.expand(power * deviation)
        power = sum(_part_of_order(power, t, m) * t ** m for m in range(1, order + 1))
        logarithm += sp.Rational((-1) ** (n + 1), n) * power
    logarithm = sp.expand(logarithm)
    return {n: sp.expand(_part_of_order(logarithm, t, n)) for n in range(1, order + 1)}


def _quartic_lagrangian_polynomial():
    """Quartic part (no trace taken) of (f_pi^2/4)(dU dU^dagger + chi U^dagger + chi^dagger U) at sigma = 0.

    Noncommutative: Phi, D1 (d_mu Phi in dU), D2 (d^mu Phi in dU^dagger), S0..S4, P0..P4 (the
    original coefficients S_0n, P_0n).  Commutative: x = a/f_a, f_pi, B0.
    """
    meson, first_derivative, second_derivative = sp.symbols('Phi D1 D2', commutative=False)
    scalar = sp.symbols('S0:5', commutative=False)
    pseudoscalar = sp.symbols('P0:5', commutative=False)
    count, angle = sp.symbols('t x', real=True)
    order = 4

    def exponential(sign):
        return sum((sign * sp.I * count * meson / f_pi) ** n / sp.factorial(n) for n in range(order + 1))

    def derivative_of_exponential(sign, derivative):
        return sum((sign * sp.I * count / f_pi) ** n / sp.factorial(n)
                   * sum(meson ** k * derivative * meson ** (n - 1 - k) for k in range(n)) for n in range(1, order + 1))

    chi = 2 * B_0 * sum((count * angle) ** n * (scalar[n] + sp.I * pseudoscalar[n]) for n in range(order + 1))
    chi_dagger = 2 * B_0 * sum((count * angle) ** n * (scalar[n] - sp.I * pseudoscalar[n]) for n in range(order + 1))
    kinetic = _part_of_order(sp.expand(derivative_of_exponential(1, first_derivative)
                                       * derivative_of_exponential(-1, second_derivative)), count, order)
    potential = _part_of_order(sp.expand(chi * exponential(-1) + chi_dagger * exponential(1)), count, order)
    return sp.expand(f_pi ** 2 / 4 * (kinetic + potential))


def quartic_lagrangian_stated():
    """L4 = (1/(48 f_pi^2)) <[Phi, D1][Phi, D2]> + B0 f_pi^2 [<S0 Phi^4>/(24 f_pi^4) - x <P1 Phi^3>/(6 f_pi^3)
            - x^2 <S2 Phi^2>/(2 f_pi^2) + x^3 <P3 Phi>/f_pi + x^4 <S4>]  (x = a/f_a)."""
    meson, first_derivative, second_derivative = sp.symbols('Phi D1 D2', commutative=False)
    scalar = sp.symbols('S0:5', commutative=False)
    pseudoscalar = sp.symbols('P0:5', commutative=False)
    angle = sp.Symbol('x', real=True)
    commutator = lambda first, second: first * second - second * first
    return sp.expand(commutator(meson, first_derivative) * commutator(meson, second_derivative) / (48 * f_pi ** 2)
                     + B_0 * f_pi ** 2 * (scalar[0] * meson ** 4 / (24 * f_pi ** 4)
                                          - angle * pseudoscalar[1] * meson ** 3 / (6 * f_pi ** 3)
                                          - angle ** 2 * scalar[2] * meson ** 2 / (2 * f_pi ** 2)
                                          + angle ** 3 * pseudoscalar[3] * meson / f_pi
                                          + angle ** 4 * scalar[4]))


def derive_quartic_lagrangian():
    """Trace classes of (derived quartic Lagrangian - stated L4): empty if the statement holds."""
    return cyclic_trace_classes(_quartic_lagrangian_polynomial() - quartic_lagrangian_stated())


def _leg_symbols():
    matrices = sp.symbols('F_a F_i F_j F_k', commutative=False)
    alp = sp.symbols('alpha_a alpha_i alpha_j alpha_k', real=True)
    products = {}
    for x, y in itertools.combinations(range(4), 2):
        products[(x, y)] = products[(y, x)] = sp.Symbol(f'k{LEG_LABELS[x]}k{LEG_LABELS[y]}', real=True)
    return matrices, alp, products


def _commutator(first, second):
    return first * second - second * first


def quartic_vertex_stated(matrices, alp, products, scalar, pseudoscalar):
    """Fully symmetric four-leg vertex of the stated L4 on legs (F_x, alpha_x) and momenta k_x.

    V4 = B0 f_pi^2 [ (1/(24 f_pi^4)) sum_{S4} <S0 F F F F> - (1/(6 f_a f_pi^3)) sum_x alpha_x sum_{S3} <P1 F F F>
         - (1/(f_a^2 f_pi^2)) sum_{x<y} alpha_x alpha_y <S2 {F_u, F_v}> + (6/(f_a^3 f_pi)) sum_x prod_{y!=x} alpha_y <P3 F_x>
         + (24/f_a^4) prod alpha <S4> ]
         - (1/(24 f_pi^2)) sum_{x<y} (k_x.k_y) (<[F_u,F_x][F_v,F_y]> + <[F_v,F_x][F_u,F_y]>),
    {u, v} the complement of {x, y}.  Works with sympy symbols or numpy matrices (traces taken by
    the caller for sympy, inside for numpy).
    """
    legs = range(4)
    total = 0
    for permutation in itertools.permutations(legs):
        total += scalar[0] * _product([matrices[n] for n in permutation]) / (24 * f_pi ** 4)
    for x in legs:
        others = [y for y in legs if y != x]
        for permutation in itertools.permutations(others):
            total += -alp[x] * pseudoscalar[1] * _product([matrices[n] for n in permutation]) / (6 * f_a * f_pi ** 3)
    for x, y in itertools.combinations(legs, 2):
        u, v = [z for z in legs if z not in (x, y)]
        total += -alp[x] * alp[y] * scalar[2] * (matrices[u] * matrices[v] + matrices[v] * matrices[u]) / (
            f_a ** 2 * f_pi ** 2)
    for x in legs:
        total += 6 * sp.prod([alp[y] for y in legs if y != x]) * pseudoscalar[3] * matrices[x] / (f_a ** 3 * f_pi)
    total += 24 * sp.prod(alp) * scalar[4] / f_a ** 4
    total = B_0 * f_pi ** 2 * total
    for x, y in itertools.combinations(legs, 2):
        u, v = [z for z in legs if z not in (x, y)]
        total += -products[(x, y)] / (24 * f_pi ** 2) * (
            _commutator(matrices[u], matrices[x]) * _commutator(matrices[v], matrices[y])
            + _commutator(matrices[v], matrices[x]) * _commutator(matrices[u], matrices[y]))
    return sp.expand(total)


def _product(factors):
    """Ordered product of sympy (noncommutative) factors."""
    result = factors[0]
    for factor in factors[1:]:
        result = result * factor
    return result


def _matrix_product(factors):
    """Ordered matrix product of numpy arrays."""
    result = factors[0]
    for factor in factors[1:]:
        result = result @ factor
    return result


def derive_quartic_vertex():
    """Brute-force fourth derivative of the derived quartic Lagrangian against quartic_vertex_stated.

    Phi -> sum_x c_x F_x, x -> sum_x c_x alpha_x/f_a, D1 -> sum_x c_x (i k_x) F_x,
    D2 -> sum_x c_x (i k_x) F_x (plane waves e^{i k.x}, so D1 D2 -> -(k_x.k_y)); the coefficient of
    c_a c_i c_j c_k is the sum over all assignments of the four legs to the four field slots of
    each word.  Returns the trace classes of the difference (empty if the stated vertex is right).
    """
    meson, first_derivative, second_derivative = sp.symbols('Phi D1 D2', commutative=False)
    scalar = sp.symbols('S0:5', commutative=False)
    pseudoscalar = sp.symbols('P0:5', commutative=False)
    angle = sp.Symbol('x', real=True)
    matrices, alp, products = _leg_symbols()
    brute_force = 0
    for coefficient, word in _words(_quartic_lagrangian_polynomial()):
        angular_power = sp.degree(coefficient, angle)
        rest = coefficient.coeff(angle, angular_power) if angular_power else coefficient
        field_positions = [n for n, symbol in enumerate(word) if symbol in (meson, first_derivative, second_derivative)]
        assert angular_power + len(field_positions) == 4
        for permutation in itertools.permutations(range(4)):
            value = rest * sp.prod([alp[leg] / f_a for leg in permutation[:angular_power]])
            assigned = list(word)
            derivative_legs = {}
            for position, leg in zip(field_positions, permutation[angular_power:]):
                if word[position] != meson:
                    derivative_legs[word[position]] = leg
                assigned[position] = matrices[leg]
            if derivative_legs:
                value *= -products[(derivative_legs[first_derivative], derivative_legs[second_derivative])]
            brute_force += value * sp.Mul(*assigned)
    stated = quartic_vertex_stated(matrices, alp, products, scalar, pseudoscalar)
    return cyclic_trace_classes(sp.expand(brute_force) - stated)


def derive_kinetic_channel_identity():
    """Momentum part of the vertex.  Returns (trace classes of sum_r C_r, residual of
    k_a.k_r + k_s.k_t - (s_st - (1/2) sum_x k_x^2) for momentum conservation), where
    C_r = <[F_s,F_a][F_t,F_r]> + <[F_t,F_a][F_s,F_r]> for r in {i, j, k}, {s, t} the other daughters."""
    matrices, _, _ = _leg_symbols()
    total = 0
    for r, (s, t) in ((1, (2, 3)), (2, (1, 3)), (3, (1, 2))):
        total += (_commutator(matrices[s], matrices[0]) * _commutator(matrices[t], matrices[r])
                  + _commutator(matrices[t], matrices[0]) * _commutator(matrices[s], matrices[r]))
    momenta = sp.Matrix(3, 3, lambda row, column: sp.Symbol(f'k{row}{column}', real=True))  # daughters x 3 components
    metric = sp.diag(1, -1, -1)
    daughters = [momenta[n, :] for n in range(3)]
    parent = -(daughters[0] + daughters[1] + daughters[2])
    dot = lambda first, second: (first * metric * second.T)[0]
    legs = [parent] + daughters
    squares = sum(dot(leg, leg) for leg in legs)
    r, s, t = 1, 2, 3
    residual = sp.expand(dot(legs[0], legs[r]) + dot(legs[s], legs[t])
                         - (dot(legs[s] + legs[t], legs[s] + legs[t]) - squares / 2))
    return cyclic_trace_classes(sp.expand(total)), residual


def quartic_vertex_on_pole_legs(inverse_images, scalar, pseudoscalar, invariants):
    """Braces of V4 = (2 B0 f_pi^2/(3 f_a^4)) prod_x N(m_x^2)^{-1/2} {...} for four legs with an ALP
    component, R_x = M^{-1}_{m_x^2}(P_01):
      B0^4 sum_{S4} <S0 R R R R> - 2 B0^3 sum_x sum_{S3} <P1 R R R> - 6 B0^2 sum_{x<y} <S2 {R_u, R_v}>
      + 18 B0 sum_x <P3 R_x> + 36 <S4> - B0^3 sum_{(r|st)} s_st (<[R_s,R_a][R_t,R_r]> + <[R_t,R_a][R_s,R_r]>).
    invariants: {(s, t): s_st} for the daughter pairs (1, 2), (1, 3), (2, 3)."""
    legs = range(4)
    total = 0
    for permutation in itertools.permutations(legs):
        total += B_0 ** 4 * scalar[0] * _product([inverse_images[n] for n in permutation])
    for x in legs:
        for permutation in itertools.permutations([y for y in legs if y != x]):
            total += -2 * B_0 ** 3 * pseudoscalar[1] * _product([inverse_images[n] for n in permutation])
    for x, y in itertools.combinations(legs, 2):
        u, v = [z for z in legs if z not in (x, y)]
        total += -6 * B_0 ** 2 * scalar[2] * (inverse_images[u] * inverse_images[v] + inverse_images[v] * inverse_images[u])
    for x in legs:
        total += 18 * B_0 * pseudoscalar[3] * inverse_images[x]
    total += 36 * scalar[4]
    for r, (s, t) in ((1, (2, 3)), (2, (1, 3)), (3, (1, 2))):
        total += -B_0 ** 3 * invariants[(s, t)] * (
            _commutator(inverse_images[s], inverse_images[0]) * _commutator(inverse_images[t], inverse_images[r])
            + _commutator(inverse_images[t], inverse_images[0]) * _commutator(inverse_images[s], inverse_images[r]))
    return total


def derive_pole_leg_form():
    """quartic_vertex_stated with F_x = (2 B0 f_pi/f_a) R_x, alpha_x = 1 and on-shell momenta
    against (2 B0 f_pi^2/(3 f_a^4)) quartic_vertex_on_pole_legs.  Returns the difference classes."""
    inverse_images = sp.symbols('R_a R_i R_j R_k', commutative=False)
    scalar = sp.symbols('S0:5', commutative=False)
    pseudoscalar = sp.symbols('P0:5', commutative=False)
    masses = sp.symbols('m_a2 m_i2 m_j2 m_k2', real=True)
    invariants = {(1, 2): sp.Symbol('s_ij', real=True), (1, 3): sp.Symbol('s_ik', real=True),
                  (2, 3): sp.Symbol('s_jk', real=True)}
    channel = {1: (2, 3), 2: (1, 3), 3: (1, 2)}
    products = {}
    for r in (1, 2, 3):
        s, t = channel[r]
        parent_dot = -(masses[0] + masses[r] - invariants[(s, t)]) / 2  # k_a.k_r = -p.p_r
        products[(0, r)] = products[(r, 0)] = parent_dot
        products[(s, t)] = products[(t, s)] = (invariants[(s, t)] - masses[s] - masses[t]) / 2
    scale = 2 * B_0 * f_pi / f_a
    general = quartic_vertex_stated([scale * inverse_image for inverse_image in inverse_images], [1, 1, 1, 1], products,
                                    scalar, pseudoscalar)
    stated = 2 * B_0 * f_pi ** 2 / (3 * f_a ** 4) * quartic_vertex_on_pole_legs(inverse_images, scalar, pseudoscalar,
                                                                                invariants)
    return cyclic_trace_classes(sp.expand(general - stated))


def _diagonal_sources(masses, charges):
    """S_0n, P_0n of diag[m_q exp(-i n_q x)] as sympy matrices, n = 0..4."""
    entries = lambda n: [masses[q] * (-sp.I * charges[q]) ** n / sp.factorial(n) for q in range(3)]
    scalar = [sp.diag(*[sp.re(value) for value in entries(n)]) for n in range(5)]
    pseudoscalar = [sp.diag(*[sp.im(value) for value in entries(n)]) for n in range(5)]
    return scalar, pseudoscalar


def _matrix_vertex(matrices, alp, products, scalar, pseudoscalar):
    """quartic_vertex_stated evaluated with explicit 3x3 sympy matrices (trace taken)."""
    legs = range(4)
    total = 0
    for permutation in itertools.permutations(legs):
        total += (scalar[0] * _product([matrices[n] for n in permutation])).trace() / (24 * f_pi ** 4)
    for x in legs:
        for permutation in itertools.permutations([y for y in legs if y != x]):
            total += -alp[x] * (pseudoscalar[1] * _product([matrices[n] for n in permutation])).trace() / (
                6 * f_a * f_pi ** 3)
    for x, y in itertools.combinations(legs, 2):
        u, v = [z for z in legs if z not in (x, y)]
        total += -alp[x] * alp[y] * (scalar[2] * (matrices[u] * matrices[v] + matrices[v] * matrices[u])).trace() / (
            f_a ** 2 * f_pi ** 2)
    for x in legs:
        total += 6 * sp.prod([alp[y] for y in legs if y != x]) * (pseudoscalar[3] * matrices[x]).trace() / (
            f_a ** 3 * f_pi)
    total += 24 * sp.prod(alp) * scalar[4].trace() / f_a ** 4
    total = B_0 * f_pi ** 2 * total
    for x, y in itertools.combinations(legs, 2):
        u, v = [z for z in legs if z not in (x, y)]
        total += -products[(x, y)] / (24 * f_pi ** 2) * (
            _commutator(matrices[u], matrices[x]) * _commutator(matrices[v], matrices[y])
            + _commutator(matrices[v], matrices[x]) * _commutator(matrices[u], matrices[y])).trace()
    return sp.expand(total)


def charged_pair_form(masses, charges, parent_diagonal, daughter_diagonal, parent_alp, daughter_alp, invariants,
                      pair=(0, 1)):
    """V4 for a parent and a neutral daughter i (diagonal meson matrices, ALP components) and the
    pair j = sqrt2 |f><f'|, k = sqrt2 |f'><f| (pair = (f, f')), diagonal PQ-charge sources:

    B0 { (m_f + m_f')/(6 f_pi^2) (A_f + A_f')(I_f + I_f') + (m_f A_f I_f + m_f' A_f' I_f')/(3 f_pi^2)
         + (1/(3 f_a f_pi)) sum_{(q,q')} m_q n_q [alpha_i (2 A_q + A_q') + alpha_a (2 I_q + I_q')]
         + (alpha_a alpha_i/f_a^2) sum_q m_q n_q^2 }
    + (A_f - A_f')(I_f - I_f') (2 s_jk - s_ij - s_ik)/(12 f_pi^2),
    A = parent_diagonal, I = daughter_diagonal; (q, q') runs over (f, f'), (f', f).
    invariants = (s_ij, s_ik, s_jk)."""
    f, g = pair
    A, I = parent_diagonal, daughter_diagonal
    s_ij, s_ik, s_jk = invariants
    value = ((masses[f] + masses[g]) / (6 * f_pi ** 2) * (A[f] + A[g]) * (I[f] + I[g])
             + (masses[f] * A[f] * I[f] + masses[g] * A[g] * I[g]) / (3 * f_pi ** 2)
             + sum(masses[q] * charges[q] * (daughter_alp * (2 * A[q] + A[p]) + parent_alp * (2 * I[q] + I[p]))
                   for q, p in ((f, g), (g, f))) / (3 * f_a * f_pi)
             + parent_alp * daughter_alp / f_a ** 2 * sum(masses[q] * charges[q] ** 2 for q in (f, g)))
    return B_0 * value + (A[f] - A[g]) * (I[f] - I[g]) * (2 * s_jk - s_ij - s_ik) / (12 * f_pi ** 2)


def derive_diagonal_forms():
    """Generic vertex with diagonal PQ-charge sources (symbolic n_q).  Returns residuals of
    (a) four diagonal legs against B0 f_pi^2 sum_q m_q prod_x ((F_x)_qq/f_pi + n_q alpha_x/f_a);
    (b) charged pair (u, d) and (d, s) against charged_pair_form (momenta k from s_uv on shell)."""
    masses = sp.symbols('m_u m_d m_s', positive=True)
    charges = sp.symbols('n_u n_d n_s', real=True)
    scalar, pseudoscalar = _diagonal_sources(masses, charges)
    diagonals = {label: sp.symbols(f'F_{label}u F_{label}d F_{label}s', real=True) for label in LEG_LABELS}
    alp = sp.symbols('alpha_a alpha_i alpha_j alpha_k', real=True)
    _, _, products = _leg_symbols()
    neutral = _matrix_vertex([sp.diag(*diagonals[label]) for label in LEG_LABELS], alp, products, scalar,
                             pseudoscalar)
    expected = B_0 * f_pi ** 2 * sum(masses[q] * sp.prod([diagonals[label][q] / f_pi + charges[q] * alp[n] / f_a
                                                          for n, label in enumerate(LEG_LABELS)]) for q in range(3))
    residuals = [sp.simplify(neutral - sp.expand(expected))]
    unit = lambda row, column: sp.Matrix(3, 3, lambda x, y: 1 if (x, y) == (row, column) else 0)
    masses_squared = sp.symbols('m_a2 m_i2 m_j2 m_k2', real=True)
    s_ij, s_ik, s_jk = sp.symbols('s_ij s_ik s_jk', real=True)
    pair_products = {}
    for (x, y), value in {(0, 1): -(masses_squared[0] + masses_squared[1] - s_jk) / 2,
                          (0, 2): -(masses_squared[0] + masses_squared[2] - s_ik) / 2,
                          (0, 3): -(masses_squared[0] + masses_squared[3] - s_ij) / 2,
                          (2, 3): (s_jk - masses_squared[2] - masses_squared[3]) / 2,
                          (1, 3): (s_ik - masses_squared[1] - masses_squared[3]) / 2,
                          (1, 2): (s_ij - masses_squared[1] - masses_squared[2]) / 2}.items():
        pair_products[(x, y)] = pair_products[(y, x)] = value
    for pair in ((0, 1), (1, 2)):
        f, g = pair
        legs = [sp.diag(*diagonals['a']), sp.diag(*diagonals['i']), sp.sqrt(2) * unit(f, g), sp.sqrt(2) * unit(g, f)]
        direct = _matrix_vertex(legs, [alp[0], alp[1], 0, 0], pair_products, scalar, pseudoscalar)
        stated = charged_pair_form(masses, charges, diagonals['a'], diagonals['i'], alp[0], alp[1],
                                   (s_ij, s_ik, s_jk), pair)
        residuals.append(sp.simplify(direct - sp.expand(stated)))
    return residuals


def _neutral_sector_symbols():
    masses = sp.symbols('m_u m_d m_s', positive=True)
    up_charge = sp.Symbol('n_u', real=True)
    singlet = sp.Symbol('m_1^2', positive=True)
    return masses, (up_charge, 0, 0), singlet


def cosine_direction_masses_squared_symbolic(masses, charges):
    """2 B0 m_q (1 + n_q^2 f_pi^2/(2 f_a^2)): squared masses of the cosine directions without the anomaly."""
    return [2 * B_0 * masses[q] * (1 + charges[q] ** 2 * f_pi ** 2 / (2 * f_a ** 2)) for q in range(3)]


def derive_neutral_pole_relation():
    """Diagonal neutral sector at sigma = 0, original fields (Phi_uu, Phi_dd, Phi_ss, a), only u charged:
    L2 = (1/4) sum_q (d Phi_qq)^2 + (1/2)(d a)^2 - (B0 f_pi^2/2) sum_q m_q (Phi_qq/f_pi + n_q a/f_a)^2
         - (m_1^2/12)(sum_q Phi_qq)^2.
    With E_q = dL2/dPhi_qq and E_a = dL2/da at p^2 = z, checks for q = u, d, s that
    (2/f_pi) E_q + (n_q/f_a) E_a = (z - 2 B0 m_q (1 + n_q^2 f_pi^2/(2 f_a^2))) (Phi_qq/f_pi + n_q a/f_a)
                                   - (m_1^2/(3 f_pi)) sum_p Phi_pp,
    and that E_a = 0 gives Phi_uu/f_pi + n_u a/f_a = z f_a a/(B0 m_u n_u f_pi^2).  Returns residuals."""
    masses, charges, singlet = _neutral_sector_symbols()
    fields = sp.symbols('Phi_uu Phi_dd Phi_ss', real=True)
    angular, z = sp.symbols('a z', real=True)
    argument = [fields[q] / f_pi + charges[q] * angular / f_a for q in range(3)]
    lagrangian = (z / 4 * sum(field ** 2 for field in fields) + z / 2 * angular ** 2
                  - B_0 * f_pi ** 2 / 2 * sum(masses[q] * argument[q] ** 2 for q in range(3))
                  - singlet / 12 * sum(fields) ** 2)
    equations = [sp.diff(lagrangian, field) for field in fields]
    alp_equation = sp.diff(lagrangian, angular)
    nodes = cosine_direction_masses_squared_symbolic(masses, charges)
    residuals = [sp.simplify(2 / f_pi * equations[q] + charges[q] / f_a * alp_equation
                             - ((z - nodes[q]) * argument[q] - singlet / (3 * f_pi) * sum(fields))) for q in range(3)]
    solved = sp.solve(alp_equation, angular)[0]
    up_bracket = sp.simplify(argument[0].subs(angular, solved) - z * f_a * solved / (B_0 * masses[0] * charges[0] * f_pi ** 2))
    return residuals, up_bracket


def derive_neutral_secular_product():
    """prod over the four poles of (z - m_x^2) = det(z K - H)/det K for the sector of
    derive_neutral_pole_relation (K = diag(1/2, 1/2, 1/2, 1)).  Returns, for q = u, d, s, the residual
    of prod_x (M_q - m_x^2) = -(2 B0 m_1^2/3) m_q prod_{p != q} (M_q - M_p) at M_q = 2 B0 m_q (1 + n_q^2 f_pi^2/(2 f_a^2))."""
    masses, charges, singlet = _neutral_sector_symbols()
    fields = sp.symbols('Phi_uu Phi_dd Phi_ss', real=True)
    angular, z = sp.symbols('a z', real=True)
    coordinates = list(fields) + [angular]
    argument = [fields[q] / f_pi + charges[q] * angular / f_a for q in range(3)]
    potential = B_0 * f_pi ** 2 / 2 * sum(masses[q] * argument[q] ** 2 for q in range(3)) + singlet / 12 * sum(fields) ** 2
    hessian = sp.hessian(potential, coordinates)
    kinetic = sp.diag(sp.Rational(1, 2), sp.Rational(1, 2), sp.Rational(1, 2), 1)
    characteristic = sp.expand((z * kinetic - hessian).det() / kinetic.det())
    nodes = cosine_direction_masses_squared_symbolic(masses, charges)
    return [sp.simplify(characteristic.subs(z, nodes[q]) + 2 * B_0 * singlet / 3 * masses[q]
                        * sp.prod([nodes[q] - nodes[p] for p in range(3) if p != q])) for q in range(3)]


def derive_three_term_identity():
    """sum over cyclic (x, y, z) of 1/((x - y)(x - z)), and the contact sum
    sum_q m_q/prod_x(M_q - m_x^2) = -(3/(2 B0 m_1^2)) sum_q 1/prod_{p != q}(M_q - M_p) with the three
    explicit squared masses.  Returns both (zero)."""
    first, second, third = sp.symbols('x y z', real=True)
    identity = sp.simplify(1 / ((first - second) * (first - third)) + 1 / ((second - first) * (second - third))
                           + 1 / ((third - first) * (third - second)))
    masses, charges, singlet = _neutral_sector_symbols()
    nodes = cosine_direction_masses_squared_symbolic(masses, charges)
    contact = sp.simplify(sum(-3 / (2 * B_0 * singlet) / sp.prod([nodes[q] - nodes[p] for p in range(3) if p != q])
                              for q in range(3)))
    return identity, contact


def derive_only_up_radial_vertex():
    """Only-up f_a g_xy(q^2) sqrt(N_x N_y) for neutral legs with an ALP component.

    Leg meson matrix (2 B0 f_pi/f_a) R_x/sqrt(N_x) with R_x = M^{-1}_{m_x^2}(P_01) and (R_x)_11 = -m_u E_x,
    E_x = <e_1 M^{-1}_{m_x^2}(e_1)>; ALP component 1/sqrt(N_x).  Pole equation
    m_x^2 - B0 m_u f_pi^2/f_a^2 = -(2 B0^2 m_u^2 f_pi^2/f_a^2) E_x.
    Returns (closed form, residual against m_x^2 m_y^2 f_a^2/(B0 m_u f_pi^2) + q^2 - m_x^2 - m_y^2).
    """
    up, momentum, mass_x, mass_y, norm_x, norm_y = sp.symbols('m_u q2 m_x2 m_y2 N_x N_y', positive=True)
    inverse_x, inverse_y = sp.symbols('E_x E_y', real=True)
    leg = {}
    for label, inverse, norm in (('x', inverse_x, norm_x), ('y', inverse_y, norm_y)):
        leg[label] = (-(2 * B_0 * f_pi * up / f_a) * inverse / sp.sqrt(norm), 1 / sp.sqrt(norm))
    factor = lambda c: c[0] + f_pi / f_a * c[1]
    f_a_g = (B_0 * up * factor(leg['x']) * factor(leg['y'])
             + leg['x'][1] * leg['y'][1] * (momentum - mass_x - mass_y))
    pole = {inverse_x: sp.solve(sp.Eq(mass_x - B_0 * up * f_pi ** 2 / f_a ** 2,
                                      -2 * B_0 ** 2 * up ** 2 * f_pi ** 2 / f_a ** 2 * inverse_x), inverse_x)[0],
            inverse_y: sp.solve(sp.Eq(mass_y - B_0 * up * f_pi ** 2 / f_a ** 2,
                                      -2 * B_0 ** 2 * up ** 2 * f_pi ** 2 / f_a ** 2 * inverse_y), inverse_y)[0]}
    closed = sp.simplify(sp.expand(f_a_g.subs(pole) * sp.sqrt(norm_x * norm_y)))
    expected = mass_x * mass_y * f_a ** 2 / (B_0 * up * f_pi ** 2) + momentum - mass_x - mass_y
    return closed, sp.simplify(closed - expected)


# ---------------------------------------------------------------------------
# (2) Parameter sets
# ---------------------------------------------------------------------------

def charge_model_parameters(charges, quark_masses=(0.0022, 0.0047, 0.095), decay_constant=0.4,
                            quartic_value=0.8, B0=2.6, pion_decay_constant=0.0924, singlet_mass_squared=0.75):
    """M = diag[m_q (1 + sigma/f_a)^{n_q} exp(-i n_q a/f_a)], coefficients with r + s <= 4."""
    sources = {'S': {}, 'P': {}}
    for (r, s) in itertools.product(range(MAXIMUM_SOURCE_ORDER + 1), repeat=2):
        if r + s > MAXIMUM_SOURCE_ORDER:
            continue
        entries = [complex(mass * np.prod([n - j for j in range(r)]) / math.factorial(r)
                           * (-1j * n) ** s / math.factorial(s)) for mass, n in zip(quark_masses, charges)]
        sources['S'][(r, s)] = np.diag([e.real for e in entries])
        sources['P'][(r, s)] = np.diag([e.imag for e in entries])
    parameters = {'B_0': B0, 'f_pi': pion_decay_constant, 'f_a': decay_constant,
                  'm_1_squared': singlet_mass_squared, 'lambda': quartic_value, 'sources': sources}
    return impose_stationarity(parameters)


def axial_rotation_parameters(axial_coupling, quark_masses=(0.0022, 0.0047, 0.095), decay_constant=0.4,
                              quartic_value=0.8):
    """M = (1 + sigma/f_a) exp(i x c) m exp(i x c) with a Hermitian (flavour-violating) c.

    Satisfies M(-x) = M(x)^dagger for any Hermitian c: parity conserving but flavour violating.
    """
    mass = np.diag(quark_masses).astype(complex)
    rotation = 1j * axial_coupling
    taylor = []
    for s in range(MAXIMUM_SOURCE_ORDER + 1):
        total = np.zeros((3, 3), dtype=complex)
        for r in range(s + 1):
            total += (np.linalg.matrix_power(rotation, r) @ mass @ np.linalg.matrix_power(rotation, s - r)
                      / (math.factorial(r) * math.factorial(s - r)))
        taylor.append(total)
    sources = {'S': {}, 'P': {}}
    for r in range(2):
        for s in range(MAXIMUM_SOURCE_ORDER + 1 - r):
            sources['S'][(r, s)] = ((taylor[s] + taylor[s].conj().T) / 2)
            sources['P'][(r, s)] = ((taylor[s] - taylor[s].conj().T) / 2j)
    parameters = {'B_0': 2.6, 'f_pi': 0.0924, 'f_a': decay_constant, 'm_1_squared': 0.75,
                  'lambda': quartic_value, 'sources': sources}
    return impose_stationarity(parameters)


def with_extra_sources(parameters, extra):
    """Copy with additional source coefficients, extra = {('S' or 'P', (r, s)): matrix}."""
    result = dict(parameters)
    result['sources'] = {'S': dict(parameters['sources']['S']), 'P': dict(parameters['sources']['P'])}
    for (kind, order), matrix in extra.items():
        result['sources'][kind][order] = result['sources'][kind].get(order, np.zeros((3, 3))) + matrix
    return impose_stationarity(result)


# ---------------------------------------------------------------------------
# (2) Numerical checks on the exact Lagrangian
# ---------------------------------------------------------------------------

def flavour_direction(matrix, angular_component=0.0):
    """Eleven-field direction (sigma, a, P_0..P_8) of a meson matrix, components <lambda_k X>/2."""
    direction = np.zeros(NUMBER_OF_FIELDS, dtype=complex)
    direction[1] = angular_component
    direction[2:] = [0.5 * np.trace(LAMBDA[k] @ matrix) for k in range(9)]
    return direction


def elementary(row, column):
    matrix = np.zeros((3, 3), dtype=complex)
    matrix[row, column] = 1.0
    return matrix


CHARGED_MATRICES = {  # sqrt2 |u><d| etc. (u, d, s = 0, 1, 2); one convention for each pair
    'pi+': np.sqrt(2) * elementary(1, 0), 'pi-': np.sqrt(2) * elementary(0, 1),
    'K+': np.sqrt(2) * elementary(2, 0), 'K-': np.sqrt(2) * elementary(0, 2),
    'K0': np.sqrt(2) * elementary(2, 1), 'K0bar': np.sqrt(2) * elementary(1, 2)}


def rayleigh_mass_squared(data, direction):
    """Squared mass of a direction inside an eigenspace of the polar Hessian."""
    return float((direction.conj() @ data['mass_squared'] @ direction).real / (direction.conj() @ direction).real)


def multilinear_complex(function, directions):
    """Evaluate a function multilinear in real directions on complex directions (dict of values)."""
    options = []
    for direction in directions:
        choices = [(1.0, direction.real.copy())]
        if np.any(direction.imag != 0):
            choices.append((1j, direction.imag.copy()))
        options.append(choices)
    total = {}
    for combination in itertools.product(*options):
        weight = np.prod([w for w, _ in combination])
        for key, value in function([v for _, v in combination]).items():
            total[key] = total.get(key, 0) + weight * np.asarray(value)
    return total


def radial_channels(parameters, data, directions, gram, invariants):
    """sigma exchange of each channel (k|ij) from the exact Lagrangian, real directions."""
    unit = np.eye(4)
    radial = np.eye(NUMBER_OF_FIELDS)[0]
    values = []
    for attached, (first, second) in PAIRINGS:
        internal = -(unit[0] + unit[attached])
        left = vertex(parameters, 'polar', [directions[0], directions[attached], radial],
                      [unit[0], unit[attached], internal], gram)
        right = vertex(parameters, 'polar', [radial, directions[first], directions[second]],
                       [-internal, unit[first], unit[second]], gram)
        momentum_squared = invariants[frozenset({first, second})]
        propagator = np.linalg.inv(data['mass_squared'] - momentum_squared * data['kinetic'])[0, 0]
        values.append(left * propagator * right)
    return np.array(values)


def direct_pieces(parameters, data, directions, masses_squared, s12, s13):
    """Direct polar amplitude for complex leg directions (parent first).

    Returns dict with 'total', 'contact', 'exchange', 'radial_exchange', 'radial_channels' and
    'non_radial_exchange' (meson and ALP exchange, i.e. V3 G V3).
    """
    invariants = dalitz_invariants(masses_squared[0], masses_squared[1:], s12, s13)
    gram = decay_gram(masses_squared[0], masses_squared[1:], invariants)

    def real_pieces(real_directions):
        vectors = np.column_stack(real_directions)
        result = direct_amplitude(parameters, 'polar', data['mass_squared'], data['kinetic'],
                                  np.array(masses_squared), vectors, (0, 1, 2, 3), invariants)
        return {'total': result['total'], 'contact': result['contact'], 'exchange': result['exchange'],
                'radial_exchange': result['radial_exchange'],
                'radial_channels': radial_channels(parameters, data, real_directions, gram, invariants)}

    pieces = multilinear_complex(real_pieces, directions)
    pieces['non_radial_exchange'] = pieces['exchange'] - pieces['radial_exchange']
    pieces['invariants'] = invariants
    pieces['gram'] = gram
    return pieces


def maximum_cubic_vertex(parameters, seed=1):
    """Largest |V3| among all cubic vertices of (a, P_0..P_8) at random momenta, and the largest
    sigma-P-P vertex for scale."""
    generator = np.random.default_rng(seed)
    momenta = generator.normal(size=(3, 3))
    gram = momenta @ np.diag([1.0, -1.0, -1.0]) @ momenta.T
    identity = np.eye(NUMBER_OF_FIELDS)
    largest, radial = 0.0, 0.0
    for x, y, z in itertools.combinations_with_replacement(range(NUMBER_OF_FIELDS), 3):
        value = abs(vertex(parameters, 'polar', [identity[x], identity[y], identity[z]], np.eye(3), gram))
        if x == 0 and y > 0:
            radial = max(radial, value)
        elif x > 0:
            largest = max(largest, value)
    return largest, radial


# --- neutral-sector cosine model ------------------------------------------

def neutral_sector_model(charges, quark_masses, B0, pion_decay_constant, decay_constant, singlet_mass_squared,
                         bare_alp_mass_squared=0.0):
    """Diagonal neutral sector at sigma = 0 in canonical fields (Phi_uu, Phi_dd, Phi_ss, a)/(sqrt2, sqrt2, sqrt2, 1).

    V = -B0 f_pi^2 sum_q m_q cos(Phi_qq/f_pi + n_q a/f_a) + (m_1^2/12)(sum_q Phi_qq)^2 + mu^2 a^2/2.
    Returns (pole masses squared, pole vectors as columns, contact(list of four vectors)).
    The contact is V4 = B0 f_pi^2 sum_q m_q prod_x (cosine argument of flavour q on v_x); cubic
    vertices vanish.
    """
    arguments = np.zeros((3, 4))  # cosine argument of flavour q in canonical coordinates
    for q in range(3):
        arguments[q, q] = np.sqrt(2) / pion_decay_constant
        arguments[q, 3] = charges[q] / decay_constant
    weights = B0 * pion_decay_constant ** 2 * np.asarray(quark_masses)
    singlet = np.array([1.0, 1.0, 1.0, 0.0])
    hessian = (sum(weights[q] * np.outer(arguments[q], arguments[q]) for q in range(3))
               + singlet_mass_squared / 3 * np.outer(singlet, singlet))
    hessian[3, 3] += bare_alp_mass_squared
    masses_squared, vectors = np.linalg.eigh(hessian)

    def contact(legs):
        return float(sum(weights[q] * np.prod([arguments[q] @ leg for leg in legs]) for q in range(3)))

    return masses_squared, vectors, contact


def neutral_coordinates(direction):
    """Canonical (Phi_uu/sqrt2, Phi_dd/sqrt2, Phi_ss/sqrt2, a) of an eleven-field direction (diagonal mesons)."""
    diagonal = np.einsum('k,kii->i', direction[2:], LAMBDA).real / np.sqrt(2)
    return np.concatenate([diagonal, [direction[1].real]])


def cosine_direction_masses_squared(parameters, charges):
    """2 B0 m_q (1 + n_q^2 f_pi^2/(2 f_a^2)) for q = u, d, s (GeV^2)."""
    masses = np.diag(parameters['sources']['S'][(0, 0)]).real
    return 2 * parameters['B_0'] * masses * (1 + np.asarray(charges, dtype=float) ** 2
                                             * parameters['f_pi'] ** 2 / (2 * parameters['f_a'] ** 2))


# --- pole legs and radial vertex ------------------------------------------

def column_from_direction(direction):
    """(meson matrix, ALP component) of an eleven-field direction in the original polar fields."""
    return np.einsum('k,kij->ij', direction[2:], LAMBDA), direction[1]


def inverse_mass_operator_of_p01(parameters, mass_squared):
    """M^{-1}_{m^2}(P_01) (GeV^-1)."""
    return inverse_mass_operator(parameters, mass_squared, parameters['sources']['P'][(0, 1)])


def pole_leg(parameters, mass_squared):
    """Original-field meson matrix and ALP component of a pole with an ALP component:
    ((2 B0 f_pi/f_a) M^{-1}_{m^2}(P_01)/sqrt(N(m^2)), 1/sqrt(N(m^2)))."""
    root = np.sqrt(pole_normalization(parameters, mass_squared))
    scale = 2 * parameters['B_0'] * parameters['f_pi'] / parameters['f_a']
    return scale * inverse_mass_operator_of_p01(parameters, mass_squared) / root, 1 / root


def radial_vertex_from_direction(parameters, first, second, momentum_squared):
    """g_xy(q^2) of the radial-exchange record on (direction, m^2) pairs."""
    (direction_x, mass_x), (direction_y, mass_y) = first, second
    return radial_vertex_trace_form(parameters, (column_from_direction(direction_x), mass_x),
                                    (column_from_direction(direction_y), mass_y), momentum_squared)


def only_up_radial_vertex_closed_form(parameters, mass_x, mass_y, momentum_squared):
    """Only-up, neutral legs with an ALP component:
    g_xy(q^2) = [q^2 - m_x^2 - m_y^2 + m_x^2 m_y^2 f_a^2/(B0 m_u f_pi^2)]/(f_a sqrt(N(m_x^2) N(m_y^2)))."""
    B0, fpi, fa = parameters['B_0'], parameters['f_pi'], parameters['f_a']
    up = parameters['sources']['S'][(1, 0)][0, 0]
    return ((momentum_squared - mass_x - mass_y + mass_x * mass_y * fa ** 2 / (B0 * up * fpi ** 2))
            / (fa * np.sqrt(pole_normalization(parameters, mass_x) * pole_normalization(parameters, mass_y))))


def radial_coupling_invariant_form(parameters, mass_x, mass_y, momentum_squared):
    """g_xy(q^2) for two legs with an ALP component, written with invariant quantities only:
    f_a sqrt(N(m_x^2) N(m_y^2)) g_xy(q^2)
      = (2 B0^3 f_pi^2/f_a^2) <S_10 {M^{-1}_{m_x^2}(P_01), M^{-1}_{m_y^2}(P_01)}>
        - (2 B0^2 f_pi^2/f_a^2) <P_11 (M^{-1}_{m_x^2}(P_01) + M^{-1}_{m_y^2}(P_01))>
        - (2 B0 f_pi^2/f_a^2) <S_12> + q^2 - m_x^2 - m_y^2  (GeV^-1)."""
    B0, fpi, fa = parameters['B_0'], parameters['f_pi'], parameters['f_a']
    S, P = parameters['sources']['S'], parameters['sources']['P']
    inverse_x = inverse_mass_operator_of_p01(parameters, mass_x)
    inverse_y = inverse_mass_operator_of_p01(parameters, mass_y)
    braces = (2 * B0 ** 3 * fpi ** 2 / fa ** 2 * np.trace(S[(1, 0)] @ (inverse_x @ inverse_y + inverse_y @ inverse_x))
              - 2 * B0 ** 2 * fpi ** 2 / fa ** 2 * np.trace(P[(1, 1)] @ (inverse_x + inverse_y))
              - 2 * B0 * fpi ** 2 / fa ** 2 * np.trace(S[(1, 2)])
              + momentum_squared - mass_x - mass_y)
    return braces.real / (fa * np.sqrt(pole_normalization(parameters, mass_x) * pole_normalization(parameters, mass_y)))


def generator_components(matrix):
    """Components <lambda_k X>/2, k = 0..8, of a Hermitian flavour matrix on the U(3) generators."""
    return np.array([0.5 * np.trace(LAMBDA[k] @ matrix).real for k in range(9)])


def radial_coupling_generator_component_form(parameters, mass_x, mass_y, momentum_squared):
    """The earlier form of g_xy(q^2) in the radial-decay record, used only as a cross-check: the
    meson matrix of each pole rebuilt from its components c_k on the generators,
    X = sum_k c_k lambda_k with c_k = <lambda_k (2 B0 f_pi/f_a) M^{-1}_{m_x^2}(P_01)>/2, N = 1 + sum_k c_k^2,
    g_xy = {(B0/2) <S_10 {X, Y}> - (B0 f_pi/f_a) <P_11 (X + Y)> - (2 B0 f_pi^2/f_a^2) <S_12>
            + q^2 - m_x^2 - m_y^2}/(f_a sqrt(N_x N_y))."""
    B0, fpi, fa = parameters['B_0'], parameters['f_pi'], parameters['f_a']
    S, P = parameters['sources']['S'], parameters['sources']['P']
    rebuilt, norms = [], []
    for mass_squared in (mass_x, mass_y):
        components = generator_components(alp_shift_matrix(parameters, mass_squared))
        rebuilt.append(np.einsum('k,kij->ij', components, LAMBDA))
        norms.append(1 + components @ components)
    first, second = rebuilt
    braces = (B0 / 2 * np.trace(S[(1, 0)] @ (first @ second + second @ first))
              - B0 * fpi / fa * np.trace(P[(1, 1)] @ (first + second))
              - 2 * B0 * fpi ** 2 / fa ** 2 * np.trace(S[(1, 2)])
              + momentum_squared - mass_x - mass_y)
    return braces.real / (fa * np.sqrt(norms[0] * norms[1]))


def no_mixing_field_change_invariant_form(parameters, parent_mass_squared):
    """Linear field change to the basis without ALP-meson mixing (before the meson-only rotation),
    as written in the radial-decay record with invariant quantities:
      Phi = Phi'' - (2 B0 f_pi/f_a) M^{-1}_{m_a^2}(P_01) (B0 f_pi/f_a) <M^{-1}_{m_a^2}(P_01) Phi''>/(sqrt(N)(sqrt(N) + 1))
            + (2 B0 f_pi/f_a) M^{-1}_{m_a^2}(P_01) a_phys/sqrt(N),
      a   = -(B0 f_pi/f_a) <M^{-1}_{m_a^2}(P_01) Phi''>/sqrt(N) + a_phys/sqrt(N),   N = N(m_a^2).
    Returns the 10 x 10 matrix mapping (P''_0..P''_8, a_phys) to (P_0..P_8, a) (Phi = sum_k P_k lambda_k)."""
    B0, fpi, fa = parameters['B_0'], parameters['f_pi'], parameters['f_a']
    inverse = inverse_mass_operator_of_p01(parameters, parent_mass_squared)
    shift = alp_shift_matrix(parameters, parent_mass_squared)
    root = np.sqrt(pole_normalization(parameters, parent_mass_squared))
    matrix = np.zeros((10, 10))
    for j in range(9):
        overlap = B0 * fpi / fa * np.trace(inverse @ LAMBDA[j]).real
        meson = LAMBDA[j] - shift * overlap / (root * (root + 1))
        matrix[:9, j] = generator_components(meson)
        matrix[9, j] = -overlap / root
    matrix[:9, 9] = generator_components(shift / root)
    matrix[9, 9] = 1 / root
    return matrix


def no_mixing_field_change_component_form(parameters, parent_mass_squared):
    """The earlier component form of the same matrix, used only as a cross-check:
    [[delta_ij - c_i c_j/(sqrt(N)(sqrt(N) + 1)), c_i/sqrt(N)], [-c_j/sqrt(N), 1/sqrt(N)]],
    c_k = <lambda_k (2 B0 f_pi/f_a) M^{-1}_{m_a^2}(P_01)>/2, N = 1 + sum_k c_k^2."""
    components = generator_components(alp_shift_matrix(parameters, parent_mass_squared))
    root = np.sqrt(1 + components @ components)
    matrix = np.zeros((10, 10))
    matrix[:9, :9] = np.eye(9) - np.outer(components, components) / (root * (root + 1))
    matrix[:9, 9] = components / root
    matrix[9, :9] = -components / root
    matrix[9, 9] = 1 / root
    return matrix


# --- quartic vertex in invariant form ------------------------------------------

def quartic_vertex_from_legs(parameters, legs, gram):
    """Fully symmetric quartic vertex of the original polar fields (quartic_vertex_stated).

    legs = [(meson matrix F_x, ALP component alpha_x)] x 4 (parent first), gram = k_x.k_y of the
    all-outgoing momenta (valid on and off shell).  Uses only S_00, P_01, S_02, P_03, S_04.
    """
    B0, fpi, fa = parameters['B_0'], parameters['f_pi'], parameters['f_a']
    zero = np.zeros((3, 3))
    S, P = parameters['sources']['S'], parameters['sources']['P']
    scalar = [S.get((0, n), zero) for n in range(5)]
    pseudoscalar = [P.get((0, n), zero) for n in range(5)]
    matrices = [np.asarray(leg[0], dtype=complex) for leg in legs]
    alp = [complex(leg[1]) for leg in legs]
    trace = np.trace
    commutator = lambda first, second: first @ second - second @ first
    total = 0.0
    for permutation in itertools.permutations(range(4)):
        total += trace(scalar[0] @ _matrix_product([matrices[n] for n in permutation])) / (24 * fpi ** 4)
    for x in range(4):
        for permutation in itertools.permutations([y for y in range(4) if y != x]):
            total -= alp[x] * trace(pseudoscalar[1] @ _matrix_product([matrices[n] for n in permutation])) / (6 * fa * fpi ** 3)
    for x, y in itertools.combinations(range(4), 2):
        u, v = [z for z in range(4) if z not in (x, y)]
        total -= alp[x] * alp[y] * trace(scalar[2] @ (matrices[u] @ matrices[v] + matrices[v] @ matrices[u])) / (
            fa ** 2 * fpi ** 2)
    for x in range(4):
        total += 6 * np.prod([alp[y] for y in range(4) if y != x]) * trace(pseudoscalar[3] @ matrices[x]) / (fa ** 3 * fpi)
    total += 24 * np.prod(alp) * trace(scalar[4]) / fa ** 4
    total *= B0 * fpi ** 2
    for x, y in itertools.combinations(range(4), 2):
        u, v = [z for z in range(4) if z not in (x, y)]
        total -= gram[x, y] / (24 * fpi ** 2) * trace(
            commutator(matrices[u], matrices[x]) @ commutator(matrices[v], matrices[y])
            + commutator(matrices[v], matrices[x]) @ commutator(matrices[u], matrices[y]))
    return complex(total)


def original_basis_direct_quartic(parameters, directions, gram):
    """Quartic vertex of the exact polar Lagrangian in the original fields on complex directions."""
    return multilinear_complex(lambda real: {'contact': vertex(parameters, 'polar', real, np.eye(4), gram)},
                               directions)['contact']


def alp_shift_matrix(parameters, parent_mass_squared):
    """(2 B0 f_pi/f_a) M^{-1}_{m_a^2}(P_01): the meson matrix multiplying a' in the linear part
    Phi = Phi' + a' (2 B0 f_pi/f_a) M^{-1}_{m_a^2}(P_01) of the chiral shift eq:murayama-chiral-shift."""
    return (2 * parameters['B_0'] * parameters['f_pi'] / parameters['f_a']
            * inverse_mass_operator_of_p01(parameters, parent_mass_squared))


def unmixed_basis_lagrangian(fields, parameters, parent_mass_squared, mask_dot):
    """Polar Lagrangian in the fields (sigma, a', Phi') of the chiral shift eq:murayama-chiral-shift:
    U = g U' g, g = exp(i a' B0 M^{-1}_{m_a^2}(P_01)/f_a), a = a'.
    <Phi> = <Phi'> + a' (2 B0 f_pi/f_a) <M^{-1}_{m_a^2}(P_01)>."""
    rotation_generator = alp_shift_matrix(parameters, parent_mass_squared)
    decay_constant, pion_decay_constant = parameters['f_a'], parameters['f_pi']
    radial, angular = fields[0], fields[1]
    kinetic = (derivative_product(radial, radial, mask_dot) * 0.5
               + (1 + radial / decay_constant) * (1 + radial / decay_constant)
               * derivative_product(angular, angular, mask_dot) * 0.5)
    meson_matrix = radial.constant_like(np.zeros((3, 3)))
    for k in range(9):
        meson_matrix = meson_matrix + fields[2 + k].times_matrix(LAMBDA[k])
    half = angular.times_matrix(rotation_generator * (1j / (2 * pion_decay_constant)))
    rotation, rotation_dagger = matrix_exponential(half), matrix_exponential(-half)
    chiral_field = rotation.matmul(matrix_exponential(meson_matrix * (1j / pion_decay_constant))).matmul(rotation)
    chiral_field_dagger = rotation_dagger.matmul(
        matrix_exponential(meson_matrix * (-1j / pion_decay_constant))).matmul(rotation_dagger)
    kinetic = kinetic + derivative_product(chiral_field, chiral_field_dagger, mask_dot, matrices=True).trace() * (
        pion_decay_constant ** 2 / 4)
    scalar, pseudoscalar = source_matrices(parameters['sources'], radial, angular, decay_constant)
    chi = (scalar + pseudoscalar * 1j) * (2 * parameters['B_0'])
    chi_plus_trace = (chi.matmul(chiral_field_dagger) + chi.dagger().matmul(chiral_field)).trace()
    singlet_trace = meson_matrix.trace() + angular * np.trace(rotation_generator).real
    rho_squared = (radial + decay_constant) * (radial + decay_constant)
    radial_potential = rho_squared * (parameters['m_phi_squared'] / 2) + rho_squared * rho_squared * (
        parameters['lambda'] / 4)
    return (kinetic + chi_plus_trace * (pion_decay_constant ** 2 / 4)
            - singlet_trace * singlet_trace * (parameters['m_1_squared'] / 12) - radial_potential)


def unmixed_basis_vertex(parameters, parent_mass_squared, directions, gram):
    """n-point vertex of unmixed_basis_lagrangian on real directions (leg l has momentum k_l)."""
    size = 1 << len(directions)
    mask_momenta = np.zeros((size, len(directions)))
    for mask in range(size):
        for leg in range(len(directions)):
            if mask >> leg & 1:
                mask_momenta[mask, leg] = 1.0
    mask_dot = mask_momenta @ gram @ mask_momenta.T
    return unmixed_basis_lagrangian(field_point(directions), parameters, parent_mass_squared,
                                    mask_dot).coefficients[-1].real


def unmixed_basis_direct_quartic(parameters, parent_mass_squared, directions, gram):
    """Internal cross-check: quartic vertex in the basis of the chiral shift on the same poles.

    directions are eleven-field directions in the original fields; the linear part of the shift,
    Phi' = Phi - a (2 B0 f_pi/f_a) M^{-1}_{m_a^2}(P_01), gives the directions in the new fields
    (the meson matrix of each direction is shifted by its ALP component times alp_shift_matrix)."""
    shift = flavour_direction(alp_shift_matrix(parameters, parent_mass_squared))
    transformed = [np.array(direction, dtype=complex) - shift * direction[1] for direction in directions]
    return multilinear_complex(
        lambda real: {'contact': unmixed_basis_vertex(parameters, parent_mass_squared, real, gram)},
        transformed)['contact']


def off_shell_gram(seed=3):
    """Gram matrix of four random all-outgoing momenta with k_0 = -(k_1 + k_2 + k_3) (GeV^2)."""
    generator = np.random.default_rng(seed)
    daughters = 0.3 * generator.normal(size=(3, 4))
    momenta = np.vstack([-daughters.sum(axis=0), daughters])
    return momenta @ np.diag([1.0, -1.0, -1.0, -1.0]) @ momenta.T


def diagonal_neutral_quartic_vertex(parameters, charges, masses_squared):
    """Diagonal PQ-charge sources, four legs with an ALP component:
    V4 = (B0 f_pi^2/f_a^4) prod_x N(m_x^2)^{-1/2} sum_q m_q prod_x [2 B0 (M^{-1}_{m_x^2}(P_01))_qq + n_q]."""
    B0, fpi, fa = parameters['B_0'], parameters['f_pi'], parameters['f_a']
    masses = np.diag(parameters['sources']['S'][(0, 0)]).real
    total = 0.0
    for q in range(3):
        total += masses[q] * np.prod([2 * B0 * inverse_mass_operator_of_p01(parameters, m2)[q, q].real + charges[q] for m2 in masses_squared])
    return B0 * fpi ** 2 / fa ** 4 * total / np.prod([np.sqrt(pole_normalization(parameters, m2)) for m2 in masses_squared])


def single_flavour_neutral_quartic_vertex(parameters, charges, masses_squared):
    """One PQ-charged flavour, four legs with an ALP component:
    V4 = (B0 f_pi^2/f_a^4)(2 B0 m_1^2/3)^4 prod_x <M^{-1}_{m_x^2}(P_01)>/sqrt(N(m_x^2))
         * sum_q m_q / prod_x (m_x^2 - 2 B0 m_q (1 + n_q^2 f_pi^2/(2 f_a^2)))."""
    B0, fpi, fa = parameters['B_0'], parameters['f_pi'], parameters['f_a']
    masses = np.diag(parameters['sources']['S'][(0, 0)]).real
    nodes = cosine_direction_masses_squared(parameters, charges)
    traces = np.prod([np.trace(inverse_mass_operator_of_p01(parameters, m2)).real / np.sqrt(pole_normalization(parameters, m2))
                      for m2 in masses_squared])
    total = sum(masses[q] / np.prod([m2 - nodes[q] for m2 in masses_squared]) for q in range(3))
    return B0 * fpi ** 2 / fa ** 4 * (2 * B0 * parameters['m_1_squared'] / 3) ** 4 * traces * total


def diagonal_charged_pair_quartic_vertex(parameters, charges, pair, parent_mass_squared, daughter_mass_squared,
                                         invariants):
    """charged_pair_form with A = diag of the parent meson matrix, I = diag of the neutral daughter's,
    alpha = 1/sqrt(N); invariants = (s_ij, s_ik, s_jk) (GeV^2)."""
    masses = np.diag(parameters['sources']['S'][(0, 0)]).real
    parent_matrix, parent_alp = pole_leg(parameters, parent_mass_squared)
    daughter_matrix, daughter_alp = pole_leg(parameters, daughter_mass_squared)
    expression = charged_pair_form(list(masses), list(charges), list(np.diag(parent_matrix).real),
                                   list(np.diag(daughter_matrix).real), parent_alp, daughter_alp, invariants, pair)
    return float(sp.N(sp.sympify(expression).subs({f_pi: parameters['f_pi'], f_a: parameters['f_a'],
                                                   B_0: parameters['B_0']})))


# ---------------------------------------------------------------------------
# Drivers
# ---------------------------------------------------------------------------

def neutral_legs(data):
    """Parent and the other neutral poles with an ALP component, as (direction, m^2)."""
    parent = data['parent']
    others = [n for n in data['mixed'] if n != parent]
    as_leg = lambda n: (data['vectors'][:, n].astype(complex), data['masses_squared'][n])
    return as_leg(parent), [as_leg(n) for n in others]


def charged_leg(data, name):
    direction = flavour_direction(CHARGED_MATRICES[name])
    return direction, rayleigh_mass_squared(data, direction)


def final_states(data):
    """Named daughter sets for the only-up model (neutral poles ordered by mass)."""
    _, neutral = neutral_legs(data)
    neutral = sorted(neutral, key=lambda leg: leg[1])
    pion, eta, eta_prime = neutral
    charged = {name: charged_leg(data, name) for name in CHARGED_MATRICES}
    return {
        'pi0 eta eta_prime': [pion, eta, eta_prime],
        'pi0 pi0 eta': [pion, pion, eta],
        'pi0 pi0 pi0': [pion, pion, pion],
        'pi0 pi+ pi-': [pion, charged['pi+'], charged['pi-']],
        'eta pi+ pi-': [eta, charged['pi+'], charged['pi-']],
        'pi0 K+ K-': [pion, charged['K+'], charged['K-']],
        'eta K+ K-': [eta, charged['K+'], charged['K-']],
        'pi0 K0 K0bar': [pion, charged['K0'], charged['K0bar']],
    }


def invariant_legs(parameters, legs):
    """(meson matrix, ALP component) for each (direction, m^2): pole_leg if it has an ALP component,
    otherwise its meson generator and 0."""
    result = []
    for direction, mass_squared in legs:
        if abs(direction[1]) > 0:
            result.append(pole_leg(parameters, mass_squared))
        else:
            result.append((column_from_direction(direction)[0], 0.0))
    return result


def assemble(parameters, data, daughters, s12, s13):
    """Direct pieces, the invariant quartic vertex, the vertex in the basis without ALP-meson mixing,
    and the assembled formula (invariant V4 + trace-form sigma exchange)."""
    parent, _ = neutral_legs(data)
    legs = [parent] + daughters
    masses_squared = [leg[1] for leg in legs]
    pieces = direct_pieces(parameters, data, [leg[0] for leg in legs], masses_squared, s12, s13)
    invariants, gram = pieces['invariants'], pieces['gram']
    pieces['invariant_quartic'] = quartic_vertex_from_legs(parameters, invariant_legs(parameters, legs), gram)
    pieces['unmixed_basis_contact'] = unmixed_basis_direct_quartic(parameters, parent[1], [leg[0] for leg in legs], gram)
    channels = []
    radial_mass_squared = radial_mass_squared_trace_form(parameters)
    for attached, (first, second) in PAIRINGS:
        momentum_squared = invariants[frozenset({first, second})]
        channels.append(radial_vertex_from_direction(parameters, legs[0], legs[attached], momentum_squared)
                        * radial_vertex_from_direction(parameters, legs[first], legs[second], momentum_squared)
                        / (radial_mass_squared - momentum_squared))
    pieces['radial_channels_trace_form'] = np.array(channels)
    pieces['assembled'] = pieces['invariant_quartic'] + sum(channels)
    return pieces


def main():
    print('Units: GeV. Tree level, LO U(3), polar coordinates, no expansion in f_pi/f_a.\n')
    _, violations = derive_charge_source_parity()
    print('PQ-charge sources m(1+rho)^n e^{-inx}: parity violations', violations)
    general, parity, kinetic = derive_cubic_trace_parity()
    print('cubic trace words, general sources:', general)
    print('cubic trace words, S_odd = P_even = 0:', parity, '; kinetic:', kinetic)
    chiral_map = derive_symmetric_chiral_map()
    print('log(e^X e^Y e^X) by order:', chiral_map)
    print('quartic Lagrangian (derived - stated):', derive_quartic_lagrangian())
    print('quartic vertex (brute force - stated):', derive_quartic_vertex())
    print('kinetic channel identity:', derive_kinetic_channel_identity())
    print('pole-leg form (difference):', derive_pole_leg_form())
    print('diagonal forms (neutral, pair ud, pair ds):', derive_diagonal_forms())
    print('neutral pole relation:', derive_neutral_pole_relation())
    print('secular product residuals:', derive_neutral_secular_product())
    print('three-term identity, contact sum:', derive_three_term_identity())
    print('only-up sigma vertex closed form:', derive_only_up_radial_vertex(), '\n')

    parameters = charge_model_parameters((1, 0, 0))
    data = spectrum(parameters)
    largest, radial = maximum_cubic_vertex(parameters)
    print(f'only-up: max |V3| = {largest:.2e}, max |sigma P P| = {radial:.2e}')
    for name, daughters in final_states(data).items():
        pieces = assemble(parameters, data, daughters, 0.3, 0.5)
        scale = abs(pieces['total'])
        print(f'  {name:18s} total {complex(pieces["total"]).real:+.6e}  V3GV3 {abs(pieces["non_radial_exchange"]):.1e}'
              f'  contact {complex(pieces["contact"]).real:+.3e}  sigma channels '
              f'{np.round(pieces["radial_channels"].real, 12)}  |assembled-direct|/|M| '
              f'{abs(pieces["assembled"] - pieces["total"]) / scale:.1e}  |V4 unmixed - V4 original|/|M| '
              f'{abs(pieces["unmixed_basis_contact"] - pieces["contact"]) / scale:.1e}')
    for label, sources in (('only-up', parameters), ('random sources (seed 2026)', random_parameters(2026))):
        data = spectrum(sources)
        parent, others = neutral_legs(data)
        directions = [leg[0] for leg in [parent] + others[:3]]
        masses_squared = [leg[1] for leg in [parent] + others[:3]]
        on_shell = decay_gram(masses_squared[0], masses_squared[1:],
                              dalitz_invariants(masses_squared[0], masses_squared[1:], 0.3, 0.5))
        for kind, gram in (('on shell', on_shell), ('off shell', off_shell_gram())):
            original = original_basis_direct_quartic(sources, directions, gram)
            unmixed = unmixed_basis_direct_quartic(sources, masses_squared[0], directions, gram)
            print(f'{label}, {kind}: V4 original fields {complex(original).real:+.6e}, '
                  f'basis without mixing {complex(unmixed).real:+.6e}')


if __name__ == '__main__':
    main()
