"""Exact checks for the selection of the observable A -> 3 pi0 (OpenSpec change
estimate-loop-tree-crossing-finite-mixing, task 1.2).

Checks
    S1  the neutral kinetic term is free of interactions (closed form in the fields), so every
        neutral vertex of the dressed frame is momentum independent;
    S2  on-shell structure of contact amplitudes polynomial in the scalar products;
    S3  the bubble function Bbar on the Dalitz plane (Taylor coefficients, coefficient of e3 at
        the symmetric point, dependence on the subtraction scale);
    S4  coefficients of Bbar(s; m_r, m_t) with an ALP-like internal line ("A-threshold bubbles")
        in the one-loop A -> 3 pi0 amplitude, and their behaviour when the leading tree
        coefficient vanishes and in the isospin limit;
    S5  equality of these coefficients in the dressed frame and in the kappa frame.

Conventions (Alda:2026nrd, as in verify_alp_loop_observable_candidates.py, which is imported and
not modified).  Units GeV.  U = exp(i Phi/f_pi), tr(lambda_i lambda_j) = 2 delta_ij, N_f = 3,
x = a/f_a, epsilon = f_pi/f_a.  Dressed frame, LO U(3):
    L = f^2/4 <dU dU^dag> + f^2 B0/2 <M_hat(x) U^dag + U M_hat(x)^dag> - f^2 m0^2/(4 N_f) X_TC^2
        + (da)^2/2 - m_a0^2 a^2/2,     X_TC = -tr Phi/f - 2 Theta_hat(x).
Neutral fields N_q = Phi_qq.  The invariants S_hat_n, P_hat_n, Theta_hat_n are sympy symbols;
no value is assumed for any of them.  Meson-sector rationals (meson_sector_inputs of the imported
module, m_u = m_d) are test inputs used only in S4 and S5.

Loop functions (MSbar, dimensional regularization), the convention of
derive_alp_invariant_loop_crossing.py and memory/alp/loop-corrections/finite-mixing/power-counting-and-validity/crossing-and-validity/alp-invariant-loop-crossing-meson-poles.md:
    A0(m^2) = m^2/(16 pi^2) [1 - ln(m^2/mu^2)],
    Bbar(s; m1, m2) = -int_0^1 dz ln[(z m1^2 + (1-z) m2^2 - z(1-z) s - i0)/mu^2].
In this convention Bbar carries the whole mu dependence of the finite part; the full scalar
bubble is B0(s) = 1/eps_hat + Bbar(s), with the s-independent pole term 1/eps_hat =
2/(4-d) - gamma_E + ln(4 pi) removed by MSbar.  The amplitude bubble coefficients of S4/S5 are
those of amplitude_one_loop: sym V4(r,t,A,pi) V4(r,t,pi,pi)/(n0_r n0_t), units 1/(16 pi^2),
sym = 1/2 for r = t, one channel.

Mandelstam variables s_ij = (p_i + p_j)^2 of the three pions; on shell
s12 + s13 + s23 = m_A^2 + 3 m_pi^2.  e1, e2, e3 are the elementary symmetric polynomials of
(s12, s13, s23); e2', e3' those of the deviations s_ij - s0 from the symmetric point
s0 = (m_A^2 + 3 m_pi^2)/3.

Approximations: none beyond those of the imported module (LO U(3), dressed sources flavour
diagonal, CP-conserving sources in S4/S5, neutral internal lines only, truncation in epsilon at
AMPLITUDE_ORDER = 5).

Minimal invocation:
    .venv/bin/python scripts/verify_alp_loop_observable_selection.py
Tests:
    .venv/bin/python -m unittest discover -s scripts -p test_verify_alp_loop_observable_selection.py
"""

from __future__ import annotations

import itertools
import sys
import time
from functools import lru_cache
from pathlib import Path

import sympy as sp
from sympy.polys.polyfuncs import symmetrize

sys.path.insert(0, str(Path(__file__).resolve().parent))

import verify_alp_loop_observable_candidates as candidates  # noqa: E402

ALP_THRESHOLD_PAIRS = (('pi0', 'alp'), ('eta', 'alp'), ('eta_prime', 'alp'), ('alp', 'alp'))
ORDER = candidates.AMPLITUDE_ORDER

# ============================================================== kinematic symbols
alp_mass = sp.Symbol('m_A', positive=True)
pion_mass = sp.Symbol('m_pi', positive=True)
meson_mass = sp.Symbol('m_P', positive=True)
meson_mass_excess = sp.Symbol('delta_P', nonnegative=True)      # m_P = m_pi + delta_P
alp_mass_excess = sp.Symbol('delta_A', positive=True)           # m_A = 3 m_pi + delta_A (decay open)
first_mass, second_mass = sp.symbols('m_1 m_2', positive=True)
mandelstam = sp.symbols('s_12 s_13 s_23', real=True)
squared_momentum = sp.Symbol('s', real=True)
feynman_parameter = sp.Symbol('z', positive=True)
positive_parameter = sp.Symbol('t', positive=True)               # z = t/(1+t)
log_scale = sp.Symbol('L_mu', real=True)                         # ln(mu^2)
plane_sum = alp_mass**2 + 3*pion_mass**2                         # s12 + s13 + s23 on shell
symmetric_point = plane_sum/3


# ============================================================== S1: neutral kinetic term
def neutral_kinetic_term_check():
    """f^2/4 tr(d_mu U d_nu U^dag) - 1/4 tr(d_mu Phi d_nu Phi) for U = exp(i Phi/f),
    Phi = diag(phi_u, phi_d, phi_s) functions of two coordinates (x0, x1); exact in the fields.
    Returns the 2x2 matrix of differences (mu, nu in {0, 1}); bilinearity in the derivatives makes
    this sufficient for any metric contraction."""
    f = sp.Symbol('f', positive=True)
    coordinates = sp.symbols('x0 x1', real=True)
    fields = [sp.Function(f'phi_{q}', real=True)(*coordinates) for q in candidates.FLAVOURS]
    phi_matrix = sp.diag(*fields)
    goldstone = sp.diag(*[sp.exp(sp.I*field/f) for field in fields])   # exp of a diagonal matrix
    conjugate = sp.diag(*[sp.exp(-sp.I*field/f) for field in fields])
    assert sp.simplify(goldstone*conjugate - sp.eye(3)) == sp.zeros(3, 3)
    differences = sp.zeros(2, 2)
    for mu, nu in itertools.product(range(2), repeat=2):
        kinetic = f**2/4*(goldstone.diff(coordinates[mu])*conjugate.diff(coordinates[nu])).trace()
        free = (phi_matrix.diff(coordinates[mu])*phi_matrix.diff(coordinates[nu])).trace()/4
        differences[mu, nu] = sp.simplify(sp.expand(kinetic - free))
    return differences


@lru_cache(maxsize=None)
def builder_neutral_derivative_vertices():
    """Builder (dressed frame, neutral fields only, field degree <= 4): derivative vertex flags,
    number of derivative quartic vertices, and whether the kinetic matrix is epsilon independent;
    for CP-conserving and general sources."""
    result = {}
    for cp_conserving in (True, False):
        vertices = candidates.frame_vertices('dressed', None, 'example', cp_conserving)
        size = candidates.NEUTRAL_SIZE
        kinetic_constant = all(candidates.leading_power(candidates.truncate(
            vertices.kinetic[i][j] - candidates.epsilon_coefficient(vertices.kinetic[i][j], 0), 10)) is None
            for i in range(size) for j in range(size))
        result[cp_conserving] = {'flags': {k: list(v) for k, v in vertices.flags.items()},
                                 'derivative_quartic': len(vertices.derivative_quartic),
                                 'kinetic_constant': kinetic_constant}
    return result


# ============================================================== S2: contact amplitudes on shell
def scalar_product_symbols():
    """Symbols P_ij = p_i.p_j (1 <= i <= j <= 3) of the pion momenta."""
    return {(i, j): sp.Symbol(f'pp_{i}{j}', real=True) for i in range(1, 4) for j in range(i, 4)}


def scalar_product(first, second, products):
    """p_I.p_J for I, J in {0 (ALP, p = p1 + p2 + p3), 1, 2, 3} in terms of the pion products."""
    left = [1, 2, 3] if first == 0 else [first]
    right = [1, 2, 3] if second == 0 else [second]
    return sum(products[tuple(sorted((i, j)))] for i in left for j in right)


def on_shell_replacement(products):
    """p_i^2 = m_pi^2, p_i.p_j = (s_ij - 2 m_pi^2)/2, s23 eliminated with the on-shell plane."""
    s12, s13, _ = mandelstam
    s23 = plane_sum - s12 - s13
    values = {(1, 2): s12, (1, 3): s13, (2, 3): s23}
    replacement = {products[(i, i)]: pion_mass**2 for i in range(1, 4)}
    replacement.update({products[key]: (value - 2*pion_mass**2)/2 for key, value in values.items()})
    return replacement


def general_symmetric_contact_amplitude(degree):
    """Most general polynomial of degree <= `degree` in the ten products p_I.p_J
    (I <= J in {0, 1, 2, 3}) with symbolic coefficients c_k, symmetrized over the pion labels."""
    products = scalar_product_symbols()
    pairs = [(i, j) for i in range(4) for j in range(i, 4)]
    monomials = [()]
    for d in range(1, degree + 1):
        monomials += list(itertools.combinations_with_replacement(pairs, d))
    coefficients = sp.symbols(f'c_0:{len(monomials)}', real=True)
    total = 0
    for permutation in itertools.permutations((1, 2, 3)):
        relabel = {0: 0, 1: permutation[0], 2: permutation[1], 3: permutation[2]}
        for coefficient, monomial in zip(coefficients, monomials):
            term = coefficient
            for first, second in monomial:
                term *= scalar_product(relabel[first], relabel[second], products)
            total += term
    return sp.expand(total/6), products


def contact_structure_check():
    """(a) degree <= 1: derivatives of the on-shell amplitude w.r.t. s12, s13;
    (b) degree <= 2: on-shell amplitude minus alpha - beta sum s_ij^2 (beta from the s12^2
        coefficient), must be free of s12, s13;
    (c) on the plane: sum s_ij^2 and sum s_ij^3 in terms of e1, e2, e3, and the Jacobian
        d(e2, e3)/d(s12, s13) on the plane."""
    s12, s13, s23 = mandelstam
    result = {}
    linear, products = general_symmetric_contact_amplitude(1)
    linear_on_shell = sp.expand(linear.xreplace(on_shell_replacement(products)))
    result['linear_derivatives'] = (sp.diff(linear_on_shell, s12), sp.diff(linear_on_shell, s13))
    quadratic, products = general_symmetric_contact_amplitude(2)
    quadratic_on_shell = sp.expand(quadratic.xreplace(on_shell_replacement(products)))
    sum_squares = sp.expand((s12**2 + s13**2 + (plane_sum - s12 - s13)**2))
    beta = sp.Poly(quadratic_on_shell, s12, s13).coeff_monomial(s12**2)/2
    remainder = sp.expand(quadratic_on_shell - beta*sum_squares)
    result['quadratic_beta'] = sp.factor(beta)
    result['quadratic_remainder_derivatives'] = (sp.expand(sp.diff(remainder, s12)),
                                                 sp.expand(sp.diff(remainder, s13)))
    elementary = sp.symbols('e1 e2 e3')
    for power in (2, 3):
        expression, rest, definitions = symmetrize(s12**power + s13**power + s23**power,
                                                   s12, s13, s23, formal=True, symbols=elementary)
        assert rest == 0
        result[f'power_sum_{power}'] = sp.expand(expression)
    on_plane = {s23: plane_sum - s12 - s13}
    e2 = (s12*s13 + s12*s23 + s13*s23).xreplace(on_plane)
    e3 = (s12*s13*s23).xreplace(on_plane)
    result['jacobian_e2_e3'] = sp.factor(sp.Matrix([[sp.diff(e2, s12), sp.diff(e2, s13)],
                                                    [sp.diff(e3, s12), sp.diff(e3, s13)]]).det())
    result['elementary_symbols'] = elementary
    return result


# ============================================================== S3: Bbar on the Dalitz plane
def bubble_integrand(s_value, mass_one, mass_two, squared_scale=None):
    """Integrand of Bbar(s; m1, m2) in the Feynman parameter z (real part below threshold)."""
    scale = sp.exp(log_scale) if squared_scale is None else squared_scale
    z = feynman_parameter
    return -sp.log((z*mass_one**2 + (1 - z)*mass_two**2 - z*(1 - z)*s_value)/scale)


def bubble_denominator_identity():
    """z m1^2 + (1-z) m2^2 - z(1-z) s - [(z m1 - (1-z) m2)^2 + z(1-z)((m1+m2)^2 - s)] (must be 0):
    the argument of the logarithm is positive on 0 < z < 1 for s < (m1+m2)^2."""
    z, s = feynman_parameter, squared_momentum
    argument = z*first_mass**2 + (1 - z)*second_mass**2 - z*(1 - z)*s
    decomposed = (z*first_mass - (1 - z)*second_mass)**2 + z*(1 - z)*((first_mass + second_mass)**2 - s)
    return sp.expand(argument - decomposed)


def positive_parametrization(integrand):
    """integrand(z) dz with z = t/(1+t), dz = dt/(1+t)^2, t in (0, oo)."""
    t = positive_parameter
    return sp.factor(sp.together(integrand.subs(feynman_parameter, t/(1 + t))/(1 + t)**2))


def bubble_taylor_coefficients_at_zero(maximum_power=3):
    """Coefficients b_n of s^n in Bbar(s; m1, m2) about s = 0, as integrands in z (from the series
    of the integrand) and in t; equal-mass closed forms; positivity of the t integrand."""
    integrand = bubble_integrand(squared_momentum, first_mass, second_mass)
    series = sp.series(integrand, squared_momentum, 0, maximum_power + 1).removeO()
    z = feynman_parameter
    denominator = z*first_mass**2 + (1 - z)*second_mass**2
    equal = sp.Symbol('m', positive=True)
    result = {}
    for n in range(2, maximum_power + 1):
        coefficient = sp.together(series.coeff(squared_momentum, n))
        expected = (z*(1 - z))**n/(n*denominator**n)
        t_form = positive_parametrization(coefficient)
        result[n] = {'integrand_minus_expected': sp.simplify(coefficient - expected),
                     't_integrand': t_form, 't_integrand_positive': t_form.is_positive,
                     'equal_mass': sp.integrate(coefficient.subs({first_mass: equal, second_mass: equal}),
                                                (z, 0, 1))}
    return result


def generic_symmetric_expansion(maximum_degree=6):
    """sum_ij g(s_ij) about s0 with symbolic Taylor coefficients g_k (k <= maximum_degree), written
    in e2', e3' (e1' = 0 on the plane): returns (expression, g symbols, e2', e3')."""
    deviations = sp.symbols('d_12 d_13 d_23', real=True)
    taylor = sp.symbols(f'g_0:{maximum_degree + 1}', real=True)
    total = sum(taylor[k]*sum(d**k for d in deviations)/sp.factorial(k) for k in range(maximum_degree + 1))
    symbols = sp.symbols('e1p e2p e3p')
    expression, rest, _ = symmetrize(sp.expand(total), *deviations, formal=True, symbols=symbols)
    assert rest == 0
    return sp.expand(expression.subs(symbols[0], 0)), taylor, symbols[1], symbols[2]


def e3_coefficient_from_expansion(maximum_degree=6):
    """Coefficient of e3' (at e2' = 0, first order in e3') in sum_ij g(s_ij); all Taylor orders
    k <= maximum_degree are included, so the result shows that only g''' contributes."""
    expression, taylor, e2p, e3p = generic_symmetric_expansion(maximum_degree)
    return sp.expand(expression.coeff(e3p, 1).subs(e2p, 0)), taylor


def bubble_third_derivative_check():
    """Third s derivative of the Bbar integrand, positivity of its t form at s0 for the mass pairs
    (m_A, m_P) with m_P = m_pi + delta_P and (m_A, m_A) with m_A = 3 m_pi + delta_A; the distance
    Delta = (m1+m2)^2 - s0 to the threshold in both cases.  With the decomposition of
    bubble_denominator_identity and z = t/(1+t) the integrand of Bbar'''(s0) times dz/dt is
    2 t^3 / {[(t m1 - m2)^2 + t Delta]^3 (1+t)^2}, positive for t > 0 when Delta > 0, so
    Bbar'''(s0) > 0 (s0 below threshold, Bbar analytic there)."""
    integrand = bubble_integrand(squared_momentum, first_mass, second_mass)
    third = sp.together(sp.diff(integrand, squared_momentum, 3))
    z = feynman_parameter
    expected = 2*(z*(1 - z))**3/(z*first_mass**2 + (1 - z)*second_mass**2 - z*(1 - z)*squared_momentum)**3
    cases = {'alp_meson': ({first_mass: alp_mass, second_mass: pion_mass + meson_mass_excess},
                           {}),
             'alp_alp': ({first_mass: alp_mass, second_mass: alp_mass},
                         {alp_mass: 3*pion_mass + alp_mass_excess})}
    result = {'third_minus_expected': sp.simplify(third - expected)}
    for label, (masses, kinematics) in cases.items():
        threshold_distance = sp.expand(((masses[first_mass] + masses[second_mass])**2
                                        - symmetric_point).xreplace(kinematics))
        # rewrite the denominator with the decomposition of bubble_denominator_identity
        squared_difference = (z*masses[first_mass] - (1 - z)*masses[second_mass])**2
        distance_symbol = sp.Symbol('Delta', positive=True)
        positive_form = 2*(z*(1 - z))**3/(squared_difference + z*(1 - z)*distance_symbol)**3
        t = positive_parameter
        manifest_form = 2*t**3/(((t*masses[first_mass] - masses[second_mass])**2 + t*distance_symbol)**3
                                * (1 + t)**2)
        result[label] = {'threshold_distance': threshold_distance,
                         'threshold_distance_positive': threshold_distance.is_positive,
                         't_form_minus_manifest': sp.simplify(positive_parametrization(positive_form)
                                                              - manifest_form),
                         't_integrand_positive': manifest_form.is_positive,
                         'positive_form_minus_third': sp.simplify(
                             (positive_form.subs(distance_symbol, (masses[first_mass] + masses[second_mass])**2
                                                 - symmetric_point)
                              - expected.subs(masses).subs(squared_momentum, symmetric_point)))}
    return result


def subtraction_scale_dependence():
    """d Bbar/d ln mu^2 (integrand and integral over z)."""
    integrand = bubble_integrand(squared_momentum, first_mass, second_mass)
    derivative = sp.simplify(sp.expand_log(sp.diff(integrand, log_scale), force=True))
    return {'integrand': derivative, 'integral': sp.integrate(derivative, (feynman_parameter, 0, 1))}


# ============================================================== S4: A-threshold bubbles
def isovector_symbols():
    """Pairs (X_u, X_d) of every invariant with flavour components (S_hat_n, P_hat_n)."""
    pairs = []
    for n in range(2, candidates.MAXIMUM_SOURCE_ORDER + 1):
        pairs.append(candidates.scalar_invariant(n)[:2])
    for n in range(1, candidates.MAXIMUM_SOURCE_ORDER + 1):
        pairs.append(candidates.pseudoscalar_invariant(n)[:2])
    return pairs


def equate_up_and_down(series, invariant_pairs):
    """Set X_u = X_d for the given pairs (X_u, X_d)."""
    g = candidates.generator
    return series.compose([(g(up), g(down)) for up, down in invariant_pairs])


def series_symbols(polynomial):
    """Invariant symbols on which a ring element depends."""
    return sorted({candidates.RING.symbols[i] for monomial in polynomial.keys()
                   for i, exponent in enumerate(monomial) if exponent}
                  - {candidates.epsilon}, key=str)


@lru_cache(maxsize=None)
def alp_threshold_bubbles():
    """Leading power, leading coefficient and its symbols for each A-threshold bubble; the same
    after P_hat_1u = P_hat_1d, after also P_hat_3u = P_hat_3d, and in the isospin limit; the tree
    amplitude under the same conditions.  Caveat: with the test inputs m_u = m_d the tree amplitude
    is zero through epsilon^5 already for P_hat_1u = P_hat_1d (the only isospin breaking left
    comes from the invariants); with m_u != m_d it would not be."""
    one_loop = candidates.amplitude_one_loop(ORDER)
    tree = candidates.tree_amplitude(ORDER)
    p1 = candidates.pseudoscalar_invariant(1)[:2]
    p3 = candidates.pseudoscalar_invariant(3)[:2]
    conditions = {'P1_isovector_zero': [p1], 'P1_P3_isovector_zero': [p1, p3],
                  'isospin_limit': isovector_symbols()}
    result = {'tree': {'full': tree}, 'bubbles': {}}
    for name, pairs in conditions.items():
        result['tree'][name] = equate_up_and_down(tree, pairs)
    for key in ALP_THRESHOLD_PAIRS:
        series = one_loop['bubble'][key]
        power = candidates.leading_power(series)
        coefficient = candidates.epsilon_coefficient(series, power)
        entry = {'series': series, 'leading_power': power, 'leading_coefficient': coefficient,
                 'symbols': series_symbols(coefficient)}
        for name, pairs in conditions.items():
            reduced = equate_up_and_down(series, pairs)
            reduced_power = candidates.leading_power(reduced)
            entry[name] = {'series': reduced,
                           'coefficient_at_original_power': candidates.epsilon_coefficient(reduced, power),
                           'leading_power': reduced_power,
                           'leading_coefficient': None if reduced_power is None
                           else candidates.epsilon_coefficient(reduced, reduced_power)}
        result['bubbles'][key] = entry
    return result


def isospin_form(polynomial):
    """Rewrite a ring element with X_u = (X_plus + X_minus)/2, X_d = (X_plus - X_minus)/2 for each
    flavour invariant (X_minus = X_u - X_d isovector, X_plus = X_u + X_d) and factor.  Only
    meaningful for elements in which X_u and X_d were not equated beforehand."""
    expression = polynomial.as_expr()
    replacement = {}
    for up, down in isovector_symbols():
        stem = str(up)[:-2]
        plus, minus = sp.Symbol(f'{stem}_plus', real=True), sp.Symbol(f'{stem}_minus', real=True)
        replacement[up] = (plus + minus)/2
        replacement[down] = (plus - minus)/2
    return sp.factor(sp.expand(expression.xreplace(replacement)))


# ============================================================== S5: frame check
@lru_cache(maxsize=None)
def alp_threshold_bubbles_frame_check(kappa='symbolic', order=ORDER):
    """A-threshold bubble coefficients from the builder vertices and pole vectors in the dressed
    and kappa frames (CP-conserving sources), their difference through epsilon^order, and the
    difference between the dressed builder result and amplitude_one_loop."""
    unmixed = candidates.unmixed_neutral_poles(candidates.meson_sector_inputs())
    coefficients = {}
    for frame in ('dressed', 'kappa'):
        vertices = candidates.frame_vertices(frame, None, kappa, True)
        assert not vertices.cubic and not vertices.derivative_quartic
        size = candidates.NEUTRAL_SIZE
        kinetic = [row[:size] for row in vertices.kinetic[:size]]
        mass = [row[:size] for row in vertices.mass[:size]]
        poles = {}
        for name in ('pi0', 'eta', 'eta_prime', 'alp'):
            vector, value, norm = unmixed[name]
            _, series, norm_check = candidates.perturbative_pole(kinetic, mass, vector, value, order)
            assert norm_check == norm
            poles[name] = (candidates.pad(series), norm)
        alp, pion = poles['alp'][0], poles['pi0'][0]
        coefficients[frame] = {}
        for r, t in ALP_THRESHOLD_PAIRS:
            left = candidates.full_contraction(vertices.quartic, [poles[r][0], poles[t][0], alp, pion], order)
            right = candidates.full_contraction(vertices.quartic, [poles[r][0], poles[t][0], pion, pion], order)
            symmetry = sp.Rational(1, 2) if r == t else sp.Integer(1)
            coefficients[frame][(r, t)] = candidates.truncate(
                candidates.rational(symmetry/(poles[r][1]*poles[t][1]))*left*right, order)
    closed_form = candidates.amplitude_one_loop(order)['bubble']
    return {'coefficients': coefficients,
            'difference': {k: candidates.truncate(coefficients['kappa'][k] - coefficients['dressed'][k], order)
                           for k in ALP_THRESHOLD_PAIRS},
            'dressed_minus_closed_form': {k: candidates.truncate(coefficients['dressed'][k] - closed_form[k], order)
                                          for k in ALP_THRESHOLD_PAIRS}}


# ============================================================== summary
def main():
    start = time.time()
    print('S1 f^2/4 tr(dU dU^dag) - 1/4 tr(dPhi dPhi) =', list(neutral_kinetic_term_check()))
    builder = builder_neutral_derivative_vertices()
    print('   builder dressed neutral vertices (CP-even, general):',
          {k: (v['derivative_quartic'], v['kinetic_constant']) for k, v in builder.items()})
    contact = contact_structure_check()
    print('S2 (a) d/ds12, d/ds13 of degree<=1 contact amplitude:', contact['linear_derivatives'])
    print('   (b) degree<=2 minus beta sum s^2, d/ds12, d/ds13:', contact['quadratic_remainder_derivatives'])
    print('   (c) sum s^2 =', contact['power_sum_2'], '; sum s^3 =', contact['power_sum_3'],
          '; Jacobian d(e2,e3)/d(s12,s13) =', contact['jacobian_e2_e3'])
    taylor = bubble_taylor_coefficients_at_zero()
    print('S3 (a) Taylor coefficients of Bbar at s=0, positive integrand:',
          {n: v['t_integrand_positive'] for n, v in taylor.items()},
          '; equal mass:', {n: v['equal_mass'] for n, v in taylor.items()})
    e3_coefficient, _ = e3_coefficient_from_expansion()
    third = bubble_third_derivative_check()
    print('   (b) coefficient of e3\' in sum g(s_ij):', e3_coefficient, '; Bbar\'\'\'(s0) > 0:',
          {k: (third[k]['threshold_distance_positive'], third[k]['t_integrand_positive'])
           for k in ('alp_meson', 'alp_alp')})
    print('   (c) d Bbar/d ln mu^2 =', subtraction_scale_dependence()['integral'])
    bubbles = alp_threshold_bubbles()
    print('S4 tree leading power:', {k: candidates.leading_power(v) for k, v in bubbles['tree'].items()})
    for key, entry in bubbles['bubbles'].items():
        print(f'   {key}: power {entry["leading_power"]}, symbols {entry["symbols"]}')
        if key != ('alp', 'alp'):
            print('      leading coefficient (isospin form):', isospin_form(entry['leading_coefficient']))
        else:
            print('      leading coefficient: ', len(entry['leading_coefficient']), 'terms')
        for name in ('P1_isovector_zero', 'P1_P3_isovector_zero', 'isospin_limit'):
            reduced = entry[name]['coefficient_at_original_power']
            print(f'      {name}: coefficient at power {entry["leading_power"]} nonzero:', bool(reduced),
                  '; new leading power:', entry[name]['leading_power'])
            if reduced and key == ('eta', 'alp'):
                print('         =', sp.factor(reduced.as_expr()))
    frame = alp_threshold_bubbles_frame_check()
    print('S5 kappa - dressed (symbolic kappa) zero through eps^5:',
          {k: not v for k, v in frame['difference'].items()},
          '; dressed builder = closed form:', {k: not v for k, v in frame['dressed_minus_closed_form'].items()})
    print(f'runtime {time.time() - start:.1f} s')


if __name__ == '__main__':
    main()
