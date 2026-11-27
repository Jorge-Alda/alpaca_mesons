"""Constant axial redefinition of the loop-size estimate with the radial mode.

Run:
    .venv/bin/python scripts/derive_murayama_radial_frame_invariance.py
    PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_murayama_radial_frame_invariance

Conventions (memory/alp/radial-mode/couplings-and-decays/finite-mixing/murayama-sigma-two-pseudoscalar-decays.md):
partonic mass term -qbar_R M q_L - h.c., M = S + i P, expanded in the
dimensionless ratios y = sigma/f_a and x = a/f_a,
    S(x, y) = sum_{r,s} S_rs y^r x^s,   P(x, y) = sum_{r,s} P_rs y^r x^s,
first index counting sigma, second counting a.  Coefficients have mass
dimension one (GeV); x and y are dimensionless.  The polar spurion of the
radial record in this convention is m_q (1 + y)^N exp(-i N x); the paper's
M_q of DiLuzio:2026poh is its complex conjugate.

Redefinition: g_R = exp(i x kappa), g_L = g_R^dagger, with a constant
Hermitian axial generator kappa.  The redefinition depends on a only, not on
sigma.  Transport (memory/alp/chiral-transformations/source-invariants/finite-transformations/alp-general-chiral-redefinition-invariants.md):
    d/dx W_eta(x, y) = i C_eta(x, y) W_eta(x, y),   W_eta(0, y) = 1,
i.e. the transport is taken along x at fixed sigma.  Since g does not depend
on y, W'_eta = g_eta W_eta holds at every fixed y, and the dressed source
M_hat = W_R^dagger M W_L is unchanged.  Currents of d sigma, if present,
transform homogeneously under a sigma-independent g and need no transport.
Dressed topological coefficient:
    Theta_hat(x, y) = Theta(x, y) + int_0^x tr C_A(x', y) dx'.

The script does not evaluate any loop integral, loop function, loop
coefficient or loop amplitude.  It acts only on the invariant inputs and
on contractions of these inputs that enter the estimate.  All physical
parameters are symbols.
"""
import itertools

import sympy as sp

X_RATIO = sp.Symbol('x', real=True)        # a / f_a
Y_RATIO = sp.Symbol('y', real=True)        # sigma / f_a


# ---------------------------------------------------------------------------
# Truncated double series in y (radial) and x (ALP) with noncommuting
# coefficients, stored as {(r, s): coefficient of y^r x^s}.
# ---------------------------------------------------------------------------

def series_product(*factors, radial_order, alp_order):
    """Truncated product of double series, keeping the order of the factors."""
    result = factors[0]
    for factor in factors[1:]:
        product = {}
        for (r1, s1), left in result.items():
            for (r2, s2), right in factor.items():
                r, s = r1 + r2, s1 + s2
                if r <= radial_order and s <= alp_order:
                    product[(r, s)] = product.get((r, s), 0) + left * right
        result = {key: sp.expand(value) for key, value in product.items()}
    return result


def series_sum(*terms):
    total = {}
    for term in terms:
        for key, value in term.items():
            total[key] = total.get(key, 0) + value
    return {key: sp.expand(value) for key, value in total.items()}


def series_scale(series, factor):
    return {key: sp.expand(factor * value) for key, value in series.items()}


def adjoint(expression):
    """Hermitian conjugate; every noncommuting symbol used here is Hermitian."""
    expression = sp.expand(expression)
    result = 0
    for term in sp.Add.make_args(expression):
        commuting, ordered = term.args_cnc()
        result += sp.conjugate(sp.Mul(*commuting)) * sp.Mul(*reversed(ordered))
    return sp.expand(result)


def series_adjoint(series):
    return {key: adjoint(value) for key, value in series.items()}


def cyclic_trace(expression):
    """Flavour trace, with each word mapped to the smallest cyclic rotation."""
    result = 0
    for term in sp.Add.make_args(sp.expand(expression)):
        commuting, ordered = term.args_cnc()
        letters = []
        for factor in ordered:
            base, exponent = factor.as_base_exp()
            letters += [base.name] * int(exponent)
        if not letters:
            result += sp.Mul(*commuting) * sp.Symbol('N_f', positive=True)
            continue
        rotations = [tuple(letters[k:] + letters[:k]) for k in range(len(letters))]
        word = min(rotations)
        result += sp.Mul(*commuting) * sp.Symbol('tr(' + ' '.join(word) + ')')
    return sp.expand(result)


def axial_exponential(generator, sign, alp_order):
    """exp(sign i x generator) as a double series."""
    return {(0, s): sp.expand((sign * sp.I * generator) ** s / sp.factorial(s))
            for s in range(alp_order + 1)}


def transport(connection, radial_order, alp_order):
    """Solve d/dx W = i C W, W(0, y) = 1, order by order in x at fixed y."""
    transport_series = {(0, 0): sp.Integer(1)}
    for s in range(alp_order):
        for r in range(radial_order + 1):
            value = 0
            for (r1, s1), coefficient in connection.items():
                key = (r - r1, s - s1)
                if key in transport_series:
                    value += coefficient * transport_series[key]
            value = sp.expand(sp.I * value / (s + 1))
            if value != 0:
                transport_series[(r, s + 1)] = value
    return transport_series


def hermitian_split(series):
    scalar = {key: sp.expand((value + adjoint(value)) / 2) for key, value in series.items()}
    pseudoscalar = {key: sp.expand((value - adjoint(value)) / (2 * sp.I)) for key, value in series.items()}
    return scalar, pseudoscalar


def dressed_sources(mass, right_connection, left_connection, theta, radial_order, alp_order):
    """Dressed S_hat_rs, P_hat_rs and Theta_hat_rs of one frame."""
    right_transport = transport(right_connection, radial_order, alp_order)
    left_transport = transport(left_connection, radial_order, alp_order)
    dressed_mass = series_product(series_adjoint(right_transport), mass, left_transport,
                                  radial_order=radial_order, alp_order=alp_order)
    scalar, pseudoscalar = hermitian_split(dressed_mass)
    axial_connection = series_scale(series_sum(right_connection, series_scale(left_connection, -1)),
                                    sp.Rational(1, 2))
    dressed_theta = dict(theta)
    for (r, s), value in axial_connection.items():
        if s + 1 <= alp_order:
            dressed_theta[(r, s + 1)] = sp.expand(
                dressed_theta.get((r, s + 1), 0) + cyclic_trace(value) / (s + 1))
    return scalar, pseudoscalar, dressed_theta


def general_original_frame(radial_order, alp_order, with_connections):
    """Original frame with general noncommuting sources and connections."""
    mass = {}
    for r in range(radial_order + 1):
        for s in range(alp_order + 1):
            scalar = sp.Symbol(f'S_{r}{s}', commutative=False)
            pseudoscalar = sp.Symbol(f'P_{r}{s}', commutative=False)
            mass[(r, s)] = scalar + sp.I * pseudoscalar
    right_connection, left_connection = {}, {}
    if with_connections:
        for r in range(radial_order + 1):
            for s in range(alp_order):
                right_connection[(r, s)] = sp.Symbol(f'CR_{r}{s}', commutative=False)
                left_connection[(r, s)] = sp.Symbol(f'CL_{r}{s}', commutative=False)
    theta = {(r, s): sp.Symbol(f'Theta_{r}{s}', real=True)
             for r in range(radial_order + 1) for s in range(alp_order + 1)}
    return mass, right_connection, left_connection, theta


def transformed_frame(mass, right_connection, left_connection, theta, generator,
                      radial_order, alp_order):
    """Apply g_R = exp(i x kappa), g_L = g_R^dagger (constant axial generator)."""
    orders = dict(radial_order=radial_order, alp_order=alp_order)
    right_group = axial_exponential(generator, +1, alp_order)
    right_inverse = axial_exponential(generator, -1, alp_order)
    left_group, left_inverse = right_inverse, right_group
    new_mass = series_product(right_group, mass, left_inverse, **orders)
    # C' = g C g^dagger + i g d_x g^dagger; the second term is +kappa (R), -kappa (L).
    new_right = series_sum(series_product(right_group, right_connection, right_inverse, **orders)
                           if right_connection else {}, {(0, 0): generator})
    new_left = series_sum(series_product(left_group, left_connection, left_inverse, **orders)
                          if left_connection else {}, {(0, 0): -generator})
    new_right = {key: value for key, value in new_right.items() if key[1] < alp_order}
    new_left = {key: value for key, value in new_left.items() if key[1] < alp_order}
    # Anomaly: Theta' = Theta - (1/2i) log(det g_R / det g_L) = Theta - x tr kappa.
    new_theta = dict(theta)
    new_theta[(0, 1)] = sp.expand(new_theta.get((0, 1), 0) - cyclic_trace(generator))
    return new_mass, new_right, new_left, new_theta


def source_invariance(radial_order=2, alp_order=3, with_connections=True):
    """Part (a): residuals S_hat' - S_hat, P_hat' - P_hat, Theta_hat' - Theta_hat."""
    generator = sp.Symbol('kappa', commutative=False)
    original = general_original_frame(radial_order, alp_order, with_connections)
    transformed = transformed_frame(*original, generator, radial_order, alp_order)
    before = dressed_sources(*original, radial_order, alp_order)
    after = dressed_sources(*transformed, radial_order, alp_order)
    residuals = {}
    for name, old, new in zip(('S_hat', 'P_hat', 'Theta_hat'), before, after):
        for key in itertools.product(range(radial_order + 1), range(alp_order + 1)):
            residuals[(name, key)] = sp.expand(new.get(key, 0) - old.get(key, 0))
    return residuals


# ---------------------------------------------------------------------------
# Parts (b) and (c): contractions and pole data, flavour-diagonal sector.
# ---------------------------------------------------------------------------

FLAVOURS = ('u', 'd', 's')


def diagonal_frames(alp_order=2, radial_order=2):
    """Polar spurion with diagonal charges N and a diagonal axial generator.

    Returns per flavour the dressed coefficients of the original frame and the
    dressed coefficients of the transformed frame, plus the topological
    coefficients.  The sources of the transformed frame are used only as the
    input of the transport that produces its dressed coefficients.
    With commuting entries the transported dressing cancels the redefinition
    directly; the noncommuting statement is the check of source_invariance.
    This sector is used to build the pole operator inverse explicitly."""
    masses = sp.symbols('m_u m_d m_s', positive=True)
    charges = sp.symbols('n_u n_d n_s', real=True)
    kappas = sp.symbols('kappa_u kappa_d kappa_s', real=True)
    theta_one = sp.Symbol(r'\Theta_{01}', real=True)
    x, y = X_RATIO, Y_RATIO

    def truncate(expression):
        expression = sp.expand(expression)
        return sum(expression.coeff(y, r).coeff(x, s) * y**r * x**s
                   for r in range(radial_order + 1) for s in range(alp_order + 1))

    def series_exp(argument):
        return truncate(sum(argument**k / sp.factorial(k) for k in range(alp_order + 1)))

    def binomial_series(charge):
        return sum(sp.binomial(charge, r) * y**r for r in range(radial_order + 1))

    def coefficients(source):
        source = sp.expand(source)
        scalar, pseudoscalar = {}, {}
        for r in range(radial_order + 1):
            for s in range(alp_order + 1):
                value = sp.expand(sp.expand_func(source.coeff(y, r).coeff(x, s)))
                scalar[(r, s)] = sp.re(value)
                pseudoscalar[(r, s)] = sp.im(value)
        return scalar, pseudoscalar

    original, dressed_after = {}, {}
    for mass, charge, kappa, flavour in zip(masses, charges, kappas, FLAVOURS):
        source = truncate(mass * binomial_series(charge) * series_exp(-sp.I * charge * x))
        framed = truncate(series_exp(sp.I * kappa * x) * source * series_exp(sp.I * kappa * x))
        # W_R' = exp(i kappa x), W_L' = exp(-i kappa x): dressed = W_R'^dag M' W_L'.
        redressed = truncate(series_exp(-sp.I * kappa * x) * framed * series_exp(-sp.I * kappa * x))
        original[flavour] = coefficients(source)
        dressed_after[flavour] = coefficients(redressed)
    trace_kappa = sum(kappas)
    # Theta' = Theta_01 - tr kappa (anomaly) plus the transport term + tr kappa.
    theta = {'original': theta_one, 'dressed_after': theta_one - trace_kappa + trace_kappa}
    return {'original': original, 'dressed_after': dressed_after,
            'theta': theta, 'masses': masses, 'charges': charges, 'kappas': kappas}


B0, F_PI, F_A = sp.symbols('B_0 f_pi f_a', positive=True)
SINGLET_MASS_SQUARED = sp.Symbol('m_0^2', positive=True)
QUARTIC = sp.Symbol('lambda', positive=True)
POLE_I, POLE_J = sp.symbols('m_i^2 m_j^2', positive=True)
NUMBER_OF_FLAVOURS = 3


def pole_inverse(scalar_zero, source, pole_mass_squared):
    """Diagonal solution Y of M_s(Y) = source with
    M_s(Y) = B_0 {S_hat_00, Y} + (m_0^2/N_f) tr(Y) 1 - s Y."""
    unknowns = sp.symbols('Y_u Y_d Y_s')
    equations = [2 * B0 * scalar_zero[q] * unknowns[k]
                 + SINGLET_MASS_SQUARED / NUMBER_OF_FLAVOURS * sum(unknowns)
                 - pole_mass_squared * unknowns[k] - source[q]
                 for k, q in enumerate(FLAVOURS)]
    solution = sp.solve(equations, unknowns, dict=True)[0]
    return {q: sp.together(solution[unknowns[k]]) for k, q in enumerate(FLAVOURS)}


def estimate_inputs(coefficients, theta_one):
    """Pole data and contractions used by the estimate, from one set of coefficients."""
    def entry(kind, r, s):
        return {q: coefficients[q][0 if kind == 'S' else 1][(r, s)] for q in FLAVOURS}

    def trace(*factors):
        return sum(sp.Mul(*[factor[q] for factor in factors]) for q in FLAVOURS)

    scalar_zero = entry('S', 0, 0)
    pseudoscalar_one = entry('P', 0, 1)
    epsilon = F_PI / F_A
    tadpole_source = {q: B0 * pseudoscalar_one[q]
                      - SINGLET_MASS_SQUARED / NUMBER_OF_FLAVOURS * theta_one for q in FLAVOURS}
    # (c) pole data: M_s^{-1}(C) and N(s) at two pole masses.
    inverse_i = pole_inverse(scalar_zero, tadpole_source, POLE_I)
    inverse_j = pole_inverse(scalar_zero, tadpole_source, POLE_J)
    normalization_i = 1 + 2 * epsilon**2 * trace(inverse_i, inverse_i)
    # sigma-record form of the pole meson matrix, B_0 M_s^{-1}(P_hat_01).
    radial_i = {q: B0 * value for q, value in pole_inverse(scalar_zero, pseudoscalar_one, POLE_I).items()}
    radial_j = {q: B0 * value for q, value in pole_inverse(scalar_zero, pseudoscalar_one, POLE_J).items()}
    radial_normalization_i = 1 + 2 * epsilon**2 * trace(radial_i, radial_i)
    radial_normalization_j = 1 + 2 * epsilon**2 * trace(radial_j, radial_j)
    radial_mass_squared = 2 * QUARTIC * F_A**2 + B0 * epsilon**2 * (
        trace(entry('S', 1, 0)) - 2 * trace(entry('S', 2, 0)))
    anticommutator = {q: 2 * radial_i[q] * radial_j[q] for q in FLAVOURS}
    pair_sum = {q: radial_i[q] + radial_j[q] for q in FLAVOURS}

    def pair_contraction(scalar_two_meson, pseudoscalar_mixed, scalar_two_alp, kinetic):
        # Structure of f_a sqrt(N_i N_j) g_ij (radial record), for the given coefficients.
        return (2 * B0 * epsilon**2 * trace(scalar_two_meson, anticommutator)
                - 2 * B0 * epsilon**2 * trace(pseudoscalar_mixed, pair_sum)
                - 2 * B0 * epsilon**2 * trace(scalar_two_alp) + kinetic)

    radial_coupling = pair_contraction(entry('S', 1, 0), entry('P', 1, 1), entry('S', 1, 2),
                                       radial_mass_squared - POLE_I - POLE_J) / (
        F_A * sp.sqrt(radial_normalization_i * radial_normalization_j))
    return {
        'tadpole source C': tadpole_source,
        'M^{-1}(C) at m_i^2': inverse_i,
        'N(m_i^2)': normalization_i,
        'ALP pole term tr(C M^{-1}(C)) at m_i^2': trace(tadpole_source, inverse_i),
        'tr S_hat_02': trace(entry('S', 0, 2)),
        'm_sigma^2': radial_mass_squared,
        'sigma-PP coupling g_ij': radial_coupling,
        'sigma-a-P vertex P_hat_11': entry('P', 1, 1),
        'sigma-a-P contraction tr(P_hat_11 Omega_i)': trace(entry('P', 1, 1), radial_i),
        'R2: tr(S_hat_20 {Omega_i, Omega_j})': trace(entry('S', 2, 0), anticommutator),
        'R2: tr(P_hat_21 (Omega_i + Omega_j))': trace(entry('P', 2, 1), pair_sum),
        'R2: tr S_hat_22': trace(entry('S', 2, 2)),
        'R2: tr S_hat_20 (sigma^2 a^0 vertex)': trace(entry('S', 2, 0)),
    }


def sample_point(expression, shift=0):
    """Deterministic rational values for every symbol, used only to show that
    a difference is nonzero."""
    symbols = sorted(expression.free_symbols, key=lambda symbol: symbol.name)
    return {symbol: sp.Rational(3 + 2 * index + shift, 7 + index)
            for index, symbol in enumerate(symbols)}


def is_zero(value):
    """A difference is zero if it vanishes identically (exact symbolic test).
    Before that, it is evaluated at rational points with 60 significant
    digits; a value above 1e-30 shows that it is nonzero."""
    if isinstance(value, dict):
        return all(is_zero(entry) for entry in value.values())
    if value == 0:
        return True
    for shift in range(3):
        evaluated = sp.N(value.subs(sample_point(value, shift)), 60)
        if evaluated.is_finite and abs(evaluated) > sp.Float('1e-30'):
            return False
    return sp.cancel(sp.together(sp.expand(value))) == 0


def difference(first, second):
    """Difference without expansion; identical expressions cancel on construction."""
    if isinstance(first, dict):
        return {key: difference(first[key], second[key]) for key in first}
    return first - second


def contraction_checks():
    """Parts (b) and (c): dressed transformed frame minus original frame."""
    frames = diagonal_frames(alp_order=2, radial_order=2)
    original = estimate_inputs(frames['original'], frames['theta']['original'])
    dressed = estimate_inputs(frames['dressed_after'], frames['theta']['dressed_after'])
    return {name: difference(dressed[name], original[name]) for name in original}


def main():
    residuals = source_invariance(radial_order=2, alp_order=3, with_connections=True)
    nonzero = {key: value for key, value in residuals.items() if value != 0}
    print('(a) General noncommuting sources and connections, y^r x^s with r <= 2, s <= 3:')
    print(f'    {len(residuals)} residuals checked, nonzero: {nonzero if nonzero else "none"}')
    residuals = source_invariance(radial_order=2, alp_order=6, with_connections=False)
    nonzero = {key: value for key, value in residuals.items() if value != 0}
    print('(a) Connection-free polar representative, r <= 2, s <= 6:')
    print(f'    {len(residuals)} residuals checked, nonzero: {nonzero if nonzero else "none"}')
    dressed_differences = contraction_checks()
    print('\n(b), (c) Dressed transformed frame minus original frame (diagonal sector):')
    for name, value in dressed_differences.items():
        print(f'    {name}: {"0" if is_zero(value) else value}')


if __name__ == '__main__':
    main()
