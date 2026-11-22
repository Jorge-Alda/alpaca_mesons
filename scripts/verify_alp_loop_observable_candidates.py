"""Exact leading powers of epsilon = f_pi/f_a for four candidate observables (OpenSpec change
estimate-loop-tree-crossing-finite-mixing, task 1.1): tree term, each one-loop class, and the
behaviour under the constant axial redefinition g_R = g_L^dag = exp(i x kappa).

Candidates
    C1  meson pole mass m_n^2 (pi^0, eta, eta' and the charged pairs ud, us, ds);
    C2  ALP-meson mixing entry D_{pi0 a}(s) of the inverse propagator and of the one-loop
        self-energy;
    C3  ALP-like pole mass m_A^2;
    C4  on-shell amplitude A -> pi^0 pi^0 pi^0, A the ALP-like pole.

Conventions (Alda:2026nrd, and the scripts derive_alp_invariant_loop_crossing.py and
derive_alp_loop_power_counting.py).  Units GeV.  U = exp(i Phi/f_pi), tr(lambda_i lambda_j) =
2 delta_ij, N_f = 3, x = a/f_a, epsilon = f_pi/f_a (f_pi fixed, f_a = f_pi/epsilon).  Dressed
frame, LO U(3):
    L = f^2/4 <dU dU^dag> + f^2 B0/2 <M_hat(x) U^dag + U M_hat(x)^dag> - f^2 m0^2/(4 N_f) X_TC^2
        + (da)^2/2 - m_a0^2 a^2/2,       X_TC = -tr Phi/f - 2 Theta_hat(x),
    M_hat(x) = m_q + sum_{n>=1} x^n/n! (S_hat_n + i P_hat_n),  Theta_hat(x) = sum_n x^n/n! Theta_hat_n,
with flavour-diagonal invariants (S_hat_1 = 0, P_hat_0 = 0).  The diagonal entries of S_hat_n,
P_hat_n and the numbers Theta_hat_n are sympy symbols throughout: no value is assumed for any of
them.  The kappa frame is the change of variables U_hat = g^dag U g^dag of the existing builder.

Neutral fields are N_q = Phi_qq = sqrt2 P_q (as in build_lagrangian_exact), so K = diag(1/2,1/2,
1/2,1) in the dressed neutral sector.  Pole vectors u are normalized to u^T K u = n0, with n0 the
rational norm of the unmixed vector (n0 = 1 for pi^0 = lambda_3 and for the ALP); the physical
K-normalized vector is u/sqrt(n0), so every quadratic quantity is rational.  Sign conventions:
the pi^0 vector starts as lambda_3, the ALP-like vector has ALP component +1 + O(epsilon^2).

Method.  Poles and pole vectors are exact truncated power series in epsilon from Rayleigh-
Schroedinger perturbation theory of the pencil M^2(epsilon) - s K(epsilon) with the rational
reduced resolvent of the unperturbed pencil (exact rational arithmetic, polynomial ring QQ[...]).
One-loop coefficients are those of pole_shift_coefficients: per loop function, contact tadpole
plus loop-induced one-point term (coefficient of A0(m_r^2)), and bubbles -mult V^2/(32 pi^2)
(coefficient of Bbar(s_n; m_r, m_t)); bubble coefficients are given in units of 1/(32 pi^2) for
two-point functions and 1/(16 pi^2) for the amplitude.  A0(m^2) = m^2/(16 pi^2)[1 - ln(m^2/mu^2)],
<phi phi> = -A0.  Two-point vertices come from the existing builder (fourth order in the fields,
all ten fields, momentum-dependent vertices exact); the amplitude uses the closed-form neutral
dressed Lagrangian to sixth order in the fields with general Taylor coefficients.

"Leading power p": all coefficients below epsilon^p vanish identically (symbolically) and the
epsilon^p coefficient is a nonzero polynomial in the invariants.  "zero to order N": all
coefficients through epsilon^N vanish identically.

Meson-sector test inputs (rational, verification only; not physics inputs): m_u = m_d, so that
pi^0 is exactly lambda_3; m_s and m0^2 chosen so that the unmixed eta, eta' masses are rational.
The frame generator kappa is symbolic in the tree-level checks and EXAMPLE_FRAME_GENERATOR of
derive_alp_invariant_loop_crossing.py (rational) in the one-loop checks.

Minimal invocation:
    .venv/bin/python scripts/verify_alp_loop_observable_candidates.py
Tests:
    .venv/bin/python -m unittest discover -s scripts -p test_verify_alp_loop_observable_candidates.py
"""

from __future__ import annotations

import itertools
import sys
import time
from functools import lru_cache
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ
from sympy.combinatorics import Permutation
from sympy.polys.rings import ring

sys.path.insert(0, str(Path(__file__).resolve().parent))

import derive_alp_invariant_loop_crossing as crossing  # noqa: E402

FLAVOURS = ('u', 'd', 's')
NUMBER_OF_FLAVOURS = 3
NEUTRAL_SIZE = crossing.NEUTRAL_SECTOR_SIZE          # N_u, N_d, N_s, a
ALP_INDEX = crossing.ALP_INDEX
FIELD_COUNT = len(crossing.FIELD_NAMES)
CHARGED_PAIRS = crossing.CHARGED_PAIRS               # {'ud': (4, 5), 'us': (6, 7), 'ds': (8, 9)}
MAXIMUM_SOURCE_ORDER = 6
TWO_POINT_ORDER = 4                                  # epsilon order of the two-point series
AMPLITUDE_ORDER = 5                                  # epsilon order of the amplitude series

# ============================================================== symbols and polynomial ring
epsilon = sp.Symbol('epsilon', positive=True)
squared_momentum = sp.Symbol('s', real=True)
kappa_symbols = tuple(sp.Symbol(f'kappa_{q}', real=True) for q in FLAVOURS)
bare_alp_mass_symbol = sp.Symbol('m_a0_squared', positive=True)


def scalar_invariant(order):
    """Diagonal entries of S_hat_order."""
    return tuple(sp.Symbol(f'S_hat_{order}_{q}', real=True) for q in FLAVOURS)


def pseudoscalar_invariant(order):
    """Diagonal entries of P_hat_order."""
    return tuple(sp.Symbol(f'P_hat_{order}_{q}', real=True) for q in FLAVOURS)


def topological_invariant(order):
    """Theta_hat_order."""
    return sp.Symbol(f'Theta_hat_{order}', real=True)


ALL_INVARIANTS = ([v for n in range(2, MAXIMUM_SOURCE_ORDER + 1) for v in scalar_invariant(n)]
                  + [v for n in range(1, MAXIMUM_SOURCE_ORDER + 1) for v in pseudoscalar_invariant(n)]
                  + [topological_invariant(n) for n in range(1, MAXIMUM_SOURCE_ORDER + 1)])
CP_ODD_INVARIANTS = ([v for n in range(3, MAXIMUM_SOURCE_ORDER + 1, 2) for v in scalar_invariant(n)]
                     + [v for n in range(2, MAXIMUM_SOURCE_ORDER + 1, 2) for v in pseudoscalar_invariant(n)]
                     + [topological_invariant(n) for n in range(2, MAXIMUM_SOURCE_ORDER + 1, 2)])
neutral_field_symbols = (sp.Symbol('N_u', real=True), sp.Symbol('N_d', real=True),
                         sp.Symbol('N_s', real=True), crossing.meson_fields[ALP_INDEX])

RING, *_GENERATORS = ring([epsilon, squared_momentum, *kappa_symbols, bare_alp_mass_symbol,
                           *ALL_INVARIANTS, *neutral_field_symbols], QQ)
GENERATOR_INDEX = {symbol: index for index, symbol in enumerate(RING.symbols)}
FIELD_GENERATOR_INDICES = tuple(GENERATOR_INDEX[v] for v in neutral_field_symbols)
ZERO = RING(0)
ONE = RING(1)


def generator(symbol):
    return RING.gens[GENERATOR_INDEX[symbol]]


def to_ring(expression):
    return RING(sp.expand(expression))


def rational(value):
    """sympy Rational -> ring constant."""
    value = sp.Rational(value)
    return RING(QQ(int(value.p), int(value.q)))


def truncate(polynomial, order):
    """Drop the terms of order > `order` in epsilon."""
    return RING({monomial: c for monomial, c in polynomial.items() if monomial[0] <= order})


def epsilon_coefficient(polynomial, power):
    """Coefficient of epsilon^power (a ring element without epsilon)."""
    return RING({(0,) + monomial[1:]: c for monomial, c in polynomial.items() if monomial[0] == power})


def leading_power(polynomial):
    """Smallest power of epsilon with a nonzero coefficient; None if the series is zero."""
    if not polynomial:
        return None
    return min(monomial[0] for monomial in polynomial.keys())


def set_to_zero(polynomial, symbols):
    """Set the given generators to zero."""
    indices = [GENERATOR_INDEX[v] for v in symbols]
    return RING({m: c for m, c in polynomial.items() if all(m[i] == 0 for i in indices)})


def substitute(polynomial, symbol, value, order):
    """Replace a generator by a ring element and truncate in epsilon."""
    return truncate(polynomial.compose(generator(symbol), value), order)


def constant_value(polynomial):
    """Rational value of a ring constant (asserts that it has no generators)."""
    if not polynomial:
        return sp.Integer(0)
    assert list(polynomial.keys()) == [(0,)*RING.ngens], polynomial
    value = polynomial.LC
    return sp.Rational(int(value.numerator), int(value.denominator))


def series_inverse(polynomial, order):
    """1/p as a truncated series, for p with a nonzero rational epsilon^0 term."""
    leading = constant_value(epsilon_coefficient(polynomial, 0))
    assert leading != 0
    remainder = truncate(ONE - polynomial*rational(1/leading), order)
    result, power = ONE, ONE
    for _ in range(order):
        power = truncate(power*remainder, order)
        result += power
    return truncate(result*rational(1/leading), order)


def series_inverse_square_root(polynomial, order):
    """p^(-1/2) for p = 1 + O(epsilon)."""
    assert constant_value(epsilon_coefficient(polynomial, 0)) == 1
    deviation = truncate(polynomial - ONE, order)
    result, power = ONE, ONE
    for m in range(1, order + 1):
        power = truncate(power*deviation, order)
        result += rational(sp.binomial(sp.Rational(-1, 2), m))*power
    return truncate(result, order)


# ============================================================== test inputs
def meson_sector_inputs():
    """Rational meson-sector test values (verification only).  m_u = m_d; 2 B0 m_u = 1/50,
    2 B0 m_s = 199/450 and m0^2/3 = 1/4 give rational unmixed neutral masses 1/50, 109/450,
    97/100 (pi^0, eta, eta'); m_a0^2 = 3/5 is not degenerate with any meson mass."""
    b0 = sp.Rational(13, 5)
    return {'f': sp.Rational(23, 250), 'B0': b0,
            'quark_masses': (sp.Rational(1, 50)/(2*b0), sp.Rational(1, 50)/(2*b0), sp.Rational(199, 450)/(2*b0)),
            'm0_squared': sp.Rational(3, 4), 'bare_alp_mass_squared': sp.Rational(3, 5)}


def unmixed_neutral_poles(inputs):
    """Unmixed neutral poles: {name: (vector over (N_u, N_d, N_s, a), squared mass, norm n0)}.

    The meson vectors are the rational eigenvectors of (chi_q delta + m0^2/3) in the N basis
    (eigenproblem M0 u = s K0 u with K0 = 1/2), scaled to primitive integer vectors with a positive
    first nonzero entry; n0 = u^T K0 u."""
    chi = [2*inputs['B0']*m for m in inputs['quark_masses']]
    matrix = sp.diag(*chi) + inputs['m0_squared']/NUMBER_OF_FLAVOURS*sp.ones(3, 3)
    poles = []
    for value, multiplicity, vectors in matrix.eigenvects():
        assert multiplicity == 1 and value.is_Rational, (value, multiplicity)
        vector = list(vectors[0])
        denominators = sp.ilcm(*[sp.fraction(sp.nsimplify(c))[1] for c in vector])
        vector = [sp.nsimplify(c*denominators) for c in vector]
        common = sp.igcd(*[int(c) for c in vector if c != 0])
        vector = [c/common for c in vector]
        sign = 1 if next(c for c in vector if c != 0) > 0 else -1
        vector = [sign*c for c in vector]
        poles.append((value, vector))
    poles.sort(key=lambda item: item[0])
    result = {}
    for name, (value, vector) in zip(('pi0', 'eta', 'eta_prime'), poles):
        result[name] = (vector + [0], value, sum(c**2 for c in vector)/2)
    result['alp'] = ([0, 0, 0, 1], inputs['bare_alp_mass_squared'], sp.Integer(1))
    return result


def builder_parameters(inputs, cp_conserving=False, bare_alp_mass=None):
    """Parameter dict of crossing.build_lagrangian with symbolic invariants (orders <= 4)."""
    zero_set = set(CP_ODD_INVARIANTS) if cp_conserving else set()

    def keep(values):
        return tuple(0 if v in zero_set else v for v in values) if isinstance(values, tuple) \
            else (0 if values in zero_set else values)
    return {'f': inputs['f'], 'f_a': inputs['f']/epsilon, 'B0': inputs['B0'],
            'quark_masses': inputs['quark_masses'], 'm0_squared': inputs['m0_squared'],
            'bare_alp_mass_squared': inputs['bare_alp_mass_squared'] if bare_alp_mass is None else bare_alp_mass,
            'dressed_scalar': {n: keep(scalar_invariant(n)) for n in range(2, crossing.TRUNCATION_ORDER + 1)},
            'dressed_pseudoscalar': {n: keep(pseudoscalar_invariant(n)) for n in range(1, crossing.TRUNCATION_ORDER + 1)},
            'dressed_topological': {n: keep(topological_invariant(n)) for n in range(1, crossing.TRUNCATION_ORDER + 1)}}


# ============================================================== builder vertices in the N basis
def builder_lagrangian_neutral_basis(parameters, frame_generator, active_fields):
    """crossing.build_lagrangian rewritten in the variables N_q = sqrt2 P_q (rational coefficients)."""
    lagrangian = crossing.build_lagrangian(parameters, first_frame_generator=frame_generator,
                                           active_fields=active_fields)
    replacement = {}
    for q in range(3):
        replacement[crossing.meson_fields[q]] = crossing.meson_fields[q]/sp.sqrt(2)
        replacement[crossing.derivative_fields[q]] = crossing.derivative_fields[q]/sp.sqrt(2)
    return sp.expand(lagrangian.xreplace(replacement))


class FrameVertices:
    """K, M^2 (10x10), cubic V3, quartic W and derivative quartic Z (dicts over field indices) of
    one frame, with entries in RING, as returned by crossing.extract_vertices."""

    def __init__(self, lagrangian):
        vertices = crossing.extract_vertices(crossing.sympy_terms(lagrangian), numeric=False)
        self.flags = vertices['flags']
        self.kinetic = [[to_ring(vertices['K'][i, j]) for j in range(FIELD_COUNT)] for i in range(FIELD_COUNT)]
        self.mass = [[to_ring(vertices['M2'][i, j]) for j in range(FIELD_COUNT)] for i in range(FIELD_COUNT)]
        self.cubic = {k: to_ring(v) for k, v in vertices['V3'].items() if sp.expand(v) != 0}
        self.quartic = {k: to_ring(v) for k, v in vertices['W'].items() if sp.expand(v) != 0}
        self.derivative_quartic = {k: to_ring(v) for k, v in vertices['Z'].items() if sp.expand(v) != 0}


@lru_cache(maxsize=None)
def frame_vertices(frame, pair, kappa='example', cp_conserving=False):
    """Vertices from the builder: frame 'dressed' or 'kappa'; pair None or a key of CHARGED_PAIRS
    (the charged pair kept active besides the neutral fields); kappa 'example' (rational) or
    'symbolic'."""
    inputs = meson_sector_inputs()
    if frame == 'dressed':
        frame_generator = (0, 0, 0)
    else:
        frame_generator = crossing.EXAMPLE_FRAME_GENERATOR if kappa == 'example' else kappa_symbols
    active = ['P_u', 'P_d', 'P_s', 'a']
    if pair is not None:
        first, second = CHARGED_PAIRS[pair]
        active += [crossing.FIELD_NAMES[first], crossing.FIELD_NAMES[second]]
    lagrangian = builder_lagrangian_neutral_basis(builder_parameters(inputs, cp_conserving),
                                                  frame_generator, tuple(active))
    return FrameVertices(lagrangian)


# ============================================================== closed-form neutral Lagrangian
def field_degree(monomial):
    return sum(monomial[i] for i in FIELD_GENERATOR_INDICES)


def truncate_fields(polynomial, maximum_degree):
    return RING({m: c for m, c in polynomial.items() if field_degree(m) <= maximum_degree})


def homogeneous_part(polynomial, degree):
    return RING({m: c for m, c in polynomial.items() if field_degree(m) == degree})


def taylor_series(coefficients, x_powers):
    """sum_n c_n x^n/n! for a list of ring coefficients c_0, c_1, ..."""
    return sum((c*rational(1/sp.factorial(n))*x_powers[n] for n, c in enumerate(coefficients) if c), ZERO)


def trigonometric_series(angle, maximum_degree):
    """cos and sin of a ring element (degree one in the fields), truncated."""
    cosine, sine, power = ZERO, ZERO, ONE
    for k in range(maximum_degree + 1):
        if k:
            power = truncate_fields(power*angle, maximum_degree)
        if k % 2 == 0:
            cosine += rational(sp.Integer(-1)**(k//2)/sp.factorial(k))*power
        else:
            sine += rational(sp.Integer(-1)**((k - 1)//2)/sp.factorial(k))*power
    return cosine, sine


def default_source_coefficients(inputs):
    """Taylor coefficients n = 0..MAXIMUM_SOURCE_ORDER of the dressed sources as ring elements:
    S_0 = m_q, S_1 = 0, S_n = S_hat_n; P_0 = 0, P_n = P_hat_n; Theta_0 = 0, Theta_n = Theta_hat_n."""
    scalar = [[rational(m) for m in inputs['quark_masses']], [ZERO]*3]
    scalar += [[generator(v) for v in scalar_invariant(n)] for n in range(2, MAXIMUM_SOURCE_ORDER + 1)]
    pseudoscalar = [[ZERO]*3] + [[generator(v) for v in pseudoscalar_invariant(n)]
                                 for n in range(1, MAXIMUM_SOURCE_ORDER + 1)]
    topological = [ZERO] + [generator(topological_invariant(n)) for n in range(1, MAXIMUM_SOURCE_ORDER + 1)]
    return scalar, pseudoscalar, topological


def neutral_powers_of_x(inputs, maximum_degree):
    alp = generator(neutral_field_symbols[ALP_INDEX])
    x = generator(epsilon)*alp*rational(1/inputs['f'])
    powers = [ONE]
    for _ in range(maximum_degree + 2):
        powers.append(truncate_fields(powers[-1]*x, maximum_degree))
    return powers


def neutral_mass_term(inputs, scalar, pseudoscalar, maximum_degree):
    """f^2 B0 sum_q [S_q(x) cos(N_q/f) + P_q(x) sin(N_q/f)] for given Taylor coefficients."""
    f = inputs['f']
    x_powers = neutral_powers_of_x(inputs, maximum_degree)
    total = ZERO
    for q in range(3):
        cosine, sine = trigonometric_series(generator(neutral_field_symbols[q])*rational(1/f), maximum_degree)
        scalar_function = taylor_series([c[q] for c in scalar], x_powers)
        pseudoscalar_function = taylor_series([c[q] for c in pseudoscalar], x_powers)
        total += truncate_fields(scalar_function*cosine, maximum_degree)
        total += truncate_fields(pseudoscalar_function*sine, maximum_degree)
    return total*rational(f**2*inputs['B0'])


def singlet_combination(inputs, topological, maximum_degree):
    """sum_q N_q/f + 2 Theta(x) = -X_TC."""
    x_powers = neutral_powers_of_x(inputs, maximum_degree)
    fields = sum((generator(neutral_field_symbols[q]) for q in range(3)), ZERO)*rational(1/inputs['f'])
    return fields + 2*taylor_series(topological, x_powers)


def neutral_dressed_potential(inputs, maximum_degree=6, cp_conserving=False):
    """Non-derivative dressed Lagrangian of the commuting fields (N_u, N_d, N_s, a), exact in the
    Taylor coefficients, truncated at field degree `maximum_degree` (constant dropped):
        f^2 B0 sum_q [S_q(x) cos(N_q/f) + P_q(x) sin(N_q/f)] - f^2 m0^2/12 (sum N/f + 2 Theta(x))^2
        - m_a0^2 a^2/2,   x = epsilon a/f,
    S_q(x) = m_q + sum_{n>=2} x^n/n! S_hat_n,q,  P_q(x) = sum_{n>=1} x^n/n! P_hat_n,q.
    The kinetic term is (1/4) sum_q (dN_q)^2 + (1/2)(da)^2 and has no interaction."""
    scalar, pseudoscalar, topological = default_source_coefficients(inputs)
    mass_term = neutral_mass_term(inputs, scalar, pseudoscalar, maximum_degree)
    combination = singlet_combination(inputs, topological, maximum_degree)
    singlet = -rational(inputs['f']**2*inputs['m0_squared']/(4*NUMBER_OF_FLAVOURS))*truncate_fields(
        combination**2, maximum_degree)
    alp = generator(neutral_field_symbols[ALP_INDEX])
    bare = -rational(inputs['bare_alp_mass_squared']/2)*alp**2
    lagrangian = truncate_fields(mass_term + singlet + bare, maximum_degree)
    lagrangian = RING({m: c for m, c in lagrangian.items() if field_degree(m) > 0})
    if cp_conserving:
        lagrangian = set_to_zero(lagrangian, CP_ODD_INVARIANTS)
    return lagrangian


def neutral_kinetic_matrix():
    return [[rational(sp.Rational(1, 2)) if i == j and i < 3 else (ONE if i == j else ZERO)
             for j in range(NEUTRAL_SIZE)] for i in range(NEUTRAL_SIZE)]


def quadratic_mass_matrix(potential):
    """M^2_ij = -d^2 L/dF_i dF_j at zero fields (neutral closed form)."""
    quadratic = homogeneous_part(potential, 2)
    return [[-multilinear_form(quadratic, [unit_neutral(i), unit_neutral(j)], 10**6)
             for j in range(NEUTRAL_SIZE)] for i in range(NEUTRAL_SIZE)]


def unit_neutral(index):
    return [ONE if k == index else ZERO for k in range(NEUTRAL_SIZE)]


def multilinear_form(homogeneous, vectors, order, complete=True):
    """Vertex contraction V(w_1, ..., w_k): the directional derivatives sum_i w_i d/dF_i applied
    to a homogeneous polynomial in (N_u, N_d, N_s, a); vectors have 4 ring components.  With
    complete=False fewer vectors than the degree may be given and the result keeps the fields."""
    polynomial = truncate(homogeneous, order)
    for vector in vectors:
        result = ZERO
        for index, component in enumerate(vector):
            if component:
                derivative = polynomial.diff(generator(neutral_field_symbols[index]))
                if derivative:
                    result += truncate(component*derivative, order)
        polynomial = result
    if complete:
        assert all(field_degree(m) == 0 for m in polynomial.keys())
    return polynomial


# ============================================================== perturbative poles
def perturbative_pole(kinetic, mass, zeroth_vector, zeroth_mass, order):
    """Pole s(epsilon) and vector u(epsilon) of (M^2(epsilon) - s K(epsilon)) u = 0, as truncated
    series, with u^T K u = n0 = u_0^T K_0 u_0 (Rayleigh-Schroedinger, rational reduced resolvent).

    kinetic, mass: square lists of ring elements polynomial in epsilon whose epsilon^0 parts are
    rational; zeroth_vector: rational null vector of M_0 - zeroth_mass K_0."""
    size = len(zeroth_vector)
    kinetic_orders = [[[epsilon_coefficient(kinetic[i][j], k) for k in range(order + 1)]
                       for j in range(size)] for i in range(size)]
    mass_orders = [[[epsilon_coefficient(mass[i][j], k) for k in range(order + 1)]
                    for j in range(size)] for i in range(size)]
    kinetic_zero = sp.Matrix(size, size, lambda i, j: constant_value(kinetic_orders[i][j][0]))
    mass_zero = sp.Matrix(size, size, lambda i, j: constant_value(mass_orders[i][j][0]))
    vector_zero = sp.Matrix(zeroth_vector)
    pencil = mass_zero - zeroth_mass*kinetic_zero
    assert pencil*vector_zero == sp.zeros(size, 1)
    bordered = sp.zeros(size + 1, size + 1)
    bordered[:size, :size] = pencil
    bordered[:size, size] = kinetic_zero*vector_zero
    bordered[size, :size] = (kinetic_zero*vector_zero).T
    resolvent = bordered.inv()[:size, :size]
    resolvent_ring = [[rational(resolvent[i, j]) for j in range(size)] for i in range(size)]
    norm_zero = (vector_zero.T*kinetic_zero*vector_zero)[0, 0]
    vectors = [[rational(c) for c in zeroth_vector]]
    poles = [rational(zeroth_mass)]
    kinetic_zero_vector = [rational(c) for c in kinetic_zero*vector_zero]

    def apply(matrix_orders, k, vector):
        return [sum((matrix_orders[i][j][k]*vector[j] for j in range(size) if vector[j]), ZERO)
                for i in range(size)]
    for k in range(1, order + 1):
        rest = [ZERO]*size
        for j in range(1, k + 1):
            rest = [a + b for a, b in zip(rest, apply(mass_orders, j, vectors[k - j]))]
        for i in range(0, k + 1):
            for j in range(0, k - i + 1):
                l = k - i - j
                if (i, j, l) in ((0, 0, k), (k, 0, 0)) or not poles[i]:
                    continue
                product = apply(kinetic_orders, j, vectors[l])
                rest = [a - poles[i]*b for a, b in zip(rest, product)]
        pole = sum((rational(zeroth_vector[i])*rest[i] for i in range(size) if zeroth_vector[i]), ZERO)
        pole = pole*rational(1/norm_zero)
        right = [pole*kinetic_zero_vector[i] - rest[i] for i in range(size)]
        vectors.append([sum((resolvent_ring[i][j]*right[j] for j in range(size)), ZERO) for i in range(size)])
        poles.append(pole)
    eps = generator(epsilon)
    pole_series = sum((poles[k]*eps**k for k in range(order + 1)), ZERO)
    vector_series = [sum((vectors[k][i]*eps**k for k in range(order + 1)), ZERO) for i in range(size)]
    norm = ZERO
    for i in range(size):
        for j in range(size):
            if kinetic[i][j] and vector_series[i] and vector_series[j]:
                norm += truncate(truncate(kinetic[i][j]*vector_series[i], order)*vector_series[j], order)
    factor = series_inverse_square_root(truncate(norm*rational(1/norm_zero), order), order)
    vector_series = [truncate(factor*c, order) for c in vector_series]
    return pole_series, vector_series, norm_zero


def pad(vector, size=FIELD_COUNT):
    return list(vector) + [ZERO]*(size - len(vector))


@lru_cache(maxsize=None)
def frame_poles(frame, order=TWO_POINT_ORDER, kappa='example', cp_conserving=False):
    """Neutral poles of a builder frame: {name: (s series, 10-component u, n0)}."""
    vertices = frame_vertices(frame, 'ud', kappa, cp_conserving)
    kinetic = [row[:NEUTRAL_SIZE] for row in vertices.kinetic[:NEUTRAL_SIZE]]
    mass = [row[:NEUTRAL_SIZE] for row in vertices.mass[:NEUTRAL_SIZE]]
    result = {}
    for name, (vector, value, _) in unmixed_neutral_poles(meson_sector_inputs()).items():
        pole, series, norm = perturbative_pole(kinetic, mass, vector, value, order)
        result[name] = (pole, pad(series), norm)
    return result


def charged_poles(vertices, pair):
    """Charged components of a pair: exact unit vectors and constant masses (asserted)."""
    result = {}
    for index in CHARGED_PAIRS[pair]:
        for j in range(FIELD_COUNT):
            if j != index:
                assert not vertices.kinetic[index][j] and not vertices.mass[index][j]
        assert vertices.kinetic[index][index] == ONE
        value = constant_value(vertices.mass[index][index])
        result[index] = (rational(value), [ONE if k == index else ZERO for k in range(FIELD_COUNT)], sp.Integer(1))
    return result


# ============================================================== tensor contractions (builder dicts)
def contract(tensor, vectors, order):
    """Contract the last len(vectors) indices of a dict tensor with the vectors (in order).
    Returns a dict over the remaining indices (key () for a full contraction)."""
    current = tensor
    for vector in reversed(vectors):
        reduced = {}
        for key, value in current.items():
            component = vector[key[-1]]
            if component:
                product = truncate(value*component, order)
                if product:
                    reduced[key[:-1]] = reduced.get(key[:-1], ZERO) + product
        current = {k: v for k, v in reduced.items() if v}
    return current


def full_contraction(tensor, vectors, order):
    return contract(tensor, vectors, order).get((), ZERO)


def matrix_series_inverse(matrix, order):
    """Inverse of a square matrix of ring elements with invertible rational epsilon^0 part
    (Newton iteration X -> X (2 - M X))."""
    size = len(matrix)
    zero_part = sp.Matrix(size, size, lambda i, j: constant_value(epsilon_coefficient(matrix[i][j], 0)))
    inverse_zero = zero_part.inv()
    current = [[rational(inverse_zero[i, j]) for j in range(size)] for i in range(size)]

    def multiply(first, second):
        return [[truncate(sum((first[i][k]*second[k][j] for k in range(size)), ZERO), order)
                 for j in range(size)] for i in range(size)]
    precision = 0
    while precision < order:
        product = multiply(matrix, current)
        correction = [[(2*ONE if i == j else ZERO) - product[i][j] for j in range(size)] for i in range(size)]
        current = multiply(current, correction)
        precision = 2*precision + 1
    return current


@lru_cache(maxsize=None)
def neutral_mass_inverse(frame, order=TWO_POINT_ORDER, kappa='example'):
    vertices = frame_vertices(frame, 'ud', kappa)
    return matrix_series_inverse([row[:NEUTRAL_SIZE] for row in vertices.mass[:NEUTRAL_SIZE]], order)


def contraction_pieces(frame, label, order=TWO_POINT_ORDER):
    """Internal-line contraction <F F> = sum g w w^T A0, <dF.dF> = sum h w w^T A0 for a loop label:
    a neutral pole name or a charged pair."""
    if label in CHARGED_PAIRS:
        vertices = frame_vertices(frame, label)
        pieces = []
        for index, (value, vector, _) in charged_poles(vertices, label).items():
            pieces.append((vector, -ONE, -value))
        return pieces
    pole, vector, norm = frame_poles(frame, order)[label]
    g_value = -rational(1/norm)
    return [(vector, g_value, truncate(g_value*pole, order))]


def two_point_element(frame, label, left, right, s_value, order=TWO_POINT_ORDER, run_pair=None):
    """left^T [delta M^2 - s delta K] right per unit A0 of the loop label (contact tadpole) and the
    loop-induced one-point term -V3(left, right, <F>).  Returns (tadpole, one_point).

    run_pair: charged pair active in the builder run (needed when left/right are charged); by
    default the label's pair for a charged label and 'ud' otherwise (the neutral vertices are the
    same in every run)."""
    if run_pair is None:
        run_pair = label if label in CHARGED_PAIRS else 'ud'
    assert label not in CHARGED_PAIRS or label == run_pair
    vertices = frame_vertices(frame, run_pair)
    tadpole = ZERO
    source = [ZERO]*NEUTRAL_SIZE
    for vector, g_value, h_value in contraction_pieces(frame, label, order):
        tadpole += truncate(-rational(sp.Rational(1, 2))*g_value*full_contraction(
            vertices.quartic, [left, right, vector, vector], order), order)
        tadpole += truncate(-rational(sp.Rational(1, 2))*h_value*full_contraction(
            vertices.derivative_quartic, [vector, vector, left, right], order), order)
        tadpole += truncate(-rational(sp.Rational(1, 2))*truncate(s_value*g_value, order)*full_contraction(
            vertices.derivative_quartic, [left, right, vector, vector], order), order)
        reduced = contract(vertices.cubic, [vector, vector], order)
        for key, value in reduced.items():
            assert key[0] < NEUTRAL_SIZE
            source[key[0]] += truncate(rational(sp.Rational(1, 2))*g_value*value, order)
    inverse = neutral_mass_inverse(frame, order)
    expectation = pad([truncate(sum((inverse[i][j]*source[j] for j in range(NEUTRAL_SIZE)), ZERO), order)
                       for i in range(NEUTRAL_SIZE)])
    one_point = -full_contraction(vertices.cubic, [left, right, expectation], order) if vertices.cubic else ZERO
    return truncate(tadpole, order), truncate(one_point, order)


def loop_labels():
    return ('pi0', 'eta', 'eta_prime', 'alp') + tuple(CHARGED_PAIRS)


def external_state(frame, name, order=TWO_POINT_ORDER):
    """(s series, vector, n0) of an external pole: neutral pole name or a charged field index."""
    if isinstance(name, int):
        pair = next(p for p, indices in CHARGED_PAIRS.items() if name in indices)
        value, vector, norm = charged_poles(frame_vertices(frame, pair), pair)[name]
        return value, vector, norm
    return frame_poles(frame, order)[name]


@lru_cache(maxsize=None)
def pole_mass_coefficients(frame, external, order=TWO_POINT_ORDER):
    """One-loop pole-shift coefficients of an external pole (neutral name or charged field index):
    'tadpole'[label], 'one_point'[label]: series of the coefficient of A0(m_label^2) (contact
        tadpole and loop-induced one-point term; label = neutral pole name or charged pair);
    'bubble_vertex'[(r, t)]: series of the cubic vertex V3(u_n, u_r, u_t) (unnormalized vectors);
        the coefficient of Bbar(m_n^2; m_r, m_t) is -mult V3^2/(n0_n n0_r n0_t) in units of
        1/(32 pi^2), mult = 2 for r != t, so its leading power is twice that of the vertex."""
    pole, vector, norm = external_state(frame, external, order)
    if isinstance(external, int):
        pair = next(p for p, indices in CHARGED_PAIRS.items() if external in indices)
        labels = ('pi0', 'eta', 'eta_prime', 'alp', pair)
    else:
        pair = None
        labels = loop_labels()
    tadpole, one_point = {}, {}
    for label in labels:
        first, second = two_point_element(frame, label, vector, vector, pole, order, pair)
        tadpole[label] = truncate(first*rational(1/norm), order)
        one_point[label] = truncate(second*rational(1/norm), order)
    vertices = frame_vertices(frame, pair or 'ud')
    poles = frame_poles(frame, order)
    names = ('pi0', 'eta', 'eta_prime', 'alp')
    bubble_vertex = {}
    for i, r in enumerate(names):
        for t in names[i:]:
            bubble_vertex[(r, t)] = full_contraction(vertices.cubic, [vector, poles[r][1], poles[t][1]], order)
    return {'tadpole': tadpole, 'one_point': one_point, 'bubble_vertex': bubble_vertex}


def doubled(power):
    return None if power is None else 2*power


def product_power(first, second):
    """Leading power of a product of two series (the ring is an integral domain)."""
    return None if first is None or second is None else first + second


# ============================================================== exact first-order solution Y
def first_order_solution_exact(inputs, s_value):
    """Diagonal Y with (chi_q - s) Y_q + m0^2/N_f tr Y = C_q, C = B0 P_hat_1 - m0^2/N_f Theta_hat_1
    (exact version of crossing.first_order_solution; ring elements)."""
    chi = [2*inputs['B0']*m for m in inputs['quark_masses']]
    m0 = inputs['m0_squared']/NUMBER_OF_FLAVOURS
    matrix = sp.diag(*[c - s_value for c in chi]) + m0*sp.ones(3, 3)
    inverse = matrix.inv()
    source = [rational(inputs['B0'])*generator(pseudoscalar_invariant(1)[q])
              - rational(m0)*generator(topological_invariant(1)) for q in range(3)]
    return [sum((rational(inverse[q, r])*source[r] for r in range(3)), ZERO) for q in range(3)]


def quark_mass_shift_closed_form(inputs):
    """Delta_q = -(1/(2 f^2)) (S_hat_2 + 4 Y P_hat_1 - 4 Y^2 m_q)_qq, the epsilon^2 coefficient of the
    diagonal matrix whose B0 <lambda_n^2 .> gives c_{n,a} for neutral n (class A), Y at s = m_a0^2."""
    y_values = first_order_solution_exact(inputs, inputs['bare_alp_mass_squared'])
    f = inputs['f']
    return [-rational(1/(2*f**2))*(generator(scalar_invariant(2)[q])
                                   + 4*y_values[q]*generator(pseudoscalar_invariant(1)[q])
                                   - 4*rational(inputs['quark_masses'][q])*y_values[q]**2) for q in range(3)]


# ============================================================== C1, C3: pole masses
def tree_pole_check(frame='dressed', order=TWO_POINT_ORDER):
    """Series m_n^2 - m_{n,0}^2 for the neutral poles."""
    result = {}
    for name, (pole, _, _) in frame_poles(frame, order).items():
        result[name] = truncate(pole - rational(constant_value(epsilon_coefficient(pole, 0))), order)
    return result


def alp_line_coefficient(external):
    """Coefficient of A0(m_A^2) (tadpole + one-point), dressed frame."""
    coefficients = pole_mass_coefficients('dressed', external)
    return coefficients['tadpole']['alp'] + coefficients['one_point']['alp']


@lru_cache(maxsize=None)
def neutral_quark_mass_shift_check():
    """c_{n,a} at epsilon^2 for pi0, eta, eta': solve c = B0 <lambda_n^2 Delta> for a diagonal Delta
    (only Delta_u + Delta_d and Delta_s enter when m_u = m_d: two unknowns, three equations) and
    compare with the closed form of quark_mass_shift_closed_form."""
    inputs = meson_sector_inputs()
    b0 = inputs['B0']
    unmixed = unmixed_neutral_poles(inputs)
    coefficients = {name: epsilon_coefficient(alp_line_coefficient(name), 2) for name in ('pi0', 'eta', 'eta_prime')}
    rows = {}
    for name in coefficients:
        vector, _, norm = unmixed[name]
        assert vector[0]**2 == vector[1]**2
        rows[name] = (b0*vector[0]**2/norm, b0*vector[2]**2/norm)  # <lambda^2 Delta> = (u_u^2 (D_u+D_d) + u_s^2 D_s)/n0
    first, second = rows['pi0'], rows['eta_prime']
    determinant = first[0]*second[1] - first[1]*second[0]
    sum_ud = (coefficients['pi0']*rational(second[1]) - coefficients['eta_prime']*rational(first[1]))*rational(1/determinant)
    strange = (coefficients['eta_prime']*rational(first[0]) - coefficients['pi0']*rational(second[0]))*rational(1/determinant)
    eta_residual = coefficients['eta'] - rational(rows['eta'][0])*sum_ud - rational(rows['eta'][1])*strange
    closed = quark_mass_shift_closed_form(inputs)
    return {'coefficients': coefficients, 'sum_ud': sum_ud, 'strange': strange,
            'eta_residual': eta_residual, 'sum_ud_minus_closed_form': sum_ud - closed[0] - closed[1],
            'strange_minus_closed_form': strange - closed[2]}


@lru_cache(maxsize=None)
def charged_remainder_check(pair):
    """c_{qr,a} at epsilon^2 minus B0 (Delta_q + Delta_r), compared with
    (1/(2 f^2)) m_a0^2 <[Y, lambda]^2> = -(m_a0^2/f^2) (Y_q - Y_r)^2."""
    inputs = meson_sector_inputs()
    q, r = crossing.FLAVOUR_INDEX[pair[0]], crossing.FLAVOUR_INDEX[pair[1]]
    coefficient = epsilon_coefficient(alp_line_coefficient(CHARGED_PAIRS[pair][0]), 2)
    closed = quark_mass_shift_closed_form(inputs)
    y_values = first_order_solution_exact(inputs, inputs['bare_alp_mass_squared'])
    remainder = coefficient - rational(inputs['B0'])*(closed[q] + closed[r])
    expected = -rational(inputs['bare_alp_mass_squared']/inputs['f']**2)*(y_values[q] - y_values[r])**2
    return {'coefficient': coefficient, 'remainder': remainder, 'expected': expected,
            'difference': remainder - expected}


@lru_cache(maxsize=None)
def frame_comparison(external, order=TWO_POINT_ORDER):
    """Differences (kappa frame - dressed frame) of every on-shell coefficient of a pole: pole
    mass, tadpole + one-point per label, squared bubble vertex per pair (exact series)."""
    dressed = pole_mass_coefficients('dressed', external, order)
    framed = pole_mass_coefficients('kappa', external, order)
    differences = {'pole': truncate(external_state('kappa', external, order)[0]
                                    - external_state('dressed', external, order)[0], order)}
    for label in dressed['tadpole']:
        differences[label] = truncate(framed['tadpole'][label] + framed['one_point'][label]
                                      - dressed['tadpole'][label] - dressed['one_point'][label], order)
    for key in dressed['bubble_vertex']:
        differences[('bubble',) + key] = truncate(framed['bubble_vertex'][key]**2
                                                  - dressed['bubble_vertex'][key]**2, 2*order)
    return differences


def determinant_polynomial(matrix):
    """Determinant of a square list of ring elements (Leibniz formula)."""
    size = len(matrix)
    total = ZERO
    for permutation in itertools.permutations(range(size)):
        sign = Permutation(list(permutation)).signature()
        term = ONE
        for row, column in enumerate(permutation):
            term = term*matrix[row][column]
            if not term:
                break
        total += sign*term
    return total


@lru_cache(maxsize=None)
def alp_mass_absorption_check(order=TWO_POINT_ORDER):
    """d m_A^2/d m_a0^2 from the secular polynomial P(s, m_a0^2) = det(M^2 - s K), as
    -dP/dm_a0^2 / dP/ds at s = m_A^2(epsilon), compared with (u_A,a)^2/n0 (dressed frame)."""
    vertices = frame_vertices('dressed', 'ud')
    inputs = meson_sector_inputs()
    mu, s_gen = generator(bare_alp_mass_symbol), generator(squared_momentum)
    matrix = [[vertices.mass[i][j] - s_gen*vertices.kinetic[i][j] for j in range(NEUTRAL_SIZE)]
              for i in range(NEUTRAL_SIZE)]
    matrix[ALP_INDEX][ALP_INDEX] += mu - rational(inputs['bare_alp_mass_squared'])
    polynomial = determinant_polynomial(matrix)
    pole, vector, norm = frame_poles('dressed', order)['alp']
    at_pole = {}
    for key, value in (('mass', polynomial.diff(mu)), ('s', polynomial.diff(s_gen))):
        value = value.compose(mu, rational(inputs['bare_alp_mass_squared']))
        at_pole[key] = substitute(value, squared_momentum, pole, order)
    derivative = truncate(-at_pole['mass']*series_inverse(at_pole['s'], order), order)
    component_squared = truncate(vector[ALP_INDEX]**2*rational(1/norm), order)
    return {'derivative': derivative, 'alp_component_squared': component_squared,
            'difference': truncate(derivative - component_squared, order)}


# ============================================================== C2: mixing entry
def unit_vector(index):
    return [ONE if k == index else ZERO for k in range(FIELD_COUNT)]


def pion_direction():
    """lambda_3 in the N basis: (N_u, N_d) = (1, -1), n0 = 1."""
    return [ONE, -ONE] + [ZERO]*(FIELD_COUNT - 2)


@lru_cache(maxsize=None)
def mixing_entry_tree():
    """D_{pi0 a}(s) = e_pi^T (s K - M^2) e_a in both frames (symbolic kappa), exact in epsilon, the
    frame difference, its value at the unmixed pion mass and at the pion pole (series)."""
    s_gen = generator(squared_momentum)
    pion = pion_direction()
    entries = {}
    for frame in ('dressed', 'kappa'):
        vertices = frame_vertices(frame, 'ud', 'symbolic')
        entries[frame] = sum((pion[i]*(s_gen*vertices.kinetic[i][ALP_INDEX] - vertices.mass[i][ALP_INDEX])
                              for i in range(NEUTRAL_SIZE)), ZERO)
    difference = entries['kappa'] - entries['dressed']
    unmixed_mass = unmixed_neutral_poles(meson_sector_inputs())['pi0'][1]
    pole = frame_poles('dressed')['pi0'][0]
    expected = -generator(epsilon)*(generator(kappa_symbols[0]) - generator(kappa_symbols[1]))*(
        s_gen - rational(unmixed_mass))
    return {'dressed': entries['dressed'], 'kappa': entries['kappa'], 'difference': difference,
            'difference_minus_closed_form': difference - expected,
            'at_unmixed_mass': difference.compose(s_gen, rational(unmixed_mass)),
            'at_pole': substitute(difference, squared_momentum, pole, TWO_POINT_ORDER),
            'pole_shift': truncate(pole - rational(unmixed_mass), TWO_POINT_ORDER)}


@lru_cache(maxsize=None)
def mixing_entry_one_loop(order=TWO_POINT_ORDER):
    """(pi0, a) entry of the one-loop self-energy, coefficient of A0 per label (tadpole + one-point),
    as a function of s in both frames; its frame difference at s = m_pi^2(epsilon); and the
    leading powers of the bubble vertices V3(e_pi, r, t), V3(e_a, r, t) (dressed frame)."""
    pion, alp = pion_direction(), unit_vector(ALP_INDEX)
    s_gen = generator(squared_momentum)
    entries = {'dressed': {}, 'kappa': {}}
    for frame in entries:
        for label in loop_labels():
            tadpole, one_point = two_point_element(frame, label, pion, alp, s_gen, order)
            entries[frame][label] = truncate(tadpole + one_point, order)
    pole = frame_poles('dressed', order)['pi0'][0]
    on_shell_difference = {label: substitute(entries['kappa'][label] - entries['dressed'][label],
                                             squared_momentum, pole, order) for label in loop_labels()}
    vertices = frame_vertices('dressed', 'ud')
    poles = frame_poles('dressed', order)
    names = ('pi0', 'eta', 'eta_prime', 'alp')
    bubble_powers = {}
    for i, r in enumerate(names):
        for t in names[i:]:
            left = full_contraction(vertices.cubic, [pion, poles[r][1], poles[t][1]], order)
            right = full_contraction(vertices.cubic, [alp, poles[r][1], poles[t][1]], order)
            bubble_powers[(r, t)] = (leading_power(left), leading_power(right))
    return {'entries': entries, 'on_shell_difference': on_shell_difference, 'bubble_powers': bubble_powers}


# ============================================================== C4: A -> 3 pi0
@lru_cache(maxsize=None)
def amplitude_poles(order=AMPLITUDE_ORDER, cp_conserving=True):
    """Closed-form dressed potential (sixth order) and its neutral poles (4 components)."""
    inputs = meson_sector_inputs()
    potential = neutral_dressed_potential(inputs, 6, cp_conserving)
    mass = quadratic_mass_matrix(potential)
    kinetic = neutral_kinetic_matrix()
    poles = {name: perturbative_pole(kinetic, mass, vector, value, order)
             for name, (vector, value, _) in unmixed_neutral_poles(inputs).items()}
    return potential, poles


@lru_cache(maxsize=None)
def tree_amplitude(order=AMPLITUDE_ORDER):
    """M(A -> 3 pi0) = V4(u_A, u_pi, u_pi, u_pi) (n0 = 1 for both), CP-conserving sources; the
    cubic vertices vanish, so this is the full tree amplitude."""
    potential, poles = amplitude_poles(order)
    assert not homogeneous_part(potential, 3)
    alp, pion = poles['alp'][1], poles['pi0'][1]
    assert poles['alp'][2] == 1 and poles['pi0'][2] == 1
    return multilinear_form(homogeneous_part(potential, 4), [alp, pion, pion, pion], order)


@lru_cache(maxsize=None)
def tree_amplitude_conjecture_check(order=AMPLITUDE_ORDER):
    """epsilon^1 coefficient against (B0/f^2) <(2 Y m_q - P_hat_1) lambda_3> and against the
    equivalent form (B0/f^2) m_a0^2 (P_hat_1u - P_hat_1d)/(m_pi0^2 - m_a0^2) (m_u = m_d); meson
    components of the ALP-like vector at O(epsilon) against 2 Y (N basis)."""
    inputs = meson_sector_inputs()
    amplitude = tree_amplitude(order)
    alp_mass = inputs['bare_alp_mass_squared']
    y_values = first_order_solution_exact(inputs, alp_mass)
    lambda_three = (1, -1, 0)
    prefactor = rational(inputs['B0']/inputs['f']**2)
    conjecture = prefactor*sum((lambda_three[q]*(2*y_values[q]*rational(inputs['quark_masses'][q])
                                                 - generator(pseudoscalar_invariant(1)[q])) for q in range(3)), ZERO)
    pion_mass = unmixed_neutral_poles(inputs)['pi0'][1]
    reduced = prefactor*rational(alp_mass/(pion_mass - alp_mass))*(
        generator(pseudoscalar_invariant(1)[0]) - generator(pseudoscalar_invariant(1)[1]))
    _, poles = amplitude_poles(order)
    alp_vector = poles['alp'][1]
    first_order = epsilon_coefficient(amplitude, 1)
    return {'amplitude': amplitude, 'first_order': first_order,
            'first_order_minus_conjecture': first_order - conjecture,
            'first_order_minus_reduced_form': first_order - reduced,
            'alp_vector_minus_two_y': [epsilon_coefficient(alp_vector[q], 1) - 2*y_values[q] for q in range(3)]}


@lru_cache(maxsize=None)
def kappa_map_linearity_check():
    """Builder Lagrangian of the neutral fields in the kappa frame (symbolic kappa and invariants)
    minus the dressed one after N_q -> N_q - 2 epsilon kappa_q a, dN_q -> dN_q - 2 epsilon kappa_q da;
    returns the nonzero coefficients with field degree 1..4 (empty list: the map is linear)."""
    parameters = builder_parameters(meson_sector_inputs())
    active = ('P_u', 'P_d', 'P_s', 'a')
    dressed = builder_lagrangian_neutral_basis(parameters, (0, 0, 0), active)
    framed = builder_lagrangian_neutral_basis(parameters, kappa_symbols, active)
    fields = crossing.meson_fields[:NEUTRAL_SIZE]
    derivatives = crossing.derivative_fields[:NEUTRAL_SIZE]
    replacement = {}
    for q in range(3):
        replacement[fields[q]] = fields[q] - 2*epsilon*kappa_symbols[q]*fields[ALP_INDEX]
        replacement[derivatives[q]] = derivatives[q] - 2*epsilon*kappa_symbols[q]*derivatives[ALP_INDEX]
    mapped = sp.expand(dressed.xreplace(replacement))
    polynomial = sp.Poly(sp.expand(mapped - framed), *fields, *derivatives)
    return [c for monomial, c in polynomial.terms() if 0 < sum(monomial) <= crossing.TRUNCATION_ORDER]


def closed_form_against_builder():
    """Closed-form neutral potential (degree <= 4) minus the builder's non-derivative dressed
    Lagrangian of the neutral fields (symbolic invariants and epsilon)."""
    inputs = meson_sector_inputs()
    builder = builder_lagrangian_neutral_basis(builder_parameters(inputs), (0, 0, 0), ('P_u', 'P_d', 'P_s', 'a'))
    builder = builder.subs({d: 0 for d in crossing.derivative_fields})
    builder = builder.xreplace({crossing.meson_fields[q]: neutral_field_symbols[q] for q in range(3)})
    builder_ring = RING({m: c for m, c in to_ring(builder).items() if field_degree(m) > 0})
    return builder_ring - truncate_fields(neutral_dressed_potential(inputs, 4), 4)


@lru_cache(maxsize=None)
def tree_amplitude_frame_check(order=AMPLITUDE_ORDER):
    """Tree A -> 3 pi0 amplitude from the builder in the dressed and kappa (symbolic) frames,
    CP-conserving sources (no cubic and no derivative quartic vertices in either frame)."""
    amplitudes = {}
    unmixed = unmixed_neutral_poles(meson_sector_inputs())
    for frame in ('dressed', 'kappa'):
        vertices = frame_vertices(frame, None, 'symbolic', True)
        assert not vertices.cubic and not vertices.derivative_quartic
        kinetic = [row[:NEUTRAL_SIZE] for row in vertices.kinetic[:NEUTRAL_SIZE]]
        mass = [row[:NEUTRAL_SIZE] for row in vertices.mass[:NEUTRAL_SIZE]]
        _, alp, _ = perturbative_pole(kinetic, mass, *unmixed['alp'][:2], order)
        _, pion, _ = perturbative_pole(kinetic, mass, *unmixed['pi0'][:2], order)
        amplitudes[frame] = full_contraction(vertices.quartic, [pad(alp), pad(pion), pad(pion), pad(pion)], order)
    return {'amplitudes': amplitudes, 'difference': truncate(amplitudes['kappa'] - amplitudes['dressed'], order)}


@lru_cache(maxsize=None)
def amplitude_one_loop(order=AMPLITUDE_ORDER):
    """One-loop A -> 3 pi0 coefficients with neutral internal lines (dressed frame, CP-conserving):
    'tadpole'[r]: coefficient of A0(m_r^2), -(1/(2 n0_r)) V6(r, r, A, pi, pi, pi);
    'bubble'[(r, t)]: coefficient of Bbar(s_channel; m_r, m_t) for one channel (A pi_1 -> r t ->
        pi_2 pi_3), units 1/(16 pi^2): sym V4(r, t, A, pi) V4(r, t, pi, pi)/(n0_r n0_t), sym = 1/2
        for r = t and 1 otherwise (the other two channels are the same with s -> t, u);
    'external'[r]: external-leg correction, coefficient of A0(m_r^2):
        sum_legs sum_{m != leg} S^(r)_{m leg}/(s_leg - s_m) M_tree(leg -> m),
        S^(r)_{mn} = V4(r, r, m, n)/(2 n0_r) (normalized vectors).  There is no wave-function term:
        the neutral-line self-energies do not depend on s and there are no cubic vertices."""
    potential, poles = amplitude_poles(order)
    quartic = homogeneous_part(potential, 4)
    sextic = homogeneous_part(potential, 6)
    alp, pion = poles['alp'][1], poles['pi0'][1]
    names = ('pi0', 'eta', 'eta_prime', 'alp')
    tadpole = {}
    for r in names:
        vector, norm = poles[r][1], poles[r][2]
        tadpole[r] = truncate(-rational(1/(2*norm))*multilinear_form(
            sextic, [alp, pion, pion, pion, vector, vector], order), order)
    bubble = {}
    for i, r in enumerate(names):
        for t in names[i:]:
            left = multilinear_form(quartic, [alp, pion, poles[r][1], poles[t][1]], order)
            right = multilinear_form(quartic, [pion, pion, poles[r][1], poles[t][1]], order)
            symmetry = sp.Rational(1, 2) if r == t else sp.Integer(1)
            bubble[(r, t)] = truncate(rational(symmetry/(poles[r][2]*poles[t][2]))*left*right, order)
    # external legs: the three pion legs are identical, so they are counted once with weight 3
    trees, denominators = {}, {}
    for leg, weight in (('alp', 1), ('pi0', 3)):
        for m in names:
            if m == leg:
                continue
            replaced = [poles[m][1], pion, pion, pion] if leg == 'alp' else [alp, poles[m][1], pion, pion]
            trees[(leg, m)] = truncate(rational(weight)*multilinear_form(quartic, replaced, order), order)
            denominators[(leg, m)] = series_inverse(truncate(poles[leg][0] - poles[m][0], order), order)
    external = {}
    for r in names:
        vector_r, norm_r = poles[r][1], poles[r][2]
        reduced = multilinear_form(quartic, [vector_r, vector_r], order, complete=False)
        total = ZERO
        for (leg, m), tree in trees.items():
            self_energy = multilinear_form(reduced, [poles[m][1], poles[leg][1]], order)
            factor = rational(1/(2*norm_r*poles[m][2]*poles[leg][2]))
            total += truncate(truncate(self_energy*tree, order)*denominators[(leg, m)], order)*factor
        external[r] = truncate(total, order)
    return {'tadpole': tadpole, 'bubble': bubble, 'external': external}


@lru_cache(maxsize=None)
def source_shift_identity():
    """Contraction of one pure-ALP pair with <a a> = -A0, per unit A0: -(1/2) d^2 L/da^2 =
    -(1/(2 f_a^2)) d^2 L/dx^2, compared with the linear change of L under
        S_n -> S_n - S_{n+2}/(2 f_a^2),  P_n -> P_n - P_{n+2}/(2 f_a^2),
        Theta_n -> Theta_n - Theta_{n+2}/(2 f_a^2)   (n >= 0; S_0 = m_q, S_1 = P_0 = Theta_0 = 0),
    i.e. M_hat -> M_hat - M_hat''/(2 f_a^2), Theta_hat -> Theta_hat - Theta_hat''/(2 f_a^2), plus the
    pure-ALP term (f^2 m0^2/(3 f_a^2)) Theta_hat'(x)^2 from the square of X_TC.  General Taylor
    coefficients (orders <= 6); all terms of field degree 1..4 are compared, which contains the
    A -> 3 pi0 vertex n_a = 1, n_P = 3 obtained from n_a = 3, n_P = 3."""
    inputs = meson_sector_inputs()
    f = inputs['f']
    degree = 6
    potential = neutral_dressed_potential(inputs, degree)
    alp = generator(neutral_field_symbols[ALP_INDEX])
    contracted = truncate_fields(-rational(sp.Rational(1, 2))*potential.diff(alp).diff(alp), degree - 2)
    contracted = RING({m: c for m, c in contracted.items() if field_degree(m) > 0})
    half_inverse_square = generator(epsilon)**2*rational(1/(2*f**2))     # 1/(2 f_a^2)
    scalar, pseudoscalar, topological = default_source_coefficients(inputs)

    def shifted(coefficients):
        return [[-half_inverse_square*coefficients[n + 2][q] if n + 2 < len(coefficients) else ZERO
                 for q in range(3)] for n in range(len(coefficients))]
    mass_change = neutral_mass_term(inputs, shifted(scalar), shifted(pseudoscalar), degree)
    theta_shift = [-half_inverse_square*topological[n + 2] if n + 2 < len(topological) else ZERO
                   for n in range(len(topological))]
    x_powers = neutral_powers_of_x(inputs, degree)
    singlet_change = -rational(f**2*inputs['m0_squared']/3)*truncate_fields(
        singlet_combination(inputs, topological, degree)*taylor_series(theta_shift, x_powers), degree)
    theta_derivative = taylor_series(topological[1:], x_powers)
    pure_alp = rational(f**2*inputs['m0_squared']/3)*2*half_inverse_square*truncate_fields(theta_derivative**2, degree)
    total = truncate_fields(mass_change + singlet_change + pure_alp, degree - 2)
    total = RING({m: c for m, c in total.items() if field_degree(m) > 0})
    alp_three_meson = {m: c for m, c in contracted.items()
                       if m[FIELD_GENERATOR_INDICES[ALP_INDEX]] == 1 and field_degree(m) == 4}
    return {'contracted': contracted, 'residual': contracted - total,
            'alp_three_meson_terms': RING(alp_three_meson)}


# ============================================================== summary
def power_text(power):
    return 'zero' if power is None else str(power)


def lowest(powers):
    values = [p for p in powers if p is not None]
    return min(values) if values else None


def series_powers(series_dict):
    return {k: leading_power(v) for k, v in series_dict.items()}


def pole_mass_summary(external):
    """Leading powers of every class of a pole mass (dressed frame)."""
    coefficients = pole_mass_coefficients('dressed', external)
    labels = [k for k in coefficients['tadpole'] if k != 'alp']
    meson = {}
    for label in labels:
        series = coefficients['tadpole'][label] + coefficients['one_point'][label]
        meson[label] = (leading_power(series), leading_power(truncate(series - epsilon_coefficient(series, 0),
                                                                        TWO_POINT_ORDER)))
    vertices = coefficients['bubble_vertex']
    p2_zero = list(pseudoscalar_invariant(2)) + [topological_invariant(2)]
    return {
        'meson_lines': meson,
        'alp_line': leading_power(coefficients['tadpole']['alp'] + coefficients['one_point']['alp']),
        'alp_line_tadpole_cp_even': leading_power(set_to_zero(coefficients['tadpole']['alp'], CP_ODD_INVARIANTS)),
        'one_point': lowest(leading_power(v) for v in coefficients['one_point'].values()),
        'bubble_alp_alp': doubled(leading_power(vertices[('alp', 'alp')])),
        'bubble_alp_alp_without_P2_Theta2': doubled(leading_power(set_to_zero(vertices[('alp', 'alp')], p2_zero))),
        'bubble_alp_meson': lowest(doubled(leading_power(v)) for k, v in vertices.items() if 'alp' in k and k != ('alp', 'alp')),
        'bubble_meson_meson': lowest(doubled(leading_power(v)) for k, v in vertices.items() if 'alp' not in k),
        'cubic_vertices_cp_even_zero': all(not set_to_zero(v, CP_ODD_INVARIANTS) for v in vertices.values()),
    }


def table_rows():
    """Rows (candidate, quantity, tree power, per-class loop powers, ratio powers, axial)."""
    rows = []
    tree = tree_pole_check()
    for name in ('pi0', 'eta', 'eta_prime', 4, 6, 8):
        summary = pole_mass_summary(name)
        if isinstance(name, int):
            label = next(p for p, idx in CHARGED_PAIRS.items() if name in idx) + ' (charged)'
            tree_power = None
        else:
            label, tree_power = name, leading_power(tree[name])
        meson_leading = lowest(v[0] for v in summary['meson_lines'].values())
        meson_dependent = lowest(v[1] for v in summary['meson_lines'].values())
        loops = (f"A0(m_P): {power_text(meson_leading)} [eps-dependent part {power_text(meson_dependent)}]; "
                 f"A0(m_A): {power_text(summary['alp_line'])}; Bbar(A,A): {power_text(summary['bubble_alp_alp'])} "
                 f"[P2=Theta2=0: {power_text(summary['bubble_alp_alp_without_P2_Theta2'])}]; "
                 f"Bbar(A,P): {power_text(summary['bubble_alp_meson'])}; Bbar(P,P'): {power_text(summary['bubble_meson_meson'])}; "
                 f"one-point: {power_text(summary['one_point'])}")
        classes = {'A0(m_P)': meson_leading, 'A0(m_P) eps-dep.': meson_dependent, 'A0(m_A)': summary['alp_line'],
                   'Bbar(A,A)': summary['bubble_alp_alp']}
        second = 'vs (ii) m_n^2: ' + ', '.join(f'{k} {power_text(p)}' for k, p in classes.items())
        if tree_power is None:
            ratios = 'vs (i): undefined (tree ALP shift is zero); ' + second
        else:
            ratios = ('vs (i) m_n^2-m_n0^2: ' + ', '.join(f'{k} {"-" if p is None else p - tree_power}'
                                                          for k, p in classes.items()) + '; ' + second)
        axial = 'on shell: equal through eps^4' if all(not v for v in frame_comparison(name).values()) else 'on shell: DIFFER'
        rows.append(('C1', f'm_n^2 {label}', 'zero (exact)' if tree_power is None else tree_power, loops, ratios, axial))
    mixing = mixing_entry_tree()
    loop = mixing_entry_one_loop()
    meson_labels = ('pi0', 'eta', 'eta_prime', 'ud', 'us', 'ds')
    entries = loop['entries']['dressed']
    bubbles = [product_power(*v) for v in loop['bubble_powers'].values()]
    difference = loop['on_shell_difference']
    rows.append(('C2', 'D_{pi0 a}', leading_power(mixing['dressed']),
                 f"meson loops: {lowest(leading_power(entries[k]) for k in meson_labels)}; A0(m_A): "
                 f"{power_text(leading_power(entries['alp']))}; bubbles: {power_text(lowest(bubbles))}",
                 'ratio to tree entry: meson 0, A0(m_A) 2',
                 f"tree: kappa-dressed = -eps(kappa_u-kappa_d)(s-m_pi0,0^2), zero at s=m_pi0,0^2, "
                 f"power {power_text(leading_power(mixing['at_pole']))} at the pole; one loop at the pole: meson "
                 f"{power_text(lowest(leading_power(difference[k]) for k in meson_labels))}, A0(m_A) "
                 f"{power_text(leading_power(difference['alp']))} (off shell: differs)"))
    alp_summary = pole_mass_summary('alp')
    absorption = alp_mass_absorption_check()
    rows.append(('C3', 'm_A^2', leading_power(tree['alp']),
                 f"A0(m_P): {power_text(lowest(v[0] for v in alp_summary['meson_lines'].values()))}; A0(m_A): "
                 f"{power_text(alp_summary['alp_line'])} (one-point {power_text(alp_summary['one_point'])}, CP-odd); "
                 f"Bbar(A,P): {power_text(alp_summary['bubble_alp_meson'])}; Bbar(A,A): {power_text(alp_summary['bubble_alp_alp'])}; "
                 f"Bbar(P,P'): {power_text(alp_summary['bubble_meson_meson'])} (bubbles CP-odd)",
                 'vs ALP tree shift: meson 0, A0(m_A) 2',
                 ('on shell: equal through eps^4' if all(not v for v in frame_comparison('alp').values()) else 'DIFFER')
                 + ('; dm_A^2/dm_a0^2 = (v_a^A)^2 exact' if not absorption['difference'] else '; absorption FAILS')))
    amplitude = tree_amplitude()
    one_loop = amplitude_one_loop()
    tree_power = leading_power(amplitude)
    tadpoles, bubble = one_loop['tadpole'], one_loop['bubble']
    meson_names = ('pi0', 'eta', 'eta_prime')
    meson_tadpole = lowest(leading_power(tadpoles[r]) for r in meson_names)
    alp_meson_bubble = lowest(leading_power(bubble[(r, 'alp')]) for r in meson_names)
    meson_bubble = lowest(leading_power(v) for k, v in bubble.items() if 'alp' not in k)
    external = one_loop['external']
    classes = {'A0(m_A)': leading_power(tadpoles['alp']), 'A0(m_P)': meson_tadpole,
               'Bbar(A,P)': alp_meson_bubble, 'Bbar(A,A)': leading_power(bubble[('alp', 'alp')]),
               'Bbar(P,P\')': meson_bubble,
               'legs A0(m_P)': lowest(leading_power(external[r]) for r in meson_names),
               'legs A0(m_A)': leading_power(external['alp'])}
    frame = tree_amplitude_frame_check()
    rows.append(('C4', 'M(A->3pi0)', tree_power,
                 '; '.join(f'{k}: {power_text(v)}' for k, v in classes.items()),
                 'vs tree: ' + ', '.join(f'{k} {"-" if v is None else v - tree_power}' for k, v in classes.items()),
                 'tree on shell: equal exactly (neutral kappa map linear)' if not frame['difference'] else 'tree DIFFER'))
    return rows


def main():
    start = time.time()
    print('Exact leading powers of epsilon = f_pi/f_a (dressed frame, LO U(3), symbolic invariants).')
    print('columns: candidate | quantity | tree power | one-loop class powers | ratio powers | axial')
    for row in table_rows():
        print(' | '.join(str(item) for item in row))
    shift = neutral_quark_mass_shift_check()
    print('C1 quark-mass-shift form: eta residual =', shift['eta_residual'],
          '; Delta_u+Delta_d, Delta_s minus closed form =', shift['sum_ud_minus_closed_form'],
          shift['strange_minus_closed_form'])
    print('   Delta = -(eps^2/(2 f^2)) (S_hat_2 + 4 Y P_hat_1 - 4 Y^2 m_q), Y = M_{m_a^2}^{-1}(C)')
    print('C1 charged remainder minus (eps^2/(2 f^2)) m_a^2 <[Y,lambda]^2>:',
          {pair: charged_remainder_check(pair)['difference'] for pair in CHARGED_PAIRS})
    mixing = mixing_entry_tree()
    print('C2 tree D_kappa - D_dressed =', mixing['difference'].as_expr(),
          '; at the pole, eps^3 coefficient =', epsilon_coefficient(mixing['at_pole'], 3).as_expr())
    conjecture = tree_amplitude_conjecture_check()
    print('C4 tree eps^1 coefficient =', conjecture['first_order'].as_expr(),
          '; minus conjecture:', conjecture['first_order_minus_conjecture'],
          '; even powers zero:', all(not epsilon_coefficient(conjecture['amplitude'], k) for k in (0, 2, 4)))
    print('C4 source-shift identity residual:', source_shift_identity()['residual'])
    print('C4 neutral kappa map linear (residual terms):', len(kappa_map_linearity_check()))
    print('C4 charged internal lines in the amplitude: NOT computed.')
    print(f'runtime {time.time() - start:.1f} s')


if __name__ == '__main__':
    main()
