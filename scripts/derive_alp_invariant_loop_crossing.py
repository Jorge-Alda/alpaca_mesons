"""Frame-invariant one-loop ALP contributions to the meson pole masses and the
values of f_pi/f_a where they reach the size of the tree-level mass.

Conventions (Alda:2026nrd): U = u^2 = exp(i Phi/f_pi), tr(lambda_i lambda_j) = 2 delta_ij,
N_f = 3, x = a/f_a, epsilon = f_pi/f_a. Units: GeV.

The calculation is done with the dressed sources (no ALP-dependent quark currents):
    L = f^2/4 <dU dU^dag> + f^2 B0/2 <M_hat(x) U^dag + U M_hat(x)^dag>
        - f^2 m0^2/(4 N_f) X_TC^2 + (da)^2/2 - m_a0^2 a^2/2,
    X_TC = -tr Phi/f + 2 tr G(x) - 2 Theta_hat(x),
with M_hat(x) = sum_n x^n/n! (S_hat_n + i P_hat_n) and Theta_hat(x) = sum_n x^n/n! Theta_hat_n
the invariant dressed sources.  Assumption: flavour-diagonal dressed sources, which gives
S_hat_0 = m_q, P_hat_0 = 0 and S_hat_1 = -i[c^VS, m_q] = 0.  No values are assumed for any
invariant; numerical values below are test values used only for verification.

A constant axial redefinition with diagonal generator kappa (g_R = g_L^dag = exp(i x kappa)) is
applied as the change of variables U_hat = g^dag U g^dag (G(x) = kappa x, optionally plus
kappa_2 x^2/2).  The one-loop pole shift is computed in both frames.

Loop functions (MSbar, dimensional regularization):
    A0(m^2) = m^2/(16 pi^2) [1 - ln(m^2/mu^2)],   <phi phi> = -A0,
    Bbar(s; m1, m2) = -int_0^1 dx ln[(x m1^2 + (1-x) m2^2 - x(1-x) s - i0)/mu^2].
One-loop pole shift of a pole n with K-normalized tree vector v_n:
    delta s_n = v_n^T [delta M^2(s_n) - s_n delta K] v_n.

Minimal invocation:
    .venv/bin/python scripts/derive_alp_invariant_loop_crossing.py
"""

from __future__ import annotations

import itertools
import math
from fractions import Fraction
from functools import lru_cache

import numpy as np
import sympy as sp
from scipy import integrate

NUMBER_OF_FLAVOURS = 3
FIELD_NAMES = ('P_u', 'P_d', 'P_s', 'a', 'P_1', 'P_2', 'P_4', 'P_5', 'P_6', 'P_7')
NEUTRAL_SECTOR_SIZE = 4          # P_u, P_d, P_s, a: the fields that mix with the ALP
ALP_INDEX = 3
CHARGED_PAIRS = {'ud': (4, 5), 'us': (6, 7), 'ds': (8, 9)}
FLAVOUR_INDEX = {'u': 0, 'd': 1, 's': 2}

TRUNCATION_ORDER = 4
grading = sp.Symbol('grading_parameter', real=True)

meson_fields = sp.symbols(FIELD_NAMES, real=True)
derivative_fields = sp.symbols(['d' + name for name in FIELD_NAMES], real=True)

# ---------------------------------------------------------------- physical inputs (table only)
NEUTRAL_PION_MASS = 0.13498
PION_DECAY_CONSTANT = 0.09207
QUARK_MASSES = (2.2e-3, 4.7e-3, 95e-3)      # default of the repository's Murayama script
RENORMALIZATION_SCALE = 0.77
ALP_MASSES = (0.1, 0.3, 0.5, 1.0, 2.0)


# ============================================================== Lagrangian construction
def meson_generators():
    """Generators lambda_k (tr = 2 delta) in the order of FIELD_NAMES, without the ALP."""
    def unit(i, j):
        matrix = sp.zeros(3, 3)
        matrix[i, j] = 1
        return matrix
    sqrt2 = sp.sqrt(2)
    generators = [sqrt2*unit(0, 0), sqrt2*unit(1, 1), sqrt2*unit(2, 2), None]
    for i, j in ((0, 1), (0, 2), (1, 2)):
        generators.append(unit(i, j) + unit(j, i))
        generators.append(-sp.I*unit(i, j) + sp.I*unit(j, i))
    return generators


def _truncate(expression):
    expanded = sp.expand(expression)
    return sum(expanded.coeff(grading, order)*grading**order for order in range(TRUNCATION_ORDER + 1))


def _truncate_matrix(matrix):
    return matrix.applyfunc(_truncate)


def _product(*matrices):
    result = matrices[0]
    for matrix in matrices[1:]:
        result = _truncate_matrix(result*matrix)
    return result


def build_lagrangian(parameters, first_frame_generator=(0, 0, 0), second_frame_generator=(0, 0, 0),
                     active_fields=FIELD_NAMES):
    """Lagrangian up to fourth order in the fields (1D derivative trick: dF stands for d_mu F,
    and every product of two derivative factors stands for their Lorentz contraction).

    parameters: dict with f, B0, quark_masses (3), m0_squared, bare_alp_mass_squared, f_a,
    dressed_scalar {n: (3 diag)} n>=2, dressed_pseudoscalar {n: (3 diag)} n>=1,
    dressed_topological {n: value} n>=1.  Entries may be sympy symbols or rationals.
    """
    f = parameters['f']
    f_a = parameters['f_a']
    b0 = parameters['B0']
    generators = meson_generators()
    active = [name in active_fields for name in FIELD_NAMES]

    meson_matrix = sp.zeros(3, 3)
    meson_derivative = sp.zeros(3, 3)
    for index, generator in enumerate(generators):
        if generator is None or not active[index]:
            continue
        meson_matrix += grading*meson_fields[index]*generator
        meson_derivative += grading*derivative_fields[index]*generator
    alp = grading*meson_fields[ALP_INDEX] if active[ALP_INDEX] else 0
    alp_derivative = grading*derivative_fields[ALP_INDEX] if active[ALP_INDEX] else 0
    x = alp/f_a

    exponent = sp.I*meson_matrix/f
    powers = [sp.eye(3)]
    for _ in range(TRUNCATION_ORDER):
        powers.append(_product(powers[-1], exponent))
    goldstone = sum((powers[k]/sp.factorial(k) for k in range(TRUNCATION_ORDER + 1)), sp.zeros(3, 3))
    goldstone_derivative = sp.zeros(3, 3)
    for k in range(1, TRUNCATION_ORDER + 1):
        for j in range(k):
            term = _product(powers[j], sp.I*meson_derivative/f, powers[k - 1 - j])
            goldstone_derivative += term/sp.factorial(k)
    goldstone_derivative = _truncate_matrix(goldstone_derivative)

    first = sp.diag(*first_frame_generator)
    second = sp.diag(*second_frame_generator)
    frame_generator = first*x + second*x**2/2
    frame_generator_derivative = (first + second*x)*alp_derivative/f_a

    covariant = goldstone_derivative - sp.I*(frame_generator_derivative*goldstone
                                             + goldstone*frame_generator_derivative)
    covariant = _truncate_matrix(covariant)
    kinetic = f**2/4*_product(covariant, covariant.H).trace()

    scalar = parameters['dressed_scalar']
    pseudoscalar = parameters['dressed_pseudoscalar']
    diagonal = []
    for q in range(3):
        entry = parameters['quark_masses'][q]
        for order in range(1, TRUNCATION_ORDER + 1):
            s_value = scalar.get(order, (0, 0, 0))[q]
            p_value = pseudoscalar.get(order, (0, 0, 0))[q]
            entry += x**order/sp.factorial(order)*(s_value + sp.I*p_value)
        phase = sum(((2*sp.I*frame_generator[q, q])**k/sp.factorial(k)
                     for k in range(TRUNCATION_ORDER + 1)))
        diagonal.append(_truncate(phase*entry))
    mass_source = sp.diag(*diagonal)
    mass = f**2*b0/2*(_product(mass_source, goldstone.H) + _product(goldstone, mass_source.H)).trace()

    topological = parameters['dressed_topological']
    theta_hat = sum((topological.get(order, 0)*x**order/sp.factorial(order)
                     for order in range(1, TRUNCATION_ORDER + 1)))
    topological_combination = -meson_matrix.trace()/f + 2*frame_generator.trace() - 2*theta_hat
    singlet = -f**2*parameters['m0_squared']/(4*NUMBER_OF_FLAVOURS)*topological_combination**2

    alp_terms = alp_derivative**2/2 - parameters['bare_alp_mass_squared']*alp**2/2
    lagrangian = _truncate(kinetic + mass + singlet + alp_terms)
    return sp.expand(lagrangian.subs(grading, 1))


# ============================================================== exact sparse builder
# Neutral variables of the exact builder are N_q = Phi_qq = sqrt2 P_q (no sqrt2 in the
# coefficients); the change P -> N is a constant rescaling, a congruence that leaves poles and
# the coefficients of every loop function in a pole shift unchanged.
EXACT_FIELD_NAMES = ('N_u', 'N_d', 'N_s', 'a', 'P_1', 'P_2', 'P_4', 'P_5', 'P_6', 'P_7')


@lru_cache(maxsize=1)
def _exact_ring():
    from sympy.polys.domains import QQ_I
    from sympy.polys.rings import ring
    names = list(EXACT_FIELD_NAMES) + ['d' + name for name in EXACT_FIELD_NAMES]
    polynomial_ring, *generators = ring(names, QQ_I)
    return polynomial_ring, generators, QQ_I


def build_lagrangian_exact(parameters, first_frame_generator=(0, 0, 0), second_frame_generator=(0, 0, 0)):
    """Same Lagrangian as build_lagrangian, all ten fields, rational parameters only.

    Returns a list of (exponent tuple of length 20, Fraction coefficient) with the field
    variables (N_u, N_d, N_s, a, charged) followed by their derivatives."""
    polynomial_ring, generators, domain = _exact_ring()
    size = len(EXACT_FIELD_NAMES)
    fields, derivatives = generators[:size], generators[size:]

    def constant(value):
        value = sp.nsimplify(value)
        real_part, imaginary_part = sp.Rational(sp.re(value)), sp.Rational(sp.im(value))
        rational = domain.dom
        return polynomial_ring(domain(rational(real_part.p, real_part.q), rational(imaginary_part.p, imaginary_part.q)))

    one, imaginary_unit = constant(1), constant(sp.I)
    zero = polynomial_ring(0)

    def truncate(polynomial):
        return polynomial_ring({monomial: coefficient for monomial, coefficient in polynomial.items()
                                if sum(monomial) <= TRUNCATION_ORDER})

    def multiply(first, second):
        return [[truncate(sum((first[i][k]*second[k][j] for k in range(3)), zero)) for j in range(3)]
                for i in range(3)]

    def add(first, second, factor=None):
        factor = one if factor is None else factor
        return [[first[i][j] + factor*second[i][j] for j in range(3)] for i in range(3)]

    def scale(matrix, factor):
        return [[truncate(factor*matrix[i][j]) for j in range(3)] for i in range(3)]

    def dagger(matrix):
        return [[polynomial_ring({m: domain(c.x, -c.y) for m, c in matrix[j][i].items()}) for j in range(3)]
                for i in range(3)]

    def trace(matrix):
        return matrix[0][0] + matrix[1][1] + matrix[2][2]

    def zero_matrix():
        return [[zero]*3 for _ in range(3)]

    def identity():
        matrix = zero_matrix()
        for i in range(3):
            matrix[i][i] = one
        return matrix

    f = sp.Rational(parameters['f'])
    f_a = sp.Rational(parameters['f_a'])
    b0 = sp.Rational(parameters['B0'])
    meson_matrix, meson_derivative = zero_matrix(), zero_matrix()
    for q in range(3):
        meson_matrix[q][q] += fields[q]
        meson_derivative[q][q] += derivatives[q]
    for offset, (i, j) in enumerate(((0, 1), (0, 2), (1, 2))):
        real_index, imaginary_index = 4 + 2*offset, 5 + 2*offset
        meson_matrix[i][j] += fields[real_index] - imaginary_unit*fields[imaginary_index]
        meson_matrix[j][i] += fields[real_index] + imaginary_unit*fields[imaginary_index]
        meson_derivative[i][j] += derivatives[real_index] - imaginary_unit*derivatives[imaginary_index]
        meson_derivative[j][i] += derivatives[real_index] + imaginary_unit*derivatives[imaginary_index]
    x = constant(1/f_a)*fields[ALP_INDEX]
    alp_derivative = derivatives[ALP_INDEX]

    exponent = scale(meson_matrix, imaginary_unit*constant(1/f))
    exponent_derivative = scale(meson_derivative, imaginary_unit*constant(1/f))
    powers = [identity()]
    for _ in range(TRUNCATION_ORDER):
        powers.append(multiply(powers[-1], exponent))
    goldstone = zero_matrix()
    for k in range(TRUNCATION_ORDER + 1):
        goldstone = add(goldstone, scale(powers[k], constant(sp.Rational(1, math.factorial(k)))))
    goldstone_derivative = zero_matrix()
    for k in range(1, TRUNCATION_ORDER + 1):
        for j in range(k):
            term = multiply(multiply(powers[j], exponent_derivative), powers[k - 1 - j])
            goldstone_derivative = add(goldstone_derivative, scale(term, constant(sp.Rational(1, math.factorial(k)))))

    frame_generator = zero_matrix()
    frame_generator_derivative = zero_matrix()
    for q in range(3):
        first = constant(first_frame_generator[q])
        second = constant(second_frame_generator[q])
        frame_generator[q][q] = truncate(first*x + second*x**2*constant(sp.Rational(1, 2)))
        frame_generator_derivative[q][q] = truncate((first + second*x)*alp_derivative*constant(1/f_a))
    covariant = add(goldstone_derivative,
                    add(multiply(frame_generator_derivative, goldstone),
                        multiply(goldstone, frame_generator_derivative)), -imaginary_unit)
    kinetic = constant(f**2/4)*trace(multiply(covariant, dagger(covariant)))

    scalar = parameters['dressed_scalar']
    pseudoscalar = parameters['dressed_pseudoscalar']
    mass_source = zero_matrix()
    for q in range(3):
        entry = constant(parameters['quark_masses'][q])
        x_power = one
        for order in range(1, TRUNCATION_ORDER + 1):
            x_power = truncate(x_power*x)
            value = sp.Rational(scalar.get(order, (0, 0, 0))[q]) + sp.I*sp.Rational(pseudoscalar.get(order, (0, 0, 0))[q])
            entry += constant(value/math.factorial(order))*x_power
        phase, phase_power = one, one
        for k in range(1, TRUNCATION_ORDER + 1):
            phase_power = truncate(phase_power*constant(2*sp.I)*frame_generator[q][q])
            phase += constant(sp.Rational(1, math.factorial(k)))*phase_power
        mass_source[q][q] = truncate(phase*entry)
    mass = constant(f**2*b0/2)*trace(add(multiply(mass_source, dagger(goldstone)),
                                         multiply(goldstone, dagger(mass_source))))

    topological = parameters['dressed_topological']
    theta_hat, x_power = zero, one
    for order in range(1, TRUNCATION_ORDER + 1):
        x_power = truncate(x_power*x)
        theta_hat += constant(sp.Rational(topological.get(order, 0))/math.factorial(order))*x_power
    combination = -constant(1/f)*trace(meson_matrix) + constant(2)*trace(frame_generator) - constant(2)*theta_hat
    singlet = -constant(f**2*sp.Rational(parameters['m0_squared'])/(4*NUMBER_OF_FLAVOURS))*truncate(combination**2)
    alp_terms = (constant(sp.Rational(1, 2))*alp_derivative**2
                 - constant(sp.Rational(parameters['bare_alp_mass_squared'])/2)*fields[ALP_INDEX]**2)
    lagrangian = truncate(kinetic + mass + singlet + alp_terms)
    terms = []
    for monomial, coefficient in lagrangian.items():
        if sum(monomial) == 0:
            continue
        assert coefficient.y == 0, (monomial, coefficient)
        real_part = coefficient.x
        terms.append((tuple(monomial), Fraction(int(real_part.numerator), int(real_part.denominator))))
    return terms


def sympy_terms(lagrangian):
    """(exponent tuple, coefficient) pairs of a sympy Lagrangian of build_lagrangian."""
    variables = list(meson_fields) + list(derivative_fields)
    terms = []
    for term in sp.Add.make_args(lagrangian):
        coefficient, monomial = term.as_independent(*variables, as_Add=False)
        powers = monomial.as_powers_dict()
        exponents = [0]*len(variables)
        for variable, power in powers.items():
            if variable != 1:
                exponents[variables.index(variable)] = int(power)
        if sum(exponents):
            terms.append((tuple(exponents), coefficient))
    return terms


# ============================================================== vertex extraction
def extract_vertices(terms, numeric=True):
    """Return dict with K, M2 (10x10), cubic V3 (10^3), quartic W (10^4, no derivatives),
    Z (10^4: Z[I,J,K,L] = d^4 L/d dF_I d dF_J dF_K dF_L), and flags.

    terms: list of (exponent tuple over fields then derivatives, coefficient)."""
    size = len(FIELD_NAMES)
    def convert(value):
        if numeric:
            return float(value)
        if isinstance(value, Fraction):
            return sp.Rational(value.numerator, value.denominator)
        return sp.nsimplify(value)
    zero = 0.0 if numeric else sp.Integer(0)
    kinetic = [[zero]*size for _ in range(size)]
    mass = [[zero]*size for _ in range(size)]
    cubic = {}
    quartic = {}
    derivative_quartic = {}
    flags = {'linear_terms': [], 'field_derivative_quadratic': [], 'derivative_cubic': [],
             'odd_derivative_quartic': []}
    for exponents, coefficient in terms:
        field_list, derivative_list = [], []
        multiplicity = 1
        for index, power in enumerate(exponents):
            if power:
                target = field_list if index < size else derivative_list
                target.extend([index % size]*power)
                multiplicity *= math.factorial(power)
        monomial = exponents
        degree = len(field_list) + len(derivative_list)
        value = coefficient*multiplicity
        if degree == 0:
            continue
        if degree == 1:
            flags['linear_terms'].append(monomial)
        elif degree == 2:
            if len(derivative_list) == 2:
                i, j = derivative_list
                if i == j:
                    kinetic[i][i] += convert(value)
                else:
                    kinetic[i][j] += convert(value)
                    kinetic[j][i] += convert(value)
            elif len(field_list) == 2:
                i, j = field_list
                if i == j:
                    mass[i][i] -= convert(value)
                else:
                    mass[i][j] -= convert(value)
                    mass[j][i] -= convert(value)
            else:
                flags['field_derivative_quadratic'].append(monomial)
        elif degree == 3:
            if derivative_list:
                flags['derivative_cubic'].append((monomial, coefficient))
                continue
            for permutation in set(itertools.permutations(field_list)):
                cubic[permutation] = cubic.get(permutation, zero) + convert(value)
        elif degree == 4:
            if not derivative_list:
                for permutation in set(itertools.permutations(field_list)):
                    quartic[permutation] = quartic.get(permutation, zero) + convert(value)
            elif len(derivative_list) == 2:
                for first in set(itertools.permutations(derivative_list)):
                    for second in set(itertools.permutations(field_list)):
                        key = first + second
                        derivative_quartic[key] = derivative_quartic.get(key, zero) + convert(value)
            else:
                flags['odd_derivative_quartic'].append(monomial)
    if numeric:
        def dense(table, rank):
            array = np.zeros((size,)*rank)
            for key, value in table.items():
                array[key] = value
            return array
        return {'K': np.array(kinetic), 'M2': np.array(mass), 'V3': dense(cubic, 3),
                'W': dense(quartic, 4), 'Z': dense(derivative_quartic, 4), 'flags': flags}
    return {'K': sp.Matrix(kinetic), 'M2': sp.Matrix(mass), 'V3': cubic, 'W': quartic,
            'Z': derivative_quartic, 'flags': flags}


# ============================================================== tree spectrum
def tree_spectrum(vertices):
    """Poles and K-normalized pole vectors (10-component) of D(s) = sK - M^2."""
    kinetic, mass = vertices['K'], vertices['M2']
    size = len(FIELD_NAMES)
    neutral = slice(0, NEUTRAL_SECTOR_SIZE)
    assert np.allclose(kinetic[neutral, NEUTRAL_SECTOR_SIZE:], 0)
    assert np.allclose(mass[neutral, NEUTRAL_SECTOR_SIZE:], 0)
    cholesky = np.linalg.cholesky(kinetic[neutral, neutral])
    inverse = np.linalg.inv(cholesky)
    reduced = inverse @ mass[neutral, neutral] @ inverse.T
    values, vectors = np.linalg.eigh(reduced)
    pole_vectors = []
    pole_masses = []
    for k in range(NEUTRAL_SECTOR_SIZE):
        vector = np.zeros(size)
        vector[neutral] = inverse.T @ vectors[:, k]
        pole_vectors.append(vector)
        pole_masses.append(values[k])
    charged_kinetic = kinetic[NEUTRAL_SECTOR_SIZE:, NEUTRAL_SECTOR_SIZE:]
    charged_mass = mass[NEUTRAL_SECTOR_SIZE:, NEUTRAL_SECTOR_SIZE:]
    assert np.allclose(charged_kinetic, np.eye(size - NEUTRAL_SECTOR_SIZE))
    assert np.allclose(charged_mass, np.diag(np.diag(charged_mass)))
    for k in range(size - NEUTRAL_SECTOR_SIZE):
        vector = np.zeros(size)
        vector[NEUTRAL_SECTOR_SIZE + k] = 1.0
        pole_vectors.append(vector)
        pole_masses.append(charged_mass[k, k])
    return np.array(pole_masses), np.array(pole_vectors)


def alp_like_pole(pole_vectors):
    """Index of the neutral pole with the largest ALP component."""
    return int(np.argmax([abs(pole_vectors[k][ALP_INDEX]) for k in range(NEUTRAL_SECTOR_SIZE)]))


# ============================================================== one-loop pole shift
def tadpole_effective_matrices(vertices, field_contraction, derivative_contraction):
    """delta M^2 and delta K from contracting two fields of the quartic terms.

    field_contraction G_KL = <F_K F_L>, derivative_contraction H_KL = <dF_K dF_L>."""
    quartic, derivative_quartic = vertices['W'], vertices['Z']
    delta_mass = (-0.5*np.einsum('ijkl,kl->ij', quartic, field_contraction)
                  - 0.5*np.einsum('klij,kl->ij', derivative_quartic, derivative_contraction))
    delta_kinetic = 0.5*np.einsum('ijkl,kl->ij', derivative_quartic, field_contraction)
    return delta_mass, delta_kinetic


def pole_shift_coefficients(vertices, external_pole):
    """Coefficients of each loop function in the one-loop shift of the pole external_pole.

    Returns dict: 'tadpole'[r] (coefficient of A0(m_r^2), contact tadpoles),
    'one_point'[r] (coefficient of A0(m_r^2) from the loop-induced one-point function),
    'bubble'[(r, t)] r<=t (coefficient of Bbar(s_n; m_r, m_t))."""
    pole_masses, pole_vectors = tree_spectrum(vertices)
    external = pole_vectors[external_pole]
    s_value = pole_masses[external_pole]
    cubic = vertices['V3']
    size = len(pole_masses)
    tadpole = np.zeros(size)
    one_point = np.zeros(size)
    inverse_mass = np.linalg.inv(vertices['M2'])
    for r in range(size):
        projector = np.outer(pole_vectors[r], pole_vectors[r])
        delta_mass, delta_kinetic = tadpole_effective_matrices(
            vertices, -projector, -pole_masses[r]*projector)
        tadpole[r] = external @ (delta_mass - s_value*delta_kinetic) @ external
        linear_source = -0.5*np.einsum('klm,lm->k', cubic, projector)
        expectation = inverse_mass @ linear_source
        one_point[r] = -np.einsum('i,j,ijk,k->', external, external, cubic, expectation)
    bubble = {}
    for r in range(size):
        for t in range(r, size):
            vertex = np.einsum('i,j,k,ijk->', external, pole_vectors[r], pole_vectors[t], cubic)
            multiplicity = 1 if r == t else 2
            bubble[(r, t)] = -multiplicity*vertex**2/(32*np.pi**2)
    return {'tadpole': tadpole, 'one_point': one_point, 'bubble': bubble,
            'pole_masses': pole_masses, 'pole_vectors': pole_vectors}


# ============================================================== closed forms (dressed invariants)
def neutral_meson_basis(parameters):
    """Unmixed neutral diagonal mesons: flavour-basis matrix chi_q delta + m0^2/3, eigenvectors."""
    chi = [2*float(parameters['B0'])*float(m) for m in parameters['quark_masses']]
    m0_squared = float(parameters['m0_squared'])
    matrix = np.diag(chi) + m0_squared/NUMBER_OF_FLAVOURS*np.ones((3, 3))
    values, vectors = np.linalg.eigh(matrix)
    return values, vectors.T          # rows: flavour vectors e_n, (lambda_n)_qq = sqrt2 e_nq


def tadpole_source(parameters):
    """Alda:2026nrd's C = B0 P_hat_1 - m0^2/N_f Theta_hat_1 (diagonal entries)."""
    b0 = float(parameters['B0'])
    return np.array([b0*float(p) for p in parameters['dressed_pseudoscalar'][1]]) \
        - float(parameters['m0_squared'])/NUMBER_OF_FLAVOURS*float(parameters['dressed_topological'][1])


def first_order_solution(parameters, s_value):
    """Diagonal Y solving M_s(Y) = C: (chi_q - s) Y_q + m0^2/N_f tr Y = C_q."""
    chi = np.array([2*float(parameters['B0'])*float(m) for m in parameters['quark_masses']])
    m0_squared = float(parameters['m0_squared'])
    source = tadpole_source(parameters)
    trace = np.sum(source/(chi - s_value))/(1 + m0_squared/NUMBER_OF_FLAVOURS*np.sum(1/(chi - s_value)))
    return (source - m0_squared/NUMBER_OF_FLAVOURS*trace)/(chi - s_value)


def aligned_scalar_source(parameters, s_value):
    """S_hat_2 + 4 Y P_hat_1 - 4 Y^2 m_q (diagonal), Y = M_s^{-1}(C)."""
    y_value = first_order_solution(parameters, s_value)
    scalar = np.array([float(v) for v in parameters['dressed_scalar'][2]])
    pseudo = np.array([float(v) for v in parameters['dressed_pseudoscalar'][1]])
    masses = np.array([float(v) for v in parameters['quark_masses']])
    return scalar + 4*y_value*pseudo - 4*y_value**2*masses


def generator_matrix(meson):
    """3x3 generator lambda of a meson: 'ud', 'us', 'ds' give lambda_1, lambda_4, lambda_6 (one
    real component of the charged pair); a flavour vector e (3,) gives diag(sqrt2 e)."""
    if isinstance(meson, str):
        q, r = FLAVOUR_INDEX[meson[0]], FLAVOUR_INDEX[meson[1]]
        matrix = np.zeros((3, 3))
        matrix[q, r] = matrix[r, q] = 1.0
        return matrix
    return np.diag(math.sqrt(2)*np.asarray(meson, dtype=float))


def alp_line_coefficient(parameters, meson, alp_mass_squared):
    """O(epsilon^2) coefficient of A0(m_a^2) in the one-loop pole shift of a meson that does not
    mix with the ALP at leading order (external leg unmixed):
        -(epsilon^2/(2 f^2)) [ B0 <lambda^2 S'> - m_a^2 <[Y, lambda]^2> ],
    S' = S_hat_2 + 4 Y P_hat_1 - 4 Y^2 m_q,  Y = M_{m_a^2}^{-1}(C)."""
    epsilon = float(parameters['f'])/float(parameters['f_a'])
    f = float(parameters['f'])
    b0 = float(parameters['B0'])
    generator = generator_matrix(meson)
    combination = np.diag(aligned_scalar_source(parameters, alp_mass_squared))
    y_matrix = np.diag(first_order_solution(parameters, alp_mass_squared))
    commutator = y_matrix @ generator - generator @ y_matrix
    return -epsilon**2/(2*f**2)*(b0*np.trace(generator @ generator @ combination)
                                 - alp_mass_squared*np.trace(commutator @ commutator))


def alp_pair_bubble_coefficient(parameters, flavour_vector):
    """O(epsilon^4) coefficient of Bbar(m_n^2; m_a, m_a): -V_naa^2/(32 pi^2),
    V_naa = (epsilon^2/f) <lambda_n (B0 P_hat_2 - m0^2/N_f Theta_hat_2)>."""
    epsilon = float(parameters['f'])/float(parameters['f_a'])
    f = float(parameters['f'])
    b0 = float(parameters['B0'])
    source = np.array([b0*float(p) for p in parameters['dressed_pseudoscalar'][2]]) \
        - float(parameters['m0_squared'])/NUMBER_OF_FLAVOURS*float(parameters['dressed_topological'][2])
    vertex = epsilon**2/f*np.sum(math.sqrt(2)*np.asarray(flavour_vector)*source)
    return -vertex**2/(32*np.pi**2)


def meson_line_coefficients(parameters, external, alp_mass_squared, external_squared_mass):
    """Coefficients of A0(m_P^2) (meson lines) in the pole shift of the neutral diagonal meson with
    index `external` in neutral_meson_basis, through O(epsilon^2).

    alpha_j = mu_j/(m_j^2 - m_a^2) (ALP component of meson pole j), mu_j = -epsilon tr(lambda_j C),
    beta_jk = alpha_j mu_k/(m_j^2 - m_k^2) (second-order component of pole j on meson k).
    Returns ({neutral index P: coefficient}, {charged pair: coefficient summed over its two real
    components}).  external_squared_mass enters the momentum-dependent charged-line vertices."""
    epsilon = float(parameters['f'])/float(parameters['f_a'])
    f = float(parameters['f'])
    b0 = float(parameters['B0'])
    masses = np.array([float(v) for v in parameters['quark_masses']])
    pseudo = np.array([float(v) for v in parameters['dressed_pseudoscalar'][1]])
    values, vectors = neutral_meson_basis(parameters)
    diagonal = math.sqrt(2)*vectors                    # (lambda_j)_qq
    source = tadpole_source(parameters)
    mixing = -epsilon*diagonal @ source                # mu_j
    alpha = mixing/(values - alp_mass_squared)
    count = len(values)
    beta = np.zeros((count, count))
    for j in range(count):
        for k in range(count):
            if j != k:
                beta[j, k] = alpha[j]*mixing[k]/(values[j] - values[k])

    def mesons_only(i, j, k, l):
        return b0/f**2*np.sum(diagonal[i]*diagonal[j]*diagonal[k]*diagonal[l]*masses)

    def one_alp(i, j, k):
        return -epsilon*b0/f**2*np.sum(diagonal[i]*diagonal[j]*diagonal[k]*pseudo)

    n = external
    neutral = {}
    for p in range(count):
        value = mesons_only(n, n, p, p)*(1 - alpha[n]**2 - alpha[p]**2)
        value += 2*alpha[n]*one_alp(n, p, p) + 2*alpha[p]*one_alp(n, n, p)
        value += 2*sum(beta[n, k]*mesons_only(k, n, p, p) for k in range(count) if k != n)
        value += 2*sum(beta[p, k]*mesons_only(n, n, p, k) for k in range(count) if k != p)
        neutral[p] = 0.5*value
    vector = np.zeros(4)                               # (P_u, P_d, P_s, a), first order + second order
    vector[:3] = (1 - alpha[n]**2/2)*vectors[n] + sum(beta[n, k]*vectors[k] for k in range(count) if k != n)
    vector[3] = alpha[n]
    s_value = external_squared_mass
    charged = {}
    for pair in CHARGED_PAIRS:
        q, r = FLAVOUR_INDEX[pair[0]], FLAVOUR_INDEX[pair[1]]
        matrix = np.zeros((4, 4))                       # delta M^2_eff(s) per unit A0(m_qr^2)
        for first, second in ((q, r), (r, q)):
            matrix[first, first] = -(s_value - 2*b0*masses[first])/(3*f**2)
            matrix[first, ALP_INDEX] = matrix[ALP_INDEX, first] = \
                -2*math.sqrt(2)*epsilon*b0/(3*f**2)*(pseudo[first] + pseudo[second]/2)
        matrix[q, r] = matrix[r, q] = (s_value + 2*b0*(masses[q] + masses[r]))/(3*f**2)
        charged[pair] = vector @ matrix @ vector
    return neutral, charged


# ============================================================== ALP-line coefficient (symbolic)
def derive_charged_alp_line_coefficient():
    """Symbolic O(epsilon^2) coefficient of A0(m_a^2) in the pole of the ud meson (field P_1),
    from the dressed vertices, and its reduction to the pole-aligned form
        -(1/(2 f_a^2)) [B0 (S'_u + S'_d) + 2 m_a^2 (Y_u - Y_d)^2]
    after B0 (P_hat_1u - P_hat_1d) is eliminated with the defining equation of Y.

    The internal ALP-like pole vector at O(epsilon) is e = (sqrt2 eps Y_u, sqrt2 eps Y_d,
    sqrt2 eps Y_s, 1) in (P_u, P_d, P_s, a); the external leg is P_1 at s = B0(m_u + m_d)."""
    f, b0, f_a, alp_mass_squared = sp.symbols('f B0 f_a m_a2', positive=True)
    masses = sp.symbols('m_u m_d m_s', positive=True)
    scalar2 = sp.symbols('S2_u S2_d S2_s', real=True)
    pseudo1 = sp.symbols('P1_u P1_d P1_s', real=True)
    y_values = sp.symbols('Y_u Y_d Y_s', real=True)
    parameters = {'f': f, 'f_a': f_a, 'B0': b0, 'quark_masses': masses, 'm0_squared': 0,
                  'bare_alp_mass_squared': 0, 'dressed_scalar': {2: scalar2},
                  'dressed_pseudoscalar': {1: pseudo1}, 'dressed_topological': {}}
    active = ('P_u', 'P_d', 'P_s', 'a', 'P_1')
    lagrangian = build_lagrangian(parameters, active_fields=active)
    epsilon = f/f_a
    direction = [sp.sqrt(2)*epsilon*y_values[0], sp.sqrt(2)*epsilon*y_values[1],
                 sp.sqrt(2)*epsilon*y_values[2], 1]
    external, external_derivative = meson_fields[4], derivative_fields[4]
    poly = sp.Poly(lagrangian, *(list(meson_fields[:5]) + list(derivative_fields[:5])))
    quartic = sum(coefficient*sp.prod([v**p for v, p in zip(poly.gens, monomial)])
                  for monomial, coefficient in poly.terms() if sum(monomial) == 4)
    effective = 0
    for first in range(4):
        for second in range(4):
            weight = direction[first]*direction[second]
            effective += -sp.Rational(1, 2)*weight*sp.diff(quartic, meson_fields[first], meson_fields[second])
            effective += -sp.Rational(1, 2)*alp_mass_squared*weight*sp.diff(
                quartic, derivative_fields[first], derivative_fields[second])
    effective = sp.expand(effective)
    zero_fields = {v: 0 for v in list(meson_fields[:5]) + list(derivative_fields[:5])}
    mass_part = -2*effective.coeff(external, 2).subs(zero_fields)
    kinetic_part = 2*effective.coeff(external_derivative, 2).subs(zero_fields)
    s_value = b0*(masses[0] + masses[1])
    coefficient = sp.expand(mass_part - s_value*kinetic_part)
    aligned = [scalar2[q] + 4*y_values[q]*pseudo1[q] - 4*y_values[q]**2*masses[q] for q in range(2)]
    pole_aligned_form = -(b0*(aligned[0] + aligned[1]) + 2*alp_mass_squared*(y_values[0] - y_values[1])**2)/(2*f_a**2)
    # (chi_q - m_a^2) Y_q + m0^2/3 (tr Y + Theta_1) = B0 P_hat_1q  =>  difference of q = u, d:
    defining_equation = {pseudo1[0]: pseudo1[1] + ((2*b0*masses[0] - alp_mass_squared)*y_values[0]
                                                   - (2*b0*masses[1] - alp_mass_squared)*y_values[1])/b0}
    residual = sp.simplify(sp.expand((coefficient - pole_aligned_form).subs(defining_equation)))
    return {'coefficient': coefficient, 'pole_aligned_form': pole_aligned_form, 'residual': residual}



# ============================================================== loop functions
def loop_function_tadpole(squared_mass, squared_scale):
    """A0(m^2) = m^2/(16 pi^2)[1 - ln(m^2/mu^2)]."""
    return squared_mass/(16*np.pi**2)*(1 - np.log(squared_mass/squared_scale))


def loop_function_bubble(s_value, first_mass, second_mass, squared_scale):
    """Bbar(s; m1, m2) (complex above threshold), Feynman-parameter integral."""
    def real_integrand(fraction):
        argument = fraction*first_mass**2 + (1 - fraction)*second_mass**2 - fraction*(1 - fraction)*s_value
        return -np.log(abs(argument)/squared_scale)

    def imaginary_integrand(fraction):
        argument = fraction*first_mass**2 + (1 - fraction)*second_mass**2 - fraction*(1 - fraction)*s_value
        return np.pi if argument < 0 else 0.0
    real_part = integrate.quad(real_integrand, 0, 1, limit=200)[0]
    imaginary_part = integrate.quad(imaginary_integrand, 0, 1, limit=200)[0]
    return complex(real_part, imaginary_part)


def equal_mass_bubble_below_threshold(s_value, mass, squared_scale):
    """Closed form of Bbar(s; m, m) for 0 < s < 4 m^2."""
    beta = math.sqrt(4*mass**2/s_value - 1)
    return 2 - math.log(mass**2/squared_scale) - 2*beta*math.atan(1/beta)


# ============================================================== test parameter points
def test_parameters(epsilon):
    """Rational test values (verification only; not physics inputs)."""
    f = sp.Rational(23, 250)
    return {'f': f, 'f_a': f/epsilon, 'B0': sp.Rational(13, 5),
            'quark_masses': (sp.Rational(1, 400), sp.Rational(1, 200), sp.Rational(1, 10)),
            'm0_squared': sp.Rational(7, 10), 'bare_alp_mass_squared': sp.Rational(3, 10),
            'dressed_scalar': {2: (sp.Rational(-3, 400), sp.Rational(1, 300), sp.Rational(-1, 20)),
                               3: (sp.Rational(1, 500), sp.Rational(-1, 700), sp.Rational(1, 60)),
                               4: (sp.Rational(1, 90), sp.Rational(1, 80), sp.Rational(-1, 30))},
            'dressed_pseudoscalar': {1: (sp.Rational(-1, 300), sp.Rational(1, 250), sp.Rational(-3, 40)),
                                     2: (sp.Rational(1, 350), sp.Rational(-1, 450), sp.Rational(1, 25)),
                                     3: (sp.Rational(1, 200), sp.Rational(1, 600), sp.Rational(-1, 70)),
                                     4: (sp.Rational(-1, 300), sp.Rational(1, 100), sp.Rational(1, 40))},
            'dressed_topological': {1: sp.Rational(1, 3), 2: sp.Rational(-1, 5),
                                    3: sp.Rational(1, 7), 4: sp.Rational(1, 11)}}


EXAMPLE_FRAME_GENERATOR = (sp.Rational(7, 10), sp.Rational(-2, 5), sp.Rational(9, 10))


# ============================================================== invariance check (symbolic)
def symbolic_frame_difference(epsilon=sp.Rational(1, 3), generator=EXAMPLE_FRAME_GENERATOR):
    """Exact check that the tadpole self-energy in the frame with constant axial generator kappa
    differs from the dressed one by D(s) X + X^T D(s).

    Neutral internal contractions are kept as a generic symmetric matrix (G, H), transformed
    covariantly; charged contractions are symbols A_c, with H_c = m_c^2 A_c (both frames)."""
    parameters = test_parameters(epsilon)
    dressed = extract_vertices(build_lagrangian_exact(parameters), numeric=False)
    framed = extract_vertices(build_lagrangian_exact(parameters, first_frame_generator=generator), numeric=False)
    size = len(FIELD_NAMES)
    s_symbol = sp.Symbol('s')
    neutral_g = sp.Matrix(4, 4, lambda i, j: sp.Symbol(f'G_{min(i, j)}{max(i, j)}'))
    neutral_h = sp.Matrix(4, 4, lambda i, j: sp.Symbol(f'H_{min(i, j)}{max(i, j)}'))
    charged_a = sp.symbols('A_ud A_us A_ds')
    transform = sp.eye(size)          # N_hat_q = N_q - 2 epsilon kappa_q a (exact-builder variables)
    for q in range(3):
        transform[q, ALP_INDEX] = -2*epsilon*generator[q]

    def contractions(neutral_field, neutral_derivative, vertices):
        field_matrix = sp.zeros(size, size)
        derivative_matrix = sp.zeros(size, size)
        field_matrix[:4, :4] = neutral_field
        derivative_matrix[:4, :4] = neutral_derivative
        for pair_index, (first, second) in enumerate(CHARGED_PAIRS.values()):
            charged_mass = vertices['M2'][first, first]
            for index in (first, second):
                field_matrix[index, index] = -charged_a[pair_index]
                derivative_matrix[index, index] = -charged_mass*charged_a[pair_index]
        return field_matrix, derivative_matrix

    def effective(vertices, field_matrix, derivative_matrix):
        delta_mass = sp.zeros(4, 4)
        delta_kinetic = sp.zeros(4, 4)
        for (i, j, k, l), value in vertices['W'].items():
            if i < 4 and j < 4 and field_matrix[k, l] != 0:
                delta_mass[i, j] += -sp.Rational(1, 2)*value*field_matrix[k, l]
        for (k, l, i, j), value in vertices['Z'].items():
            if i < 4 and j < 4 and derivative_matrix[k, l] != 0:
                delta_mass[i, j] += -sp.Rational(1, 2)*value*derivative_matrix[k, l]
            if k < 4 and l < 4 and field_matrix[i, j] != 0:
                delta_kinetic[k, l] += sp.Rational(1, 2)*value*field_matrix[i, j]
        return (delta_mass - s_symbol*delta_kinetic).applyfunc(sp.expand)

    neutral_transform = transform[:4, :4]
    inverse_transform = neutral_transform.inv()
    framed_field = inverse_transform*neutral_g*inverse_transform.T
    framed_derivative = inverse_transform*neutral_h*inverse_transform.T
    dressed_self_energy = effective(dressed, *contractions(neutral_g, neutral_h, dressed))
    framed_self_energy = effective(framed, *contractions(framed_field, framed_derivative, framed))
    pulled_back = (inverse_transform.T*framed_self_energy*inverse_transform).applyfunc(sp.expand)
    difference = (pulled_back - dressed_self_energy).applyfunc(sp.expand)
    inverse_propagator = s_symbol*dressed['K'][:4, :4] - dressed['M2'][:4, :4]
    unknowns = sp.symbols('X0:16')
    x_matrix = sp.Matrix(4, 4, unknowns)
    residual = (difference - inverse_propagator*x_matrix - x_matrix.T*inverse_propagator).applyfunc(sp.expand)
    equations = []
    for entry in residual:
        equations.extend(sp.Poly(entry, s_symbol, *neutral_g.free_symbols, *neutral_h.free_symbols,
                                 *charged_a).coeffs())
    # X is linear in the charged A's: solve with an ansatz X = sum_c A_c X^(c)
    ansatz_symbols = []
    substitution = {}
    for index in range(16):
        components = sp.symbols(f'X{index}_ud X{index}_us X{index}_ds')
        ansatz_symbols.extend(components)
        substitution[unknowns[index]] = sum(c*a for c, a in zip(components, charged_a))
    residual_ansatz = residual.subs(substitution).applyfunc(sp.expand)
    equations = []
    for entry in residual_ansatz:
        if entry != 0:
            equations.extend(sp.Poly(entry, s_symbol, *neutral_g.free_symbols, *neutral_h.free_symbols,
                                     *charged_a).coeffs())
    solution = sp.solve(equations, ansatz_symbols, dict=True)
    return {'difference': difference, 'solution': solution, 'dressed': dressed, 'framed': framed,
            'transform': transform, 'neutral_difference_free_of_G_H':
            all(not (entry.free_symbols & (neutral_g.free_symbols | neutral_h.free_symbols))
                for entry in difference)}


def cubic_congruence_residual(dressed, framed, transform):
    """max |V3_dressed contracted with T T T - V3_frame| (exact tensors as dicts)."""
    size = len(FIELD_NAMES)
    dense_dressed = np.zeros((size,)*3)
    dense_framed = np.zeros((size,)*3)
    for key, value in dressed['V3'].items():
        dense_dressed[key] = float(value)
    for key, value in framed['V3'].items():
        dense_framed[key] = float(value)
    t = np.array(transform, dtype=float)
    pulled = np.einsum('abc,ai,bj,ck->ijk', dense_dressed, t, t, t)
    return float(np.max(np.abs(pulled - dense_framed)))


def raw_tadpole_coefficient_shift(generator=EXAMPLE_FRAME_GENERATOR):
    """Raw sources of the kappa frame built from the dressed frame, and the linear-spurion
    combination -4 m c^A (c^A - 2 c^P) compared with S_hat_2 (symbolic, diagonal entry)."""
    kappa, mass, pseudo1, scalar2 = sp.symbols('kappa m_q P1hat S2hat', real=True)
    x = sp.Symbol('x', real=True)
    dressed_source = mass + sp.I*x*pseudo1 + x**2/2*scalar2
    framed_source = sp.expand(sp.series(sp.exp(2*sp.I*kappa*x)*dressed_source, x, 0, 3).removeO())
    raw_scalar2 = sp.expand(2*sp.re(framed_source.coeff(x, 2)))
    raw_pseudo1 = sp.expand(sp.im(framed_source.coeff(x, 1)))
    axial = kappa                         # C_A of the kappa frame
    pseudo_coupling = raw_pseudo1/(2*mass)  # P_1 = {c^P, m} = 2 m c^P (diagonal)
    linear_spurion = sp.expand(-4*mass*axial*(axial - 2*pseudo_coupling))
    return {'raw_scalar2': raw_scalar2, 'linear_spurion_combination': linear_spurion,
            'sum_equals_invariant': sp.simplify(raw_scalar2 + linear_spurion - scalar2) == 0,
            'linear_spurion_minus_invariant': sp.simplify(linear_spurion - scalar2)}


def second_order_frame_vertex_check(epsilon=sp.Rational(1, 3), generator=(sp.Rational(1, 2), sp.Rational(-3, 4), sp.Rational(1, 5))):
    """Frame with G(x) = kappa_2 x^2/2: the momentum-dependent P_q a a vertex equals the dressed
    one plus [D(s) w]_q, so it is unchanged on shell.  Neutral sector only."""
    parameters = test_parameters(epsilon)
    active = ('P_u', 'P_d', 'P_s', 'a')
    dressed = build_lagrangian(parameters, active_fields=active)
    framed = build_lagrangian(parameters, second_frame_generator=generator, active_fields=active)
    s_symbol = sp.Symbol('s')
    alp, alp_derivative = meson_fields[ALP_INDEX], derivative_fields[ALP_INDEX]
    results = []
    dressed_vertices = extract_vertices(sympy_terms(dressed), numeric=False)
    inverse_propagator = s_symbol*dressed_vertices['K'][:4, :4] - dressed_vertices['M2'][:4, :4]
    vertex_difference = []
    for q in range(3):
        field, derivative = meson_fields[q], derivative_fields[q]
        dressed_vertex = 2*sp.expand(dressed).coeff(field, 1).coeff(alp, 2)
        dressed_vertex = dressed_vertex.subs({v: 0 for v in list(meson_fields) + list(derivative_fields)})
        framed_static = 2*sp.expand(framed).coeff(field, 1).coeff(alp, 2)
        framed_static = framed_static.subs({v: 0 for v in list(meson_fields) + list(derivative_fields)})
        # d P_q . a d a: legs P(p), a(k1), a(k2): -p.(k1 + k2) = p^2 = s
        framed_derivative = sp.expand(framed).coeff(derivative, 1).coeff(alp, 1).coeff(alp_derivative, 1)
        framed_derivative = framed_derivative.subs({v: 0 for v in list(meson_fields) + list(derivative_fields)})
        vertex_difference.append(sp.expand(framed_static + framed_derivative*s_symbol - dressed_vertex))
    w_symbols = sp.symbols('w_u w_d w_s')
    w_vector = sp.Matrix([*w_symbols, 0])
    target = inverse_propagator*w_vector
    solution = sp.solve([sp.Poly(vertex_difference[q] - target[q], s_symbol).coeffs()[i]
                         for q in range(3) for i in range(2)], w_symbols, dict=True)
    results.append(solution)
    return {'vertex_difference': vertex_difference, 'solution': solution}


# ============================================================== numerical table
def table_rows(alp_masses=ALP_MASSES):
    """Loop factors and crossings per unit invariant for the neutral pion (lambda_3 approximation).

    ALP line:  |delta m^2/m^2| = epsilon^2 |A0(m_a^2)| Q_A/(2 f^2), Q_A = |B0<lambda_3^2 S'>|/m_pi^2
               => epsilon_* sqrt(Q_A) = f sqrt(2/|A0|)  (size estimate |A0| -> m_a^2/16pi^2).
    ALP pair:  |delta m^2/m^2| = epsilon^4 m_pi^2 Q_B^2 |Bbar|/(32 pi^2 f^2), Q_B = |B0<lambda_3 P_hat_2>|/m_pi^2
               => epsilon_* sqrt(Q_B) = [32 pi^2 f^2/(m_pi^2 |Bbar|)]^(1/4)  (size estimate |Bbar| -> 1)."""
    squared_scale = RENORMALIZATION_SCALE**2
    b0 = NEUTRAL_PION_MASS**2/(QUARK_MASSES[0] + QUARK_MASSES[1])
    pion_mass_squared = NEUTRAL_PION_MASS**2
    f = PION_DECAY_CONSTANT
    rows = []
    for alp_mass in alp_masses:
        tadpole = loop_function_tadpole(alp_mass**2, squared_scale)
        bubble = loop_function_bubble(pion_mass_squared, alp_mass, alp_mass, squared_scale)
        rows.append({
            'alp_mass': alp_mass,
            'A0': tadpole,
            'alp_line_crossing_with_logarithm': f*math.sqrt(2/abs(tadpole)),
            'alp_line_crossing_size_estimate': 4*math.pi*f*math.sqrt(2)/alp_mass,
            'Bbar': bubble.real,
            'alp_pair_crossing_with_logarithm': (32*math.pi**2*f**2/(pion_mass_squared*abs(bubble)))**0.25,
            'alp_pair_crossing_size_estimate': (32*math.pi**2*f**2/pion_mass_squared)**0.25,
        })
    return b0, rows


def main():
    derivation = derive_charged_alp_line_coefficient()
    print('charged ALP-line coefficient minus pole-aligned form (after the Y equation):', derivation['residual'])
    result = raw_tadpole_coefficient_shift()
    print('kappa frame: raw S_2 =', result['raw_scalar2'])
    print('  linear-spurion combination -4 m c^A(c^A-2c^P) =', result['linear_spurion_combination'])
    print('  raw S_2 + linear-spurion combination == S_hat_2:', result['sum_equals_invariant'])
    b0, rows = table_rows()
    print(f'B0 = {b0:.4f} GeV, f_pi = {PION_DECAY_CONSTANT} GeV, mu = {RENORMALIZATION_SCALE} GeV, '
          f'm_pi0 = {NEUTRAL_PION_MASS} GeV')
    print('m_a[GeV]  A0(m_a^2)[GeV^2]  ALP line: L  S   Bbar(m_pi^2;m_a,m_a)  ALP pair: L  S')
    for row in rows:
        print(f"{row['alp_mass']:5.1f}  {row['A0']:+.6f}  {row['alp_line_crossing_with_logarithm']:.3f}  "
              f"{row['alp_line_crossing_size_estimate']:.3f}  {row['Bbar']:+.4f}  "
              f"{row['alp_pair_crossing_with_logarithm']:.3f}  {row['alp_pair_crossing_size_estimate']:.3f}")


if __name__ == '__main__':
    main()
