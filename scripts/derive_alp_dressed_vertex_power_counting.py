#!/usr/bin/env python3
"""Invariant form of the dressed-frame potential vertices and their powers of eps = f_pi/f_a.

Run:   .venv/bin/python scripts/derive_alp_dressed_vertex_power_counting.py
Tests: .venv/bin/python -m unittest discover -s scripts -p test_derive_alp_dressed_vertex_power_counting.py

Minimal use (run from scripts/; exact rational test inputs, not physical values):
    import sympy as sp
    from derive_alp_dressed_vertex_power_counting import potential_vertex, test_invariants, flavour_generators
    invariants = test_invariants(down_strange=True)
    generators = flavour_generators()
    print(potential_vertex(1, [generators[7], generators[7]], invariants, sp.Symbol('epsilon')))

Conventions (Alda:2026nrd): natural units, masses in GeV; U = u^2 = exp(i Phi/f),
Phi = P_i lambda_i, tr(lambda_i lambda_j) = 2 delta_ij, N_f = 3, x = a/f_a,
eps = f/f_a with f = f_pi. Dressed-frame Lagrangian
    L = f^2/4 <u_mu u^mu> + f^2/4 <chi_+> - f^2 m_0^2/(4 N_f) (<Phi>/f + 2 hat Theta(x))^2
        + (1/2)(d a)^2 - (1/2) m_{a,0}^2 a^2,
    chi = 2 B_0 hat M(x),  hat M(x) = sum_n x^n/n! (S_n + i P_n),  hat Theta(x) = sum_n x^n/n! Theta_n,
with <chi_+> = <chi U^dag + U chi^dag>. Vertices are derivatives of L with respect to the
fields at zero fields (amputated, overall factor i omitted, all momenta incoming).

Inputs (binding rules of the task): S_0 = m_q diagonal with three different entries,
P_0 = 0, Theta_0 = 0. S_1 = -i[c_q^{VS}, m_q] (first-order record
alp-dressed-mass-first-order-invariants.md), built from a test Hermitian c_q^{VS} commuting
with the quark charge: S_1 has zero diagonal (commutator with diagonal m_q), hence <S_1> = 0
follows (no ALP tadpole from the mass term); its only entries are d-s and s-d, and S_1 = 0
when c_q^{VS} is diagonal. Every other invariant is a general Hermitian matrix commuting with
the quark charge (diagonal plus a d-s entry), taken from verification_inputs of the tree
record. No relation among quark masses or invariants is used. All numbers are exact rational
test inputs used only for checks; none is a physical value. No loop integral is evaluated.

Invariant potential vertex (n_a ALP legs, n_P meson legs with generators lambda_1..lambda_nP,
Lambda = sum over the n_P! orderings sigma of lambda_s1 ... lambda_snP):
    V = eps^{n_a} f^{2 - n_P - n_a} [ B_0 (-1)^{floor(n_P/2)} (1/n_P!) <X_{n_a} Lambda> + topological ],
    X_{n_a} = S_{n_a} for even n_P,  X_{n_a} = P_{n_a} for odd n_P,
    topological: n_P = 0:            -(m_0^2/N_f) sum_{k=1}^{n_a-1} binom(n_a, k) Theta_k Theta_{n_a-k}
                 n_P = 1:            -(m_0^2/N_f) Theta_{n_a} <lambda_1>
                 n_P = 2, n_a = 0:   -(m_0^2/(2 N_f)) <lambda_1><lambda_2>
                 otherwise 0,
plus -m_{a,0}^2 for (n_a, n_P) = (2, 0). The kinetic term gives momentum-dependent vertices
with n_a = 0 only; the six-meson kinetic vertex is not covered by L1 and no formula is given.
"""
from __future__ import annotations

import itertools
import math
from functools import lru_cache

import sympy as sp
from sympy.polys.domains import QQ_I
from sympy.polys.matrices import DomainMatrix
from sympy.polys.rings import ring

from derive_alp_finite_mixing_three_meson_tree import GENERATORS as GELL_MANN_GENERATORS
from derive_alp_finite_mixing_three_meson_tree import (
    daughter_pole_series, dressed_source_arrays, inverse_pole_operator, matrix_C, normalization,
    pole_operator, unmixed_neutral_matrix, unmixed_neutral_squared_masses, verification_inputs)
from derive_alp_invariant_loop_crossing import EXACT_FIELD_NAMES, build_lagrangian_exact
from derive_alp_lo_joint_vertices import vertex as joint_vertex
from derive_alp_loop_power_counting import dressed_mass_derivative_closed_form, mass_term_vertex_closed_form

FLAVOURS = 3
ALP_LABEL = 'a'
MESON_FIELD_NAMES = tuple(name for name in EXACT_FIELD_NAMES if name != ALP_LABEL)
EPSILON = sp.Symbol('epsilon', positive=True)
ALP_POLE_SQUARED_MASS = sp.Integer(2)       # m_A^2 (GeV^2), test input as in the tree record
TEST_EPSILON = sp.Rational(1, 3)            # eps used where a rational value is needed


def unit_matrix(row, column):
    """3x3 matrix with a single 1 at (row, column)."""
    matrix = sp.zeros(3)
    matrix[row, column] = 1
    return matrix


def flavour_generators():
    """Meson generators in the field order of the exact fourth-order builder, without the ALP.

    (N_u, N_d, N_s): E_uu, E_dd, E_ss (the builder uses N_q = Phi_qq = sqrt2 P_q);
    (P_1, P_2, P_4, P_5, P_6, P_7): lambda_1, lambda_2, lambda_4, lambda_5, lambda_6, lambda_7.
    """
    generators = [unit_matrix(q, q) for q in range(3)]
    for row, column in ((0, 1), (0, 2), (1, 2)):
        generators.append(unit_matrix(row, column) + unit_matrix(column, row))
        generators.append(-sp.I * unit_matrix(row, column) + sp.I * unit_matrix(column, row))
    return generators


# ---------------------------------------------------------------------------
# Test inputs (exact rationals, checks only)
# ---------------------------------------------------------------------------


def commutator(left, right):
    """[left, right]."""
    return left * right - right * left


def anticommutator(left, right):
    """{left, right}."""
    return left * right + right * left


QUARK_CHARGE = sp.diag(sp.Rational(2, 3), sp.Rational(-1, 3), sp.Rational(-1, 3))


def first_order_scalar_source(quark_mass, down_strange=False):
    """(c_q^{VS}, S_1) with S_1 = -i[c_q^{VS}, m_q] (first-order record), exact test inputs.

    c_q^{VS} = diag(1/3, -1/5, 2/7), plus the d-s entry 3 + 2i (and its conjugate) when
    down_strange=True; it commutes with the quark charge. For diagonal m_q the commutator has
    zero diagonal, so S_1 = 0 without the d-s entry and S_1 has only d-s, s-d entries with it:
    (S_1)_ds = -i (c_q^{VS})_ds (m_s - m_d).
    """
    R = sp.Rational
    vector_scalar_invariant = sp.diag(R(1, 3), R(-1, 5), R(2, 7))
    if down_strange:
        vector_scalar_invariant[1, 2], vector_scalar_invariant[2, 1] = 3 + 2 * sp.I, 3 - 2 * sp.I
    first_order = sp.expand(-sp.I * commutator(vector_scalar_invariant, quark_mass))
    return vector_scalar_invariant, first_order


def test_invariants(down_strange=False, sixth_order=False, bare_alp_mass_squared=sp.Rational(3, 10)):
    """Invariants of verification_inputs (tree record), with S_1 replaced, padded to order 6.

    S_1 = -i[c_q^{VS}, m_q] from first_order_scalar_source (the tree record's diagonal S_1 is
    not of this form and is overridden here): zero diagonal, d-s entries only when
    down_strange=True, S_1 = 0 otherwise.
    With sixth_order=True the higher Taylor coefficients P_4, S_5, P_5, S_6, P_6, Theta_5,
    Theta_6 are switched on (needed only for n_a + n_P = 6; test inputs, flagged), with d-s
    entries when down_strange=True.
    """
    R = sp.Rational
    invariants = verification_inputs(down_strange=down_strange)
    vector_scalar_invariant, first_order = first_order_scalar_source(invariants['S'][0], down_strange)
    invariants['S'][1] = first_order
    invariants['vector_scalar_invariant'] = vector_scalar_invariant
    S = list(invariants['S']) + [sp.zeros(3), sp.zeros(3)]
    P = list(invariants['P']) + [sp.zeros(3), sp.zeros(3)]
    theta = list(invariants['theta']) + [0, 0]
    if sixth_order:
        P[4] = sp.diag(R(1, 7), R(-2, 9), R(1, 4))
        S[5] = sp.diag(R(-1, 6), R(1, 8), R(2, 7))
        P[5] = sp.diag(R(1, 9), R(1, 5), R(-1, 3))
        S[6] = sp.diag(R(1, 10), R(-1, 7), R(1, 5))
        P[6] = sp.diag(R(-1, 4), R(1, 6), R(1, 11))
        theta[5], theta[6] = R(-1, 7), R(1, 8)
        if down_strange:
            for matrix, entry in ((P[4], R(1, 5) - sp.I * R(1, 7)), (S[5], R(-1, 3) + sp.I * R(1, 9)),
                                  (P[5], R(1, 8) + sp.I * R(1, 6)), (S[6], R(1, 7)),
                                  (P[6], -sp.I * R(1, 5))):
                matrix[1, 2], matrix[2, 1] = entry, sp.conjugate(entry)
    invariants.update(S=S, P=P, theta=theta, bare_alp_mass_squared=bare_alp_mass_squared)
    return invariants


def builder_parameters(invariants, eps):
    """Parameter dict of build_lagrangian_exact (diagonal invariants only).

    With diagonal invariants S_1 = 0: the only diagonal matrix with zero diagonal, as required
    by S_1 = -i[c_q^{VS}, m_q]. The builder comparison therefore does not test S_1; the
    d-s comparisons (joint vertex code, direct expansion) do.
    """
    def diagonal(matrix):
        if any(matrix[i, j] != 0 for i in range(3) for j in range(3) if i != j):
            raise ValueError('the fourth-order builder accepts diagonal invariants only')
        return tuple(matrix[q, q] for q in range(3))
    f = invariants['f_pi']
    return {'f': f, 'f_a': f / eps, 'B0': invariants['B_0'],
            'quark_masses': diagonal(invariants['S'][0]), 'm0_squared': invariants['anomaly_mass_squared'],
            'bare_alp_mass_squared': invariants['bare_alp_mass_squared'],
            'dressed_scalar': {n: diagonal(invariants['S'][n]) for n in range(1, 5)},
            'dressed_pseudoscalar': {n: diagonal(invariants['P'][n]) for n in range(1, 5)},
            'dressed_topological': {n: invariants['theta'][n] for n in range(1, 5)}}


def _source(invariants, kind, order):
    values = invariants[kind]
    if order < len(values):
        return values[order]
    return sp.zeros(3) if kind in ('S', 'P') else sp.S.Zero


# ---------------------------------------------------------------------------
# Invariant vertex formula
# ---------------------------------------------------------------------------


def symmetrized_product(matrices):
    """sum over the n! orderings sigma of M_s1 ... M_sn (recursion over subsets of legs)."""
    count = len(matrices)
    table = {0: sp.eye(3)}
    for mask in range(1, 2**count):
        table[mask] = sum((matrices[i] * table[mask & ~(1 << i)] for i in range(count) if mask & (1 << i)),
                          sp.zeros(3))
    return table[2**count - 1]


def vertex_from_symmetrized(number_of_alps, number_of_mesons, symmetrized, generator_traces, invariants, eps,
                            include_bare_alp_mass=True):
    """Pieces {'mass', 'topological', 'bare'} of the potential vertex.

    symmetrized = sum over orderings of the meson generators, generator_traces = [<lambda_i>].
    """
    n_a, n_p = number_of_alps, number_of_mesons
    f, B_0 = invariants['f_pi'], invariants['B_0']
    m0_over_flavours = invariants['anomaly_mass_squared'] / FLAVOURS
    prefactor = eps**n_a * f**(2 - n_p - n_a)
    source = _source(invariants, 'S' if n_p % 2 == 0 else 'P', n_a)
    mass = prefactor * B_0 * (-1)**(n_p // 2) / sp.factorial(n_p) * (source * symmetrized).trace()
    topological = sp.S.Zero
    theta = [_source(invariants, 'theta', k) for k in range(n_a + 1)]
    if n_p == 0:
        topological = -m0_over_flavours * prefactor * sum(sp.binomial(n_a, k) * theta[k] * theta[n_a - k]
                                                         for k in range(n_a + 1))
    elif n_p == 1:
        topological = -m0_over_flavours * prefactor * theta[n_a] * generator_traces[0]
    elif n_p == 2 and n_a == 0:
        topological = -m0_over_flavours / 2 * generator_traces[0] * generator_traces[1]
    bare = -invariants.get('bare_alp_mass_squared', 0) if (n_a, n_p) == (2, 0) and include_bare_alp_mass else 0
    return {'mass': sp.expand(mass), 'topological': sp.expand(topological), 'bare': bare}


def potential_vertex_pieces(number_of_alps, meson_generators, invariants, eps, include_bare_alp_mass=True):
    """Pieces of the potential vertex for n_a ALP legs and the given meson generators."""
    return vertex_from_symmetrized(number_of_alps, len(meson_generators), symmetrized_product(meson_generators),
                                   [g.trace() for g in meson_generators], invariants, eps, include_bare_alp_mass)


def potential_vertex(number_of_alps, meson_generators, invariants, eps, include_bare_alp_mass=True):
    """Invariant potential vertex V_{n_a; lambda_1..lambda_nP} (see module docstring)."""
    return sp.expand(sum(potential_vertex_pieces(number_of_alps, meson_generators, invariants, eps,
                                                 include_bare_alp_mass).values()))


def single_meson_combination(number_of_alps, generator, invariants, eps):
    """eps^n f^{1-n} <lambda (B_0 P_n - (m_0^2/N_f) Theta_n 1)>  (n_P = 1)."""
    n = number_of_alps
    matrix = (invariants['B_0'] * _source(invariants, 'P', n)
              - invariants['anomaly_mass_squared'] / FLAVOURS * _source(invariants, 'theta', n) * sp.eye(3))
    return sp.expand(eps**n * invariants['f_pi']**(1 - n) * (generator * matrix).trace())


def l1_mass_vertex(number_of_alps, meson_generators, invariants, eps):
    """L1 formula (u-form with hat M_n = S_n + i P_n) of the existing power-counting code."""
    f = invariants['f_pi']
    dressed = _source(invariants, 'S', number_of_alps) + sp.I * _source(invariants, 'P', number_of_alps)
    return mass_term_vertex_closed_form(meson_generators, number_of_alps, dressed, f, f / eps, invariants['B_0'])


class SymmetrizedProductTable:
    """Sums over distinct words of fixed letter counts, exact (DomainMatrix over QQ_I).

    W(c) = sum_i G_i W(c - e_i); the sum over all n! orderings of the legs is prod(c_i!) W(c).
    """

    def __init__(self, generators):
        self.generators = [DomainMatrix.from_Matrix(sp.Matrix(g)).convert_to(QQ_I) for g in generators]
        self.table = {tuple([0] * len(generators)): DomainMatrix.eye(3, QQ_I)}

    def words(self, counts):
        counts = tuple(counts)
        if counts not in self.table:
            total = DomainMatrix.zeros((3, 3), QQ_I)
            for i, c in enumerate(counts):
                if c:
                    reduced = list(counts)
                    reduced[i] -= 1
                    total = total + self.generators[i] * self.words(reduced)
            self.table[counts] = total
        return self.table[counts]

    def symmetrized(self, counts):
        return self.words(counts).to_Matrix() * math.prod(math.factorial(c) for c in counts)


def formula_for_counts(table, generator_traces, alp_count, meson_counts, invariants, eps,
                       include_bare_alp_mass=True):
    """Potential vertex for a field monomial given by counts (uses SymmetrizedProductTable)."""
    traces = [t for t, c in zip(generator_traces, meson_counts) for _ in range(c)]
    pieces = vertex_from_symmetrized(alp_count, sum(meson_counts), table.symmetrized(meson_counts), traces,
                                     invariants, eps, include_bare_alp_mass)
    return sp.expand(sum(pieces.values()))


def field_count_tuples(number_of_fields, max_degree, min_degree=1):
    """All exponent tuples with min_degree <= total degree <= max_degree."""
    for degree in range(min_degree, max_degree + 1):
        for combination in itertools.combinations_with_replacement(range(number_of_fields), degree):
            counts = [0] * number_of_fields
            for index in combination:
                counts[index] += 1
            yield tuple(counts)


# ---------------------------------------------------------------------------
# (a), (c): comparison with the existing fourth-order builder (diagonal invariants)
# ---------------------------------------------------------------------------


@lru_cache(maxsize=4)
def _builder_terms(key):
    invariants = test_invariants(down_strange=False)
    if key == 'no_mass_term':
        invariants['B_0'] = 0
    return tuple(build_lagrangian_exact(builder_parameters(invariants, TEST_EPSILON)))


def builder_derivatives(terms, size=len(EXACT_FIELD_NAMES)):
    """{(field exponents, derivative exponents): derivative of L} from builder terms."""
    result = {}
    for exponents, coefficient in terms:
        multiplicity = math.prod(math.factorial(p) for p in exponents)
        result[(tuple(exponents[:size]), tuple(exponents[size:]))] = (
            sp.Rational(coefficient.numerator, coefficient.denominator) * multiplicity)
    return result


def compare_with_builder(max_degree=4):
    """Residuals formula - builder for every non-derivative monomial of degree <= max_degree.

    Ten fields (N_u, N_d, N_s, a, six charged components), diagonal test invariants
    (the builder accepts only diagonal invariants, so S_1 = 0 here), eps = TEST_EPSILON.
    Returns (number of monomials compared, list of (counts, residual) with residual != 0).
    """
    invariants = test_invariants(down_strange=False)
    derivatives = builder_derivatives(_builder_terms('full'))
    alp_index = EXACT_FIELD_NAMES.index(ALP_LABEL)
    generators = flavour_generators()
    table = SymmetrizedProductTable(generators)
    traces = [g.trace() for g in generators]
    zero_derivatives = tuple([0] * len(EXACT_FIELD_NAMES))
    failures, compared = [], 0
    for counts in field_count_tuples(len(EXACT_FIELD_NAMES), max_degree):
        meson_counts = counts[:alp_index] + counts[alp_index + 1:]
        formula = formula_for_counts(table, traces, counts[alp_index], meson_counts, invariants, TEST_EPSILON)
        residual = sp.simplify(formula - derivatives.get((counts, zero_derivatives), 0))
        compared += 1
        if residual != 0:
            failures.append((counts, residual))
    unmatched = [key for key in derivatives if key[1] == zero_derivatives and sum(key[0]) > max_degree]
    return compared, failures, unmatched


def single_meson_check_with_builder(max_alps=3):
    """Builder coefficient of lambda a^n versus eps^n f^{1-n} <lambda (B_0 P_n - m_0^2/N_f Theta_n)>."""
    invariants = test_invariants(down_strange=False)
    derivatives = builder_derivatives(_builder_terms('full'))
    alp_index = EXACT_FIELD_NAMES.index(ALP_LABEL)
    zero_derivatives = tuple([0] * len(EXACT_FIELD_NAMES))
    residuals = {}
    meson_indices = [i for i in range(len(EXACT_FIELD_NAMES)) if i != alp_index]
    for n_a in range(max_alps + 1):
        for position, field_index in enumerate(meson_indices):
            counts = [0] * len(EXACT_FIELD_NAMES)
            counts[field_index], counts[alp_index] = 1, n_a
            value = derivatives.get((tuple(counts), zero_derivatives), 0)
            expected = single_meson_combination(n_a, flavour_generators()[position], invariants, TEST_EPSILON)
            residuals[(n_a, EXACT_FIELD_NAMES[field_index])] = sp.simplify(value - expected)
    return residuals


def kinetic_and_topological_structure():
    """(c) From the builder: ALP fields in derivative monomials, and non-derivative monomials
    of the topological term (builder with B_0 = 0) with n_P >= 2.

    Returns (ALP-carrying derivative monomials other than (d a)^2, topological monomials with
    n_P >= 2 other than n_a = 0, n_P = 2).
    """
    size = len(EXACT_FIELD_NAMES)
    alp_index = EXACT_FIELD_NAMES.index(ALP_LABEL)
    kinetic_alp = []
    for (fields, derivatives), value in builder_derivatives(_builder_terms('full')).items():
        if sum(derivatives) == 0 or value == 0:
            continue
        alp_derivative_square = sum(fields) == 0 and derivatives[alp_index] == 2 and sum(derivatives) == 2
        if (fields[alp_index] or derivatives[alp_index]) and not alp_derivative_square:
            kinetic_alp.append((fields, derivatives, value))
    topological_extra = []
    for (fields, derivatives), value in builder_derivatives(_builder_terms('no_mass_term')).items():
        if sum(derivatives) or value == 0:
            continue
        n_a, n_p = fields[alp_index], sum(fields) - fields[alp_index]
        if n_p >= 2 and (n_a, n_p) != (0, 2):
            topological_extra.append((fields, value))
    return kinetic_alp, topological_extra


# ---------------------------------------------------------------------------
# (a) with d-s entries: existing general joint-vertex code (`vertex`)
# ---------------------------------------------------------------------------


def joint_potential_vertex(alp_count, meson_generators, invariants, eps):
    """Non-derivative vertex from derive_alp_lo_joint_vertices.vertex (all momenta set to zero)."""
    legs = list(meson_generators) + [sp.zeros(3)] * alp_count
    alps = [0] * len(meson_generators) + [1] * alp_count
    f = invariants['f_pi']
    return joint_vertex(legs, alps, [0] * len(legs), *dressed_source_arrays(invariants), f, f / eps,
                        invariants['anomaly_mass_squared'], invariants.get('bare_alp_mass_squared', 0))


# representative generator set for the d-s comparison with `vertex` (slow, about 0.35 s per vertex)
DOWN_STRANGE_REPRESENTATIVE = (1, 2, 3, 7, 8)       # E_dd, E_ss, lambda_1, lambda_6, lambda_7 (indices in flavour_generators())


def compare_with_joint_vertex(max_degree=4, generator_indices=DOWN_STRANGE_REPRESENTATIVE):
    """Residuals formula - vertex() with d-s invariants, for monomials in the given generators + a
    that contain at least one d-s generator (lambda_6 or lambda_7) or are pure ALP."""
    invariants = test_invariants(down_strange=True)
    generators = flavour_generators()
    failures, compared = [], 0
    for counts in field_count_tuples(len(generator_indices) + 1, max_degree):
        meson = [generators[generator_indices[i]] for i, c in enumerate(counts[:-1]) for _ in range(c)]
        labels = [generator_indices[i] for i, c in enumerate(counts[:-1]) for _ in range(c)]
        if meson and not ({7, 8} & set(labels)):
            continue
        formula = potential_vertex(counts[-1], meson, invariants, TEST_EPSILON)
        residual = sp.simplify(formula - joint_potential_vertex(counts[-1], meson, invariants, TEST_EPSILON))
        compared += 1
        if residual != 0:
            failures.append((counts, residual))
    return compared, failures


# ---------------------------------------------------------------------------
# (d), (b): direct Taylor expansion of f^2/4 <chi_+> and the topological term
# ---------------------------------------------------------------------------

SIXTH_ORDER_FIELDS = ('N_u', 'N_d', 'N_s', 'P_1', 'P_6', 'P_7')     # plus the ALP


def direct_potential_expansion(invariants, field_names, order, include_mass=True, include_topological=True):
    """Exact expansion of f^2 B_0/2 <hat M U^dag + U hat M^dag> - f^2 m_0^2/(4N_f)(<Phi>/f + 2 hat Theta)^2.

    x = eps a/f with eps a polynomial variable (not counted in the truncation degree).
    Returns {field counts (mesons..., a): {eps power: derivative of L}}.
    """
    meson_names = list(field_names)
    names = meson_names + ['a', 'eps']
    polynomial_ring, *variables = ring(names, QQ_I)
    mesons, alp, eps = variables[:len(meson_names)], variables[-2], variables[-1]
    size = len(meson_names) + 1
    generator_of = dict(zip(MESON_FIELD_NAMES, flavour_generators()))

    def constant(value):
        return polynomial_ring.ground_new(QQ_I.from_sympy(sp.nsimplify(value)))

    def truncate(polynomial):
        return polynomial_ring({m: c for m, c in polynomial.items() if sum(m[:size]) <= order})

    def multiply(left, right):
        return [[truncate(sum((left[i][k] * right[k][j] for k in range(3)), polynomial_ring(0)))
                 for j in range(3)] for i in range(3)]

    f = invariants['f_pi']
    phi = [[polynomial_ring(0)] * 3 for _ in range(3)]
    for field, name in zip(mesons, meson_names):
        generator = generator_of[name]
        phi = [[phi[i][j] + constant(generator[i, j]) * field for j in range(3)] for i in range(3)]
    step = [[constant(-sp.I / f) * phi[i][j] for j in range(3)] for i in range(3)]
    powers = [[[constant(1 if i == j else 0) for j in range(3)] for i in range(3)]]
    for _ in range(order):
        powers.append(multiply(powers[-1], step))
    x = constant(1 / f) * eps * alp
    x_powers = [polynomial_ring(1)]
    for _ in range(order):
        x_powers.append(truncate(x_powers[-1] * x))
    lagrangian = polynomial_ring(0)
    if include_mass:
        trace_part = polynomial_ring(0)
        for n in range(order + 1):
            source = _source(invariants, 'S', n) + sp.I * _source(invariants, 'P', n)
            for k in range(order + 1 - n):
                term = sum((constant(source[i, j]) * powers[k][j][i] for i in range(3) for j in range(3)),
                           polynomial_ring(0))
                trace_part += truncate(constant(sp.Rational(1, math.factorial(n) * math.factorial(k)))
                                       * x_powers[n] * term)
        conjugate = polynomial_ring({m: QQ_I(c.x, -c.y) for m, c in trace_part.items()})
        lagrangian += constant(f**2 * invariants['B_0'] / 2) * (trace_part + conjugate)
    if include_topological:
        theta = sum((constant(_source(invariants, 'theta', n) / math.factorial(n)) * x_powers[n]
                     for n in range(1, order + 1)), polynomial_ring(0))
        combination = constant(1 / f) * (phi[0][0] + phi[1][1] + phi[2][2]) + constant(2) * theta
        lagrangian += constant(-f**2 * invariants['anomaly_mass_squared'] / (4 * FLAVOURS)) * truncate(
            combination**2)
    result = {}
    for monomial, coefficient in truncate(lagrangian).items():
        counts = tuple(monomial[:size])
        if sum(counts) == 0:
            continue
        value = QQ_I.to_sympy(coefficient) * math.prod(math.factorial(p) for p in counts)
        result.setdefault(counts, {})[monomial[size]] = value
    return result


def compare_with_direct_expansion(invariants, field_names, order, min_degree=1):
    """Homogeneity in eps (b) and equality with the formula (a, d) for every monomial.

    Returns (number compared, non-homogeneous monomials, list of (counts, residual))."""
    expansion = direct_potential_expansion(invariants, field_names, order)
    generator_of = dict(zip(MESON_FIELD_NAMES, flavour_generators()))
    generators = [generator_of[name] for name in field_names]
    table = SymmetrizedProductTable(generators)
    traces = [g.trace() for g in generators]
    non_homogeneous, failures, compared = [], [], 0
    for counts in field_count_tuples(len(field_names) + 1, order, min_degree):
        powers = {p: v for p, v in expansion.get(counts, {}).items() if v != 0}
        if any(p != counts[-1] for p in powers):
            non_homogeneous.append((counts, powers))
        formula = formula_for_counts(table, traces, counts[-1], counts[:-1], invariants, 1,
                                     include_bare_alp_mass=False)
        residual = sp.expand(formula - powers.get(counts[-1], 0))
        compared += 1
        if residual != 0:
            failures.append((counts, residual))
    return compared, non_homogeneous, failures


def topological_monomials_with_two_or_more_mesons(invariants, field_names, order):
    """(c) at sixth order: monomials of the topological term alone with n_P >= 2, other than (0, 2)."""
    expansion = direct_potential_expansion(invariants, field_names, order, include_mass=False)
    return [(counts, powers) for counts, powers in expansion.items()
            if sum(counts[:-1]) >= 2 and (counts[-1], sum(counts[:-1])) != (0, 2)
            and any(v != 0 for v in powers.values())]


def homogeneity_with_joint_vertex(cases):
    """(b) with the existing general code: vertex(f_a = f/eps) / eps^{n_a} has no eps.

    cases: list of (n_a, meson generator indices). Returns {case: leftover eps dependence}."""
    invariants = test_invariants(down_strange=True)
    generators = flavour_generators()
    result = {}
    for alp_count, indices in cases:
        meson = [generators[i] for i in indices]
        value = sp.expand(joint_potential_vertex(alp_count, meson, invariants, EPSILON) / EPSILON**alp_count)
        result[(alp_count, tuple(indices))] = sp.diff(value, EPSILON)
    return result


# ---------------------------------------------------------------------------
# (e) cubic vertices
# ---------------------------------------------------------------------------


def cubic_vertices(invariants, eps, generator_pairs=((0, 0), (7, 7), (1, 7))):
    """V_PPa (from S_1), V_Paa (B_0 P_2 - m_0^2/N_f Theta_2), V_aaa (S_3, Theta_1 Theta_2)."""
    generators = flavour_generators()
    return {'PPa': {pair: potential_vertex(1, [generators[pair[0]], generators[pair[1]]], invariants, eps)
                    for pair in generator_pairs},
            'Paa': {index: potential_vertex(2, [generators[index]], invariants, eps) for index in (0, 1, 2, 7)},
            'aaa': potential_vertex(3, [], invariants, eps)}


def alp_two_meson_cubic_table(invariants, eps):
    """{(i, j): V_{1; lambda_i, lambda_j}} for the 45 unordered pairs i <= j of lambda_0..lambda_8.

    For (n_a, n_P) = (1, 2) the vertex is -(eps B_0/(2f)) <S_1 {lambda_i, lambda_j}> (no
    topological piece). With S_1 having only d-s, s-d entries, it is nonzero only when
    {lambda_i, lambda_j} has a d-s entry.
    """
    return {(i, j): potential_vertex(1, [GELL_MANN_GENERATORS[i], GELL_MANN_GENERATORS[j]], invariants, eps)
            for i, j in itertools.combinations_with_replacement(range(9), 2)}


def alp_two_meson_cubic_expected(first, second, invariants, eps):
    """-(eps B_0/(2f)) <S_1 {first, second}>: the ALP-two-meson vertex written out directly."""
    return sp.expand(-eps * invariants['B_0'] / (2 * invariants['f_pi'])
                     * (invariants['S'][1] * anticommutator(first, second)).trace())


def nonzero_alp_two_meson_pairs(invariants, eps=EPSILON):
    """Sorted list of the pairs (i, j), i <= j, of lambda_0..lambda_8 with V_{1; lambda_i, lambda_j} != 0."""
    return sorted(pair for pair, value in alp_two_meson_cubic_table(invariants, eps).items() if value != 0)


def joint_vertex_cubic_residuals(invariants, eps=TEST_EPSILON):
    """{(i, j): V_{1; lambda_i, lambda_j} - vertex()} for the 45 pairs (independent code)."""
    table = alp_two_meson_cubic_table(invariants, eps)
    return {pair: sp.simplify(value - joint_potential_vertex(
        1, [GELL_MANN_GENERATORS[pair[0]], GELL_MANN_GENERATORS[pair[1]]], invariants, eps))
        for pair, value in table.items()}


def representative_odd_even_entries(order_max=4):
    """In the c_q^A, c_q^P representative (L1): S_n = Re hat M_n, P_n = Im hat M_n.

    Returns {n: (S_n, P_n)} with symbolic diagonal c_q^A, c_q^P, m_q; S_n = 0 for odd n and
    P_n = 0 for even n (so no cubic vertex), Theta linear (Theta_2 = 0) by construction there.
    """
    masses = sp.symbols('m_u m_d m_s', positive=True)
    axial = sp.symbols('cA_u cA_d cA_s', real=True)
    pseudo = sp.symbols('cP_u cP_d cP_s', real=True)
    result = {}
    for n in range(1, order_max + 1):
        matrix = dressed_mass_derivative_closed_form(n, masses, axial, pseudo)
        result[n] = (matrix.applyfunc(lambda z: sp.expand(sp.re(z))), matrix.applyfunc(lambda z: sp.expand(sp.im(z))))
    return result


# ---------------------------------------------------------------------------
# (f) ALP-like pole ends
# ---------------------------------------------------------------------------


def alp_pole_overlap(invariants, alp_pole_squared_mass=ALP_POLE_SQUARED_MASS):
    """Y = M^{-1}_{m_A^2}(C) (eps independent)."""
    return inverse_pole_operator(matrix_C(invariants), invariants, alp_pole_squared_mass)


def alp_pole_leg(invariants, eps, normalized=True, alp_pole_squared_mass=ALP_POLE_SQUARED_MASS):
    """(Phi_A, v_a) = (2 eps Y, 1)/sqrt(N(m_A^2)); without 1/sqrt(N) if normalized=False."""
    Y = alp_pole_overlap(invariants, alp_pole_squared_mass)
    norm = sp.sqrt(normalization(invariants, eps, alp_pole_squared_mass)) if normalized else sp.Integer(1)
    return 2 * eps * Y / norm, 1 / norm


def contract_potential_vertex(legs, invariants, eps, include_bare_alp_mass=True):
    """Multilinear contraction of the potential vertex with legs (meson matrix, ALP component)."""
    total = sp.S.Zero
    count = len(legs)
    for alp_set in itertools.chain.from_iterable(itertools.combinations(range(count), k) for k in range(count + 1)):
        weight = sp.prod([legs[i][1] for i in alp_set])
        if weight == 0:
            continue
        mesons = [legs[i][0] for i in range(count) if i not in alp_set]
        if any(m == sp.zeros(3) for m in mesons):
            continue
        total += weight * potential_vertex(len(alp_set), mesons, invariants, eps, include_bare_alp_mass)
    return sp.expand(total)


def basis_leg(index):
    """Leg of a basis field: meson generator (index in flavour_generators()) or the ALP ('a')."""
    if index == ALP_LABEL:
        return sp.zeros(3), 1
    return flavour_generators()[index], 0


def alp_pole_end_scaling(invariants, other_leg_sets, number_of_ends):
    """For each set of other basis legs (k of them ALP), V(others, A^ends) N^{ends/2} / eps^{k + ends}.

    Returns {other legs: (contraction with unnormalized A legs / eps^{k+ends}, its eps derivative)}.
    The unnormalized legs (2 eps Y, 1) differ from the normalized ones by 1/sqrt(N) each
    (multilinearity, checked separately by normalized_contraction_residual). Only interaction
    vertices (at least three legs) are meaningful: the two-leg terms, including -m_{a,0}^2,
    belong to the exact propagator.
    """
    leg = alp_pole_leg(invariants, EPSILON, normalized=False)
    result = {}
    for others in other_leg_sets:
        legs = [basis_leg(i) for i in others] + [leg] * number_of_ends
        k = sum(1 for i in others if i == ALP_LABEL)
        value = sp.expand(contract_potential_vertex(legs, invariants, EPSILON) / EPSILON**(k + number_of_ends))
        result[tuple(others)] = (value, sp.diff(value, EPSILON))
    return result


def one_end_bracket_residual(invariants, others):
    """V(others, A) sqrt(N) - eps [V(others, a)/eps + V(others, 2Y)]  (should vanish)."""
    Y = alp_pole_overlap(invariants)
    legs = [basis_leg(i) for i in others]
    normalized = contract_potential_vertex(legs + [alp_pole_leg(invariants, EPSILON)], invariants, EPSILON)
    norm = sp.sqrt(normalization(invariants, EPSILON, ALP_POLE_SQUARED_MASS))
    bracket = (contract_potential_vertex(legs + [(sp.zeros(3), 1)], invariants, EPSILON) / EPSILON
               + contract_potential_vertex(legs + [(2 * Y, 0)], invariants, EPSILON))
    return sp.simplify(normalized * norm - EPSILON * bracket)


def normalized_contraction_residual(invariants, others, number_of_ends):
    """V(others, normalized A^ends) - V(others, unnormalized A^ends)/N^{ends/2}."""
    legs = [basis_leg(i) for i in others]
    normalized = contract_potential_vertex(legs + [alp_pole_leg(invariants, EPSILON)] * number_of_ends,
                                           invariants, EPSILON)
    plain = contract_potential_vertex(legs + [alp_pole_leg(invariants, EPSILON, normalized=False)] * number_of_ends,
                                      invariants, EPSILON)
    norm = normalization(invariants, EPSILON, ALP_POLE_SQUARED_MASS)
    return sp.simplify(normalized - plain / norm**sp.Rational(number_of_ends, 2))


def closed_alp_line_potential(first, second, invariants):
    """Potential part of V(first, second, A, A) N / eps^2 (eps independent, checked by the caller)."""
    leg = alp_pole_leg(invariants, EPSILON, normalized=False)
    value = contract_potential_vertex([(first, 0), (second, 0), leg, leg], invariants, EPSILON)
    return sp.expand(value / EPSILON**2)


def rotated_scalar_source(invariants, Y):
    """S_2' = S_2 + 2{Y,P_1} - Y^2 m_q - m_q Y^2 - 2 Y m_q Y (source rotated with g_R = g_L^dag = e^{-ixY})."""
    S, P = invariants['S'], invariants['P']
    m = S[0]
    return S[2] + 2 * anticommutator(Y, P[1]) - Y * Y * m - m * Y * Y - 2 * Y * m * Y


def rotated_pseudoscalar_source(invariants, Y):
    """P_2' = P_2 - 2{Y, S_1}."""
    return invariants['P'][2] - 2 * anticommutator(Y, invariants['S'][1])


def rotated_source_prediction(first, second, invariants, Y):
    """-(B_0/(2 f^2)) <S_2' {first, second}> (brainstormer proposal for the potential closed line)."""
    f, B_0 = invariants['f_pi'], invariants['B_0']
    return sp.expand(-B_0 / (2 * f**2) * (rotated_scalar_source(invariants, Y) * anticommutator(first, second)).trace())


def commuting_class_a_form(generator, invariants, Y):
    """-(B_0/f^2) <lambda^2 (S_2 + 4 Y P_1 - 4 Y^2 m_q)> (V(lambda, lambda, A, A) N/eps^2, commuting case)."""
    f, B_0 = invariants['f_pi'], invariants['B_0']
    S, P = invariants['S'], invariants['P']
    return sp.expand(-B_0 / f**2 * (generator * generator * (S[2] + 4 * Y * P[1] - 4 * Y * Y * S[0])).trace())


def closed_alp_line_potential_matrix_form(first, second, invariants, Y):
    """Non-commuting matrix form of the potential closed line found by the script:
      -(B_0/(2f^2)) <S_2' {a,b}>
      + (B_0/(3f^2)) <(P_1 - {Y, m_q}) ([a,[b,Y]] + [b,[a,Y]])>
      + (B_0/(6f^2)) <m_q ({[Y,[Y,a]], b} + {[Y,[Y,b]], a})>.
    """
    f, B_0 = invariants['f_pi'], invariants['B_0']
    m, P1 = invariants['S'][0], invariants['P'][1]
    a, b = first, second
    double = commutator(a, commutator(b, Y)) + commutator(b, commutator(a, Y))
    term_p = B_0 / (3 * f**2) * ((P1 - anticommutator(Y, m)) * double).trace()
    term_m = B_0 / (6 * f**2) * (m * (anticommutator(commutator(Y, commutator(Y, a)), b)
                                      + anticommutator(commutator(Y, commutator(Y, b)), a))).trace()
    return sp.expand(rotated_source_prediction(a, b, invariants, Y) + term_p + term_m)


MOMENTUM_P, MOMENTUM_K = sp.symbols('p k', real=True)


def closed_alp_line_full(first, second, invariants):
    """V(first(p), second(-p), A(k), A(-k)) N/eps^2 from the existing general code `vertex`, with
    kinetic vertices; p, k formal: p**2 = p.p, k**2 = k.k, p*k = p.k."""
    leg_matrix, leg_alp = alp_pole_leg(invariants, EPSILON, normalized=False)
    f = invariants['f_pi']
    value = joint_vertex([first, second, leg_matrix, leg_matrix], [0, 0, leg_alp, leg_alp],
                         [MOMENTUM_P, -MOMENTUM_P, MOMENTUM_K, -MOMENTUM_K], *dressed_source_arrays(invariants),
                         f, f / EPSILON, invariants['anomaly_mass_squared'], 0)
    return sp.expand(value / EPSILON**2)


def kinetic_closed_line_form(first, second, invariants, Y):
    """(1/(3 f^2)) (p^2 + k^2) <[Y,a][Y,b]>  (kinetic part of V(a,b,A,A) N/eps^2)."""
    f = invariants['f_pi']
    return sp.expand((MOMENTUM_P**2 + MOMENTUM_K**2) / (3 * f**2)
                     * (commutator(Y, first) * commutator(Y, second)).trace())


def on_shell_closed_line_form(first, second, invariants, Y, alp_pole_squared_mass=ALP_POLE_SQUARED_MASS):
    """Potential + kinetic at k^2 = m_A^2, p.k = 0, after M_{m_A^2}(Y) = C:
      -(B_0/(2f^2)) <S_2'{a,b}> + (m_A^2/f^2) <[Y,a][Y,b]>
      + (1/(6f^2)) ( <[Y,a][Y, p^2 b - B_0{m_q,b}]> + (a <-> b) ).
    """
    f, B_0, m = invariants['f_pi'], invariants['B_0'], invariants['S'][0]
    a, b = first, second
    p2 = MOMENTUM_P**2

    def piece(x, y):
        return (commutator(Y, x) * commutator(Y, p2 * y - B_0 * anticommutator(m, y))).trace()
    return sp.expand(rotated_source_prediction(a, b, invariants, Y)
                     + alp_pole_squared_mass / f**2 * (commutator(Y, a) * commutator(Y, b)).trace()
                     + (piece(a, b) + piece(b, a)) / (6 * f**2))


def tadpole_kinematics(expression, alp_pole_squared_mass=ALP_POLE_SQUARED_MASS):
    """Set k.k = m_A^2 and p.k = 0 in a polynomial of the formal momenta p, k."""
    polynomial = sp.Poly(sp.expand(expression), MOMENTUM_P, MOMENTUM_K)
    total = sp.S.Zero
    for (power_p, power_k), coefficient in polynomial.terms():
        if power_k % 2:
            continue
        total += coefficient * MOMENTUM_P**power_p * alp_pole_squared_mass**(power_k // 2)
    return sp.expand(total)


def closed_line_single_meson(generator, invariants):
    """V(lambda, A, A) N/eps^2 and the rotated form (1/f) <lambda (B_0 P_2' - m_0^2/N_f Theta_2)>."""
    leg = alp_pole_leg(invariants, EPSILON, normalized=False)
    value = sp.expand(contract_potential_vertex([(generator, 0), leg, leg], invariants, EPSILON) / EPSILON**2)
    Y = alp_pole_overlap(invariants)
    matrix = (invariants['B_0'] * rotated_pseudoscalar_source(invariants, Y)
              - invariants['anomaly_mass_squared'] / FLAVOURS * invariants['theta'][2] * sp.eye(3))
    return value, sp.expand((generator * matrix).trace() / invariants['f_pi'])


def closed_line_scalar_source_piece(generator, invariants, Y):
    """S_1 piece of V(lambda, A, A) N/eps^2, from P_2' = P_2 - 2{Y, S_1}: -(2 B_0/f) <lambda {Y, S_1}>."""
    return sp.expand(-2 * invariants['B_0'] / invariants['f_pi']
                     * (generator * anticommutator(Y, invariants['S'][1])).trace())


def closed_line_scalar_source_piece_down_strange(generator, invariants, Y):
    """Same piece for diagonal lambda and S_1 with only d-s, s-d entries, using <lambda{Y,S_1}> = <{lambda,Y}S_1>:
      -(2 B_0/f) (lambda_dd + lambda_ss) (Y_ds (S_1)_sd + Y_sd (S_1)_ds).
    Only Y_ds, Y_sd enter."""
    S_1 = invariants['S'][1]
    return sp.expand(-2 * invariants['B_0'] / invariants['f_pi'] * (generator[1, 1] + generator[2, 2])
                     * (Y[1, 2] * S_1[2, 1] + Y[2, 1] * S_1[1, 2]))


def closed_line_scalar_source_piece_down_strange_pair(generator, invariants, Y):
    """Same piece for lambda with only d-s, s-d entries (lambda_6, lambda_7), S_1 with only d-s, s-d
    entries and Y diagonal plus d-s (commutes with the quark charge):
      -(2 B_0/f) (Y_dd + Y_ss) (lambda_ds (S_1)_sd + lambda_sd (S_1)_ds).
    Only the diagonal entries Y_dd, Y_ss enter."""
    S_1 = invariants['S'][1]
    return sp.expand(-2 * invariants['B_0'] / invariants['f_pi'] * (Y[1, 1] + Y[2, 2])
                     * (generator[1, 2] * S_1[2, 1] + generator[2, 1] * S_1[1, 2]))


def without_down_strange_pseudoscalar(invariants):
    """Copy of the invariants with the d-s, s-d entries of P_1 removed (S_1 unchanged)."""
    result = dict(invariants)
    P = list(invariants['P'])
    P[1] = sp.diag(*[P[1][q, q] for q in range(3)])
    result['P'] = P
    return result


def down_strange_overlap_entry(invariants, alp_pole_squared_mass=ALP_POLE_SQUARED_MASS):
    """Y_ds = C_ds / (B_0 (m_d + m_s) - m_A^2): M_s acts entrywise on off-diagonal entries for diagonal chi_0."""
    m = invariants['S'][0]
    return matrix_C(invariants)[1, 2] / (invariants['B_0'] * (m[1, 1] + m[2, 2]) - alp_pole_squared_mass)


# ---------------------------------------------------------------------------
# Closed unmixed meson line on the ALP-two-meson vertex
# ---------------------------------------------------------------------------


def unmixed_mass_matrix(invariants):
    """eps = 0 meson mass matrix in the basis lambda_0..lambda_8: M^2_ij = <lambda_i M_0(lambda_j)>/2,
    with M_s(X) = {chi_0, X}/2 + (m_0^2/N_f) <X> 1 - s X of the tree record (B_0 m_q plus the
    m_0^2 singlet term; the generators are normalized to <lambda_i lambda_j> = 2 delta_ij)."""
    return sp.Matrix(9, 9, lambda i, j: sp.expand(
        (GELL_MANN_GENERATORS[i] * pole_operator(GELL_MANN_GENERATORS[j], invariants, 0)).trace() / 2))


def unmixed_eigen_generators(invariants):
    """[(m^2, generator)]: the three neutral diagonal eigen-generators (unmixed_neutral_matrix) and
    lambda_1, lambda_2, lambda_4, lambda_5, lambda_6, lambda_7 with m^2 = (chi_i + chi_j)/2."""
    result = [(mass, unmixed_neutral_matrix(invariants, mass)) for mass in unmixed_neutral_squared_masses(invariants)]
    chi = [2 * invariants['B_0'] * invariants['S'][0][q, q] for q in range(3)]
    for index, (row, column) in ((1, (0, 1)), (2, (0, 1)), (4, (0, 2)), (5, (0, 2)), (6, (1, 2)), (7, (1, 2))):
        result.append(((chi[row] + chi[column]) / 2, GELL_MANN_GENERATORS[index]))
    return result


def closed_meson_line_flavour_matrix(propagator):
    """sum_ij propagator_ij {lambda_i, lambda_j} (Gell-Mann basis): flavour matrix that the closed
    meson line contracts with S_1 in the ALP-two-meson vertex."""
    return sp.expand(sum((propagator[i, j] * anticommutator(GELL_MANN_GENERATORS[i], GELL_MANN_GENERATORS[j])
                          for i in range(9) for j in range(9) if propagator[i, j] != 0), sp.zeros(3)))


CLOSED_LINE_SQUARED_MOMENTUM = sp.Rational(1, 7)    # k^2 (GeV^2) for the closed meson line, test input


def closed_unmixed_meson_line(invariants, squared_momentum=CLOSED_LINE_SQUARED_MOMENTUM):
    """Closed unmixed meson line on the ALP-two-meson vertex at a fixed k^2 (no loop integral).

    Delta = (k^2 - M_unmixed^2)^{-1} (9x9, eps^0); returns (sum_ij Delta_ij {lambda_i, lambda_j},
    <S_1 sum_ij Delta_ij {lambda_i, lambda_j}>, [<S_1 lambda_n^2> for the unmixed eigen-generators]).
    """
    propagator = (squared_momentum * sp.eye(9) - unmixed_mass_matrix(invariants)).inv().applyfunc(sp.radsimp)
    flavour = closed_meson_line_flavour_matrix(propagator).applyfunc(sp.radsimp)
    S_1 = invariants['S'][1]
    trace = sp.radsimp(sp.expand((S_1 * flavour).trace()))
    eigen_traces = [sp.radsimp(sp.expand((S_1 * generator * generator).trace()))
                    for _, generator in unmixed_eigen_generators(invariants)]
    return flavour, trace, eigen_traces


# ---------------------------------------------------------------------------
# (g) ALP component of the meson-like poles
# ---------------------------------------------------------------------------


def meson_pole_alp_component_series(invariants, unmixed_squared_mass, unmixed_matrix,
                                    alp_pole_squared_mass=ALP_POLE_SQUARED_MASS, order=4):
    """Series in eps (fixed m_A^2, m_{a,0}^2 adjusted) of v_a^{(n)} for the pole continuous with
    (unmixed_matrix, unmixed_squared_mass); sign by continuity <Phi_n Phi_{n,0}> > 0.

    v_a = sign (delta/eps)/sqrt(delta^2/eps^2 + 2<Z^2>), Z = delta M^{-1}_{m_n^2}(C),
    delta = m_n^2 - m_{n,0}^2 from daughter_pole_series (tree script).
    Returns (coefficients eps^0..eps^{order-1}, expected eps^1 coefficient
    <Phi_{n,0} C>/(m_A^2 - m_{n,0}^2)).
    """
    t, shift = sp.Symbol('t'), sp.Symbol('delta')
    pole = daughter_pole_series(invariants, alp_pole_squared_mass, unmixed_squared_mass, order=3)
    delta = sp.expand((pole - unmixed_squared_mass).subs(t, EPSILON**2))
    residue_form = inverse_pole_operator(matrix_C(invariants), invariants, unmixed_squared_mass + shift).applyfunc(
        lambda entry: sp.cancel(shift * entry))
    Z = residue_form.applyfunc(lambda entry: sp.expand(sp.series(entry.subs(shift, delta), EPSILON, 0, order + 2)
                                                       .removeO()))
    sign = 1 if sp.re(sp.N((residue_form.subs(shift, 0) * unmixed_matrix).trace(), 30)) > 0 else -1
    component = sign * sp.expand(delta / EPSILON) / sp.sqrt(sp.expand(delta**2 / EPSILON**2) + 2 * (Z * Z).trace())
    series = sp.expand(sp.series(component, EPSILON, 0, order).removeO())
    expected = (unmixed_matrix * matrix_C(invariants)).trace() / (alp_pole_squared_mass - unmixed_squared_mass)
    return ([sp.radsimp(sp.expand(series.coeff(EPSILON, k))) for k in range(order)],
            sp.radsimp(sp.expand(expected)))


def kaon_like_unmixed(invariants):
    """Unmixed d-s pole: m^2 = (chi_d + chi_s)/2, matrix along the residue of M^{-1}_s(C), <X^2> = 2."""
    chi = [2 * invariants['B_0'] * invariants['S'][0][i, i] for i in range(3)]
    squared_mass = (chi[1] + chi[2]) / 2
    shift = sp.Symbol('delta')
    residue = inverse_pole_operator(matrix_C(invariants), invariants, squared_mass + shift).applyfunc(
        lambda entry: sp.limit(shift * entry, shift, 0))
    matrix = residue * sp.sqrt(2 / sp.expand((residue * residue).trace()))
    return squared_mass, matrix.applyfunc(lambda entry: sp.radsimp(sp.expand(entry)))


# ---------------------------------------------------------------------------

FORMULA_TEXT = r"""
Invariant potential vertex (dressed frame, n_a ALP legs, n_P meson legs lambda_1..lambda_nP):
  V = eps^{n_a} f^{2-n_P-n_a} { B_0 (-1)^{floor(n_P/2)} (1/n_P!) sum_{sigma in S_nP} <X_{n_a} lambda_s1...lambda_snP>
                                + topological }        X = S_{n_a} (n_P even), P_{n_a} (n_P odd)
  topological: n_P=0: -(m_0^2/N_f) sum_{k=1}^{n_a-1} binom(n_a,k) Theta_k Theta_{n_a-k}
               n_P=1: -(m_0^2/N_f) Theta_{n_a} <lambda_1>
               n_P=2, n_a=0: -(m_0^2/(2N_f)) <lambda_1><lambda_2>  (eps^0 f^0)
  plus -m_{a,0}^2 at (n_a, n_P) = (2, 0).
  n_P = 1: V = eps^{n_a} f^{1-n_a} <lambda (B_0 P_{n_a} - (m_0^2/N_f) Theta_{n_a} 1)>, = eps <lambda C> for n_a = 1.
One end on A: V(..., A) = (eps/sqrt(N)) [V(..., a)/eps + V(..., 2Y)],  Y = M^{-1}_{m_A^2}(C).
Closed A line, two meson legs a, b (potential part, times N/eps^2):
  -(B_0/(2f^2)) <S_2'{a,b}> + (B_0/(3f^2)) <(P_1 - {Y,m_q})([a,[b,Y]] + [b,[a,Y]])>
  + (B_0/(6f^2)) <m_q({[Y,[Y,a]],b} + {[Y,[Y,b]],a})>,
  S_2' = S_2 + 2{Y,P_1} - Y^2 m_q - m_q Y^2 - 2 Y m_q Y.
Kinetic part (legs a(p), b(-p), A(k), A(-k)): (1/(3f^2)) (p^2 + k^2) <[Y,a][Y,b]>, no p.k term.
"""


def main():
    print(FORMULA_TEXT)
    compared, failures, unmatched = compare_with_builder()
    print(f'(a) builder (diagonal invariants, 10 fields, degree <= 4): {compared} monomials, '
          f'{len(failures)} nonzero residuals, {len(unmatched)} unmatched')
    print('(a) n_P = 1 combination vs builder:', set(single_meson_check_with_builder().values()))
    compared, failures = compare_with_joint_vertex()
    print(f'(a) d-s invariants vs vertex(): {compared} vertices, {len(failures)} nonzero residuals')
    for label, invariants, fields, order in (('d-s, 10 fields, order 4', test_invariants(True), MESON_FIELD_NAMES, 4),
                                             ('d-s, 7 fields, order 6', test_invariants(True, True),
                                              SIXTH_ORDER_FIELDS, 6)):
        compared, non_homogeneous, failures = compare_with_direct_expansion(invariants, fields, order)
        print(f'(a,b,d) direct expansion {label}: {compared} monomials, {len(non_homogeneous)} not homogeneous, '
              f'{len(failures)} nonzero residuals')
    kinetic_alp, topological_extra = kinetic_and_topological_structure()
    print(f'(c) ALP in kinetic monomials: {len(kinetic_alp)}; topological n_P >= 2 extra: {len(topological_extra)}')
    cubic = cubic_vertices(test_invariants(True), EPSILON)
    print('(e) cubic vertices at the d-s test inputs:', cubic)
    invariants = test_invariants(True)
    print('S_1 = -i[c_q^{VS}, m_q] (d-s test input):', invariants['S'][1].tolist())
    print('(e) nonzero V_{1; lambda_i, lambda_j}, i <= j in 0..8:', nonzero_alp_two_meson_pairs(invariants))
    residuals = joint_vertex_cubic_residuals(invariants)
    print(f'    against vertex(): {len(residuals)} pairs, nonzero residuals: {[k for k, v in residuals.items() if v]}')
    flavour, trace, eigen_traces = closed_unmixed_meson_line(invariants)
    print(f'(e) closed unmixed meson line at k^2 = {CLOSED_LINE_SQUARED_MOMENTUM}: sum Delta_ij {{lambda_i, lambda_j}} =',
          flavour.tolist(), '; <S_1 ...> =', trace, '; <S_1 lambda_n^2> =', eigen_traces)
    for label, case in (('P_1, S_1 with d-s', invariants),
                        ('P_1 diagonal, S_1 with d-s', without_down_strange_pseudoscalar(invariants))):
        overlap = alp_pole_overlap(case)
        pieces = {index: closed_line_scalar_source_piece(GELL_MANN_GENERATORS[index], case, overlap)
                  for index in (3, 8, 0, 6, 7)}
        print(f'(f) S_1 piece -(2B_0/f)<lambda{{Y,S_1}}> of V(lambda, A, A) N/eps^2, {label}: {pieces}; '
              f'Y_ds = {overlap[1, 2]}, Y_ds - C_ds/(B_0(m_d+m_s) - m_A^2) = '
              f'{sp.simplify(overlap[1, 2] - down_strange_overlap_entry(case))}')
    Y = alp_pole_overlap(invariants)
    generators = flavour_generators()
    a, b = generators[1], generators[7]
    potential = closed_alp_line_potential(a, b, invariants)
    print('(f) closed line (E_dd, lambda_6), d-s inputs: potential =', potential,
          '; rotated-source residual =', sp.simplify(potential - rotated_source_prediction(a, b, invariants, Y)),
          '; matrix-form residual =', sp.simplify(potential - closed_alp_line_potential_matrix_form(a, b, invariants, Y)))
    full = closed_alp_line_full(a, b, invariants)
    kinetic = sp.expand(full - full.subs({MOMENTUM_P: 0, MOMENTUM_K: 0}))
    print('    kinetic part:', kinetic, '; residual against (p^2+k^2)<[Y,a][Y,b]>/(3f^2):',
          sp.simplify(kinetic - kinetic_closed_line_form(a, b, invariants, Y)))
    print('    on-shell form residual:',
          sp.simplify(tadpole_kinematics(full) - on_shell_closed_line_form(a, b, invariants, Y)))
    value, rotated = closed_line_single_meson(generators[7], invariants)
    print('(f) V(lambda_6, A, A) N/eps^2 - (1/f)<lambda (B_0 P_2\' - ...)> =', sp.simplify(value - rotated))
    for index, mass in enumerate(unmixed_neutral_squared_masses(invariants)):
        matrix = unmixed_neutral_matrix(invariants, mass)
        coefficients, expected = meson_pole_alp_component_series(invariants, mass, matrix)
        print(f'(g) neutral pole {index}: series {coefficients}, expected eps^1 {expected}')
    mass, matrix = kaon_like_unmixed(invariants)
    coefficients, expected = meson_pole_alp_component_series(invariants, mass, matrix)
    print(f'(g) K0-like pole: series {coefficients}, expected eps^1 {expected}')


if __name__ == '__main__':
    main()
