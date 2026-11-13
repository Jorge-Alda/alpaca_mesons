#!/usr/bin/env python3
"""Generic finite-mixing LO tree amplitude for a_phys -> P_i P_j P_k.

Run:    .venv/bin/python scripts/derive_alp_finite_mixing_three_pseudoscalars.py
Checks: .venv/bin/python scripts/derive_alp_finite_mixing_three_pseudoscalars.py --test
Tests:  .venv/bin/python -m unittest discover -s scripts -p test_derive_alp_finite_mixing_three_pseudoscalars.py

Conventions follow Alda:2026nrd: natural units, mostly-minus metric, S=1+iT,
Phi = lambda_I P_I with <lambda_I lambda_J> = 2 delta_IJ, U = exp(i Phi/f_pi),
chi_0 = 2 B_0 m_q. Squared masses, Mandelstam invariants and entries of M^2
share one unit (for example GeV^2). Vertices omit the Feynman-rule factor i.
<.> denotes the flavour trace. Pole operator at a squared mass s:
M_s(X) = B_0{X,m_q} + (m_0^2/N_f) <X> 1 - s X, with inverse M_s^{-1}(C), used at
s = m_a^2 (parent) and s = m_r^2 (daughter r); C = -B_0{m_q,c_q^AP} - (m_0^2/N_f) c_GA 1.

Coordinates are the derivative-source-free CCWZ coordinates of
memory/alp/three-pseudoscalar/alp-finite-mixing-three-pion-amplitudes.md: K=1,
the cubic vertex V3 has no momentum dependence and the quartic vertex is
V4 = V4_kin + V4_chi. The explicit V3, V4 and the matrix-form inverse of
M^2 - sK are imported from scripts/derive_alp_finite_three_pion.py.

A field direction is a pair (X, x): meson matrix X = lambda_I v_I and ALP
component x. For a(p) -> P_i(p_1) P_j(p_2) P_k(p_3), s_rs = (p_r+p_s)^2 and all
vertex momenta are incoming, (k_0,k_1,k_2,k_3) = (p,-p_1,-p_2,-p_3). The exact
LO tree amplitude is

  M = V4[a,1,2,3] + sum_{(a,r)|(u,v)} V3[a,r,I] (M^2 - s_uv K)^{-1}_{IJ} V3[J,u,v],

with the three pairings (a,1)|(2,3) at s_23, (a,2)|(1,3) at s_13 and
(a,3)|(1,2) at s_12. Pole columns are not expanded in f_pi/f_a. Identical
daughters receive no 1/n! in the amplitude; 1/S belongs to the phase space.

Route 2 (pole-aligned coordinates) applies the transformation of
memory/alp/shared/alp-lo-pole-aligned-tree-vertices.md with
w_i lambda_i = (2 f_pi/f_a) Y, Y = M_{m_a^2}^{-1}(C) at m_a^2, i.e.
U_df = exp(iaY/f_a) U' exp(iaY/f_a). The sources follow from the laws of
alp-chiral-source-transformations.md: C_V' = 0, C_A' = -Y, M' = exp(-ixY)
widehat M exp(-ixY) (Taylor law of the general-source record), Theta -> Theta + x <Y>,
so d_a^n chi_+|_0 = 4 B_0 S_n'/f_a^n and d_a^n chi_-|_0 = 4 i B_0 P_n'/f_a^n. The script checks these against the
joint CCWZ expansion (derive_alp_lo_joint_vertices.py), the congruence and
pole condition of the single-pole record, and route 1 at LO poles.
"""
import itertools
import sys
import unittest

import sympy as sp

from derive_alp_finite_three_pion import (
    B_0, anticommutator, ccwz_vacuum_derivatives, commutator, cubic_pair_source, cubic_vertex, f_a, f_pi,
    inverse_meson_mass_operator, inverse_mixed_mass_operator,
    linear_source_ccwz_derivatives, quartic_kinetic_channel,
    quartic_kinetic_vertex, quartic_mass_vertex)
from derive_alp_lo_joint_vertices import quadratic_blocks as joint_quadratic_blocks
from derive_alp_lo_joint_vertices import vertex as joint_vertex
from symbolic_utils import matrix_unit

# Pairings: (daughter attached to the parent vertex, the two daughters at the other vertex).
PAIRINGS = ((1, (2, 3)), (2, (1, 3)), (3, (1, 2)))


# ---------------------------------------------------------------------------
# Kinematics
# ---------------------------------------------------------------------------
def decay_scalar_products(parent_mass_sq, daughter_masses_sq, pair_invariants):
    """All-incoming scalar products k_r.k_s for legs (parent, 1, 2, 3).

    daughter_masses_sq maps 1,2,3 to squared masses; pair_invariants maps
    frozenset({u, v}) to s_uv = (p_u+p_v)^2. Uses p.p_r = (m_a^2+m_r^2-s_uv)/2
    and p_r.p_u = (s_ru-m_r^2-m_u^2)/2, i.e. momentum conservation.
    """
    products = {}
    for daughter in (1, 2, 3):
        others = frozenset({1, 2, 3} - {daughter})
        products[(0, daughter)] = -(parent_mass_sq + daughter_masses_sq[daughter]
                                    - pair_invariants[others])/2
    for first, second in itertools.combinations((1, 2, 3), 2):
        products[(first, second)] = (pair_invariants[frozenset({first, second})]
                                     - daughter_masses_sq[first] - daughter_masses_sq[second])/2
    return products


# ---------------------------------------------------------------------------
# Generic amplitude (component contraction over an explicit field basis)
# ---------------------------------------------------------------------------
def contact_vertex(directions, scalar_products, plus, minus):
    """V4 = V4_kin + V4_chi on four directions; plus/minus are chi_+-' vacuum derivatives."""
    kinetic = quartic_kinetic_vertex(directions, scalar_products)
    mass = quartic_mass_vertex(directions, plus[0], minus[1], plus[2], minus[3], plus[4])
    return sp.expand((kinetic + mass).doit())


def component_quadratic_blocks(generators, quark_masses, anomaly_mass_sq,
                               invariant_coupling, alp_mass_entry):
    """K and M^2 in the field basis (generators..., ALP) of the three-pion record, section 2.

    For generators normalized as <lambda_I lambda_J> = 2 delta_IJ, K = 1 as in
    the record; for another basis, K_IJ = <E_I E_J>/2 and the same bilinear
    forms give M^2.
    """
    size = len(generators)
    flavours = quark_masses.rows
    chi_0 = 2*B_0*quark_masses
    mass_squared = sp.zeros(size + 1)
    for row, first in enumerate(generators):
        for column, second in enumerate(generators):
            mass_squared[row, column] = (sp.trace(anticommutator(first, second)*chi_0)/4
                                         + anomaly_mass_sq*sp.trace(first)*sp.trace(second)/(2*flavours))
        mass_squared[row, size] = -f_pi*sp.trace(invariant_coupling*first)/f_a
        mass_squared[size, row] = mass_squared[row, size]
    mass_squared[size, size] = alp_mass_entry
    kinetic = sp.eye(size + 1)
    kinetic[:size, :size] = sp.Matrix(size, size, lambda row, column:
                                      sp.trace(generators[row]*generators[column])/2)
    return kinetic, sp.expand(mass_squared)


def hermitian_unit_basis(dimension=3):
    """Complete Hermitian basis without square roots: E_aa, E_ab+E_ba, i(E_ab-E_ba).

    It is not normalized as <lambda lambda>=2 delta; the kinetic block
    K_IJ = <E_I E_J>/2 of component_quadratic_blocks accounts for that.
    The contraction is basis independent (constant congruence of K, M^2, V).
    """
    unit = lambda row, column: sp.Matrix(matrix_unit(row, column, dimension))
    basis = [unit(a, a) for a in range(dimension)]
    for a, b in itertools.combinations(range(dimension), 2):
        basis += [unit(a, b) + unit(b, a), sp.I*(unit(a, b) - unit(b, a))]
    return basis


def field_basis_directions(generators):
    """Unit directions (lambda_I, 0) for every generator and (0, 1) for the ALP."""
    zero = sp.zeros(generators[0].rows)
    return [(sp.Matrix(generator), 0) for generator in generators] + [(zero, 1)]


def generic_amplitude(parent, daughters, parent_mass_sq, daughter_masses_sq,
                      pair_invariants, plus, minus, field_basis, kinetic, mass_squared):
    """Exact generic a_phys -> P_i P_j P_k amplitude; returns (total, contact, exchanges).

    parent and daughters[0..2] are pole directions (X, x) of the parent and of
    daughters 1, 2, 3. The exchange of pairing (a,r)|(u,v) is
    sum_IJ V3[a,r,e_I] (M^2 - s_uv K)^{-1}_IJ V3[e_J,u,v], summed over the
    complete field basis e_I, with the full (not pole-summed) inverse.
    """
    legs = dict(zip((1, 2, 3), daughters))
    products = decay_scalar_products(parent_mass_sq, daughter_masses_sq, pair_invariants)
    contact = contact_vertex([parent] + list(daughters), products, plus, minus)
    exchanges = {}
    for attached, (first, second) in PAIRINGS:
        invariant = pair_invariants[frozenset({first, second})]
        left = sp.Matrix([cubic_vertex(parent, legs[attached], direction,
                                       plus[1], minus[2], plus[3]).doit()
                          for direction in field_basis])
        right = sp.Matrix([cubic_vertex(direction, legs[first], legs[second],
                                        plus[1], minus[2], plus[3]).doit()
                           for direction in field_basis])
        solved = (mass_squared - invariant*kinetic).LUsolve(right)
        exchanges[(attached, first, second)] = sp.together((left.T*solved)[0])
    return contact + sum(exchanges.values()), contact, exchanges


def permuted_inputs(daughters, daughter_masses_sq, pair_invariants, order):
    """Relabel daughters with their momenta: new daughter n is old daughter order[n-1]."""
    new_daughters = [daughters[old - 1] for old in order]
    new_masses = {new: daughter_masses_sq[old] for new, old in zip((1, 2, 3), order)}
    new_invariants = {frozenset({first, second}):
                      pair_invariants[frozenset({order[first - 1], order[second - 1]})]
                      for first, second in itertools.combinations((1, 2, 3), 2)}
    return new_daughters, new_masses, new_invariants


# ---------------------------------------------------------------------------
# Record-side three-pion amplitudes (memory record section 4, matrix form)
# ---------------------------------------------------------------------------
def three_pion_record_amplitude(parent, daughters, mandelstam_s, mandelstam_t, mandelstam_u,
                                plus, minus, quark_masses, anomaly_mass_sq,
                                invariant_coupling, alp_mass_entry, numbers=None):
    """Boxed M_000 or M_0+- of the three-pion record with the record's s,t,u.

    daughters = (pi0, pi0, pi0) or (pi0, pi+, pi-). Contact uses the channel
    form of V4_kin; exchanges use the matrix-form inverse of M^2 - qK applied
    to the internal-leg source of the second cubic vertex. numbers, if
    given, substitutes values of f_pi, f_a, B_0 before simplification.
    """
    numbers = numbers or {}
    poles = [parent] + list(daughters)
    contact = (quartic_kinetic_channel(poles, mandelstam_s, mandelstam_t, mandelstam_u)
               + quartic_mass_vertex(poles, plus[0], minus[1], plus[2], minus[3], plus[4]).doit())
    total = sp.expand(contact.subs(numbers))
    # Record pairings: s with (a,0)|(+,-), t with (a,+)|(0,-), u with (a,-)|(0,+).
    for invariant, attached, (first, second) in ((mandelstam_s, 0, (1, 2)),
                                                 (mandelstam_t, 1, (0, 2)),
                                                 (mandelstam_u, 2, (0, 1))):
        source = cubic_pair_source(daughters[first], daughters[second],
                                   plus[1], minus[2], plus[3])
        source = (sp.Matrix(source[0].doit()), sp.sympify(source[1]).doit())

        def meson_inverse(matrix, invariant=invariant):
            return inverse_meson_mass_operator(matrix, quark_masses, anomaly_mass_sq, invariant)
        field = inverse_mixed_mass_operator(source, invariant_coupling, alp_mass_entry,
                                            meson_inverse, invariant)
        field = (sp.Matrix(field[0]).doit(), sp.sympify(field[1]).doit())
        total += sp.cancel(cubic_vertex(parent, daughters[attached], field,
                                        plus[1], minus[2], plus[3]).doit().subs(numbers))
    return total


def three_pion_substitution_residuals():
    """Generic formula minus the record's M_000 and M_0+-, exact in s,t,u.

    Returns channel -> (residual, {pairing: exchange is nonzero}). For the
    linear source convention the (a,0)|(+,-) exchange of pi0 pi+ pi- vanishes
    identically: V3[J,pi+,pi-] only involves the diagonal of
    d_a chi_+' ~ [c_q^{VS}, m_q], which is zero for diagonal m_q.

    Inputs are exact rational sources; the neutral columns are the record's
    P(s) = (2 f_pi M^{-1}(C)/f_a, 1)/sqrt(N) at fixed squared masses, without
    the factor 1/sqrt(N): both expressions are linear in every external
    column, so each normalization is a common factor. The identity is
    algebraic in the directions, so on-shell roots are not needed.
    """
    numbers = {f_pi: 1, f_a: 4, B_0: 2}
    quark_masses = sp.diag(1, 2, 7)
    axial = sp.Matrix([[1, 0, 1], [0, -1, 2], [1, 2, 1]])
    pseudoscalar = sp.Matrix([[2, 1, 0], [1, 1, -1], [0, -1, 3]])
    vector = sp.Matrix([[0, 1, 0], [1, 2, 1], [0, 1, -1]])
    scalar = sp.Matrix([[1, 0, 2], [0, 0, 1], [2, 1, 0]])
    anomaly_mass_sq, gluon_coupling, alp_mass_entry = sp.Rational(9, 2), sp.Rational(1, 3), sp.Rational(23, 2)
    plus, minus = linear_source_ccwz_derivatives(axial, pseudoscalar, vector, scalar, quark_masses)
    plus = [sp.Matrix(matrix).subs(numbers) for matrix in plus]
    minus = [sp.Matrix(matrix).subs(numbers) for matrix in minus]
    axial_pseudoscalar = axial - pseudoscalar
    gluon_axial = gluon_coupling + sp.trace(axial)
    invariant_coupling = (-B_0*anticommutator(quark_masses, axial_pseudoscalar)
                          - anomaly_mass_sq*gluon_axial*sp.eye(3)/3).subs(numbers)

    def pole_column(squared_mass):
        inverse = inverse_meson_mass_operator(invariant_coupling, quark_masses,
                                              anomaly_mass_sq, squared_mass).subs(numbers)
        return (2*numbers[f_pi]*inverse/numbers[f_a], 1)

    parent_mass_sq, neutral_mass_sq, charged_mass_sq = sp.Rational(31, 2), sp.Rational(5, 2), sp.Rational(3, 1)
    parent = pole_column(parent_mass_sq)
    neutral = pole_column(neutral_mass_sq)
    charged_plus = (sp.sqrt(2)*sp.Matrix(matrix_unit(0, 1)), 0)
    charged_minus = (sp.sqrt(2)*sp.Matrix(matrix_unit(1, 0)), 0)
    generators = hermitian_unit_basis()
    kinetic, mass_squared = component_quadratic_blocks(generators, quark_masses, anomaly_mass_sq,
                                                       invariant_coupling, alp_mass_entry)
    mass_squared = mass_squared.subs(numbers)
    basis = field_basis_directions(generators)
    mandelstam_s, mandelstam_t, mandelstam_u = sp.symbols('s t u', real=True)
    residuals = {}
    for channel, daughters, masses in (
            ('000', (neutral, neutral, neutral), (neutral_mass_sq,)*3),
            ('0+-', (neutral, charged_plus, charged_minus),
             (neutral_mass_sq, charged_mass_sq, charged_mass_sq))):
        # Generic labels: daughter 1,2,3 = record legs; s=s_23, t=s_13, u=s_12.
        invariants = {frozenset({2, 3}): mandelstam_s, frozenset({1, 3}): mandelstam_t,
                      frozenset({1, 2}): mandelstam_u}
        generic, _, exchanges = generic_amplitude(parent, daughters, parent_mass_sq,
                                          dict(zip((1, 2, 3), masses)), invariants,
                                          plus, minus, basis, kinetic, mass_squared)
        record = three_pion_record_amplitude(parent, daughters, mandelstam_s, mandelstam_t,
                                             mandelstam_u, plus, minus, quark_masses,
                                             anomaly_mass_sq, invariant_coupling, alp_mass_entry,
                                             numbers)
        residuals[channel] = (sp.cancel(sp.together((generic - record).subs(numbers))),
                              {pairing: sp.cancel(value.subs(numbers)) != 0
                               for pairing, value in exchanges.items()})
    return residuals


# ---------------------------------------------------------------------------
# Kinetic contact: reduction of the S_4 sum to pair invariants
# ---------------------------------------------------------------------------
def kinetic_contact_reduction():
    """Reduce V4_kin = -(1/48 f_pi^2) sum_S4 (k_s2.k_s4) <[X_s1,X_s2][X_s3,X_s4]>.

    Trace symbols T_{a1|23} = <[X_a,X_1][X_2,X_3]> etc. are canonicalized
    with antisymmetry in each commutator and <AB>=<BA>. Returns the S_4
    sum after substituting the decay scalar products, and the reduced form.
    """
    names = ('a', '1', '2', '3')
    structures = {frozenset({0, 1}): sp.Symbol(r'T_{a1|23}'),
                  frozenset({0, 2}): sp.Symbol(r'T_{a2|13}'),
                  frozenset({0, 3}): sp.Symbol(r'T_{a3|12}')}

    def canonical(first, second, third, fourth):
        sign = 1
        if 0 in (third, fourth):
            first, second, third, fourth = third, fourth, first, second
        if second == 0:
            first, second, sign = second, first, -sign
        if third > fourth:
            sign = -sign
        return sign*structures[frozenset({first, second})]

    parent_mass_sq = sp.Symbol('m_a^2')
    masses = {r: sp.Symbol(f'm_{r}^2') for r in (1, 2, 3)}
    invariants = {frozenset({u, v}): sp.Symbol(f's_{u}{v}') for u, v in itertools.combinations((1, 2, 3), 2)}
    products = decay_scalar_products(parent_mass_sq, masses, invariants)
    summed = -sum(products[tuple(sorted((order[1], order[3])))]*canonical(*order)
                  for order in itertools.permutations(range(4)))/(48*f_pi**2)
    s12, s13, s23 = (invariants[frozenset(pair)] for pair in ((1, 2), (1, 3), (2, 3)))
    reduced = -(structures[frozenset({0, 1})]*(s13 - s12)
                + structures[frozenset({0, 2})]*(s23 - s12)
                + structures[frozenset({0, 3})]*(s23 - s13))/(24*f_pi**2)
    return sp.expand(summed), reduced, names


# ---------------------------------------------------------------------------
# Index-level contraction, symmetry factors and weak-mixing power counting
# ---------------------------------------------------------------------------
def index_level_amplitude():
    """Generic contraction with two meson fields and one ALP, and its O(epsilon) term.

    Vertex tensors are derivatives of a generic Lagrangian polynomial (no
    factorials removed). Power counting, epsilon = f_pi/f_a with the source
    Taylor derivatives in x = a/f_a held fixed: a vertex with n ALP legs is
    O(epsilon^n) (in units of f_pi), the pure-meson cubic vertex is zero, the
    meson-ALP mass mixing is O(epsilon), parent meson components and daughter
    ALP components are O(epsilon), and the other components are 1 + O(epsilon^2).
    """
    epsilon = sp.Symbol('epsilon', positive=True)
    fields = sp.symbols('P_1 P_2 a')
    alp = 2

    def lagrangian(order, name):
        return sum(sp.Symbol(f'{name}_' + ''.join(map(str, indices)))
                   * epsilon**indices.count(alp)*sp.prod(fields[i] for i in indices)
                   for indices in itertools.combinations_with_replacement(range(3), order)
                   if not (order == 3 and alp not in indices))
    cubic = lagrangian(3, 'g')
    quartic = lagrangian(4, 'h')
    cubic_tensor = {indices: sp.diff(cubic, *[fields[i] for i in indices])
                    for indices in itertools.product(range(3), repeat=3)}
    quartic_tensor = {indices: sp.diff(quartic, *[fields[i] for i in indices])
                      for indices in itertools.product(range(3), repeat=4)}
    parent = sp.Matrix([epsilon*sp.Symbol('u_1'), epsilon*sp.Symbol('u_2'),
                        1 + epsilon**2*sp.Symbol('r_a')])
    daughters = [sp.Matrix([sp.Symbol(f'w_{d}1') + epsilon**2*sp.Symbol(f'd_{d}1'),
                            sp.Symbol(f'w_{d}2') + epsilon**2*sp.Symbol(f'd_{d}2'),
                            epsilon*sp.Symbol(f'xi_{d}')]) for d in (1, 2, 3)]
    mixing = sp.Matrix(sp.symbols('kappa_1 kappa_2'))

    def contact(vectors):
        return sp.expand(sum(quartic_tensor[idx]*sp.prod(vectors[n][idx[n]] for n in range(4))
                             for idx in itertools.product(range(3), repeat=4)))

    def cubic_component(vectors):
        return sp.Matrix([sum(cubic_tensor[(i, j, k)]*vectors[0][i]*vectors[1][j]
                              for i in range(3) for j in range(3)) for k in range(3)])

    invariants = {pair: sp.Symbol('s_' + ''.join(map(str, sorted(pair))))
                  for pair in map(frozenset, itertools.combinations((1, 2, 3), 2))}
    leading_propagator = {}
    exchange_terms = []
    for attached, (first, second) in PAIRINGS:
        invariant = invariants[frozenset({first, second})]
        # (M^2 - s)^{-1} with M^2 - s = A + epsilon B, A block diagonal (meson block and
        # ALP entry), B the meson-ALP mass mixing. The inverse of A at s_uv is kept as
        # generic symmetric symbols; Neumann series to O(epsilon^2).
        label = ''.join(map(str, sorted((first, second))))
        meson_inverse = sp.Matrix(2, 2, lambda i, j: sp.Symbol(f'G^{label}_{min(i, j)+1}{max(i, j)+1}'))
        inverse_block = sp.diag(meson_inverse, sp.Symbol(f'G^{label}_aa'))
        off_diagonal = sp.zeros(3)
        off_diagonal[0:2, 2] = mixing
        off_diagonal[2, 0:2] = mixing.T
        propagator = (inverse_block - epsilon*inverse_block*off_diagonal*inverse_block
                      + epsilon**2*inverse_block*off_diagonal*inverse_block*off_diagonal*inverse_block)
        leading_propagator[invariant] = inverse_block
        left = cubic_component([parent, daughters[attached - 1]])
        right = cubic_component([daughters[first - 1], daughters[second - 1]])
        exchange_terms.append(sp.expand((left.T*propagator*right)[0]))
    contact_term = contact([parent] + daughters)
    exchange = sum(exchange_terms)
    amplitude = contact_term + exchange
    series = {order: sp.expand(amplitude).coeff(epsilon, order) for order in range(3)}
    exchange_series = {order: sp.expand(exchange).coeff(epsilon, order) for order in range(3)}
    zeroth = [daughter.subs(epsilon, 0) for daughter in daughters]
    expected_leading = sp.expand(contact([sp.Matrix([0, 0, 1])] + zeroth).coeff(epsilon, 1)
                                 + contact([sp.Matrix([sp.Symbol('u_1'), sp.Symbol('u_2'), 0])] + zeroth)
                                 .coeff(epsilon, 0))
    return {'epsilon': epsilon, 'fields': fields, 'cubic_tensor': cubic_tensor,
            'quartic_tensor': quartic_tensor, 'parent': parent, 'daughters': daughters,
            'contact': contact_term, 'exchange': exchange, 'amplitude': amplitude,
            'series': series, 'exchange_series': exchange_series,
            'expected_leading': expected_leading, 'contact_function': contact}


def symmetry_factor_multiplicities():
    """Vertex multiplicities of a P_i P_j P_k Lagrangian monomial and width factors 1/S.

    Each Lagrangian coefficient g multiplies a*P_i*P_j*P_k with the fields
    repeated as in the channel; the amplitude is the derivative with respect
    to the four external legs (no 1/n!), the width carries 1/S, S = prod n_r!.
    """
    coupling = sp.Symbol('g')
    alp, first, second, third = sp.symbols('a P_1 P_2 P_3')
    channels = {'iii': (first, first, first), 'iij': (first, first, second),
                'ijk': (first, second, third)}
    result = {}
    for name, legs in channels.items():
        monomial = coupling*alp*sp.prod(legs)
        multiplicity = sp.diff(monomial, alp, *legs)/coupling
        identical = sp.prod(sp.factorial(legs.count(field)) for field in set(legs))
        result[name] = (multiplicity, identical)
    return result


# ---------------------------------------------------------------------------
# Weak-mixing limit in flavour form and the small-mixing records
# ---------------------------------------------------------------------------
def weak_mixing_leading_term(inverse_operator_matrix, daughters, scalar_products,
                             quark_masses, first_minus_derivative):
    """O(f_pi/f_a) term: V4[(2 f_pi Y/f_a, 1), (lambda_1,0), (lambda_2,0), (lambda_3,0)].

    Y = M_{m_a^2}^{-1}(C) at the parent mass. Only V4_kin, the chi_+'|_0 X^4 term
    with one parent meson leg and the d_a chi_-' term with the parent ALP leg
    survive at this order; the other V4_chi terms need two or more ALP legs.
    """
    zero = sp.zeros(quark_masses.rows)
    parent = (2*f_pi*inverse_operator_matrix/f_a, 1)
    plus = [4*B_0*quark_masses, zero, zero, zero, zero]
    minus = [zero, first_minus_derivative, zero, zero, zero]
    return contact_vertex([parent] + [(matrix, 0) for matrix in daughters],
                          scalar_products, plus, minus)


def small_mixing_ijk(axial_hat, pseudoscalar_hat, quark_masses, daughters, parent_dot_products):
    """K[T_{12;3} p.p_3 + T_{13;2} p.p_2 + T_{23;1} p.p_1] + M T^P_{123} (generic-generator record)."""
    coupling_k = 1/(6*f_a*f_pi)
    coupling_m = -B_0/(6*f_a*f_pi)
    source_p = anticommutator(pseudoscalar_hat, quark_masses)

    def axial_trace(first, second, differentiated):
        return sp.trace(axial_hat*(commutator(first, commutator(second, differentiated))
                                   + commutator(second, commutator(first, differentiated))))
    one, two, three = daughters
    ordered = sum((sp.prod(order, start=sp.eye(quark_masses.rows))
                   for order in itertools.permutations(daughters)), sp.zeros(quark_masses.rows))
    return sp.expand(coupling_k*(axial_trace(one, two, three)*parent_dot_products[3]
                                 + axial_trace(one, three, two)*parent_dot_products[2]
                                 + axial_trace(two, three, one)*parent_dot_products[1])
                     + coupling_m*sp.trace(ordered*source_p))


def weak_mixing_residuals():
    """Weak limit of the generic amplitude minus the small-mixing forms.

    Symbolic generic Y = M_{m_a^2}^{-1}(C) and c_q^{AP} (not Hermitian, not
    diagonal), symbolic diagonal m_q. Daughters are eigen-directions of the LO
    meson mass operator, m_r^2 lambda_r = B_0{m_q,lambda_r} + (m_0^2/N_f) <lambda_r>,
    which is their LO on-shell condition at zero mixing.
    """
    quark_masses = sp.diag(*sp.symbols('m_u m_d m_s', positive=True))
    m_u, m_d, m_s = (quark_masses[i, i] for i in range(3))
    inverse_matrix = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f'y_{i}{j}'))
    axial_pseudoscalar = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f'c_{i}{j}'))
    first_minus = -4*sp.I*B_0*anticommutator(axial_pseudoscalar, quark_masses)/f_a
    axial_hat = -inverse_matrix
    pseudoscalar_hat = -axial_pseudoscalar - inverse_matrix
    parent_mass_sq = sp.Symbol('m_a^2', positive=True)
    s12, s13, s23 = sp.symbols('s_12 s_13 s_23', real=True)
    unit = lambda row, column: sp.Matrix(matrix_unit(row, column))
    neutral_pion = sp.diag(1, -1, 0)
    pion_plus, pion_minus = sp.sqrt(2)*unit(0, 1), sp.sqrt(2)*unit(1, 0)
    flavour_ud, flavour_ds, flavour_us = unit(0, 1) + unit(1, 0), unit(1, 2) + unit(2, 1), unit(0, 2) + unit(2, 0)
    isospin = {m_d: m_u}
    cases = {
        'iii (ds-ds-ds)': ((flavour_ds,)*3, (B_0*(m_d + m_s),)*3, {}),
        'iij (ds-ds-us)': ((flavour_ds, flavour_ds, flavour_us),
                           (B_0*(m_d + m_s), B_0*(m_d + m_s), B_0*(m_u + m_s)), {}),
        'ijk (ud-ds-us)': ((flavour_ud, flavour_ds, flavour_us),
                           (B_0*(m_u + m_d), B_0*(m_d + m_s), B_0*(m_u + m_s)), {}),
        '000': ((neutral_pion,)*3, (2*B_0*m_u,)*3, isospin),
        '0+-': ((neutral_pion, pion_plus, pion_minus), (2*B_0*m_u,)*3, isospin),
    }
    residuals = {}
    for name, (daughters, masses, constraint) in cases.items():
        for daughter, mass in zip(daughters, masses):  # LO on-shell eigen-condition
            assert sp.expand((B_0*anticommutator(quark_masses, daughter) - mass*daughter).subs(constraint)) == sp.zeros(3)
            assert sp.trace(daughter) == 0
        mass_map = dict(zip((1, 2, 3), masses))
        invariants = {frozenset({1, 2}): s12, frozenset({1, 3}): s13, frozenset({2, 3}): s23}
        sum_rule = {s13: parent_mass_sq + sum(masses) - s12 - s23}
        products = decay_scalar_products(parent_mass_sq, mass_map, invariants)
        weak = weak_mixing_leading_term(inverse_matrix, daughters, products, quark_masses, first_minus)
        parent_dots = {r: -products[(0, r)] for r in (1, 2, 3)}
        if name.startswith('iii') or name == '000':
            coupling_m = -B_0/(6*f_a*f_pi)
            reference = 6*coupling_m*sp.trace(daughters[0]**3*anticommutator(pseudoscalar_hat, quark_masses))
        elif name.startswith('iij'):
            lam_i, lam_j = daughters[0], daughters[2]
            t_axial = sp.trace(axial_hat*commutator(lam_i, commutator(lam_i, lam_j)))
            t_pseudo = sp.trace((lam_i*lam_i*lam_j + lam_i*lam_j*lam_i + lam_j*lam_i*lam_i)
                                * anticommutator(pseudoscalar_hat, quark_masses))
            reference = (t_axial*(parent_mass_sq + 3*masses[2] - 3*s12)/(12*f_a*f_pi)
                         - 2*B_0*t_pseudo/(6*f_a*f_pi))
        else:  # ijk and the record's 0+- (same T^A, T^P structure)
            reference = small_mixing_ijk(axial_hat, pseudoscalar_hat, quark_masses, daughters, parent_dots)
        difference = sp.expand((weak - reference).subs(sum_rule).subs(constraint))
        residuals[name] = sp.simplify(difference)
        residuals[name + ' leading term nonzero'] = sp.expand(weak.subs(sum_rule).subs(constraint)) != 0
        explicit = weak_mixing_explicit_form(inverse_matrix, axial_pseudoscalar, quark_masses,
                                             daughters, invariants)
        residuals[name + ' explicit flavour form'] = sp.expand(weak - explicit)
    # The equality needs the daughter on-shell condition: with an arbitrary
    # mass for the us daughter in the iij case the difference does not vanish.
    daughters = (flavour_ds, flavour_ds, flavour_us)
    masses = (B_0*(m_d + m_s), B_0*(m_d + m_s), sp.Symbol('m_j^2', positive=True))
    invariants = {frozenset({1, 2}): s12, frozenset({1, 3}): s13, frozenset({2, 3}): s23}
    products = decay_scalar_products(parent_mass_sq, dict(zip((1, 2, 3), masses)), invariants)
    weak = weak_mixing_leading_term(inverse_matrix, daughters, products, quark_masses, first_minus)
    t_axial = sp.trace(axial_hat*commutator(flavour_ds, commutator(flavour_ds, flavour_us)))
    t_pseudo = sp.trace((flavour_ds**2*flavour_us + flavour_ds*flavour_us*flavour_ds + flavour_us*flavour_ds**2)
                        * anticommutator(pseudoscalar_hat, quark_masses))
    reference = (t_axial*(parent_mass_sq + 3*masses[2] - 3*s12)/(12*f_a*f_pi) - 2*B_0*t_pseudo/(6*f_a*f_pi))
    off_shell = sp.expand((weak - reference).subs(s13, parent_mass_sq + sum(masses) - s12 - s23))
    residuals['iij needs on-shell daughters (difference nonzero off shell)'] = off_shell != 0
    return residuals


def weak_mixing_explicit_form(inverse_matrix, axial_pseudoscalar, quark_masses, daughters, invariants):
    """Explicit O(1/f_a) flavour form written in the memory record (linear sources).

    - (1/(12 f_pi f_a)) [<[Y,l1][l2,l3]>(s13-s12) + <[Y,l2][l1,l3]>(s23-s12) + <[Y,l3][l1,l2]>(s23-s13)]
    + (B_0/(12 f_pi f_a)) sum_{S_4} <m_q . ordered product of Y, l1, l2, l3>
    + (B_0/(6 f_pi f_a)) sum_{S_3} <{c_q^AP, m_q} l l l>.
    """
    one, two, three = daughters
    s12, s13, s23 = (invariants[frozenset(pair)] for pair in ((1, 2), (1, 3), (2, 3)))
    kinetic = -(sp.trace(commutator(inverse_matrix, one)*commutator(two, three))*(s13 - s12)
                + sp.trace(commutator(inverse_matrix, two)*commutator(one, three))*(s23 - s12)
                + sp.trace(commutator(inverse_matrix, three)*commutator(one, two))*(s23 - s13))/(12*f_pi*f_a)
    ordered_four = sum((sp.prod(order, start=sp.eye(3))
                        for order in itertools.permutations((inverse_matrix, one, two, three))), sp.zeros(3))
    ordered_three = sum((sp.prod(order, start=sp.eye(3))
                         for order in itertools.permutations(daughters)), sp.zeros(3))
    mass = B_0*sp.trace(quark_masses*ordered_four)/(12*f_pi*f_a)
    source = B_0*sp.trace(anticommutator(axial_pseudoscalar, quark_masses)*ordered_three)/(6*f_pi*f_a)
    return sp.expand(kinetic + mass + source)


# ---------------------------------------------------------------------------
# Pole-aligned coordinates (route 2): transformation taken from the notes
# ---------------------------------------------------------------------------
# memory/alp/shared/alp-lo-pole-aligned-tree-vertices.md, step 1, with
# w_i lambda_i = (2 f_pi/f_a) M_{m_a^2}^{-1}(C) (three-pion record, section 2):
#   U_df = exp(i a Y/f_a) U' exp(i a Y/f_a),  Y = M_{m_a^2}^{-1}(C).
# In the convention U' = g_R U g_L^dagger of alp-chiral-source-transformations.md
# this is g_R = exp(-i x Y), g_L = exp(+i x Y) (x = a/f_a), i.e. G_R = -G_L = -Y
# of alp-pole-alignment-chiral-generator.md. The source laws of
# alp-chiral-source-transformations.md then give, from zero derivative sources,
#   r' = -(d x) Y,  l' = +(d x) Y   (C_V' = 0, C_A' = -Y exactly),
#   chi' -> g_R chi' g_L^dagger,  Theta' -> Theta' + x <Y>.
# The functions below only implement these stated laws; they are checked
# against the joint CCWZ expansion and against route 1.

def aligned_ccwz_derivatives(inverse_matrix, plus, minus):
    """a-derivatives of chi_+- at Phi=0 after the alignment rotation g_R=exp(-iaY/f_a), g_L=g_R^dagger."""
    return ccwz_vacuum_derivatives(inverse_matrix, sp.zeros(inverse_matrix.rows), plus, minus)


def dressed_coefficients(plus, minus, order=4):
    """Dressed Taylor coefficients (x-derivatives) from the derivative-free blocks.

    In the derivative-free coordinates the spurion is chi = 2 B_0 widehat M(x)
    (general-source record), so chi_+|_0 = 4 B_0 widehat S and chi_-|_0 = 4 i B_0 widehat P:
    widehat S_n = f_a^n d_a^n chi_+'|_0/(4 B_0), widehat P_n = f_a^n d_a^n chi_-'|_0/(4 i B_0).
    """
    return ([f_a**n*plus[n]/(4*B_0) for n in range(order + 1)],
            [f_a**n*minus[n]/(4*sp.I*B_0) for n in range(order + 1)])


def aligned_source_coefficients(hat_axial, dressed_scalar, dressed_pseudoscalar, order=4):
    """Taylor coefficients S_n', P_n' in pole-aligned coordinates.

    Exact Taylor law of the general-source record, M_n' = sum_{r+s+t=n}
    n!/(r!s!t!) g_{R,r} M_t g^dagger_{L,s}, applied to the derivative-free frame
    (M_t = widehat M_t, no currents) with g_R = exp(i x c_A_hat), g_L = g_R^dagger,
    c_A_hat = -M_{m_a^2}^{-1}(C); then S' and P' are its Hermitian and anti-Hermitian parts.
    """
    dressed = [dressed_scalar[n] + sp.I*dressed_pseudoscalar[n] for n in range(order + 1)]
    rotation = sp.I*hat_axial
    transformed = []
    for n in range(order + 1):
        total = sp.zeros(hat_axial.rows)
        for r in range(n + 1):
            for t in range(n - r + 1):
                s_power = n - r - t
                total += (sp.factorial(n)/(sp.factorial(r)*sp.factorial(s_power)*sp.factorial(t))
                          * rotation**r*dressed[t]*rotation**s_power)
        transformed.append(sp.expand(total))
    scalar = [sp.expand((m + m.H)/2) for m in transformed]
    pseudoscalar = [sp.expand((m - m.H)/(2*sp.I)) for m in transformed]
    return scalar, pseudoscalar


def joint_vertex_sources(plus, minus, right, left, theta_slope, numbers, anomaly_mass_sq,
                         bare_alp_mass_sq):
    """Source tuple of derive_alp_lo_joint_vertices (x-derivatives of chi = (chi_+ + chi_-)/2)."""
    fa = numbers[f_a]
    chi = [sp.Matrix(fa**n*(plus[n] + minus[n])/2).subs(numbers) for n in range(len(plus))]
    theta = [0, theta_slope] + [0]*(len(plus) - 2)
    return (chi, right, left, theta, sp.sympify(numbers[f_pi]), sp.sympify(fa),
            sp.sympify(anomaly_mass_sq), sp.sympify(bare_alp_mass_sq))


def derivative_free_and_aligned_sources(inverse_matrix, plus, minus, gluon_axial, numbers,
                                        anomaly_mass_sq, bare_alp_mass_sq):
    """Joint-vertex sources in derivative-free and in pole-aligned coordinates."""
    size = inverse_matrix.rows
    zero = [sp.zeros(size)]*len(plus)
    free = joint_vertex_sources(plus, minus, zero, zero, gluon_axial, numbers,
                                anomaly_mass_sq, bare_alp_mass_sq)
    aligned_plus, aligned_minus = aligned_ccwz_derivatives(inverse_matrix, plus, minus)
    aligned_plus = [sp.Matrix(matrix).subs(numbers) for matrix in aligned_plus]
    aligned_minus = [sp.Matrix(matrix).subs(numbers) for matrix in aligned_minus]
    right = [-inverse_matrix] + [sp.zeros(size)]*(len(plus) - 1)
    left = [inverse_matrix] + [sp.zeros(size)]*(len(plus) - 1)
    aligned = joint_vertex_sources(aligned_plus, aligned_minus, right, left,
                                   gluon_axial + sp.trace(inverse_matrix), numbers,
                                   anomaly_mass_sq, bare_alp_mass_sq)
    return free, aligned, aligned_plus, aligned_minus


def basis_components(matrix, generators):
    """Components c_I of a matrix in a complete basis: c = K^{-1} (<E_J X>/2), K_IJ = <E_I E_J>/2."""
    gram = sp.Matrix(len(generators), len(generators),
                     lambda i, j: sp.trace(generators[i]*generators[j])/2)
    projections = sp.Matrix([sp.trace(generator*matrix)/2 for generator in generators])
    return gram.LUsolve(projections)


def pole_operator_inverse(source, quark_masses, b0_value, anomaly_mass_sq, squared_mass):
    """Solve M_s(X) = B_0{X,m_q} + (m_0^2/N_f) <X> 1 - s X = source with numerical B_0."""
    size = quark_masses.rows
    identity = sp.eye(size)
    half = b0_value*quark_masses - squared_mass*identity/2
    operator = (sp.kronecker_product(identity, half) + sp.kronecker_product(half.T, identity)
                + anomaly_mass_sq*identity.vec()*identity.vec().T/size)
    solution = operator.LUsolve(source.vec())
    return sp.Matrix(size, size, lambda row, column: solution[column*size + row])


def joint_contact(directions, products, sources):
    """Four-leg vertex of the joint CCWZ expansion, momenta replaced by k_r.k_s."""
    momenta = sp.symbols('_k0:4', real=True)
    expression = joint_vertex(
        [direction[0] for direction in directions], [direction[1] for direction in directions],
        list(momenta), *sources)
    polynomial = sp.Poly(sp.expand(expression), *momenta)
    total = sp.S.Zero
    for powers, coefficient in polynomial.terms():
        legs = [leg for leg, power in enumerate(powers) for _ in range(power)]
        if not legs:
            total += coefficient
        else:
            first, second = legs
            total += coefficient*(sp.Symbol('_m2_%d' % first) if first == second
                                  else products[tuple(sorted((first, second)))])
    return sp.expand(total)


def aligned_kinetic_contact(inverse_matrix, parent_alp, daughters, products):
    """Kinetic part of the pole-aligned contact for a parent column (0, parent_alp).

    From (f_pi^2/4) <u_mu u^mu> with axial source C_A' = -Y:
      (1/(6 f_pi f_a)) d_mu a <c_A' [Phi,[Phi,d^mu Phi]]>   (one ALP leg, three mesons)
      (1/(4 f_a^2)) (d a)^2 <[c_A', Phi]^2>                 (two ALP legs, two mesons)
    with c_A' = -Y. Only the parent ALP leg can carry the single-ALP term
    because the parent has no meson component in these coordinates.
    """
    axial = -inverse_matrix
    legs = dict(zip((1, 2, 3), daughters))

    def axial_trace(first, second, differentiated):
        return sp.trace(axial*(commutator(first, commutator(second, differentiated))
                               + commutator(second, commutator(first, differentiated))))
    meson = {r: legs[r][0] for r in legs}
    parent_dot = {r: -products[(0, r)] for r in legs}
    one_alp = (axial_trace(meson[1], meson[2], meson[3])*parent_dot[3]
               + axial_trace(meson[1], meson[3], meson[2])*parent_dot[2]
               + axial_trace(meson[2], meson[3], meson[1])*parent_dot[1])/(6*f_pi*f_a)
    two_alp = sum(legs[r][1]*parent_dot[r]
                  * sp.trace(commutator(axial, meson[u])*commutator(axial, meson[v]))
                  for r, (u, v) in PAIRINGS)/f_a**2
    return sp.expand(parent_alp*(one_alp + two_alp))


def aligned_contact(inverse_matrix, parent_alp, daughters, products, aligned_plus, aligned_minus):
    """Full pole-aligned contact: kinetic part plus V4_chi with the aligned chi derivatives."""
    zero = sp.zeros(inverse_matrix.rows)
    mass = quartic_mass_vertex([(zero, parent_alp)] + list(daughters), aligned_plus[0],
                               aligned_minus[1], aligned_plus[2], aligned_minus[3],
                               aligned_plus[4]).doit()
    return sp.expand(aligned_kinetic_contact(inverse_matrix, parent_alp, daughters, products) + mass)


def transformation_checks():
    """Verify the notes' alignment transformation with exact rational data.

    Returns a dict of booleans: derivative-free blocks of the joint expansion
    equal the three-pion record's blocks; aligned blocks equal the congruence
    T^T (K, M^2) T with Phi = Phi' + a W, W = 2 f_pi Y / f_a; the meson-ALP
    entry of m_a^2 K' - M'^2 vanishes for Y = M_{m_a^2}^{-1}(C) at m_a^2; the first
    aligned chi derivatives are the hat couplings of the generic-generator
    record; cubic vertices transform with the linear shift only; the aligned
    contact (kinetic and chi_+ parts) equals the joint CCWZ expansion.
    """
    numbers = {f_pi: 1, f_a: 3, B_0: 1}
    quark_masses = sp.diag(1, 2, 4)
    axial = sp.Matrix([[1, 1, 0], [1, 0, 1], [0, 1, 2]])
    pseudoscalar = sp.Matrix([[0, 1, 1], [1, 1, 0], [1, 0, -1]])
    vector = sp.Matrix([[0, 1, 0], [1, 1, 2], [0, 2, 0]])
    scalar_input = sp.Matrix([[1, 0, 1], [0, 0, 1], [1, 1, 0]])
    anomaly_mass_sq, gluon, bare_alp_mass_sq, parent_mass_sq = 3, sp.Rational(1, 2), 13, sp.Rational(27, 2)
    plus, minus = linear_source_ccwz_derivatives(axial, pseudoscalar, vector, scalar_input, quark_masses)
    plus = [sp.Matrix(m).subs(numbers) for m in plus]
    minus = [sp.Matrix(m).subs(numbers) for m in minus]
    axial_pseudoscalar = axial - pseudoscalar
    gluon_axial = gluon + sp.trace(axial)
    coupling = (-B_0*anticommutator(quark_masses, axial_pseudoscalar)
                - anomaly_mass_sq*gluon_axial*sp.eye(3)/3).subs(numbers)
    inverse_matrix = inverse_meson_mass_operator(coupling, quark_masses, anomaly_mass_sq,
                                                 parent_mass_sq).subs(numbers)
    free, aligned, aligned_plus, aligned_minus = derivative_free_and_aligned_sources(
        inverse_matrix, plus, minus, gluon_axial, numbers, anomaly_mass_sq, bare_alp_mass_sq)
    generators = hermitian_unit_basis()
    results = {}
    kinetic, mass_squared = joint_quadratic_blocks(generators, *free)
    record_kinetic, record_mass = component_quadratic_blocks(
        generators, quark_masses, anomaly_mass_sq, coupling,
        bare_alp_mass_sq - f_pi**2*sp.trace(plus[2])/4
        + 2*f_pi**2*anomaly_mass_sq*gluon_axial**2/(3*f_a**2))
    results['derivative-free blocks = three-pion record section 2'] = (
        sp.simplify(kinetic - record_kinetic) == sp.zeros(10)
        and sp.simplify(mass_squared - record_mass.subs(numbers)) == sp.zeros(10))
    shift = 2*numbers[f_pi]*inverse_matrix/numbers[f_a]
    conversion = sp.eye(10)
    conversion[0:9, 9] = basis_components(shift, generators)
    aligned_kinetic, aligned_mass = joint_quadratic_blocks(generators, *aligned)
    results['aligned blocks = T^T (K, M^2) T (congruence of the single-pole record)'] = (
        sp.simplify(aligned_kinetic - conversion.T*kinetic*conversion) == sp.zeros(10)
        and sp.simplify(aligned_mass - conversion.T*mass_squared*conversion) == sp.zeros(10))
    mixing = (parent_mass_sq*aligned_kinetic - aligned_mass)[0:9, 9]
    results['meson-ALP entry of m_a^2 K\' - M\'^2 vanishes at the parent pole'] = (
        sp.simplify(mixing) == sp.zeros(9, 1) and aligned_mass[0:9, 9] != sp.zeros(9, 1))
    results['K\'_aa = N_a = 1 + 2 (f_pi/f_a)^2 <Y^2>'] = sp.simplify(
        aligned_kinetic[9, 9] - (1 + 2*sp.trace(inverse_matrix**2)*numbers[f_pi]**2/numbers[f_a]**2)) == 0
    hat_pseudoscalar = -axial_pseudoscalar - inverse_matrix
    results['first aligned derivatives: d chi_+ unchanged, d chi_- = 4iB_0{c_P_hat,m_q}/f_a'] = (
        sp.simplify(aligned_plus[1] - plus[1]) == sp.zeros(3)
        and sp.simplify(aligned_minus[1] - (4*sp.I*B_0*anticommutator(hat_pseudoscalar, quark_masses)
                                            / f_a).subs(numbers)) == sp.zeros(3))
    # Aligned derivatives from the notes' invariants: chi_+ = 4B_0 S', chi_- = 4iB_0 P'.
    dressed_scalar, dressed_pseudoscalar = dressed_coefficients(plus, minus)
    dressed_scalar = [sp.Matrix(m).subs(numbers) for m in dressed_scalar]
    dressed_pseudoscalar = [sp.Matrix(m).subs(numbers) for m in dressed_pseudoscalar]
    hat_axial = -inverse_matrix
    aligned_scalar, aligned_pseudoscalar = aligned_source_coefficients(hat_axial, dressed_scalar, dressed_pseudoscalar)
    taylor_ok = all(
        sp.simplify((4*B_0*aligned_scalar[n]/f_a**n).subs(numbers) - aligned_plus[n]) == sp.zeros(3)
        and sp.simplify((4*sp.I*B_0*aligned_pseudoscalar[n]/f_a**n).subs(numbers) - aligned_minus[n]) == sp.zeros(3)
        for n in range(5))
    results["aligned d_a^n chi_+- = 4B_0 S_n'/f_a^n, 4iB_0 P_n'/f_a^n (Taylor law), orders 0-4"] = taylor_ok
    masses_num = quark_masses
    results['dressed first order = first-order record: S_1 = -i[c_VS,m_q], P_1 = -{c_AP,m_q}'] = (
        sp.simplify(dressed_scalar[1] + sp.I*commutator(vector - scalar_input, masses_num)) == sp.zeros(3)
        and sp.simplify(dressed_pseudoscalar[1] + anticommutator(axial_pseudoscalar, masses_num)) == sp.zeros(3))
    # Second-order record, Alda specialization with S_2 = P_2 = C_V1 = C_A1 = 0 (linear sources).
    second_scalar = (commutator(vector, commutator(2*scalar_input - vector, masses_num))
                     + anticommutator(axial, anticommutator(2*pseudoscalar - axial, masses_num)))
    second_pseudoscalar = (-sp.I*commutator(vector, anticommutator(2*pseudoscalar - axial, masses_num))
                           - sp.I*anticommutator(axial, commutator(2*scalar_input - vector, masses_num)))
    results['dressed second order = second-order record (linear sources)'] = (
        sp.simplify(dressed_scalar[2] - second_scalar) == sp.zeros(3)
        and sp.simplify(dressed_pseudoscalar[2] - second_pseudoscalar) == sp.zeros(3))
    explicit_first = (sp.simplify(aligned_scalar[1] - dressed_scalar[1]) == sp.zeros(3)
                      and sp.simplify(aligned_pseudoscalar[1] - dressed_pseudoscalar[1]
                                      - anticommutator(hat_axial, masses_num)) == sp.zeros(3))
    explicit_second = (
        sp.simplify(aligned_scalar[2] - (dressed_scalar[2] - 2*anticommutator(hat_axial, dressed_pseudoscalar[1])
                                 - anticommutator(hat_axial, anticommutator(hat_axial, masses_num)))) == sp.zeros(3)
        and sp.simplify(aligned_pseudoscalar[2] - (dressed_pseudoscalar[2]
                                           + 2*anticommutator(hat_axial, dressed_scalar[1]))) == sp.zeros(3))
    # Inverse direction: real-mass closed form of the second-order record with C_V = 0,
    # C_A0 = c_A_hat, C_A1 = 0 maps (S_n', P_n') back to widehat S_2, widehat P_2.
    closed_back = (
        sp.simplify(aligned_scalar[2] + 2*anticommutator(hat_axial, aligned_pseudoscalar[1])
                    - anticommutator(hat_axial, anticommutator(hat_axial, masses_num)) - dressed_scalar[2]) == sp.zeros(3)
        and sp.simplify(aligned_pseudoscalar[2] - 2*anticommutator(hat_axial, aligned_scalar[1]) - dressed_pseudoscalar[2]) == sp.zeros(3))
    results["S_1' = widehat S_1, P_1' = widehat P_1 + {c_A_hat,m_q}"] = explicit_first
    results["S_2', P_2' explicit forms; second-order closed form maps them back to widehat S_2, P_2"] = (
        explicit_second and closed_back)
    # Cubic vertices: linear shift only (no quadratic term in the coordinate map).
    directions = [(sp.Matrix([[1, 0, 1], [0, 0, 1], [1, 1, -1]]), sp.Rational(1, 2)),
                  (sp.Matrix([[0, 1, 0], [1, 2, 0], [0, 0, 1]]), sp.Rational(-1, 3)),
                  (sp.Matrix([[2, 0, 0], [0, -1, 1], [0, 1, 0]]), 1)]
    shifted = [(matrix + alp*shift, alp) for matrix, alp in directions]
    aligned_cubic = cubic_vertex(*directions, aligned_plus[1], aligned_minus[2], aligned_plus[3]).doit()
    free_cubic = cubic_vertex(*shifted, plus[1], minus[2], plus[3]).doit()
    results['aligned cubic vertex = derivative-free cubic on shifted directions'] = (
        sp.simplify((aligned_cubic - free_cubic).subs(numbers)) == 0)
    # Contact with parent (0, 1): explicit kinetic and chi_+ parts against the joint expansion.
    masses = sp.symbols('_m2_0:4')
    invariants = dict(zip(map(frozenset, ((1, 2), (1, 3), (2, 3))), sp.symbols('s_12 s_13 s_23')))
    products = decay_scalar_products(masses[0], dict(zip((1, 2, 3), masses[1:])), invariants)
    joint = joint_contact([(sp.zeros(3), 1)] + directions, products, aligned)
    explicit = aligned_contact(inverse_matrix, 1, directions, products, aligned_plus, aligned_minus).subs(numbers)
    explicit_mass = quartic_mass_vertex([(sp.zeros(3), 1)] + directions, aligned_plus[0], aligned_minus[1],
                                        aligned_plus[2], aligned_minus[3], aligned_plus[4]).doit().subs(numbers)
    momentum_free = joint.subs({symbol: 0 for symbol in invariants.values()}).subs({m: 0 for m in masses})
    results['chi_+ contact (momentum-independent part) = V4_chi with aligned derivatives'] = (
        sp.simplify(momentum_free - explicit_mass) == 0 and explicit_mass != 0)
    results['full aligned contact = joint CCWZ expansion'] = sp.simplify(sp.expand(joint - explicit)) == 0
    return results


def numerical_route_comparison(off_shell_shift=0):
    """Route 1 (derivative-free contraction) versus route 2 (pole-aligned) at LO poles.

    Double-precision poles of the same exact sources (spectrum solver). Returns
    (route 1, route 2, contact difference, parent and daughter column checks).
    off_shell_shift moves daughter 1's squared mass in the kinematics only.
    """
    from solve_alp_meson_spectrum import solve_spectrum_numerical
    numbers = {f_pi: 1, f_a: 3, B_0: 1}
    quark_masses = sp.diag(1, 2, 4)
    axial = sp.Matrix([[1, 1, 0], [1, 0, 1], [0, 1, 2]])
    pseudoscalar = sp.Matrix([[0, 1, 1], [1, 1, 0], [1, 0, -1]])
    vector = sp.Matrix([[0, 1, 0], [1, 1, 2], [0, 2, 0]])
    scalar = sp.Matrix([[1, 0, 1], [0, 0, 1], [1, 1, 0]])
    anomaly_mass_sq, gluon, bare_alp_mass_sq = 3, sp.Rational(1, 2), 13
    plus, minus = linear_source_ccwz_derivatives(axial, pseudoscalar, vector, scalar, quark_masses)
    plus = [sp.Matrix(m).subs(numbers) for m in plus]
    minus = [sp.Matrix(m).subs(numbers) for m in minus]
    gluon_axial = gluon + sp.trace(axial)
    coupling = (-B_0*anticommutator(quark_masses, axial - pseudoscalar)
                - anomaly_mass_sq*gluon_axial*sp.eye(3)/3).subs(numbers)
    generators = hermitian_unit_basis()
    alp_entry = (bare_alp_mass_sq - f_pi**2*sp.trace(plus[2])/4
                 + 2*f_pi**2*anomaly_mass_sq*gluon_axial**2/(3*f_a**2))
    kinetic, mass_squared = component_quadratic_blocks(generators, quark_masses, anomaly_mass_sq,
                                                       coupling, alp_entry)
    mass_squared = mass_squared.subs(numbers)
    spectrum = solve_spectrum_numerical(sp.matrix2numpy(kinetic, dtype=float),
                                        sp.matrix2numpy(mass_squared, dtype=float))
    vectors = spectrum.pole_vectors
    poles = spectrum.pole_squared_masses
    parent_index = max(range(10), key=lambda n: abs(vectors[9, n]))
    meson_like = sorted((n for n in range(10) if n != parent_index), key=lambda n: poles[n])

    def column(n):
        return (sum((sp.Float(vectors[i, n])*generators[i] for i in range(9)), sp.zeros(3)),
                sp.Float(vectors[9, n]))
    parent = column(parent_index)
    chosen = [meson_like[0], meson_like[1], meson_like[0]]
    daughters = [column(n) for n in chosen]
    parent_mass_sq = sp.Float(poles[parent_index])
    masses = {r: sp.Float(poles[n]) for r, n in zip((1, 2, 3), chosen)}
    s12, s13 = sp.Float(3.1), sp.Float(4.7)
    kinematic_masses = dict(masses)
    kinematic_masses[1] = masses[1] + off_shell_shift
    invariants = {frozenset({1, 2}): s12, frozenset({1, 3}): s13,
                  frozenset({2, 3}): parent_mass_sq + sum(kinematic_masses.values()) - s12 - s13}
    basis = field_basis_directions(generators)
    route_one, contact_one, _ = generic_amplitude(parent, daughters, parent_mass_sq, kinematic_masses,
                                                  invariants, plus, minus, basis, kinetic, mass_squared)
    route_one = route_one.subs(numbers)
    # Route 2: W from the pole column; check it against the operator inverse.
    shift = parent[0]/parent[1]
    inverse_matrix = numbers[f_a]*shift/(2*numbers[f_pi])
    operator_inverse = pole_operator_inverse(coupling, quark_masses, numbers[B_0], anomaly_mass_sq,
                                             parent_mass_sq)
    norm_parent = 1 + 2*sp.trace(operator_inverse**2)*numbers[f_pi]**2/numbers[f_a]**2
    checks = {'parent column = (2 f_pi M_{m_a^2}^{-1}(C)/f_a, 1)/sqrt(N_a)':
              max(abs(x) for x in (shift - 2*numbers[f_pi]*operator_inverse/numbers[f_a])) < 1e-8
              and abs(parent[1]**2*norm_parent - 1) < 1e-8}
    aligned_daughters = [(matrix - alp*shift, alp) for matrix, alp in daughters]
    for r, ((matrix, alp), n) in enumerate(zip(aligned_daughters, chosen), start=1):
        daughter_inverse = pole_operator_inverse(coupling, quark_masses, numbers[B_0], anomaly_mass_sq,
                                                 sp.Float(poles[n]))
        expected = 2*numbers[f_pi]*(daughter_inverse - operator_inverse)*alp/numbers[f_a]
        checks[f'aligned daughter {r} = (2 f_pi (M_{{m_r^2}}^{{-1}} - M_{{m_a^2}}^{{-1}})(C)/f_a, 1)/sqrt(N_r)'] = (
            max(abs(x) for x in (matrix - expected)) < 1e-7)
        # Same column through c_A_hat: (M_{m_r^2}^{-1} - M_{m_a^2}^{-1})(C) = (m_a^2 - m_r^2) M_{m_r^2}^{-1}(c_A_hat).
        through_axial = pole_operator_inverse(-operator_inverse, quark_masses, numbers[B_0], anomaly_mass_sq,
                                              sp.Float(poles[n]))
        expected_axial = 2*numbers[f_pi]*(parent_mass_sq - sp.Float(poles[n]))*through_axial*alp/numbers[f_a]
        checks[f'aligned daughter {r} = 2 f_pi (m_a^2 - m_r^2) M_{{m_r^2}}^{{-1}}(c_A_hat)/(f_a sqrt(N_r))'] = (
            max(abs(x) for x in (matrix - expected_axial)) < 1e-7)
    aligned_plus, aligned_minus = aligned_ccwz_derivatives(inverse_matrix, plus, minus)
    aligned_plus = [sp.Matrix(m).subs(numbers) for m in aligned_plus]
    aligned_minus = [sp.Matrix(m).subs(numbers) for m in aligned_minus]
    products = decay_scalar_products(parent_mass_sq, kinematic_masses, invariants)
    contact_two = aligned_contact(inverse_matrix, parent[1], aligned_daughters,
                                  products, aligned_plus, aligned_minus).subs(numbers)
    conversion = sp.eye(10)
    conversion[0:9, 9] = basis_components(shift, generators)
    aligned_kinetic = conversion.T*kinetic*conversion
    aligned_mass = conversion.T*mass_squared*conversion
    route_two = contact_two
    legs = dict(zip((1, 2, 3), aligned_daughters))
    parent_aligned = (sp.zeros(3), parent[1])
    for attached, (first, second) in PAIRINGS:
        invariant = invariants[frozenset({first, second})]
        left = sp.Matrix([cubic_vertex(parent_aligned, legs[attached], direction, aligned_plus[1],
                                       aligned_minus[2], aligned_plus[3]).doit().subs(numbers)
                          for direction in basis])
        right = sp.Matrix([cubic_vertex(direction, legs[first], legs[second], aligned_plus[1],
                                        aligned_minus[2], aligned_plus[3]).doit().subs(numbers)
                           for direction in basis])
        route_two += (left.T*(aligned_mass - invariant*aligned_kinetic).LUsolve(right))[0]
    # Step-4 compact amplitude with the step-2 columns built from M_{m_r^2}^{-1}(c_A_hat).
    hat_axial = -inverse_matrix
    record_daughters = {}
    for label, n, (_, alp) in zip('ijk', chosen, aligned_daughters):
        inverse_axial_at = pole_operator_inverse(hat_axial, quark_masses, numbers[B_0], anomaly_mass_sq,
                                                 sp.Float(poles[n]))
        record_daughters[label] = (daughter_column(inverse_axial_at, parent_mass_sq - sp.Float(poles[n]),
                                                   alp)[0].subs(numbers), alp)
    fa_value = numbers[f_a]
    coefficients_s = [fa_value**n*aligned_plus[n]/(4*numbers[B_0]) for n in range(5)]
    coefficients_p = [fa_value**n*aligned_minus[n]/(4*sp.I*numbers[B_0]) for n in range(5)]
    labels = dict(zip('ijk', (1, 2, 3)))
    record_dot = {r: -products[(0, labels[r])] for r in 'ijk'}
    record_invariants = {'jk': invariants[frozenset({2, 3})], 'ik': invariants[frozenset({1, 3})],
                         'ij': invariants[frozenset({1, 2})]}
    compact = compact_amplitude(hat_axial, coefficients_p[1], sp.I*coefficients_s[1], coefficients_s,
                                coefficients_p, parent_aligned, record_daughters, record_dot, record_invariants,
                                basis, aligned_kinetic, aligned_mass).subs(numbers)
    checks['step-4 compact amplitude = generic_amplitude (route 1) at LO poles'] = (
        abs(complex(sp.N(compact)) - complex(sp.N(route_one))) < 1e-10*abs(complex(sp.N(route_one))))
    return (complex(sp.N(route_one)), complex(sp.N(route_two)),
            complex(sp.N((contact_one.subs(numbers) - contact_two))), checks)


# ---------------------------------------------------------------------------
# Compact record formulas (step 2 columns, step 3 vertices, step 4 amplitude),
# typed independently and compared with the vertex functions above.
# ---------------------------------------------------------------------------
# A column is (lambda_x, xi_x): meson matrix and ALP component, step 2 table.
# Daughter r with ALP component:
#   lambda_r = 2 (f_pi/f_a) (m_a^2 - m_r^2) M_{m_r^2}^{-1}(c_A_hat)/sqrt(N(m_r^2)),
#   xi_r = 1/sqrt(N(m_r^2)); without ALP component (lambda_r, 0); parent (0, 1/sqrt(N(m_a^2))).
# vs_commutator = [c_q^VS, m_q] = i S_1'; hat_p = {c_P_hat, m_q} = P_1'.

def daughter_column(inverse_axial_at_daughter, mass_gap, inverse_root):
    """Step-2 column of a daughter with ALP component."""
    return (2*f_pi/f_a*mass_gap*inverse_root*inverse_axial_at_daughter, inverse_root)


def tensor_axial(hat_axial, first, second, third):
    """T^A_{ij;k}(Y_i,Y_j,Y_k) = <c_A_hat([Y_i,[Y_j,Y_k]] + [Y_j,[Y_i,Y_k]])>."""
    c = commutator
    return sp.trace(hat_axial*(c(first, c(second, third)) + c(second, c(first, third))))


def tensor_pseudoscalar(hat_p, first, second, third):
    """T^P_{ijk}(Y_i,Y_j,Y_k) = <sum over the six orderings Y Y Y  P_hat>."""
    ordered = sum((sp.prod(order, start=sp.eye(first.rows))
                   for order in itertools.permutations((first, second, third))), sp.zeros(first.rows))
    return sp.trace(ordered*hat_p)


def compact_quartic(hat_axial, hat_p, s2, p3, s4, parent, daughters, parent_dot):
    """Boxed V^(4)_{aijk} of the record; daughters = {'i','j','k'} columns, parent_dot = p.p_r."""
    xi_a = parent[1]
    lam = {r: daughters[r][0] for r in 'ijk'}
    xi = {r: daughters[r][1] for r in 'ijk'}
    others = {'i': ('j', 'k'), 'j': ('i', 'k'), 'k': ('i', 'j')}
    one_alp = (sp.Rational(1, 6)/(f_pi*f_a)*(tensor_axial(hat_axial, lam['j'], lam['k'], lam['i'])*parent_dot['i']
                                             + tensor_axial(hat_axial, lam['i'], lam['k'], lam['j'])*parent_dot['j']
                                             + tensor_axial(hat_axial, lam['i'], lam['j'], lam['k'])*parent_dot['k'])
               - B_0/(6*f_pi*f_a)*tensor_pseudoscalar(hat_p, lam['i'], lam['j'], lam['k']))
    more_alp = sum(xi[r]*parent_dot[r]/f_a**2*sp.trace(commutator(hat_axial, lam[s])*commutator(hat_axial, lam[t]))
                   - B_0/(2*f_a**2)*xi[r]*sp.trace(s2*anticommutator(lam[s], lam[t]))
                   + B_0*f_pi/f_a**3*xi[s]*xi[t]*sp.trace(p3*lam[r])
                   for r, (s, t) in others.items())
    four_alp = B_0*f_pi**2/f_a**4*xi_a*xi['i']*xi['j']*xi['k']*sp.trace(s4)
    return sp.expand(xi_a*(one_alp + more_alp) + four_alp)


def compact_cubic(vs_commutator, p2, s3, first, second, third):
    """Boxed V^(3)_{xyz} of the record for three columns."""
    (lx, xx), (ly, xy), (lz, xz) = first, second, third
    ac = anticommutator
    return sp.expand(B_0/f_a*(sp.I/2*sp.trace(vs_commutator*(xx*ac(ly, lz) + xy*ac(lz, lx) + xz*ac(lx, ly)))
                              + f_pi/f_a*sp.trace(p2*(xx*xy*lz + xx*xz*ly + xy*xz*lx))
                              + f_pi**2/f_a**2*xx*xy*xz*sp.trace(s3)))


def compact_amplitude(hat_axial, hat_p, vs_commutator, coefficients_s, coefficients_p, parent, daughters,
                      parent_dot, invariants, field_basis, kinetic, mass_squared):
    """Boxed amplitude of step 4: compact quartic plus the three exchanges with G = (M'^2 - s K')^{-1}.

    invariants maps 'jk', 'ik', 'ij' to s_jk, s_ik, s_ij; field_basis lists internal columns.
    """
    total = compact_quartic(hat_axial, hat_p, coefficients_s[2], coefficients_p[3], coefficients_s[4],
                            parent, daughters, parent_dot)
    for attached, pair in (('i', 'jk'), ('j', 'ik'), ('k', 'ij')):
        left = sp.Matrix([compact_cubic(vs_commutator, coefficients_p[2], coefficients_s[3],
                                        parent, daughters[attached], internal) for internal in field_basis])
        right = sp.Matrix([compact_cubic(vs_commutator, coefficients_p[2], coefficients_s[3],
                                         internal, daughters[pair[0]], daughters[pair[1]]) for internal in field_basis])
        total += (left.T*(mass_squared - invariants[pair]*kinetic).LUsolve(right))[0]
    return total


def record_route_one_kinetic_explicit(case, inverse_coupling, inverse_coupling_at, inverse_root,
                                      parent_root, flavour, invariants):
    """Route-1 kinetic contact of the record for cases 'A' and 'B' (derivative-free columns)."""
    eps = f_pi/f_a
    c = commutator
    s12, s13, s23 = (invariants[frozenset(pair)] for pair in ((1, 2), (1, 3), (2, 3)))
    y = inverse_coupling
    if case == 'A':
        b1, b2, b3 = (inverse_coupling_at[r] for r in (1, 2, 3))
        prefactor = -2*eps**4/(3*f_pi**2)*parent_root*inverse_root[1]*inverse_root[2]*inverse_root[3]
    else:
        b1, b2, b3 = inverse_coupling_at[1], flavour[2], flavour[3]
        prefactor = -eps**2/(6*f_pi**2)*parent_root*inverse_root[1]
    return sp.expand(prefactor*(sp.trace(c(y, b1)*c(b2, b3))*(s13 - s12)
                                + sp.trace(c(y, b2)*c(b1, b3))*(s23 - s12)
                                + sp.trace(c(y, b3)*c(b1, b2))*(s23 - s13)))


def record_formula_checks():
    """Compare the record's compact formulas with the vertex functions, exact rational inputs.

    The identities are algebraic in the column entries, so arbitrary rational
    values are used for M_{m_r^2}^{-1}(c_A_hat), m_a^2 - m_r^2 and 1/sqrt(N).
    Checks: quartic for 0, 1, 2, 3 daughters with ALP component (kinetic plus
    chi_+ parts against aligned_kinetic_contact + quartic_mass_vertex); cubic
    for parent-daughter-internal and internal-daughter-daughter legs, internal
    meson and ALP; amplitude of step 4 against the vertex-function assembly
    with exact congruent K', M'^2; route-1 kinetic contacts.
    """
    numbers = {f_pi: 1, f_a: 3, B_0: 2}
    rational = lambda *rows: sp.Matrix(rows)
    inverse_matrix = rational([1, 2, 0], [2, -1, 1], [0, 1, 3])       # Y = M_{m_a^2}^{-1}(C)
    hat_axial = -inverse_matrix
    inverse_axial = {'i': rational([2, 1, 0], [1, 0, 1], [0, 1, -1]), 'j': rational([0, 1, 1], [1, 3, 0], [1, 0, 2]),
                     'k': rational([1, 0, 2], [0, 2, 1], [2, 1, 0])}
    mass_gap = {'i': sp.Rational(3, 2), 'j': sp.Rational(-1, 3), 'k': sp.Rational(5, 4)}
    inverse_root = {'i': sp.Rational(1, 3), 'j': sp.Rational(2, 5), 'k': sp.Rational(3, 7)}
    parent = (sp.zeros(3), sp.Rational(5, 6))
    flavour = {'i': rational([1, 0, 0], [0, -1, 0], [0, 0, 0]), 'j': sp.sqrt(2)*sp.Matrix(matrix_unit(0, 2)),
               'k': sp.sqrt(2)*sp.Matrix(matrix_unit(2, 0))}
    scalar_coefficients = [rational([n, 1, 0], [1, 2*n, 1], [0, 1, -n]) for n in range(5)]      # S_n'
    pseudoscalar_coefficients = [sp.zeros(3)] + [rational([0, n, 1], [n, 0, 2], [1, 2, 1]) for n in range(1, 5)]
    plus = [4*B_0*scalar_coefficients[n]/f_a**n for n in range(5)]
    minus = [4*sp.I*B_0*pseudoscalar_coefficients[n]/f_a**n for n in range(5)]
    vs_commutator = sp.I*scalar_coefficients[1]
    hat_p = pseudoscalar_coefficients[1]
    masses = sp.symbols('_m2_0:4')
    labels = dict(zip('ijk', (1, 2, 3)))
    invariants = dict(zip(map(frozenset, ((1, 2), (1, 3), (2, 3))), sp.symbols('s_ij s_ik s_jk')))
    products = decay_scalar_products(masses[0], dict(zip((1, 2, 3), masses[1:])), invariants)
    parent_dot = {r: -products[(0, labels[r])] for r in 'ijk'}
    results = {}
    patterns = {'0': set(), '1': {'i'}, '2': {'i', 'j'}, '3': {'i', 'j', 'k'}}
    for name, neutral in patterns.items():
        daughters = {r: (daughter_column(inverse_axial[r], mass_gap[r], inverse_root[r]) if r in neutral
                         else (flavour[r], 0)) for r in 'ijk'}
        ordered = [daughters[r] for r in 'ijk']
        general = (aligned_kinetic_contact(inverse_matrix, parent[1], ordered, products)
                   + quartic_mass_vertex([parent] + ordered, plus[0], minus[1], plus[2], minus[3], plus[4]).doit())
        compact = compact_quartic(hat_axial, hat_p, scalar_coefficients[2], pseudoscalar_coefficients[3],
                                  scalar_coefficients[4], parent, daughters, parent_dot)
        results[f'compact quartic, {name} daughters with ALP component'] = sp.simplify(
            sp.expand((general - compact).subs(numbers))) == 0
    basis = [(sp.Matrix(generator), 0) for generator in hermitian_unit_basis()] + [(sp.zeros(3), 1)]
    neutral_column = daughter_column(inverse_axial['i'], mass_gap['i'], inverse_root['i'])
    other_neutral = daughter_column(inverse_axial['j'], mass_gap['j'], inverse_root['j'])
    triples = []
    for internal in basis:
        triples += [(parent, neutral_column, internal), (parent, (flavour['j'], 0), internal),
                    (internal, neutral_column, other_neutral), (internal, neutral_column, (flavour['k'], 0)),
                    (internal, (flavour['j'], 0), (flavour['k'], 0))]
    triples.append(((rational([1, 0, 1], [0, 2, 0], [1, 0, 0]), sp.Rational(1, 2)), neutral_column, other_neutral))
    cubic_ok = all(sp.simplify(sp.expand((cubic_vertex(*triple, plus[1], minus[2], plus[3]).doit()
                                          - compact_cubic(vs_commutator, pseudoscalar_coefficients[2],
                                                          scalar_coefficients[3], *triple)).subs(numbers))) == 0
                   for triple in triples)
    results['compact cubic: parent/internal legs, internal meson and ALP, 0-2 daughters with ALP component'] = cubic_ok
    # Amplitude of step 4 against the vertex-function assembly, exact rational K', M'^2.
    daughters = {r: daughter_column(inverse_axial[r], mass_gap[r], inverse_root[r]) if r != 'k'
                 else (flavour['k'], 0) for r in 'ijk'}
    kinetic_matrix = sp.eye(10)
    kinetic_matrix[0:9, 9] = sp.Matrix([sp.Rational(n, 7) for n in range(1, 10)])
    kinetic_matrix[9, 0:9] = kinetic_matrix[0:9, 9].T
    kinetic_matrix[9, 9] = 3
    mass_matrix = sp.Matrix(10, 10, lambda r, c: sp.Rational(1 + (r*c) % 5, 1 + (r + c) % 3)) + 20*sp.eye(10)
    mass_matrix = (mass_matrix + mass_matrix.T)/2
    numeric_invariants = {'ij': sp.Rational(7, 2), 'ik': sp.Rational(9, 4), 'jk': sp.Rational(11, 3)}
    numeric_products = decay_scalar_products(
        sp.Rational(23, 2), {1: 2, 2: 3, 3: 3},
        {frozenset({1, 2}): numeric_invariants['ij'], frozenset({1, 3}): numeric_invariants['ik'],
         frozenset({2, 3}): numeric_invariants['jk']})
    numeric_dot = {r: -numeric_products[(0, labels[r])] for r in 'ijk'}
    compact = compact_amplitude(hat_axial, hat_p, vs_commutator, scalar_coefficients, pseudoscalar_coefficients,
                                parent, daughters, numeric_dot, numeric_invariants, basis, kinetic_matrix, mass_matrix)
    ordered = [daughters[r] for r in 'ijk']
    general = (aligned_kinetic_contact(inverse_matrix, parent[1], ordered, numeric_products)
               + quartic_mass_vertex([parent] + ordered, plus[0], minus[1], plus[2], minus[3], plus[4]).doit())
    for attached, (first, second) in PAIRINGS:
        left = sp.Matrix([cubic_vertex(parent, ordered[attached - 1], d, plus[1], minus[2], plus[3]).doit()
                          for d in basis])
        right = sp.Matrix([cubic_vertex(d, ordered[first - 1], ordered[second - 1], plus[1], minus[2], plus[3]).doit()
                           for d in basis])
        pair = 'ijk'[first - 1] + 'ijk'[second - 1]
        general += (left.T*(mass_matrix - numeric_invariants[pair]*kinetic_matrix).LUsolve(right))[0]
    results['compact amplitude of step 4 = vertex-function assembly'] = sp.simplify(
        sp.expand((compact - general).subs(numbers))) == 0
    # Route-1 kinetic contact for three daughters with ALP component and for one.
    inverse_coupling_at = {r: rational([r, 1, 1], [1, 0, r], [1, r, 2]) for r in (1, 2, 3)}
    flavour_numbered = {1: flavour['i'], 2: flavour['j'], 3: flavour['k']}
    root_numbered = {1: inverse_root['i'], 2: inverse_root['j'], 3: inverse_root['k']}
    for case, neutral in (('A', {1, 2, 3}), ('B', {1})):
        df_parent = (2*f_pi/f_a*parent[1]*inverse_matrix, parent[1])
        df_daughters = [((2*f_pi/f_a*root_numbered[r]*inverse_coupling_at[r], root_numbered[r]) if r in neutral
                         else (flavour_numbered[r], 0)) for r in (1, 2, 3)]
        generic = quartic_kinetic_vertex([df_parent] + df_daughters, products).doit()
        explicit = record_route_one_kinetic_explicit(case, inverse_matrix, inverse_coupling_at, root_numbered,
                                                     parent[1], flavour_numbered, invariants)
        results[f'record route-1 kinetic contact, {"three daughters" if case == "A" else "one daughter"} with ALP component'] = (
            sp.simplify(sp.expand((generic - explicit).subs(numbers))) == 0)
    return results


def identical_daughter_residuals():
    """Record's M_aaa and M_aab(s_jk) against the generic small-mixing forms.

    Symbolic Y = M_{m_a^2}^{-1}(C), c_q^AP, diagonal m_q and arbitrary generators
    (species alpha, beta); the record puts the identical pair on j, k and beta on i.
    """
    quark_masses = sp.diag(*sp.symbols('m_u m_d m_s', positive=True))
    inverse_matrix = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f'y_{i}{j}'))
    axial_pseudoscalar = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f'c_{i}{j}'))
    hat_axial = -inverse_matrix
    hat_pseudoscalar = -axial_pseudoscalar - inverse_matrix
    hat_p = anticommutator(hat_pseudoscalar, quark_masses)
    alpha = sp.Matrix([[1, 2, 0], [0, -1, 1], [3, 0, 0]])
    beta = sp.Matrix([[0, 1, 1], [1, 0, 2], [0, 1, 1]])
    parent_mass_sq, alpha_mass_sq, beta_mass_sq = sp.symbols('m_a^2 m_alpha^2 m_beta^2', positive=True)
    s_ij, s_ik, s_jk = sp.symbols('s_ij s_ik s_jk')
    results = {}
    reference = small_mixing_ijk(hat_axial, hat_pseudoscalar, quark_masses, [alpha]*3, {1: 0, 2: 0, 3: 0})
    record = -B_0/(6*f_a*f_pi)*tensor_pseudoscalar(hat_p, alpha, alpha, alpha)
    results['M_aaa = -B_0/(6 f_a f_pi) T^P_aaa = small-mixing form'] = sp.expand(record - reference) == 0
    results['T^P_aaa = 6 <lambda^3 P_hat>, T^A_aa;a = 0'] = (
        sp.expand(tensor_pseudoscalar(hat_p, alpha, alpha, alpha) - 6*sp.trace(alpha**3*hat_p)) == 0
        and sp.expand(tensor_axial(hat_axial, alpha, alpha, alpha)) == 0)
    invariants = {frozenset({1, 2}): s_ij, frozenset({1, 3}): s_ik, frozenset({2, 3}): s_jk}
    products = decay_scalar_products(parent_mass_sq, {1: beta_mass_sq, 2: alpha_mass_sq, 3: alpha_mass_sq}, invariants)
    reference = small_mixing_ijk(hat_axial, hat_pseudoscalar, quark_masses, [beta, alpha, alpha],
                                 {r: -products[(0, r)] for r in (1, 2, 3)})
    record = (tensor_axial(hat_axial, alpha, alpha, beta)/(24*f_a*f_pi)*(parent_mass_sq + 3*beta_mass_sq - 3*s_jk)
              - B_0/(6*f_a*f_pi)*tensor_pseudoscalar(hat_p, alpha, alpha, beta))
    sum_rule = {s_ij: parent_mass_sq + beta_mass_sq + 2*alpha_mass_sq - s_ik - s_jk}
    results['M_aab(s_jk) = small-mixing form (identical j,k; species beta on i)'] = (
        sp.expand((record - reference).subs(sum_rule)) == 0)
    results['T^A_aa;b = 2<c_A_hat[a,[a,b]]>, T^P_aab = 2<(a^2 b + a b a + b a^2) P_hat>'] = (
        sp.expand(tensor_axial(hat_axial, alpha, alpha, beta)
                  - 2*sp.trace(hat_axial*commutator(alpha, commutator(alpha, beta)))) == 0
        and sp.expand(tensor_pseudoscalar(hat_p, alpha, alpha, beta)
                      - 2*sp.trace((alpha**2*beta + alpha*beta*alpha + beta*alpha**2)*hat_p)) == 0)
    return results


def aligned_weak_residuals():
    """Route 2 at O(f_pi/f_a): contact on unmixed daughters equals the small-mixing ijk form.

    Symbolic Y = M_{m_a^2}^{-1}(C), c_q^{AP} and diagonal m_q; daughters are
    arbitrary (not eigen-directions, not on shell): in pole-aligned
    coordinates no on-shell condition is needed.
    """
    quark_masses = sp.diag(*sp.symbols('m_u m_d m_s', positive=True))
    inverse_matrix = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f'y_{i}{j}'))
    axial_pseudoscalar = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f'c_{i}{j}'))
    zero = sp.zeros(3)
    plus = [4*B_0*quark_masses] + [zero]*4
    minus = [zero, -4*sp.I*B_0*anticommutator(axial_pseudoscalar, quark_masses)/f_a] + [zero]*3
    aligned_plus, aligned_minus = aligned_ccwz_derivatives(inverse_matrix, plus, minus)
    daughters = [sp.Matrix([[1, 2, 0], [0, -1, 1], [3, 0, 0]]), sp.Matrix([[0, 1, 1], [1, 0, 2], [0, 1, 1]]),
                 sp.Matrix([[2, 0, 1], [1, 1, 0], [0, 3, -1]])]
    parent_mass_sq = sp.Symbol('m_a^2', positive=True)
    masses = dict(zip((1, 2, 3), sp.symbols('m_1^2 m_2^2 m_3^2', positive=True)))
    invariants = dict(zip(map(frozenset, ((1, 2), (1, 3), (2, 3))), sp.symbols('s_12 s_13 s_23')))
    products = decay_scalar_products(parent_mass_sq, masses, invariants)
    aligned = aligned_contact(inverse_matrix, 1, [(matrix, 0) for matrix in daughters], products,
                              aligned_plus, aligned_minus)
    reference = small_mixing_ijk(-inverse_matrix, -axial_pseudoscalar - inverse_matrix, quark_masses,
                                 daughters, {r: -products[(0, r)] for r in (1, 2, 3)})
    return sp.expand(aligned - reference)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def derive():
    """Return the exact generic amplitude pieces, its weak-mixing term and the checks' inputs."""
    kinetic_summed, kinetic_reduced, _ = kinetic_contact_reduction()
    index_level = index_level_amplitude()
    epsilon = index_level['epsilon']
    return {
        'kinetic_contact_summed': kinetic_summed,
        'kinetic_contact_reduced': kinetic_reduced,
        'exact_amplitude': index_level['amplitude'],
        'weak_mixing_leading': epsilon*index_level['series'][1],
        'index_level': index_level,
        'symmetry_factors': symmetry_factor_multiplicities(),
    }


def run_checks(verbose=True):
    """Run all deterministic checks; return a dict of name -> bool."""
    results = {}
    derived = derive()
    results['kinetic S4 sum equals reduced pair-invariant form'] = (
        sp.expand(derived['kinetic_contact_summed'] - derived['kinetic_contact_reduced']) == 0)
    index_level = derived['index_level']
    results['no O(epsilon^0) or O(epsilon^2) term'] = (index_level['series'][0] == 0
                                                       and index_level['series'][2] == 0)
    results['exchanges start beyond O(epsilon^2)'] = all(index_level['exchange_series'][n] == 0 for n in range(3))
    results['O(epsilon) term = parent ALP leg + parent meson leg on unmixed daughters'] = (
        sp.expand(index_level['series'][1] - index_level['expected_leading']) == 0)
    factors = derived['symmetry_factors']
    results['vertex multiplicities 3!, 2!, 1 and S = 3!, 2!, 1'] = (
        factors == {'iii': (6, 6), 'iij': (2, 2), 'ijk': (1, 1)})
    for channel, (residual, exchanges_nonzero) in three_pion_substitution_residuals().items():
        results[f'generic = record M_{channel} (exact in s,t,u)'] = residual == 0
        expected = {pairing: not (channel == '0+-' and pairing == (1, 2, 3))
                    for pairing in exchanges_nonzero}
        results[f'exchange pattern in M_{channel}: all nonzero, except (a,0)|(+,-) for 0+-'] = exchanges_nonzero == expected
    for name, residual in weak_mixing_residuals().items():
        if 'nonzero' in name:
            results[f'weak limit: {name}'] = residual is True
        elif 'explicit' in name:
            results[f'weak limit {name} = V4 on unmixed daughters'] = residual == 0
        else:
            results[f'weak limit = small-mixing form {name}'] = residual == 0
    for name, passed in transformation_checks().items():
        results[f'alignment: {name}'] = passed
    for name, passed in record_formula_checks().items():
        results[f'record formula: {name}'] = passed
    for name, passed in identical_daughter_residuals().items():
        results[f'record identical daughters: {name}'] = passed
    results['route 2 at O(f_pi/f_a) = small-mixing ijk form for any daughter generators'] = (
        aligned_weak_residuals() == 0)
    on_shell = numerical_route_comparison()
    off_shell = numerical_route_comparison(off_shell_shift=0.37)
    results['route 1 = route 2 on shell (relative 1e-10)'] = (
        abs(on_shell[0] - on_shell[1]) < 1e-10*abs(on_shell[0]))
    results['contact terms of the two routes agree on shell, differ off shell'] = (
        abs(on_shell[2]) < 1e-10*abs(on_shell[0]) and abs(off_shell[2]) > 1e-6)
    for name, passed in on_shell[3].items():
        results[f'invariant columns: {name}'] = passed
    results = {name: bool(passed) for name, passed in results.items()}
    if verbose:
        for name, passed in results.items():
            print(f"{'PASS' if passed else 'FAIL'}: {name}")
        print(f'route 1 = {on_shell[0].real:.15g}, route 2 = {on_shell[1].real:.15g} (on shell)')
        print(f'contact difference off shell (daughter 1 mass^2 shifted by 0.37): {off_shell[2].real:.6g}')
    return results


def main():
    derived = derive()
    print('Pairings: (a,1)|(2,3) at s_23, (a,2)|(1,3) at s_13, (a,3)|(1,2) at s_12.')
    print('Exact amplitude: M = V4[a,1,2,3] + sum_pairings V3[a,r,I](M^2-s_uv K)^{-1}_IJ V3[J,u,v].')
    print('\nV4_kin as the S_4 sum, after the decay scalar products:')
    print(sp.latex(sp.collect(derived['kinetic_contact_summed'], [sp.Symbol(r'T_{a1|23}'),
                                                                   sp.Symbol(r'T_{a2|13}'),
                                                                   sp.Symbol(r'T_{a3|12}')])))
    print('Reduced form (masses cancel):')
    print(sp.latex(derived['kinetic_contact_reduced']))
    print('\nIndex-level exact amplitude (two meson fields and one ALP, generic tensors):')
    print(f"  number of terms: {len(sp.Add.make_args(derived['exact_amplitude']))}")
    print('Leading weak-mixing term:')
    print(sp.latex(derived['weak_mixing_leading']))
    print('\nVertex multiplicity and phase-space factor S per channel:',
          derived['symmetry_factors'])
    print('\nRoute 2 (pole-aligned), sources from the general-source record Taylor law with')
    print("g_R = exp(i x c_A_hat), g_L = g_R^dagger, c_A_hat = -M_{m_a^2}^{-1}(C), starting from widehat M:")
    print("  M_n' = sum_{r+s+t=n} n!/(r!s!t!) (i c_A_hat)^r widehat M_t (i c_A_hat)^s,  S', P' = Herm./anti-Herm. parts")
    print("  d_a^n chi_+|_0 = 4 B_0 S_n'/f_a^n,  d_a^n chi_-|_0 = 4 i B_0 P_n'/f_a^n")
    print("  S_1' = widehat S_1 = -i[c_VS, m_q],  P_1' = widehat P_1 + {c_A_hat, m_q} = {c_P_hat, m_q}")
    print("  S_2' = widehat S_2 - 2{c_A_hat, widehat P_1} - {c_A_hat,{c_A_hat, m_q}},  P_2' = widehat P_2 + 2{c_A_hat, widehat S_1}")
    print()
    results = run_checks(verbose=True)
    print(f"\n{sum(results.values())}/{len(results)} checks passed")
    return results


class ThreePseudoscalarSelfChecks(unittest.TestCase):
    def test_all_checks(self):
        results = run_checks(verbose=False)
        self.assertTrue(all(results.values()), [name for name, ok in results.items() if not ok])


if __name__ == '__main__':
    if '--test' in sys.argv:
        unittest.main(argv=[sys.argv[0]])
    else:
        main()
