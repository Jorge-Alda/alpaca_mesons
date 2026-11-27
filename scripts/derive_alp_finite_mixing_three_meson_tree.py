#!/usr/bin/env python3
"""Exact tree amplitude A -> P_{n1} P_{n2} P_{n3} of the ALP-like pole at finite epsilon.

Run:   .venv/bin/python scripts/derive_alp_finite_mixing_three_meson_tree.py
Tests: .venv/bin/python -m unittest discover -s scripts -p test_derive_alp_finite_mixing_three_meson_tree.py

Minimal use (run from scripts/; exact rational inputs of `verification_inputs`):
    import sympy as sp
    from derive_alp_finite_mixing_three_meson_tree import compare_amplitude, verification_inputs
    result = compare_amplitude(verification_inputs(), sp.Integer(2), sp.Rational(17, 300),
                               sp.Rational(1, 2), sp.Rational(3, 5))
    print(result.eps_squared, sp.N(result.closed), sp.N(result.solver))

Conventions (Alda:2026nrd): natural units, masses in GeV, squared masses and
Mandelstam invariants in GeV^2; U = exp(i Phi/f_pi), Phi = P_i lambda_i,
tr(lambda_i lambda_j) = 2 delta_ij, N_f = 3, x = a/f_a, eps = f_pi/f_a.
Dressed frame (memory/alp/chiral-transformations/source-invariants/dressed-sources/finite-mixing/alp-dressed-frame-derivative-coupling-removal.md):
K = 1, no ALP currents, mass source chi = 2 B_0 hat M(x),
hat M(x) = sum_n x^n/n! (hat S_n + i hat P_n), topological coefficient
hat Theta(x) = sum_n x^n/n! hat Theta_n, LO Lagrangian
    L = f_pi^2/4 <u_mu u^mu> + f_pi^2/4 <chi_+>
        - f_pi^2 m_0^2/(4 N_f) (<Phi>/f_pi + 2 hat Theta(x))^2
        + (d a)^2/2 - m_{a,0}^2 a^2/2 .
Vertices are derivatives of L (overall factor i omitted), all momenta
incoming; the four-leg reduced amplitude is contact minus the three exchanges
-(V3 D(q^2)^{-1} V3) (memory/alp/decay-amplitudes/vertices-and-kinematics/tree-vertices/finite-mixing/alp-lo-pole-projected-tree-amplitudes.md),
no symmetry factor for identical daughters.

Inputs (decisions of the OpenSpec task 2.1): hat S_0 = m_q diagonal with three
different entries, hat P_0 = 0, hat Theta_0 = 0, <hat S_1> = 0 (no ALP tadpole,
an assumption), every other invariant a general Hermitian matrix commuting with
the quark charge (diagonal plus d-s entries); V_ext = m_{a,0}^2 a^2/2 only.

Closed formula ingredients (matrix notation):
    C = B_0 hat P_1 - (m_0^2/N_f) hat Theta_1 1,
    M_s(X) = 1/2 {chi_0, X} + (m_0^2/N_f) <X> 1 - s X,  chi_0 = 2 B_0 m_q,
    D_a(s) = s - M^2_aa + 2 eps^2 <C M_s^{-1}(C)>,
    M^2_aa = m_{a,0}^2 - eps^2 B_0 <hat S_2> + 2 eps^2 (m_0^2/N_f) hat Theta_1^2,
    N(s) = 1 + 2 eps^2 <M_s^{-1}(C)^2> = D_a'(s).
Pole vectors: A has meson matrix 2 eps Y/sqrt(N(m_A^2)), Y = M^{-1}_{m_A^2}(C),
ALP component 1/sqrt(N(m_A^2)); a mixed meson-like pole n has
sign_n (2 eps M^{-1}_{m_n^2}(C), 1)/sqrt(N(m_n^2)) with sign_n fixed so that the
meson matrix has positive trace product with the unmixed matrix (continuity);
poles without ALP component use their unmixed matrix.

Expansion convention (stated choice): fixed invariants, fixed m_A^2 (m_{a,0}^2
adjusts), fixed s_12 and s_23; s_13 = m_A^2 + sum_i m_{n_i}^2 - s_12 - s_23
follows the exact daughter poles.
"""
from __future__ import annotations

import itertools
from dataclasses import dataclass

import sympy as sp

from derive_alp_dressed_frame_unitarity import on_shell_gram
from derive_alp_lo_joint_vertices import quadratic_blocks, vertex
from derive_alp_pole_projected_amplitudes import projected_amplitude
from solve_alp_meson_spectrum import PoleSpectrum, solve_spectrum_exact
from symbolic_utils import u3_generators

FLAVOURS = 3
GENERATORS = [sp.Matrix(g) for g in u3_generators()]


def unit_matrix(row, column):
    """3x3 matrix with a single 1 at (row, column)."""
    matrix = sp.zeros(3)
    matrix[row, column] = 1
    return matrix


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------


def verification_inputs(down_strange=False):
    """Exact rational invariants (GeV units) used only for the checks.

    chi_0 = 2 B_0 m_q = diag(1, 6, 27)/50 GeV^2 and m_0^2 = 12/25 GeV^2 give
    rational unmixed neutral masses (found by a small search). All CP-odd
    invariants (hat S_1, hat S_3, hat P_2, hat Theta_2) are switched on, so every
    contact class and every exchange contributes. With down_strange=True the
    invariants hat S_1..hat S_4, hat P_1..hat P_3 also get complex d-s entries.
    """
    R = sp.Rational
    B_0 = sp.Integer(2)
    quark_mass = sp.diag(R(1, 50), R(6, 50), R(27, 50)) / (2 * B_0)
    S = [quark_mass,
         sp.diag(R(1, 3), -R(1, 5), -R(2, 15)),
         sp.diag(R(1, 2), R(-1, 4), R(1, 7)),
         sp.diag(R(-1, 3), R(2, 5), R(1, 6)),
         sp.diag(R(1, 4), R(1, 9), R(-1, 2))]
    P = [sp.zeros(3),
         sp.diag(R(1, 5), R(-1, 3), R(1, 2)),
         sp.diag(R(1, 6), R(1, 4), R(-1, 5)),
         sp.diag(R(-1, 2), R(1, 3), R(1, 8)),
         sp.zeros(3)]
    if down_strange:
        for matrix, entry in ((S[1], R(1, 4) + sp.I * R(1, 6)), (S[2], R(-1, 5) + sp.I * R(1, 3)),
                              (S[3], R(1, 7) - sp.I * R(1, 4)), (S[4], R(1, 3)),
                              (P[1], R(1, 3) - sp.I * R(1, 5)), (P[2], R(1, 4) + sp.I * R(1, 2)),
                              (P[3], R(-1, 6) + sp.I * R(1, 7))):
            matrix[1, 2] = entry
            matrix[2, 1] = sp.conjugate(entry)
    theta = [0, R(1, 3), R(-1, 4), R(1, 5), R(1, 6)]
    return dict(B_0=B_0, f_pi=R(1, 10), anomaly_mass_squared=R(12, 25), S=S, P=P, theta=theta)


ALP_POLE_SQUARED_MASS = sp.Integer(2)
# shifts m_pi^2 - m_{pi,0}^2 (GeV^2) of the three verification points; they fix eps^2
OFFSETS = (-sp.Rational(1, 1000), -sp.Rational(1, 300), -sp.Rational(1, 100))
KINEMATIC_POINT = (sp.Rational(1, 2), sp.Rational(3, 5))  # (s_12, s_23) in GeV^2


# ---------------------------------------------------------------------------
# Flavour-space operator, its inverse and the ALP-like pole
# ---------------------------------------------------------------------------


def chi_values(invariants):
    """Diagonal entries chi_i of chi_0 = 2 B_0 m_q."""
    return [2 * invariants['B_0'] * invariants['S'][0][i, i] for i in range(3)]


def pole_operator(flavour_matrix, invariants, squared_momentum):
    """M_s(X) = 1/2 {chi_0, X} + (m_0^2/N_f) <X> 1 - s X."""
    chi_0 = 2 * invariants['B_0'] * invariants['S'][0]
    return ((chi_0 * flavour_matrix + flavour_matrix * chi_0) / 2
            + invariants['anomaly_mass_squared'] / FLAVOURS * flavour_matrix.trace() * sp.eye(3)
            - squared_momentum * flavour_matrix)


def inverse_pole_operator(flavour_matrix, invariants, squared_momentum):
    """M_s^{-1}(Z) for diagonal chi_0 (explicit inverse of the operator-inverse record)."""
    chi = chi_values(invariants)
    topological = invariants['anomaly_mass_squared'] / FLAVOURS
    s = squared_momentum
    trace = (sum(flavour_matrix[k, k] / (chi[k] - s) for k in range(3))
             / (1 + topological * sum(1 / (chi[k] - s) for k in range(3))))
    result = sp.zeros(3)
    for i in range(3):
        for j in range(3):
            if i == j:
                result[i, i] = (flavour_matrix[i, i] - topological * trace) / (chi[i] - s)
            else:
                result[i, j] = flavour_matrix[i, j] / ((chi[i] + chi[j]) / 2 - s)
    return result


def matrix_C(invariants):
    """C = B_0 hat P_1 - (m_0^2/N_f) hat Theta_1 1."""
    return (invariants['B_0'] * invariants['P'][1]
            - invariants['anomaly_mass_squared'] / FLAVOURS * invariants['theta'][1] * sp.eye(3))


def alp_mass_entry(invariants, eps, bare_alp_mass_squared):
    """Dressed-frame M^2_aa = m_{a,0}^2 - eps^2 B_0 <hat S_2> + 2 eps^2 (m_0^2/N_f) hat Theta_1^2."""
    return (bare_alp_mass_squared - eps**2 * invariants['B_0'] * invariants['S'][2].trace()
            + 2 * eps**2 * invariants['anomaly_mass_squared'] / FLAVOURS * invariants['theta'][1]**2)


def alp_inverse_propagator(invariants, eps, bare_alp_mass_squared, squared_momentum):
    """D_a(s) = s - M^2_aa + 2 eps^2 <C M_s^{-1}(C)> (Schur complement of D(s))."""
    source = matrix_C(invariants)
    return (squared_momentum - alp_mass_entry(invariants, eps, bare_alp_mass_squared)
            + 2 * eps**2 * (source * inverse_pole_operator(source, invariants, squared_momentum)).trace())


def normalization(invariants, eps, squared_momentum):
    """N(s) = 1 + 2 eps^2 <M_s^{-1}(C)^2>."""
    overlap = inverse_pole_operator(matrix_C(invariants), invariants, squared_momentum)
    return 1 + 2 * eps**2 * (overlap * overlap).trace()


def bare_alp_mass_for_pole(invariants, eps, alp_pole_squared_mass):
    """m_{a,0}^2 such that D_a(m_A^2) = 0 (D_a is linear in m_{a,0}^2)."""
    bare = sp.Symbol('bare_alp_mass_squared')
    return sp.solve(alp_inverse_propagator(invariants, eps, bare, alp_pole_squared_mass), bare)[0]


def tuned_point(invariants, alp_pole_squared_mass, daughter_squared_mass):
    """(eps^2, m_{a,0}^2) with D_a(m_A^2) = D_a(m_P^2) = 0 (both rational).

    D_a(m_A^2) - D_a(m_P^2) does not contain M^2_aa and is linear in eps^2.
    """
    eps_squared = sp.Symbol('eps_squared', positive=True)
    eps = sp.sqrt(eps_squared)
    difference = (alp_inverse_propagator(invariants, eps, 0, alp_pole_squared_mass)
                  - alp_inverse_propagator(invariants, eps, 0, daughter_squared_mass))
    value = sp.solve(sp.expand(difference), eps_squared)[0]
    return value, sp.simplify(bare_alp_mass_for_pole(invariants, sp.sqrt(value), alp_pole_squared_mass))


# ---------------------------------------------------------------------------
# Unmixed meson matrices and pole legs; a leg is (meson matrix, ALP component)
# ---------------------------------------------------------------------------


def unmixed_neutral_squared_masses(invariants):
    """Roots of 1 + (m_0^2/N_f) sum_i 1/(chi_i - s) = 0 (diagonal neutral mesons, eps = 0)."""
    s = sp.Symbol('s')
    chi = chi_values(invariants)
    polynomial = sp.expand(sp.prod([c - s for c in chi]) + invariants['anomaly_mass_squared'] / FLAVOURS * sum(
        sp.prod([chi[i] - s for i in range(3) if i != j]) for j in range(3)))
    return sorted(sp.roots(sp.Poly(polynomial, s)), key=lambda r: sp.N(r))


def unmixed_neutral_matrix(invariants, unmixed_squared_mass):
    """Unmixed diagonal meson matrix with <X^2> = 2, X_ii proportional to 1/(m_{n,0}^2 - chi_i).

    Sign (a choice, identical to the pole solver's rule at eps = 0 in the basis
    (lambda_3, lambda_8, lambda_0)): the component of largest modulus is positive.
    """
    chi = chi_values(invariants)
    matrix = sp.diag(*[1 / (unmixed_squared_mass - c) for c in chi])
    matrix = matrix * sp.sqrt(2 / (matrix * matrix).trace())
    components = [(matrix * GENERATORS[k]).trace() / 2 for k in (3, 8, 0)]
    largest = max(components, key=lambda c: abs(sp.N(c)))
    return (-matrix if sp.N(largest) < 0 else matrix).applyfunc(sp.radsimp)


def alp_like_leg(invariants, eps, alp_pole_squared_mass, normalized=True):
    """(2 eps Y, 1)/sqrt(N(m_A^2)), Y = M^{-1}_{m_A^2}(C); without 1/sqrt(N) if normalized=False."""
    overlap = inverse_pole_operator(matrix_C(invariants), invariants, alp_pole_squared_mass)
    norm = sp.sqrt(normalization(invariants, eps, alp_pole_squared_mass)) if normalized else sp.Integer(1)
    return 2 * eps * overlap / norm, 1 / norm


def mixed_meson_leg(invariants, eps, pole_squared_mass, unmixed_matrix, normalized=True):
    """sign (2 eps M^{-1}_{m_n^2}(C), 1)/sqrt(N(m_n^2)), sign by continuity with the unmixed matrix."""
    overlap = inverse_pole_operator(matrix_C(invariants), invariants, pole_squared_mass)
    norm = sp.sqrt(normalization(invariants, eps, pole_squared_mass)) if normalized else sp.Integer(1)
    sign = 1 if sp.N(eps * (overlap * unmixed_matrix).trace(), 30) > 0 else -1
    return sign * 2 * eps * overlap / norm, sign / norm


# ---------------------------------------------------------------------------
# Vertices of the dressed Lagrangian, as multilinear polarizations
# ---------------------------------------------------------------------------


def commutator(left, right):
    return left * right - right * left


def anticommutator(left, right):
    return left * right + right * left


def cubic_current(invariants, eps, first_leg, second_leg):
    """(J, J_a) with V3(first, second, (X, x)) = <J X> + J_a x.

    From L_3 = -(eps B_0/(2 f_pi)) a <hat S_1 Phi^2>
               + (eps^2/(2 f_pi)) a^2 <(B_0 hat P_2 - (m_0^2/N_f) hat Theta_2) Phi>
               + (eps^3/(6 f_pi)) (B_0 <hat S_3> - 6 (m_0^2/N_f) hat Theta_1 hat Theta_2) a^3.
    No PPP vertex exists (hat P_0 = 0, the kinetic term is even in Phi).
    """
    B_0, f_pi, m0 = invariants['B_0'], invariants['f_pi'], invariants['anomaly_mass_squared']
    S, P, theta = invariants['S'], invariants['P'], invariants['theta']
    (phi_1, alp_1), (phi_2, alp_2) = first_leg, second_leg
    second_order = B_0 * P[2] - m0 / FLAVOURS * theta[2] * sp.eye(3)
    third_order = B_0 * S[3].trace() - 6 * m0 / FLAVOURS * theta[1] * theta[2]
    current = (-eps * B_0 / (2 * f_pi) * (alp_1 * anticommutator(S[1], phi_2) + alp_2 * anticommutator(S[1], phi_1))
               + eps**2 / f_pi * alp_1 * alp_2 * second_order)
    alp_current = (-eps * B_0 / (2 * f_pi) * (S[1] * anticommutator(phi_1, phi_2)).trace()
                   + eps**2 / f_pi * (alp_1 * (second_order * phi_2).trace() + alp_2 * (second_order * phi_1).trace())
                   + eps**3 / f_pi * third_order * alp_1 * alp_2)
    return current, alp_current


def cubic_vertex(invariants, eps, legs):
    """Three-leg vertex contracted with three legs."""
    current, alp_current = cubic_current(invariants, eps, legs[0], legs[1])
    return sp.expand((current * legs[2][0]).trace() + alp_current * legs[2][1])


def quartic_contact_classes(invariants, eps, legs, products):
    """Four-leg contact vertex split by class; products[i][j] = p_i . p_j (incoming).

    kinetic : -(1/(24 f^2)) sum over the 3 partitions {i,j}|{k,l} of
              (p_i.p_j + p_k.p_l)(<[Phi_k,Phi_i][Phi_l,Phi_j]> + <[Phi_l,Phi_i][Phi_k,Phi_j]>)
    PPPP    : (B_0/(24 f^2)) sum_{sigma in S_4} <hat S_0 Phi_s1 Phi_s2 Phi_s3 Phi_s4>
    PPPa    : -(eps B_0/(6 f^2)) sum_l alpha_l sum_{S_3} <hat P_1 Phi Phi Phi>
    PPaa    : -(eps^2 B_0/(2 f^2)) sum_{l<m} alpha_l alpha_m <hat S_2 {Phi_k, Phi_n}>
    Paaa    : (eps^3/f^2) sum_l prod_{m!=l} alpha_m <(B_0 hat P_3 - (m_0^2/N_f) hat Theta_3) Phi_l>
    aaaa    : (eps^4/f^2)(B_0 <hat S_4> - (m_0^2/N_f)(8 hat Theta_1 hat Theta_3 + 6 hat Theta_2^2)) prod alpha
    """
    B_0, f_pi, m0 = invariants['B_0'], invariants['f_pi'], invariants['anomaly_mass_squared']
    S, P, theta = invariants['S'], invariants['P'], invariants['theta']
    phi = [leg[0] for leg in legs]
    alp = [leg[1] for leg in legs]
    kinetic = 0
    for (i, j), (k, l) in (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))):
        structure = ((commutator(phi[k], phi[i]) * commutator(phi[l], phi[j])).trace()
                     + (commutator(phi[l], phi[i]) * commutator(phi[k], phi[j])).trace())
        kinetic += (products[i][j] + products[k][l]) * structure
    kinetic = -kinetic / (24 * f_pi**2)
    mass = B_0 / (24 * f_pi**2) * sum(
        (S[0] * phi[a] * phi[b] * phi[c] * phi[d]).trace() for a, b, c, d in itertools.permutations(range(4)))
    first = -eps * B_0 / (6 * f_pi**2) * sum(
        alp[l] * sum((P[1] * phi[a] * phi[b] * phi[c]).trace()
                     for a, b, c in itertools.permutations([m for m in range(4) if m != l]))
        for l in range(4))
    second = -eps**2 * B_0 / (2 * f_pi**2) * sum(
        alp[l] * alp[m] * (S[2] * anticommutator(*[phi[k] for k in range(4) if k not in (l, m)])).trace()
        for l, m in itertools.combinations(range(4), 2))
    third_matrix = B_0 * P[3] - m0 / FLAVOURS * theta[3] * sp.eye(3)
    third = eps**3 / f_pi**2 * sum(
        sp.prod([alp[m] for m in range(4) if m != l]) * (third_matrix * phi[l]).trace() for l in range(4))
    fourth = (eps**4 / f_pi**2 * (B_0 * S[4].trace() - m0 / FLAVOURS * (8 * theta[1] * theta[3] + 6 * theta[2]**2))
              * sp.prod(alp))
    return dict(kinetic=kinetic, PPPP=mass, PPPa=first, PPaa=second, Paaa=third, aaaa=fourth)


def exchange(invariants, eps, bare_alp_mass_squared, left_legs, right_legs, squared_momentum):
    """V3(left, .) D(q^2)^{-1} V3(., right) in matrix form (enters the amplitude with a minus sign).

    = -2 <J_L M_q^{-1}(J_R)>
      + (2 eps <J_L M_q^{-1}(C)> + J_{L,a})(2 eps <J_R M_q^{-1}(C)> + J_{R,a}) / D_a(q^2).
    """
    left, left_alp = cubic_current(invariants, eps, *left_legs)
    right, right_alp = cubic_current(invariants, eps, *right_legs)
    source = matrix_C(invariants)
    s = squared_momentum
    meson_block = -2 * (left * inverse_pole_operator(right, invariants, s)).trace()
    overlap = inverse_pole_operator(source, invariants, s)
    left_direction = 2 * eps * (left * overlap).trace() + left_alp
    right_direction = 2 * eps * (right * overlap).trace() + right_alp
    return meson_block + left_direction * right_direction / alp_inverse_propagator(
        invariants, eps, bare_alp_mass_squared, s)


def decay_products(masses, mandelstam):
    """p_i.p_j for incoming (A, n1, n2, n3); mandelstam[(i, j)] = (p_{n_i}+p_{n_j})^2, i<j in 1..3."""
    products = [[0] * 4 for _ in range(4)]
    for i in range(4):
        products[i][i] = masses[i]
    for i, j in itertools.combinations((1, 2, 3), 2):
        products[i][j] = products[j][i] = (mandelstam[(i, j)] - masses[i] - masses[j]) / 2
    for i in (1, 2, 3):
        j, k = [m for m in (1, 2, 3) if m != i]
        products[0][i] = products[i][0] = (mandelstam[(j, k)] - masses[0] - masses[i]) / 2
    return products


def three_meson_amplitude(invariants, eps, bare_alp_mass_squared, legs, masses, mandelstam):
    """Closed-formula M(A -> n1 n2 n3) = contact - sum of the three exchanges; returns the pieces.

    legs = [A, n1, n2, n3] as (meson matrix, ALP component); masses are the pole
    masses squared; mandelstam[(i, j)] = (p_{n_i} + p_{n_j})^2.
    """
    pieces = quartic_contact_classes(invariants, eps, legs, decay_products(masses, mandelstam))
    for i in (1, 2, 3):
        j, k = [m for m in (1, 2, 3) if m != i]
        pieces[f'exchange A n{i} | n{j} n{k}'] = -exchange(
            invariants, eps, bare_alp_mass_squared, (legs[0], legs[i]), (legs[j], legs[k]), mandelstam[(j, k)])
    pieces['total'] = sum(pieces.values())
    return pieces


def leading_coefficient(invariants, alp_pole_squared_mass, unmixed_legs, unmixed_masses, mandelstam,
                        overlap_argument=None):
    """Coefficient of eps^1: four-meson contact with A -> (2Y, 0) plus the hat P_1 contact with A -> (0, 1).

    unmixed_legs are the unmixed daughter matrices; Y is evaluated at
    overlap_argument (default m_A^2).
    """
    argument = alp_pole_squared_mass if overlap_argument is None else overlap_argument
    overlap = inverse_pole_operator(matrix_C(invariants), invariants, argument)
    masses = [alp_pole_squared_mass] + list(unmixed_masses)
    products = decay_products(masses, mandelstam)
    daughters = [(matrix, 0) for matrix in unmixed_legs]
    meson = quartic_contact_classes(invariants, 1, [(2 * overlap, 0)] + daughters, products)
    alp = quartic_contact_classes(invariants, 1, [(sp.zeros(3), 1)] + daughters, products)
    return meson['kinetic'] + meson['PPPP'] + alp['PPPa']


# ---------------------------------------------------------------------------
# Independent evaluation: quadratic blocks, pole solver, CCWZ vertices, D^{-1}
# ---------------------------------------------------------------------------


def dressed_source_arrays(invariants):
    """Taylor arrays (chi_n = 2 B_0 (hat S_n + i hat P_n), C_R = C_L = 0, hat Theta_n), n <= 4."""
    zero = sp.zeros(3)
    chi = [2 * invariants['B_0'] * (invariants['S'][n] + sp.I * invariants['P'][n]) for n in range(5)]
    return chi, [zero] * 5, [zero] * 5, list(invariants['theta'])


class ExactVertexContraction:
    """Exact multilinear contraction of CCWZ vertices (`vertex`) of the basis fields.

    Basis-field vertices are computed once at f_a = 1 and scaled by f_a^{-n_a};
    the call signature matches the vertex functions of projected_amplitude.
    """

    def __init__(self, meson_basis, invariants):
        self.basis = meson_basis
        self.size = len(meson_basis) + 1
        self.sources = dressed_source_arrays(invariants)
        self.f_pi = invariants['f_pi']
        self.anomaly = invariants['anomaly_mass_squared']
        self.cache = {}

    def basis_terms(self, labels):
        if labels not in self.cache:
            symbols = sp.symbols(f'_r0:{len(labels)}', real=True)
            zero = sp.zeros(3)
            expression = vertex([self.basis[i] if i < self.size - 1 else zero for i in labels],
                                [1 if i == self.size - 1 else 0 for i in labels], symbols, *self.sources,
                                self.f_pi, sp.Integer(1), self.anomaly, 0)
            terms = []
            if expression != 0:
                for powers, coefficient in sp.Poly(expression, *symbols).terms():
                    terms.append((coefficient, tuple(i for i, p in enumerate(powers) for _ in range(p))))
            self.cache[labels] = terms
        return self.cache[labels]

    def function(self, f_a):
        def contraction(directions, momenta, gram):
            count = len(directions)
            products = [[(momenta[i].T * gram * momenta[j])[0] for j in range(count)] for i in range(count)]
            total = sp.S.Zero
            for labels in itertools.product(range(self.size), repeat=count):
                weight = sp.prod([directions[k][labels[k]] for k in range(count)])
                if weight == 0:
                    continue
                order = sorted(range(count), key=lambda k: labels[k])
                alp_legs = sum(1 for label in labels if label == self.size - 1)
                for coefficient, legs in self.basis_terms(tuple(labels[k] for k in order)):
                    factor = products[order[legs[0]]][order[legs[1]]] if legs else 1
                    total += weight * coefficient * factor / f_a**alp_legs
            return sp.expand(total)
        return contraction


def solver_theory(invariants, eps, bare_alp_mass_squared, meson_basis):
    """(K, M^2) from quadratic_blocks of the dressed sources at f_a = f_pi/eps."""
    sources = dressed_source_arrays(invariants)
    return quadratic_blocks(meson_basis, sources[0][:3], sources[1][:3], sources[2][:3], sources[3][:3],
                            invariants['f_pi'], invariants['f_pi'] / eps, invariants['anomaly_mass_squared'],
                            bare_alp_mass_squared)


def null_space_pole_vector(kinetic, mass_squared, pole):
    """K-normalized null vector of D(m^2) = m^2 K - M^2 at a simple pole (the solver's algorithm)."""
    null = (pole * kinetic - mass_squared).nullspace(simplify=True)
    if len(null) != 1:
        raise ValueError(f'pole {pole} is not simple')
    vector = null[0]
    return (vector / sp.sqrt(sp.simplify((vector.T * kinetic * vector)[0]))).applyfunc(sp.radsimp)


def leg_from_vector(vector, meson_basis):
    """(sum_k v_k basis_k, v_a)."""
    return sum((vector[k] * meson_basis[k] for k in range(len(meson_basis))), sp.zeros(3)), vector[-1]


def align_with_unmixed(vector, meson_basis, unmixed_matrix):
    """Flip the vector so that its meson matrix has positive trace product with the unmixed matrix."""
    matrix, _ = leg_from_vector(vector, meson_basis)
    return -vector if sp.N((matrix * unmixed_matrix).trace(), 30) < 0 else vector


def solver_amplitude(kinetic, mass_squared, vectors, masses, mandelstam, contraction):
    """Contact plus exchanges with the direct inverse D(q^2)^{-1} (projected_amplitude)."""
    spectrum = PoleSpectrum(field_labels=tuple(f'F{i}' for i in range(kinetic.rows)), kinetic=kinetic,
                            mass_squared=mass_squared, pole_squared_masses=[], pole_vectors=sp.zeros(kinetic.rows, 1),
                            residues=[], exact=True)
    gram = on_shell_gram(masses, mandelstam[(2, 3)], mandelstam[(1, 3)])
    return projected_amplitude(spectrum, [0, 1, 2, 3], gram, contraction, pole_vectors=vectors)


def is_exactly_zero(expression):
    """Exact zero test for expressions with nested square roots."""
    value = sp.radsimp(sp.expand(expression))
    if value == 0:
        return True
    return sp.simplify(value) == 0


NEUTRAL_BASIS = [unit_matrix(0, 0), unit_matrix(1, 1), unit_matrix(2, 2)]
CHARGED_BASIS = NEUTRAL_BASIS + [GENERATORS[1], GENERATORS[2]]
DOWN_STRANGE_BASIS = NEUTRAL_BASIS + [GENERATORS[6], GENERATORS[7]]


@dataclass
class ComparisonResult:
    """Closed formula and independent evaluation at one tuned point.

    rational_form = M with unnormalized legs (2 eps M^{-1}(C), 1), which is rational
    at a tuned point (only even powers of eps = sqrt(eps^2) survive, because the
    meson part of every mixed leg carries one factor eps); the normalized amplitude
    is closed = rational_form / sqrt(normalization_product).
    """
    eps_squared: object
    bare_alp_mass_squared: object
    daughter_squared_mass: object
    closed: object
    solver: object
    pieces: dict
    rational_form: object
    normalization_product: object
    spectrum: list


CHANNEL_BASES = {'neutral': NEUTRAL_BASIS, 'charged': CHARGED_BASIS, 'down_strange': DOWN_STRANGE_BASIS}


def compare_amplitude(invariants, alp_pole_squared_mass, daughter_squared_mass, s12, s23, channel='neutral'):
    """Closed formula versus quadratic blocks + pole vectors + CCWZ vertices + direct D(q^2)^{-1}.

    channel 'neutral': A -> 3 identical pi0-like daughters (solve_spectrum_exact);
    'charged': A -> pi+ pi- pi0-like (n1 = pi+, n2 = pi-, n3 = pi0-like). Outgoing
    pi+ uses the meson matrix sqrt(2) E_21 = (lambda_1 - i lambda_2)/sqrt(2) and
    outgoing pi- uses sqrt(2) E_12 (phase convention); they have no ALP component;
    'down_strange': as 'neutral' with d-s entries in the invariants; the A and pi0-like
    vectors are the K-normalized null vectors of D at the two rational poles (the
    remaining poles are not needed and are not found in closed form).
    The pi0-like daughter is the pole continuously connected to the lowest unmixed
    neutral meson; m_A^2 and its squared mass are imposed exactly by tuned_point.
    """
    basis = CHANNEL_BASES[channel]
    eps_squared, bare = tuned_point(invariants, alp_pole_squared_mass, daughter_squared_mass)
    eps = sp.sqrt(eps_squared)
    unmixed = unmixed_neutral_matrix(invariants, unmixed_neutral_squared_masses(invariants)[0])
    kinetic, mass_squared = solver_theory(invariants, eps, bare, basis)
    if channel == 'down_strange':
        poles = []
        alp_vector = null_space_pole_vector(kinetic, mass_squared, alp_pole_squared_mass)
        neutral_vector = null_space_pole_vector(kinetic, mass_squared, daughter_squared_mass)
    else:
        spectrum = solve_spectrum_exact(kinetic, mass_squared)
        poles = spectrum.pole_squared_masses
        alp_vector = spectrum.pole_vector(poles.index(alp_pole_squared_mass))
        neutral_vector = spectrum.pole_vector(poles.index(daughter_squared_mass))
    if sp.N(alp_vector[-1]) < 0:
        alp_vector = -alp_vector
    neutral_vector = align_with_unmixed(neutral_vector, basis, unmixed)
    neutral_leg = mixed_meson_leg(invariants, eps, daughter_squared_mass, unmixed)
    neutral_bare_leg = mixed_meson_leg(invariants, eps, daughter_squared_mass, unmixed, normalized=False)
    alp_leg = alp_like_leg(invariants, eps, alp_pole_squared_mass)
    alp_bare_leg = alp_like_leg(invariants, eps, alp_pole_squared_mass, normalized=False)
    if channel == 'charged':
        charged_mass = sum(chi_values(invariants)[:2]) / 2
        plus = sp.Matrix([0, 0, 0, 1, -sp.I, 0]) / sp.sqrt(2)
        minus = sp.Matrix([0, 0, 0, 1, sp.I, 0]) / sp.sqrt(2)
        residue = sum((spectrum.residues[n] for n, p in enumerate(poles) if p == charged_mass), sp.zeros(6))
        if residue != plus * minus.T + minus * plus.T:
            raise ValueError('charged pole residue is not the (P_1, P_2) projector')
        masses = [alp_pole_squared_mass, charged_mass, charged_mass, daughter_squared_mass]
        vectors = [alp_vector, plus, minus, neutral_vector]
        charged_legs = [(sp.sqrt(2) * unit_matrix(1, 0), 0), (sp.sqrt(2) * unit_matrix(0, 1), 0)]
        legs = [alp_leg] + charged_legs + [neutral_leg]
        bare_legs = [alp_bare_leg] + charged_legs + [neutral_bare_leg]
        norms = normalization(invariants, eps, alp_pole_squared_mass) * normalization(
            invariants, eps, daughter_squared_mass)
    else:
        masses = [alp_pole_squared_mass] + [daughter_squared_mass] * 3
        vectors = [alp_vector] + [neutral_vector] * 3
        legs = [alp_leg] + [neutral_leg] * 3
        bare_legs = [alp_bare_leg] + [neutral_bare_leg] * 3
        norms = normalization(invariants, eps, alp_pole_squared_mass) * normalization(
            invariants, eps, daughter_squared_mass)**3
    mandelstam = {(1, 2): s12, (2, 3): s23, (1, 3): sum(masses) - s12 - s23}
    contraction = ExactVertexContraction(basis, invariants).function(invariants['f_pi'] / eps)
    solver = solver_amplitude(kinetic, mass_squared, vectors, masses, mandelstam, contraction)
    pieces = three_meson_amplitude(invariants, eps, bare, legs, masses, mandelstam)
    rational_form = sp.radsimp(sp.expand(
        three_meson_amplitude(invariants, eps, bare, bare_legs, masses, mandelstam)['total']))
    return ComparisonResult(eps_squared, bare, daughter_squared_mass, pieces['total'], solver, pieces,
                            rational_form, norms, poles)


def amplitudes_at_both_signs(invariants, alp_pole_squared_mass, daughter_squared_mass, s12, s23):
    """Closed formula for A -> 3 pi0-like at +eps and -eps (same eps^2, poles and m_{a,0}^2).

    The legs follow the sign convention (A: positive ALP component; daughter:
    continuity with the unmixed matrix), so the two values must be opposite.
    """
    eps_squared, bare = tuned_point(invariants, alp_pole_squared_mass, daughter_squared_mass)
    unmixed = unmixed_neutral_matrix(invariants, unmixed_neutral_squared_masses(invariants)[0])
    masses = [alp_pole_squared_mass] + [daughter_squared_mass] * 3
    mandelstam = {(1, 2): s12, (2, 3): s23, (1, 3): sum(masses) - s12 - s23}
    values = []
    for eps in (sp.sqrt(eps_squared), -sp.sqrt(eps_squared)):
        legs = [alp_like_leg(invariants, eps, alp_pole_squared_mass)] + [
            mixed_meson_leg(invariants, eps, daughter_squared_mass, unmixed)] * 3
        values.append(three_meson_amplitude(invariants, eps, bare, legs, masses, mandelstam)['total'])
    return values


# ---------------------------------------------------------------------------
# Expansion in eps at fixed invariants, fixed m_A^2, fixed s_12 and s_23
# ---------------------------------------------------------------------------


def daughter_pole_series(invariants, alp_pole_squared_mass, unmixed_squared_mass, order=3):
    """m_n^2 = m_{n,0}^2 + sum_{k=1}^{order} c_k t^k (t = eps^2) from D_a(m_n^2) = 0 with D_a(m_A^2) = 0."""
    t, shift = sp.symbols('t delta')
    eps = sp.sqrt(t)
    bare = bare_alp_mass_for_pole(invariants, eps, alp_pole_squared_mass)
    numerator = sp.numer(sp.together(alp_inverse_propagator(invariants, eps, bare, unmixed_squared_mass + shift)))
    numerator = sp.expand(numerator)
    coefficients = sp.symbols(f'c1:{order + 1}')
    trial = sum(c * t**(k + 1) for k, c in enumerate(coefficients))
    expanded = sp.expand(numerator.subs(shift, trial))
    solution = {}
    for k in range(1, order + 1):
        equation = sp.expand(expanded.coeff(t, k).subs(solution))
        solution[coefficients[k - 1]] = sp.solve(equation, coefficients[k - 1])[0]
    return sp.expand(unmixed_squared_mass + trial.subs(solution))


def _truncate(expression, symbol, order):
    return sp.expand(sp.series(expression, symbol, 0, order).removeO())


def amplitude_series(invariants, alp_pole_squared_mass, unmixed_index, s12, s23, order=5):
    """Series in eps of M(A -> 3 identical neutral daughters) from the closed formula.

    Returns (coefficients of eps^0..eps^{order-1}, eps symbol). The daughter leg
    is written in the form regular at eps = 0: with delta = m_n^2 - m_{n,0}^2 and
    Z = delta M^{-1}_{m_n^2}(C), meson matrix 2 sign Z/sqrt(delta^2/eps^2 + 2<Z^2>)
    and ALP component sign (delta/eps)/sqrt(delta^2/eps^2 + 2<Z^2>).
    """
    e = sp.Symbol('e', positive=True)
    t, shift = sp.symbols('t delta')
    unmixed_mass = unmixed_neutral_squared_masses(invariants)[unmixed_index]
    unmixed = unmixed_neutral_matrix(invariants, unmixed_mass)
    pole = daughter_pole_series(invariants, alp_pole_squared_mass, unmixed_mass, order=3).subs(t, e**2)
    delta = sp.expand(pole - unmixed_mass)
    residue_form = inverse_pole_operator(matrix_C(invariants), invariants, unmixed_mass + shift).applyfunc(
        lambda entry: sp.cancel(shift * entry))
    Z = residue_form.applyfunc(lambda entry: _truncate(entry.subs(shift, delta), e, order + 2))
    root = sp.sqrt(sp.expand(delta**2 / e**2) + 2 * (Z * Z).trace())
    sign = 1 if sp.N((Z.subs(e, 0) * unmixed).trace()) > 0 else -1
    daughter = (Z.applyfunc(lambda entry: _truncate(2 * sign * entry / root, e, order + 2)),
                _truncate(sign * sp.expand(delta / e) / root, e, order + 2))
    bare = sp.expand(bare_alp_mass_for_pole(invariants, e, alp_pole_squared_mass))
    alp_leg = alp_like_leg(invariants, e, alp_pole_squared_mass)
    alp_leg = (alp_leg[0].applyfunc(lambda entry: _truncate(entry, e, order + 2)),
               _truncate(alp_leg[1], e, order + 2))
    masses = [alp_pole_squared_mass] + [pole] * 3
    mandelstam = {(1, 2): s12, (2, 3): s23, (1, 3): sp.expand(sum(masses) - s12 - s23)}
    legs = [alp_leg] + [daughter] * 3
    series = sum(_drop(value, e, order) for value in quartic_contact_classes(
        invariants, e, legs, decay_products(masses, mandelstam)).values())
    for i in (1, 2, 3):
        j, k = [m for m in (1, 2, 3) if m != i]
        series -= _exchange_series(invariants, e, bare, (legs[0], legs[i]), (legs[j], legs[k]),
                                   mandelstam[(j, k)], order)
    series = _drop(series, e, order)
    return [series.coeff(e, k) for k in range(order)], e


def _drop(expression, symbol, order):
    """Polynomial in symbol truncated below symbol^order."""
    expanded = sp.expand(expression)
    return sum((expanded.coeff(symbol, k) * symbol**k for k in range(order)), sp.S.Zero)


def _exchange_series(invariants, e, bare, left_legs, right_legs, squared_momentum, order):
    """Truncated series of `exchange` for polynomial legs (same formula, term by term)."""
    def polynomial_current(legs):
        current, alp_current = cubic_current(invariants, e, *legs)
        return current.applyfunc(lambda y: _drop(y, e, order + 2)), _drop(alp_current, e, order + 2)

    left, left_alp = polynomial_current(left_legs)
    right, right_alp = polynomial_current(right_legs)

    def series(entry):
        return _truncate(entry, e, order + 2)

    meson_block = -2 * (left * inverse_pole_operator(right, invariants, squared_momentum).applyfunc(series)).trace()
    overlap = inverse_pole_operator(matrix_C(invariants), invariants, squared_momentum).applyfunc(series)
    left_direction = _drop(2 * e * (left * overlap).trace() + left_alp, e, order + 2)
    right_direction = _drop(2 * e * (right * overlap).trace() + right_alp, e, order + 2)
    propagator = series(1 / alp_inverse_propagator(invariants, e, bare, squared_momentum))
    return _drop(meson_block + left_direction * right_direction * propagator, e, order)


def overlap_argument_difference_series(invariants, alp_pole_squared_mass, unmixed_index, s12, s23, order=5):
    """Series of eps [leading(Y at m_{a,0}^2) - leading(Y at m_A^2)] (expected O(eps^3))."""
    e = sp.Symbol('e', positive=True)
    unmixed_mass = unmixed_neutral_squared_masses(invariants)[unmixed_index]
    unmixed = unmixed_neutral_matrix(invariants, unmixed_mass)
    masses = [unmixed_mass] * 3
    mandelstam = {(1, 2): s12, (2, 3): s23, (1, 3): alp_pole_squared_mass + sum(masses) - s12 - s23}
    bare = bare_alp_mass_for_pole(invariants, e, alp_pole_squared_mass)
    difference = e * (leading_coefficient(invariants, alp_pole_squared_mass, [unmixed] * 3, masses, mandelstam,
                                          overlap_argument=bare)
                      - leading_coefficient(invariants, alp_pole_squared_mass, [unmixed] * 3, masses, mandelstam))
    series = _truncate(difference, e, order)
    return [series.coeff(e, k) for k in range(order)]


# ---------------------------------------------------------------------------


FORMULA_TEXT = r"""
Closed formula (dressed frame, all momenta incoming, legs l = A, n1, n2, n3 with
meson matrix Phi_l = sum_k v^(l)_k lambda_k and ALP component alpha_l = v^(l)_a):
  Phi_A = 2 eps Y / sqrt(N(m_A^2)),  alpha_A = 1/sqrt(N(m_A^2)),  Y = M^{-1}_{m_A^2}(C)
  Phi_n = sgn_n 2 eps M^{-1}_{m_n^2}(C)/sqrt(N(m_n^2)), alpha_n = sgn_n/sqrt(N(m_n^2)) (mixed poles)
  M = V4_kin + V4_PPPP + V4_PPPa + V4_PPaa + V4_Paaa + V4_aaaa - sum_{i} E_i
  V4_kin  = -1/(24 f^2) sum_i [(p_nj+p_nk)^2 - (m_A^2+m_n1^2+m_n2^2+m_n3^2)/2]
              (<[Phi_nj,Phi_A][Phi_nk,Phi_ni]> + <[Phi_nk,Phi_A][Phi_nj,Phi_ni]>)
  V4_PPPP = B_0/(24 f^2) sum_{S_4} <m_q Phi Phi Phi Phi>
  V4_PPPa = -eps B_0/(6 f^2) sum_l alpha_l sum_{S_3} <hat P_1 Phi Phi Phi>
  V4_PPaa = -eps^2 B_0/(2 f^2) sum_{l<m} alpha_l alpha_m <hat S_2 {Phi_k, Phi_n}>
  V4_Paaa = eps^3/f^2 sum_l (prod_{m!=l} alpha_m) <(B_0 hat P_3 - m_0^2/N_f hat Theta_3) Phi_l>
  V4_aaaa = eps^4/f^2 [B_0 <hat S_4> - m_0^2/N_f (8 hat Theta_1 hat Theta_3 + 6 hat Theta_2^2)] prod alpha
  E_i (partition A n_i | n_j n_k, q^2 = (p_nj+p_nk)^2), with sum_k V_{l m k} lambda_k = 2 J(l,m):
     J(l,m)  = -eps B_0/(2f)(alpha_l {hat S_1,Phi_m} + alpha_m {hat S_1,Phi_l})
               + eps^2/f alpha_l alpha_m (B_0 hat P_2 - m_0^2/N_f hat Theta_2 1)
     J_a(l,m)= -eps B_0/(2f) <hat S_1 {Phi_l,Phi_m}> + eps^2/f (alpha_l <(..P_2..)Phi_m> + alpha_m <(..P_2..)Phi_l>)
               + eps^3/f (B_0 <hat S_3> - 6 m_0^2/N_f hat Theta_1 hat Theta_2) alpha_l alpha_m
     E_i = -2 <J(A,ni) M^{-1}_{q^2}(J(nj,nk))>
           + (2 eps <J(A,ni) M^{-1}_{q^2}(C)> + J_a(A,ni)) (2 eps <J(nj,nk) M^{-1}_{q^2}(C)> + J_a(nj,nk)) / D_a(q^2)
  eps^1 coefficient: V4_kin + V4_PPPP at (Phi_A -> 2Y, Phi_n -> unmixed, m_n -> m_n0)
                     - B_0/(6 f^2) sum_{S_3} <hat P_1 Phi_n1,0 Phi_n2,0 Phi_n3,0>
"""


def main():
    print(FORMULA_TEXT)
    invariants = verification_inputs()
    print('Unmixed neutral squared masses (GeV^2):', unmixed_neutral_squared_masses(invariants))
    m_A2 = ALP_POLE_SQUARED_MASS
    eps_symbol, bare_symbol = sp.symbols('epsilon m_a0sq', positive=True)
    kinetic, mass_squared = solver_theory(invariants, eps_symbol, bare_symbol, NEUTRAL_BASIS)
    print('M^2_aa from quadratic_blocks:', sp.factor(mass_squared[-1, -1]),
          '; expected:', sp.factor(alp_mass_entry(invariants, eps_symbol, bare_symbol)))
    unmixed_mass = unmixed_neutral_squared_masses(invariants)[0]
    s12, s23 = KINEMATIC_POINT
    print(f'\nA -> 3 pi0-like, m_A^2 = {m_A2}, s12 = {s12}, s23 = {s23} (GeV^2)')
    coefficients, _ = amplitude_series(invariants, m_A2, 0, s12, s23)
    unmixed = unmixed_neutral_matrix(invariants, unmixed_mass)
    mandelstam0 = {(1, 2): s12, (2, 3): s23, (1, 3): m_A2 + 3 * unmixed_mass - s12 - s23}
    leading = sp.radsimp(leading_coefficient(invariants, m_A2, [unmixed] * 3, [unmixed_mass] * 3, mandelstam0))
    print('series coefficients eps^0..eps^4:', [sp.N(c, 12) for c in coefficients])
    print('closed eps^1 coefficient:', leading, '=', sp.N(leading, 15))
    print('eps^3 coefficient (series):', sp.radsimp(coefficients[3]), '=', sp.N(coefficients[3], 15))
    shift = overlap_argument_difference_series(invariants, m_A2, 0, s12, s23)
    print('eps [leading(Y(m_a0^2)) - leading(Y(m_A^2))] series eps^0..eps^4:', [sp.N(c, 10) for c in shift])
    print('\n eps^2 | m_a0^2 | m_pi^2 | M exact | M(unnormalized legs) | eps^1 term | (M - eps^1 term)/eps^3 | closed == solver')
    for offset in OFFSETS:
        result = compare_amplitude(invariants, m_A2, unmixed_mass + offset, s12, s23)
        eps = sp.sqrt(result.eps_squared)
        difference = (result.closed - eps * leading) / eps**3
        print(f' {result.eps_squared} | {result.bare_alp_mass_squared} | {result.daughter_squared_mass} | '
              f'{sp.N(result.closed, 15)} | {result.rational_form} | {sp.N(eps * leading, 15)} | '
              f'{sp.N(difference, 15)} | {is_exactly_zero(result.closed - result.solver)}')
    plus, minus = amplitudes_at_both_signs(invariants, m_A2, unmixed_mass + OFFSETS[1], s12, s23)
    print('M(+eps) + M(-eps) = 0 exactly:', is_exactly_zero(plus + minus))
    for channel, data in (('charged', invariants), ('down_strange', verification_inputs(down_strange=True))):
        result = compare_amplitude(data, m_A2, unmixed_mass + OFFSETS[1], s12, s23, channel=channel)
        print(f'\n{channel}: eps^2 = {result.eps_squared}, m_a0^2 = {result.bare_alp_mass_squared}, '
              f'M = {sp.N(result.closed, 15).as_real_imag()[0]}, M(unnormalized legs) = {result.rational_form}, '
              f'closed == solver: {is_exactly_zero(result.closed - result.solver)}')
        print('  pieces (real parts; imaginary parts vanish exactly in the sum):',
              {key: sp.N(value, 8).as_real_imag()[0] for key, value in result.pieces.items()})

if __name__ == '__main__':
    main()
