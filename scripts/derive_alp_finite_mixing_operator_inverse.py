#!/usr/bin/env python3
"""LO finite ALP-meson mixing written with the matrix operator of Alda:2026nrd.

Run:   .venv/bin/python scripts/derive_alp_finite_mixing_operator_inverse.py
Tests: .venv/bin/python -m unittest discover -s scripts -p test_derive_alp_finite_mixing_operator_inverse.py

Conventions (Alda:2026nrd; memory/alp/shared/finite-alp-meson-mixing.md;
memory/nlo/lsz-operator-inverse.md). N_f = 3, complete Hermitian U(3) basis
lambda_i (i = 0..8) with tr(lambda_i lambda_j) = 2 delta_ij, fields
Phi = (P_i, a), L_2 = 1/2 dPhi^T K dPhi - 1/2 Phi^T M^2 Phi, D(s) = s K - M^2.
epsilon = f_pi/f_a, chi_0 = diag(chi_u, chi_d, chi_s), Hermitian c^A, c^P (all
entries, including off-diagonal ones), real c_G, c^AP = c^A - c^P,
c_GA = c_G + tr(c^A). LO blocks:
  K_ij = delta_ij, K_ia = -v_i, K_aa = 1 + v.v, v_i = epsilon tr(c^A lambda_i),
  M2_ij = tr({lambda_i, lambda_j} chi_0)/4 + m_0^2 tr(lambda_i) tr(lambda_j)/(2 N_f),
  M2_ia = -(epsilon/2) tr(lambda_i {c^P, chi_0}) + epsilon (m_0^2/N_f) c_G tr(lambda_i),
  M2_aa = m_a0^2 + 2 epsilon^2 (m_0^2/N_f) c_G^2.
The blocks are taken from scripts/derive_alp_lo_neutral_blocks.py
(closed_form_blocks with f_pi -> epsilon, f_a -> 1).
Meson overlap: w_i(s) = [(M2 - s K)^{-1}]_ij (s K_ja - M2_ja), W = w_i lambda_i.
LO pole operator: M_a[X] = {chi_0, X}/2 + (m_0^2/N_f) tr(X) 1 - s X.
Invariant matrix of Alda:2026nrd: calC = -{chi_0, c^AP}/2 - (m_0^2/N_f) c_GA 1.
Flavor-matrix components: x_i = tr(lambda_i X)/2.

Derived and checked here (step numbers follow main()):
  1. M_a[W] = C(s), C = epsilon({c^P, chi_0} - 2 s c^A - 2 (m_0^2/N_f) c_G 1).
  2. C(s) - M_a[2 epsilon c^A] = 2 epsilon calC (independent of s), hence
     w_i lambda_i = 2 epsilon [c^A + M_a^{-1}(calC)].
  3. Explicit inverse:
     M_a^{-1}(calC)_ij = calC_ij/((chi_i + chi_j)/2 - s)
                       = -((chi_i + chi_j)/2) c^AP_ij/((chi_i + chi_j)/2 - s)   (i != j),
     M_a^{-1}(calC)_ii = (calC_ii - (m_0^2/N_f) tr M_a^{-1}(calC))/(chi_i - s),
     calC_ii = -chi_i c^AP_ii - (m_0^2/N_f) c_GA,
     tr M_a^{-1}(calC) = [sum_k calC_kk/(chi_k - s)]/[1 + (m_0^2/N_f) sum_k 1/(chi_k - s)].
  4. N(s) = K_aa + 2 K_ia w_i + w_i w_i = 1 + 2 epsilon^2 tr[M_a^{-1}(calC)^2].
  5. M^2_{aa,can} = M2_aa + 2 v.M2_Pa + v.M2.v
           = m_a0^2 + 2 epsilon^2 (m_0^2/N_f) c_GA^2 - 2 epsilon^2 [tr(c^A {c^P, chi_0}) - tr(chi_0 c^A c^A)].
  6. The ALP inverse propagator with the mesons integrated out,
     D_a(s) = s - M^2_{aa,can} + 2 epsilon^2 tr[calC M_a^{-1}(calC)], is the Schur
     complement of D_PP in D(s): det D(s) = det(D_PP(s)) D_a(s) = -det(M2_PP - s) D_a(s).
     Explicitly,
     tr[calC M_a^{-1}(calC)] = sum_{i != j} ((chi_i + chi_j)/2)^2 |c^AP_ij|^2/((chi_i + chi_j)/2 - s)
                               + sum_i calC_ii^2/(chi_i - s)
                               - (m_0^2/N_f) [sum_i calC_ii/(chi_i - s)]^2/[1 + (m_0^2/N_f) sum_i 1/(chi_i - s)].
  7. N(s) = dD_a/ds.
  8. Cross-checks: numerical pole solver, epsilon -> 0 limit, shift of the
     couplings that keeps c^AP, c_GA and calC fixed, quartic polynomial of the
     flavor-diagonal sector (P_3, P_8, P_0, a).
  9. One-ALP amplitude: M = (V_a + V_i w_i)/sqrt(N)
                          = (V_a + epsilon tr{V [c^A + M_a^{-1}(calC)]})/sqrt(N), V = V_i lambda_i.

Checks. Steps 1, 2, 3, 5, 6 (Schur complement), 8b, 8c and 9 are checked
with generic symbols for all inputs. Steps 4 and 7 are checked with generic
symbolic sources and rational chi_i, m_0^2 (and, for 7, with symbolic chi_i,
m_0^2 and rational sources), because the fully symbolic form is too slow in
SymPy. The determinant identity, the sign of the meson term below the meson
poles and the comparison with scripts/solve_alp_meson_spectrum.py use exact
rational random Hermitian data for the seeds in SEEDS. For generic sources
each pair i < j gives one real meson combination orthogonal to c^AP_ij that
decouples exactly at s = (chi_i + chi_j)/2 (v_a = 0, w_i undefined); det D
keeps a simple zero there.

Units: s, chi_i, m_0^2, m_a0^2, M^2, calC, D_a(s) in GeV^2; epsilon, c^A, c^P,
c_G, K, w_i, N dimensionless. Rational inputs are illustrative, not fits.
Assumptions as in the finite-mixing record: LO only, canonical bare ALP
term, quark spurion linear in a (no a^2 spurion term), no loops.
"""

from __future__ import annotations

import random

import numpy as np
import sympy as sp

from derive_alp_lo_neutral_blocks import ALL_GENERATOR_INDICES, closed_form_blocks
from derive_nlo_lsz_operator_inverse import inverse_from_trace_system
from solve_alp_meson_spectrum import solve_spectrum_numerical
from symbolic_utils import anticommutator, u3_generators

FLAVOR_COUNT = 3
GENERATORS = u3_generators()
SQUARED_MOMENTUM = sp.Symbol("s", real=True)
SEEDS = (1, 2, 3)


# ------------------------------------------------------------ basic objects


def generator_components(flavor_matrix):
    """Components x_i = tr(lambda_i X)/2 of a flavor matrix X = x_i lambda_i."""
    return sp.Matrix([sp.expand(sp.trace(generator * flavor_matrix) / 2) for generator in GENERATORS])


def matrix_from_components(components):
    """Flavor matrix sum_i x_i lambda_i."""
    return sum((component * generator for component, generator in zip(components, GENERATORS)), sp.zeros(3))


def symbolic_hermitian_matrix(name):
    """Generic 3x3 Hermitian matrix built from nine real symbols."""
    matrix = sp.zeros(3)
    for row in range(3):
        matrix[row, row] = sp.Symbol(f"{name}_{row}{row}", real=True)
        for column in range(row + 1, 3):
            real_part = sp.Symbol(f"{name}_re{row}{column}", real=True)
            imaginary_part = sp.Symbol(f"{name}_im{row}{column}", real=True)
            matrix[row, column] = real_part + sp.I * imaginary_part
            matrix[column, row] = real_part - sp.I * imaginary_part
    return matrix


def random_rational_hermitian(generator, denominator=6):
    """Hermitian matrix with exact rational random entries (real and imaginary parts)."""
    def entry():
        return sp.Rational(generator.choice([-1, 1]) * generator.randint(1, denominator), denominator)
    matrix = sp.zeros(3)
    for row in range(3):
        matrix[row, row] = entry()
        for column in range(row + 1, 3):
            matrix[row, column] = entry() + sp.I * entry()
            matrix[column, row] = sp.conjugate(matrix[row, column])
    return matrix


def symbolic_inputs():
    """Generic symbols for all inputs (c^A, c^P generic Hermitian)."""
    return {
        "axial_source": symbolic_hermitian_matrix("cA"),
        "pseudoscalar_source": symbolic_hermitian_matrix("cP"),
        "chi_values": sp.symbols("chi_u chi_d chi_s", positive=True),
        "gluon_coupling": sp.Symbol("c_G", real=True),
        "epsilon": sp.Symbol("epsilon", positive=True),
        "anomaly_mass_squared": sp.Symbol("m_0^2", positive=True),
        "bare_alp_mass_squared": sp.Symbol("m_a0_sq", positive=True),
    }


def rational_inputs(seed):
    """Exact rational illustrative inputs (GeV units) with generic Hermitian sources."""
    generator = random.Random(seed)
    chi_values = sorted(sp.Rational(generator.randint(1, 60), 100) for _ in range(3))
    while len(set(chi_values)) < 3:
        chi_values = sorted(sp.Rational(generator.randint(1, 60), 100) for _ in range(3))
    return {
        "axial_source": random_rational_hermitian(generator),
        "pseudoscalar_source": random_rational_hermitian(generator),
        "chi_values": tuple(chi_values),
        "gluon_coupling": sp.Rational(generator.randint(-6, 6), 6),
        "epsilon": sp.Rational(generator.randint(1, 5), 10),
        "anomaly_mass_squared": FLAVOR_COUNT * sp.Rational(generator.randint(20, 35), 100),
        "bare_alp_mass_squared": sp.Rational(generator.randint(5, 80), 100),
    }


def lo_blocks(inputs):
    """10x10 K and M^2 on (P_0..P_8, a) from the closed forms of derive_alp_lo_neutral_blocks.py."""
    return closed_form_blocks(
        axial_source=inputs["axial_source"],
        pseudoscalar_source=inputs["pseudoscalar_source"],
        quark_spurion=sp.diag(*inputs["chi_values"]),
        gluon_coupling=inputs["gluon_coupling"],
        pion_decay_constant=inputs["epsilon"],
        alp_decay_constant=1,
        anomaly_mass_squared=inputs["anomaly_mass_squared"],
        bare_alp_mass_squared=inputs["bare_alp_mass_squared"],
        generator_indices=ALL_GENERATOR_INDICES,
    )


def lo_pole_operator(flavor_matrix, chi_values, anomaly_mass_squared, squared_momentum):
    """M_a[X] = {chi_0, X}/2 + (m_0^2/N_f) tr(X) 1 - s X."""
    quark_spurion = sp.diag(*chi_values)
    return (anticommutator(quark_spurion, flavor_matrix) / 2
            + anomaly_mass_squared / FLAVOR_COUNT * sp.trace(flavor_matrix) * sp.eye(3)
            - squared_momentum * flavor_matrix)


def pole_source(inputs, squared_momentum):
    """C(s) = epsilon ({c^P, chi_0} - 2 s c^A - 2 (m_0^2/N_f) c_G 1)."""
    return inputs["epsilon"] * (
        anticommutator(inputs["pseudoscalar_source"], sp.diag(*inputs["chi_values"]))
        - 2 * squared_momentum * inputs["axial_source"]
        - 2 * inputs["anomaly_mass_squared"] / FLAVOR_COUNT * inputs["gluon_coupling"] * sp.eye(3))


def axial_minus_pseudoscalar_coupling(inputs):
    """c^AP = c^A - c^P (notation of Alda:2026nrd)."""
    return inputs["axial_source"] - inputs["pseudoscalar_source"]


def gluon_plus_axial_trace_coupling(inputs):
    """c_GA = c_G + tr(c^A) (notation of Alda:2026nrd)."""
    return inputs["gluon_coupling"] + sp.trace(inputs["axial_source"])


def invariant_matrix(inputs):
    """calC = -{chi_0, c^AP}/2 - (m_0^2/N_f) c_GA 1 (Alda:2026nrd); independent of s."""
    return (-anticommutator(sp.diag(*inputs["chi_values"]), axial_minus_pseudoscalar_coupling(inputs)) / 2
            - inputs["anomaly_mass_squared"] / FLAVOR_COUNT * gluon_plus_axial_trace_coupling(inputs) * sp.eye(3))


def inverse_lo_pole_operator(source, chi_values, anomaly_mass_squared, squared_momentum):
    """M_a^{-1}(Y): Y_ij/((chi_i + chi_j)/2 - s) off the diagonal, (Y_ii - (m_0^2/N_f) tr M_a^{-1}(Y))/(chi_i - s) on it."""
    m0_squared_over_flavors = anomaly_mass_squared / FLAVOR_COUNT
    chi_minus_s = [chi - squared_momentum for chi in chi_values]
    trace_of_inverse = (sum(source[index, index] / chi_minus_s[index] for index in range(3))
                        / (1 + m0_squared_over_flavors * sum(1 / value for value in chi_minus_s)))
    solution = sp.zeros(3)
    for row in range(3):
        for column in range(3):
            if row == column:
                solution[row, row] = (source[row, row] - m0_squared_over_flavors * trace_of_inverse) / chi_minus_s[row]
            else:
                average_chi = (chi_values[row] + chi_values[column]) / 2
                solution[row, column] = source[row, column] / (average_chi - squared_momentum)
    return solution


def _diagonal_ingredients(inputs, squared_momentum):
    """calC_ii, chi_i - s and tr M_a^{-1}(calC), written with c^AP and c_GA."""
    chi_values = inputs["chi_values"]
    m0_squared_over_flavors = inputs["anomaly_mass_squared"] / FLAVOR_COUNT
    cAP, cGA = axial_minus_pseudoscalar_coupling(inputs), gluon_plus_axial_trace_coupling(inputs)
    chi_minus_s = [chi - squared_momentum for chi in chi_values]
    calC_diagonal = [-chi_values[index] * cAP[index, index] - m0_squared_over_flavors * cGA for index in range(3)]
    trace_of_inverse = (sum(calC_diagonal[index] / chi_minus_s[index] for index in range(3))
                        / (1 + m0_squared_over_flavors * sum(1 / value for value in chi_minus_s)))
    return m0_squared_over_flavors, cAP, chi_minus_s, calC_diagonal, trace_of_inverse


OFF_DIAGONAL_PAIRS = [(row, column) for row in range(3) for column in range(3) if row != column]


def closed_form_inverse_of_invariant_matrix(inputs, squared_momentum):
    """M_a^{-1}(calC) of step 3, written explicitly with c^AP, c_GA, chi_i and s."""
    chi_values = inputs["chi_values"]
    m0_squared_over_flavors, cAP, chi_minus_s, calC_diagonal, trace_of_inverse = _diagonal_ingredients(inputs, squared_momentum)
    solution = sp.zeros(3)
    for row in range(3):
        solution[row, row] = (calC_diagonal[row] - m0_squared_over_flavors * trace_of_inverse) / chi_minus_s[row]
    for row, column in OFF_DIAGONAL_PAIRS:
        average_chi = (chi_values[row] + chi_values[column]) / 2
        solution[row, column] = -average_chi * cAP[row, column] / (average_chi - squared_momentum)
    return solution


def closed_form_overlap(inputs, squared_momentum):
    """W of step 3: off-diagonal 2 epsilon ((chi_i+chi_j)/2 c^P_ij - s c^A_ij)/((chi_i+chi_j)/2 - s); diagonal from 2 epsilon [c^A + M_a^{-1}(calC)]."""
    inverse = closed_form_inverse_of_invariant_matrix(inputs, squared_momentum)
    epsilon = inputs["epsilon"]
    overlap = sp.zeros(3)
    for row in range(3):
        overlap[row, row] = 2 * epsilon * (inputs["axial_source"][row, row] + inverse[row, row])
    for row, column in OFF_DIAGONAL_PAIRS:
        average_chi = (inputs["chi_values"][row] + inputs["chi_values"][column]) / 2
        overlap[row, column] = 2 * epsilon * (
            average_chi * inputs["pseudoscalar_source"][row, column]
            - squared_momentum * inputs["axial_source"][row, column]) / (average_chi - squared_momentum)
    return overlap


def closed_form_normalization(inputs, squared_momentum):
    """N(s) of step 4 as explicit sums: 1 + 2 epsilon^2 tr[M_a^{-1}(calC)^2]."""
    chi_values = inputs["chi_values"]
    m0_squared_over_flavors, cAP, chi_minus_s, calC_diagonal, trace_of_inverse = _diagonal_ingredients(inputs, squared_momentum)
    off_diagonal = sum(((chi_values[i] + chi_values[j]) / 2)**2 * cAP[i, j] * cAP[j, i]
                       / ((chi_values[i] + chi_values[j]) / 2 - squared_momentum)**2 for i, j in OFF_DIAGONAL_PAIRS)
    diagonal = sum((calC_diagonal[i] - m0_squared_over_flavors * trace_of_inverse)**2 / chi_minus_s[i]**2 for i in range(3))
    return 1 + 2 * inputs["epsilon"]**2 * (off_diagonal + diagonal)


def closed_form_M2_aa_can(inputs):
    """M^2_{aa,can} of step 5."""
    epsilon, axial = inputs["epsilon"], inputs["axial_source"]
    quark_spurion = sp.diag(*inputs["chi_values"])
    return (inputs["bare_alp_mass_squared"]
            + 2 * epsilon**2 * inputs["anomaly_mass_squared"] / FLAVOR_COUNT * gluon_plus_axial_trace_coupling(inputs)**2
            - 2 * epsilon**2 * (sp.trace(axial * anticommutator(inputs["pseudoscalar_source"], quark_spurion))
                                - sp.trace(quark_spurion * axial * axial)))


def closed_form_alp_inverse_propagator(inputs, squared_momentum):
    """D_a(s) of step 6 as explicit sums: s - M^2_{aa,can} + 2 epsilon^2 tr[calC M_a^{-1}(calC)]."""
    chi_values = inputs["chi_values"]
    m0_squared_over_flavors, cAP, chi_minus_s, calC_diagonal, _ = _diagonal_ingredients(inputs, squared_momentum)
    off_diagonal = sum(((chi_values[i] + chi_values[j]) / 2)**2 * cAP[i, j] * cAP[j, i]
                       / ((chi_values[i] + chi_values[j]) / 2 - squared_momentum) for i, j in OFF_DIAGONAL_PAIRS)
    diagonal_ratio_sum = sum(calC_diagonal[i] / chi_minus_s[i] for i in range(3))
    diagonal = (sum(calC_diagonal[i]**2 / chi_minus_s[i] for i in range(3))
                - m0_squared_over_flavors * diagonal_ratio_sum**2
                / (1 + m0_squared_over_flavors * sum(1 / value for value in chi_minus_s)))
    return squared_momentum - closed_form_M2_aa_can(inputs) + 2 * inputs["epsilon"]**2 * (off_diagonal + diagonal)


def meson_block_determinant(chi_values, anomaly_mass_squared, squared_momentum):
    """det(M2_PP - s) = prod_{i<j} ((chi_i+chi_j)/2 - s)^2 prod_i (chi_i - s) [1 + (m_0^2/N_f) sum_i 1/(chi_i - s)]."""
    chi_minus_s = [chi - squared_momentum for chi in chi_values]
    off_diagonal = sp.prod([((chi_values[i] + chi_values[j]) / 2 - squared_momentum)**2
                            for i in range(3) for j in range(i + 1, 3)])
    return off_diagonal * sp.prod(chi_minus_s) * (1 + anomaly_mass_squared / FLAVOR_COUNT * sum(1 / value for value in chi_minus_s))


def is_zero(expression):
    """Exact zero test for rational functions with algebraic coefficients."""
    if sp.expand(expression) == 0:
        return True
    numerator = sp.expand(sp.numer(sp.together(expression)))
    return numerator == 0 or sp.expand(sp.radsimp(numerator)) == 0


def matrix_is_zero(matrix):
    return all(is_zero(entry) for entry in matrix)


# ------------------------------------------------------------- the steps


def step_pole_equation_matrix_form(inputs):
    """Step 1: components of M_a[X] and C(s) reproduce (M2 - sK)_PP x and sK_Pa - M2_Pa."""
    kinetic, mass_squared = lo_blocks(inputs)
    components = sp.Matrix(sp.symbols("x_0:9", real=True))
    operator_components = generator_components(lo_pole_operator(
        matrix_from_components(components), inputs["chi_values"], inputs["anomaly_mass_squared"], SQUARED_MOMENTUM))
    block_product = (mass_squared[:9, :9] - SQUARED_MOMENTUM * kinetic[:9, :9]) * components
    source_components = generator_components(pole_source(inputs, SQUARED_MOMENTUM))
    block_source = SQUARED_MOMENTUM * kinetic[:9, 9] - mass_squared[:9, 9]
    return {
        "(1/2) tr(lambda_i M_a[X]) = (M2 - sK)_ij x_j": matrix_is_zero(operator_components - block_product),
        "(1/2) tr(lambda_i C(s)) = s K_ia - M2_ia": matrix_is_zero(source_components - block_source),
    }


def step_canonical_source(inputs):
    """Step 2: C(s) - M_a[2 epsilon c^A] = 2 epsilon calC and v_i are the components of 2 epsilon c^A."""
    kinetic, _ = lo_blocks(inputs)
    axial_overlap = 2 * inputs["epsilon"] * inputs["axial_source"]
    operator_on_axial = lo_pole_operator(axial_overlap, inputs["chi_values"], inputs["anomaly_mass_squared"], SQUARED_MOMENTUM)
    return {
        "v_i = (1/2) tr(lambda_i 2 epsilon c^A)": matrix_is_zero(generator_components(axial_overlap) + kinetic[:9, 9]),
        "C(s) - M_a[2 epsilon c^A] = 2 epsilon calC": matrix_is_zero(
            pole_source(inputs, SQUARED_MOMENTUM) - operator_on_axial - 2 * inputs["epsilon"] * invariant_matrix(inputs)),
        "calC has no s dependence": all(not entry.has(SQUARED_MOMENTUM) for entry in invariant_matrix(inputs)),
    }


def step_explicit_inverse(inputs):
    """Step 3: the closed-form M_a^{-1}(calC) and W solve the operator equations."""
    inverse = closed_form_inverse_of_invariant_matrix(inputs, SQUARED_MOMENTUM)
    operator_on_inverse = lo_pole_operator(inverse, inputs["chi_values"], inputs["anomaly_mass_squared"], SQUARED_MOMENTUM)
    overlap = closed_form_overlap(inputs, SQUARED_MOMENTUM)
    operator_on_overlap = lo_pole_operator(overlap, inputs["chi_values"], inputs["anomaly_mass_squared"], SQUARED_MOMENTUM)
    general_inverse = inverse_lo_pole_operator(invariant_matrix(inputs), inputs["chi_values"],
                                               inputs["anomaly_mass_squared"], SQUARED_MOMENTUM)
    _, _, _, _, trace_of_inverse = _diagonal_ingredients(inputs, SQUARED_MOMENTUM)
    return {
        "M_a[M_a^{-1}(calC)] = calC": matrix_is_zero(operator_on_inverse - invariant_matrix(inputs)),
        "M_a[W] = C(s) for the closed-form W": matrix_is_zero(operator_on_overlap - pole_source(inputs, SQUARED_MOMENTUM)),
        "W = 2 epsilon [c^A + M_a^{-1}(calC)] (all entries)": matrix_is_zero(
            overlap - 2 * inputs["epsilon"] * (inputs["axial_source"] + inverse)),
        "closed form = general trace inverse applied to calC": matrix_is_zero(inverse - general_inverse),
        "tr M_a^{-1}(calC) equals its closed form": is_zero(sp.trace(inverse) - trace_of_inverse),
        "M_a^{-1}(calC) is Hermitian": matrix_is_zero(inverse - inverse.H),
    }


def step_explicit_inverse_against_trace_system(inputs, squared_momentum):
    """Step 3 cross-check: the NLO-record trace-system inverse at LO (rational data)."""
    reference, _, _ = inverse_from_trace_system(invariant_matrix(inputs), list(inputs["chi_values"]),
                                                squared_momentum, inputs["anomaly_mass_squared"] / FLAVOR_COUNT)
    inverse = closed_form_inverse_of_invariant_matrix(inputs, squared_momentum)
    return {"M_a^{-1}(calC) equals inverse_from_trace_system(calC) at LO": matrix_is_zero(inverse - reference)}


def step_component_completeness():
    """Completeness of the lambda_i: sum_i a_i b_i = tr(A B)/2 for A = a_i lambda_i, B = b_i lambda_i (generic a_i, b_i)."""
    first, second = sp.Matrix(sp.symbols("a_0:9", real=True)), sp.Matrix(sp.symbols("b_0:9", real=True))
    trace_form = sp.trace(matrix_from_components(first) * matrix_from_components(second)) / 2
    return {"sum_i a_i b_i = tr(A B)/2": is_zero((first.T * second)[0] - trace_form)}


def step_normalization(inputs):
    """Step 4: N(s) = K_aa + 2 K_ia w_i + w_i K_ij w_j.

    (a) For generic w_i (9 symbols): N = 1 + sum_i (w_i - v_i)^2.
    (b) With completeness and w_i lambda_i = 2 epsilon [c^A + M_a^{-1}(calC)]
        this is 1 + 2 epsilon^2 tr[M_a^{-1}(calC)^2]; the flavor-basis
        expression with the closed-form W is compared with it and with the
        explicit sums.
    """
    kinetic, _ = lo_blocks(inputs)
    generic_overlap = sp.Matrix(sp.symbols("w_0:9", real=True))
    normalization_generic = (kinetic[9, 9] + 2 * (kinetic[:9, 9].T * generic_overlap)[0]
                             + (generic_overlap.T * kinetic[:9, :9] * generic_overlap)[0])
    generic_difference = generic_overlap + kinetic[:9, 9]
    overlap = closed_form_overlap(inputs, SQUARED_MOMENTUM)
    axial_overlap = 2 * inputs["epsilon"] * inputs["axial_source"]
    normalization_flavor_basis = (kinetic[9, 9] - sp.trace(axial_overlap * overlap) + sp.trace(overlap * overlap) / 2)
    inverse = closed_form_inverse_of_invariant_matrix(inputs, SQUARED_MOMENTUM)
    trace_form = 1 + 2 * inputs["epsilon"]**2 * sp.trace(inverse * inverse)
    return {
        "N = 1 + sum_i (w_i - v_i)^2 for generic w": is_zero(normalization_generic - 1 - (generic_difference.T * generic_difference)[0]),
        "K_aa - tr(2 epsilon c^A W) + tr(W^2)/2 = 1 + 2 epsilon^2 tr[M_a^{-1}(calC)^2]": is_zero(normalization_flavor_basis - trace_form),
        "1 + 2 epsilon^2 tr[M_a^{-1}(calC)^2] equals the explicit sums": is_zero(trace_form - closed_form_normalization(inputs, SQUARED_MOMENTUM)),
    }


def step_M2_aa_can(inputs):
    """Step 5: M^2_{aa,can} = M2_aa + 2 v.M2_Pa + v.M2.v."""
    kinetic, mass_squared = lo_blocks(inputs)
    shift_vector = -kinetic[:9, 9]
    from_blocks = (mass_squared[9, 9] + 2 * (shift_vector.T * mass_squared[:9, 9])[0]
                   + (shift_vector.T * mass_squared[:9, :9] * shift_vector)[0])
    quark_spurion = sp.diag(*inputs["chi_values"])
    axial, epsilon = inputs["axial_source"], inputs["epsilon"]
    written_with_cAP = (inputs["bare_alp_mass_squared"]
                        + 2 * epsilon**2 * inputs["anomaly_mass_squared"] / FLAVOR_COUNT * gluon_plus_axial_trace_coupling(inputs)**2
                        - 2 * epsilon**2 * (sp.trace(quark_spurion * axial * axial)
                                            - sp.trace(axial * anticommutator(axial_minus_pseudoscalar_coupling(inputs), quark_spurion))))
    canonical_kinetic = kinetic[9, 9] - (shift_vector.T * shift_vector)[0]
    return {
        "M^2_{aa,can} equals the closed form": is_zero(from_blocks - closed_form_M2_aa_can(inputs)),
        "M^2_{aa,can} = m_a0^2 + 2eps^2 (m_0^2/N_f) c_GA^2 - 2eps^2[tr(chi c^A c^A) - tr(c^A {c^AP, chi})]": is_zero(from_blocks - written_with_cAP),
        "K_{aa,can} = K_aa + 2 v.K_Pa + v.v = 1": is_zero(canonical_kinetic - 1),
    }


def step_alp_inverse_propagator(inputs):
    """Step 6: the Schur complement D_aa - D_aP D_PP^{-1} D_Pa (blocks of D(s)) equals D_a(s).

    D_PP^{-1} D_Pa = -w, so the Schur complement is s K_aa - M2_aa + C_i w_i
    = s K_aa - M2_aa + tr(C W)/2 (completeness). It is compared, in the flavor
    basis, with s - M^2_{aa,can} + 2 epsilon^2 tr[calC M_a^{-1}(calC)] and with the
    explicit sums of D_a(s).
    """
    kinetic, mass_squared = lo_blocks(inputs)
    overlap = closed_form_overlap(inputs, SQUARED_MOMENTUM)
    schur_complement = (SQUARED_MOMENTUM * kinetic[9, 9] - mass_squared[9, 9]
                        + sp.trace(pole_source(inputs, SQUARED_MOMENTUM) * overlap) / 2)
    inverse = closed_form_inverse_of_invariant_matrix(inputs, SQUARED_MOMENTUM)
    trace_form = (SQUARED_MOMENTUM - closed_form_M2_aa_can(inputs)
                  + 2 * inputs["epsilon"]**2 * sp.trace(invariant_matrix(inputs) * inverse))
    return {
        "Schur complement = s - M^2_{aa,can} + 2 epsilon^2 tr[calC M_a^{-1}(calC)]": is_zero(schur_complement - trace_form),
        "trace form = explicit sums of D_a(s)": is_zero(trace_form - closed_form_alp_inverse_propagator(inputs, SQUARED_MOMENTUM)),
    }


def rational_flavor_basis():
    """Basis E_kk, E_kl + E_lk, i(E_kl - E_lk) (k < l) of Hermitian 3x3 matrices, with rational entries."""
    basis = [sp.Matrix(3, 3, lambda row, column: int(row == column == index)) for index in range(3)]
    for row in range(3):
        for column in range(row + 1, 3):
            symmetric, antisymmetric = sp.zeros(3), sp.zeros(3)
            symmetric[row, column] = symmetric[column, row] = 1
            antisymmetric[row, column], antisymmetric[column, row] = sp.I, -sp.I
            basis += [symmetric, antisymmetric]
    return basis


def step_determinant_identity(inputs):
    """Step 6 (rational data): det D(s) = -det(M2_PP - s) D_a(s) as polynomials in s, monic of degree 10.

    To avoid the square roots of lambda_0 and lambda_8, D(s) is first written
    in the rational basis of rational_flavor_basis(): with J_ik = tr(lambda_i B_k)/2
    and P = J x, D_x = diag(J, 1)^T D diag(J, 1) has rational entries and
    det D = det D_x / det(J)^2 (same for the meson block).
    """
    kinetic, mass_squared = lo_blocks(inputs)
    change_of_basis = sp.diag(sp.Matrix.hstack(*[generator_components(matrix) for matrix in rational_flavor_basis()]), 1)
    jacobian_squared = sp.radsimp(sp.expand(change_of_basis.det()**2))
    transformed = (change_of_basis.T * (SQUARED_MOMENTUM * kinetic - mass_squared) * change_of_basis).applyfunc(
        lambda entry: sp.expand(sp.radsimp(sp.expand(entry))))
    full_determinant = sp.expand(transformed.det(method="berkowitz") / jacobian_squared)
    meson_determinant = sp.expand(-transformed[:9, :9].det(method="berkowitz") / jacobian_squared)
    closed_meson_determinant = meson_block_determinant(inputs["chi_values"], inputs["anomaly_mass_squared"], SQUARED_MOMENTUM)
    polynomial_from_alp_inverse_propagator = sp.cancel(-closed_meson_determinant * closed_form_alp_inverse_propagator(inputs, SQUARED_MOMENTUM))
    full_polynomial = sp.Poly(full_determinant, SQUARED_MOMENTUM)
    return {
        "change of basis has rational entries in D_x": all(entry.is_rational or entry.is_polynomial(SQUARED_MOMENTUM)
                                                          and not entry.has(sp.sqrt(2), sp.sqrt(3)) for entry in transformed),
        "det(M2_PP - s) = prod ((chi_i+chi_j)/2 - s)^2 prod (chi_i - s) [1 + (m_0^2/N_f) sum 1/(chi_i - s)]":
            is_zero(meson_determinant - closed_meson_determinant),
        "det D(s) = -det(M2_PP - s) D_a(s)": is_zero(full_determinant - polynomial_from_alp_inverse_propagator),
        "det D(s) is monic of degree 10 (det K = 1)": full_polynomial.degree() == 10 and full_polynomial.LC() == 1,
    }


def step_sign_below_meson_poles(inputs):
    """Step 6 sign: for 0 <= s < min_i chi_i the operator M_a is positive, so the meson term of D_a is positive.

    Then D_a(s) > s - M^2_{aa,can} there, i.e. mixing with the heavier mesons lowers the
    ALP-like root. Checked at s = 0 and s = min(chi)/2 with rational data,
    together with the positivity of M_a (all eigenvalues of M2_PP - s K_PP).
    """
    _, mass_squared = lo_blocks(inputs)
    results = {}
    for squared_momentum in (sp.Integer(0), min(inputs["chi_values"]) / 2):
        meson_term = (closed_form_alp_inverse_propagator(inputs, squared_momentum) - squared_momentum
                      + closed_form_M2_aa_can(inputs))
        smallest_eigenvalue = min(np.linalg.eigvalsh(np.array((mass_squared[:9, :9] - squared_momentum * sp.eye(9)).evalf(30), dtype=float)))
        results[f"M_a positive and meson term of D_a positive at s = {squared_momentum}"] = (
            smallest_eigenvalue > 0 and sp.cancel(sp.expand(meson_term)) > 0)
    return results


def step_normalization_is_derivative(inputs):
    """Step 7: N(s) = dD_a/ds."""
    return {"N(s) = dD_a/ds": is_zero(sp.diff(closed_form_alp_inverse_propagator(inputs, SQUARED_MOMENTUM), SQUARED_MOMENTUM)
                                     - closed_form_normalization(inputs, SQUARED_MOMENTUM))}


def step_dependence_on_couplings(inputs):
    """Step 2/8c: M_a^{-1}(calC), N, D_a - s + M^2_{aa,can} depend on the couplings only via c^AP and c_GA; M^2_{aa,can} does not.

    The shift (c^A, c^P, c_G) -> (c^A + delta, c^P + delta, c_G - tr delta),
    for any Hermitian delta, keeps c^AP, c_GA and calC fixed. Under it
    M^2_{aa,can} -> M^2_{aa,can} - 2 epsilon^2 [tr(delta {c^P, chi_0}) + tr(chi_0 delta^2)].
    This change is compensated by the a^2 term that the transformed quark
    spurion acquires (see memory/alp/shared/alp-lo-bare-alp-mass-bound.md);
    the frame used in this script has no a^2 spurion term, so M^2_{aa,can} changes.
    """
    delta = symbolic_hermitian_matrix("delta")
    shifted = dict(inputs)
    shifted["axial_source"] = inputs["axial_source"] + delta
    shifted["pseudoscalar_source"] = inputs["pseudoscalar_source"] + delta
    shifted["gluon_coupling"] = inputs["gluon_coupling"] - sp.trace(delta)
    quark_spurion = sp.diag(*inputs["chi_values"])
    mass_change = closed_form_M2_aa_can(shifted) - closed_form_M2_aa_can(inputs)
    expected_change = -2 * inputs["epsilon"]**2 * (sp.trace(delta * anticommutator(inputs["pseudoscalar_source"], quark_spurion))
                                                   + sp.trace(quark_spurion * delta * delta))
    unit_shift = sp.Symbol("t", real=True)
    unit_change = mass_change.subs({entry: 0 for entry in delta.free_symbols if "_re" in entry.name or "_im" in entry.name})
    unit_change = unit_change.subs({sp.Symbol(f"delta_{index}{index}", real=True): unit_shift for index in range(3)})
    return {
        "c^AP, c_GA and calC unchanged": matrix_is_zero(axial_minus_pseudoscalar_coupling(shifted) - axial_minus_pseudoscalar_coupling(inputs))
        and is_zero(gluon_plus_axial_trace_coupling(shifted) - gluon_plus_axial_trace_coupling(inputs))
        and matrix_is_zero(invariant_matrix(shifted) - invariant_matrix(inputs)),
        "M_a^{-1}(calC) unchanged": matrix_is_zero(closed_form_inverse_of_invariant_matrix(shifted, SQUARED_MOMENTUM)
                                                   - closed_form_inverse_of_invariant_matrix(inputs, SQUARED_MOMENTUM)),
        "N unchanged": is_zero(closed_form_normalization(shifted, SQUARED_MOMENTUM) - closed_form_normalization(inputs, SQUARED_MOMENTUM)),
        "D_a(s) + M^2_{aa,can} unchanged": is_zero(closed_form_alp_inverse_propagator(shifted, SQUARED_MOMENTUM) + closed_form_M2_aa_can(shifted)
                                           - closed_form_alp_inverse_propagator(inputs, SQUARED_MOMENTUM) - closed_form_M2_aa_can(inputs)),
        "M^2_{aa,can} change = -2eps^2[tr(delta{c^P,chi}) + tr(chi delta^2)]": is_zero(mass_change - expected_change),
        "M^2_{aa,can} changes for delta = t 1, c_G -> c_G - 3t (no a^2 spurion term in this frame)": not is_zero(unit_change),
        "delta = t 1: change = -2eps^2[2t tr(c^P chi) + t^2 tr chi]": is_zero(
            unit_change + 2 * inputs["epsilon"]**2 * (2 * unit_shift * sp.trace(inputs["pseudoscalar_source"] * quark_spurion)
                                                      + unit_shift**2 * sp.trace(quark_spurion))),
        "_mass_change": sp.factor(sp.expand(unit_change)),
    }


def step_weak_mixing_limit(inputs):
    """Step 8b: W/epsilon is epsilon-independent and D_a(s)|_{epsilon=0} = s - m_a0^2.

    Hence the ALP-like pole tends to m_a0^2 and W/epsilon tends to
    M_a^{-1}(C/epsilon) at s = m_a0^2, the first-order solution of Alda:2026nrd.
    """
    epsilon = inputs["epsilon"]
    overlap = closed_form_overlap(inputs, SQUARED_MOMENTUM)
    first_order = inverse_lo_pole_operator(pole_source(inputs, SQUARED_MOMENTUM) / epsilon, inputs["chi_values"],
                                           inputs["anomaly_mass_squared"], SQUARED_MOMENTUM)
    alp_inverse_propagator_at_zero = closed_form_alp_inverse_propagator(inputs, SQUARED_MOMENTUM).subs(epsilon, 0)
    return {
        "W/epsilon has no epsilon dependence": all(not sp.cancel(entry / epsilon).has(epsilon) for entry in overlap),
        "W/epsilon = M_a^{-1}(C(s)/epsilon)": matrix_is_zero(overlap / epsilon - first_order),
        "D_a(s) at epsilon = 0 is s - m_a0^2": is_zero(alp_inverse_propagator_at_zero - SQUARED_MOMENTUM + inputs["bare_alp_mass_squared"]),
        "N(s) at epsilon = 0 is 1": is_zero(closed_form_normalization(inputs, SQUARED_MOMENTUM).subs(epsilon, 0) - 1),
    }


def step_flavor_diagonal_quartic():
    """Step 8d: flavor-diagonal sources, sector (P_3, P_8, P_0, a).

    det D(s) = -D_a(s) prod_i (chi_i - s) [1 + (m_0^2/N_f) sum_i 1/(chi_i - s)], and,
    with (i, j, k) a permutation of (u, d, s),
      det D(s) = -(s - M^2_{aa,can}) [prod_i (chi_i - s) + (m_0^2/N_f) sum_{i<j} (chi_i - s)(chi_j - s)]
                 - 2 epsilon^2 [sum_i calC_ii^2 (chi_j - s)(chi_k - s)
                                + (m_0^2/N_f) sum_{i<j} (calC_ii - calC_jj)^2 (chi_k - s)],
    a monic quartic in s.
    """
    inputs = symbolic_inputs()
    inputs["axial_source"] = sp.diag(*sp.symbols("cA_u cA_d cA_s", real=True))
    inputs["pseudoscalar_source"] = sp.diag(*sp.symbols("cP_u cP_d cP_s", real=True))
    kinetic, mass_squared = lo_blocks(inputs)
    positions = [3, 8, 0, 9]
    sector_kinetic, sector_mass = kinetic.extract(positions, positions), mass_squared.extract(positions, positions)
    determinant = sp.expand(sp.radsimp((SQUARED_MOMENTUM * sector_kinetic - sector_mass).det(method="berkowitz")))
    m0_squared_over_flavors, epsilon = inputs["anomaly_mass_squared"] / FLAVOR_COUNT, inputs["epsilon"]
    chi_minus_s = [chi - SQUARED_MOMENTUM for chi in inputs["chi_values"]]
    trace_factor = sp.expand(sp.prod(chi_minus_s) + m0_squared_over_flavors * (chi_minus_s[0] * chi_minus_s[1]
                                                                                + chi_minus_s[0] * chi_minus_s[2]
                                                                                + chi_minus_s[1] * chi_minus_s[2]))
    numerator, denominator = sp.fraction(sp.together(closed_form_alp_inverse_propagator(inputs, SQUARED_MOMENTUM)))
    quotient, remainder = sp.div(sp.expand(-numerator * trace_factor), sp.expand(denominator), SQUARED_MOMENTUM)
    calC_diagonal = [invariant_matrix(inputs)[index, index] for index in range(3)]
    triples = ((0, 1, 2), (1, 0, 2), (2, 0, 1))
    explicit_form = (-(SQUARED_MOMENTUM - closed_form_M2_aa_can(inputs)) * trace_factor
                     - 2 * epsilon**2 * (sum(calC_diagonal[i]**2 * chi_minus_s[j] * chi_minus_s[k] for i, j, k in triples)
                                         + m0_squared_over_flavors * sum((calC_diagonal[j] - calC_diagonal[k])**2 * chi_minus_s[i]
                                                                         for i, j, k in triples)))
    determinant_polynomial = sp.Poly(determinant, SQUARED_MOMENTUM)
    return {
        "-D_a(s) prod (chi_i - s) [1 + (m_0^2/N_f) sum 1/(chi_i - s)] is a polynomial": remainder == 0,
        "det D = -D_a(s) prod (chi_i - s) [1 + (m_0^2/N_f) sum 1/(chi_i - s)]": sp.expand(determinant - quotient) == 0,
        "det D = explicit quartic form": sp.expand(determinant - explicit_form) == 0,
        "degree 4": determinant_polynomial.degree() == 4,
        "leading coefficient 1": sp.expand(determinant_polynomial.LC() - 1) == 0,
    }


def step_numerical_pole_solver(inputs, tolerance=1e-12):
    """Step 8a: roots of D_a, w_i and N against solve_spectrum_numerical (K-orthonormal pole vectors).

    The roots of the polynomial -det(M2_PP - s) D_a(s) = det D(s) are compared
    with all ten solver poles. For each pair i < j only the real direction of
    c^AP_ij in the (Re, Im) plane couples to the ALP; the orthogonal
    direction is an exactly decoupled meson at s = (chi_i + chi_j)/2 with
    v_a = 0, where w_i is undefined. These poles are counted and excluded from
    the w_i and N comparison. w_i is compared relative to max_i |w_i|, N
    relative to N.
    """
    kinetic, mass_squared = lo_blocks(inputs)
    spectrum = solve_spectrum_numerical(np.array(kinetic.evalf(30), dtype=float),
                                        np.array(mass_squared.evalf(30), dtype=float))
    polynomial_from_alp_inverse_propagator = sp.Poly(sp.cancel(-meson_block_determinant(inputs["chi_values"], inputs["anomaly_mass_squared"], SQUARED_MOMENTUM)
                                           * closed_form_alp_inverse_propagator(inputs, SQUARED_MOMENTUM)), SQUARED_MOMENTUM)
    roots = sorted(sp.re(root) for root in polynomial_from_alp_inverse_propagator.nroots(n=40))
    overlap_components = generator_components(closed_form_overlap(inputs, SQUARED_MOMENTUM))
    normalization = closed_form_normalization(inputs, SQUARED_MOMENTUM)
    average_chi_values = sorted((inputs["chi_values"][i] + inputs["chi_values"][j]) / 2 for i in range(3) for j in range(i + 1, 3))
    root_errors, overlap_errors, normalization_errors, decoupled_poles = [], [], [], []
    for index, (pole, root) in enumerate(zip(spectrum.pole_squared_masses, roots)):
        root_errors.append(abs(float(root) - pole) / max(1.0, abs(pole)))
        vector = spectrum.pole_vector(index)
        if abs(vector[9]) < 1e-10:
            decoupled_poles.append(pole)
            continue
        overlap_at_pole = np.array([complex(entry.subs(SQUARED_MOMENTUM, root).evalf(30)) for entry in overlap_components])
        expected_overlap = vector[:9] / vector[9]
        # Floating-point solver: the error of v_i/v_a grows like 1/|v_a| for nearly decoupled poles.
        allowed_error = max(tolerance, 1e-14 / abs(vector[9]))
        overlap_errors.append(float(np.max(np.abs(overlap_at_pole - expected_overlap)) / np.max(np.abs(expected_overlap))) / allowed_error)
        normalization_at_pole = complex(normalization.subs(SQUARED_MOMENTUM, root).evalf(30))
        normalization_errors.append(abs(normalization_at_pole - 1 / vector[9]**2) * vector[9]**2 / allowed_error)
    return {
        "ten poles from the solver and ten roots of det D": len(roots) == 10 and len(spectrum.pole_squared_masses) == 10,
        "roots of det D = -det(M2_PP - s) D_a(s) agree with the solver poles": max(root_errors) < tolerance,
        "three decoupled poles (v_a = 0) at s = (chi_i + chi_j)/2": len(decoupled_poles) == 3
        and max(abs(pole - float(value)) for pole, value in zip(sorted(decoupled_poles), average_chi_values)) < tolerance,
        "w_i(m_n^2) agree with v_i/v_a (relative error < max(1e-12, 1e-14/|v_a|))": max(overlap_errors) < 1,
        "N(m_n^2) agrees with 1/v_a^2 (relative error < max(1e-12, 1e-14/|v_a|))": max(normalization_errors) < 1,
        "_largest pole error; largest w_i and N errors in units of the allowed error":
            (max(root_errors), max(overlap_errors), max(normalization_errors)),
        "_smallest |v_a| among compared poles": min(abs(spectrum.pole_vector(index)[9]) for index in range(10)
                                                     if abs(spectrum.pole_vector(index)[9]) >= 1e-10),
    }


def step_one_alp_amplitude():
    """Step 9: V_i w_i = epsilon tr{V [c^A + M_a^{-1}(calC)]} for w_i lambda_i = 2 epsilon [c^A + M_a^{-1}(calC)].

    Identity of components with generic V_i and generic components x_i of the
    matrix c^A + M_a^{-1}(calC).
    """
    epsilon = sp.Symbol("epsilon", positive=True)
    vertex_components = sp.Matrix(sp.symbols("V_0:9", real=True))
    bracket_components = sp.Matrix(sp.symbols("x_0:9", real=True))
    vertex_matrix, bracket_matrix = matrix_from_components(vertex_components), matrix_from_components(bracket_components)
    overlap_components = 2 * epsilon * bracket_components
    return {"V_i w_i = epsilon tr{V [c^A + M_a^{-1}(calC)]}": is_zero(
        (vertex_components.T * overlap_components)[0] - epsilon * sp.trace(vertex_matrix * bracket_matrix))}


# ------------------------------------------------------------------- driver


def symbolic_checks():
    """Steps with generic symbols for every input (c^A, c^P generic Hermitian; chi_i, m_0^2, c_G, epsilon, m_a0^2, s)."""
    inputs = symbolic_inputs()
    results = {"0 completeness": step_component_completeness()}
    for name, step in (("1 matrix form", step_pole_equation_matrix_form),
                       ("2 invariant matrix", step_canonical_source),
                       ("3 explicit inverse", step_explicit_inverse),
                       ("5 M^2_{aa,can}", step_M2_aa_can),
                       ("6 ALP inverse propagator D_a", step_alp_inverse_propagator),
                       ("2/8c coupling dependence", step_dependence_on_couplings),
                       ("8b weak mixing", step_weak_mixing_limit)):
        results[name] = step(inputs)
    return results


def partially_symbolic_checks(seed=1):
    """Steps 4 and 7, whose fully symbolic form is too slow for SymPy.

    (a) generic Hermitian c^A, c^P and symbolic c_G, epsilon, m_a0^2, s, with
        rational chi_i and m_0^2; (b) symbolic chi_i, m_0^2, c_G, epsilon,
        m_a0^2, s, with rational c^A, c^P (step 7 only).
    """
    rational = rational_inputs(seed)
    symbolic_sources = symbolic_inputs()
    symbolic_sources.update({key: rational[key] for key in ("chi_values", "anomaly_mass_squared")})
    symbolic_masses = symbolic_inputs()
    symbolic_masses.update({key: rational[key] for key in ("axial_source", "pseudoscalar_source")})
    return {
        "4 normalization (symbolic sources)": step_normalization(symbolic_sources),
        "7 N = dD_a/ds (symbolic sources)": step_normalization_is_derivative(symbolic_sources),
        "7 N = dD_a/ds (symbolic chi_i, m_0^2)": step_normalization_is_derivative(symbolic_masses),
    }


def rational_checks(seed):
    """Exact rational data: trace-system inverse, steps 4, 6, 7, determinant, numerical solver."""
    inputs = rational_inputs(seed)
    test_momentum = sp.Rational(7, 50)
    return {
        "3 trace-system inverse": step_explicit_inverse_against_trace_system(inputs, test_momentum),
        "4 normalization": step_normalization(inputs),
        "6 ALP inverse propagator D_a": step_alp_inverse_propagator(inputs),
        "6 determinant": step_determinant_identity(inputs),
        "6 sign below the meson poles": step_sign_below_meson_poles(inputs),
        "7 N = dD_a/ds": step_normalization_is_derivative(inputs),
        "8a numerical solver": step_numerical_pole_solver(inputs),
    }


def print_results(results):
    for step, checks in results.items():
        print(f"  step {step}:")
        for label, value in checks.items():
            if label.startswith("_"):
                print(f"    {label[1:]}: {value}")
            else:
                print(f"    [{'ok' if value else 'FAILED'}] {label}")


def main():
    print("LO finite ALP-meson mixing via M_a[X] = {chi_0, X}/2 + (m_0^2/N_f) tr(X) 1 - s X (N_f = 3, GeV^2 units)")
    print("calC = -{chi_0, c^AP}/2 - (m_0^2/N_f) c_GA 1, c^AP = c^A - c^P, c_GA = c_G + tr(c^A)")
    print("Step 1: M_a[W] = C(s), C = epsilon({c^P, chi_0} - 2 s c^A - 2 (m_0^2/N_f) c_G 1)")
    print("Step 2: C(s) - M_a[2 epsilon c^A] = 2 epsilon calC, so w_i lambda_i = 2 epsilon [c^A + M_a^{-1}(calC)]")
    print("Step 3: M_a^{-1}(calC)_ij = calC_ij/((chi_i + chi_j)/2 - s) (i != j),")
    print("        M_a^{-1}(calC)_ii = (calC_ii - (m_0^2/N_f) tr M_a^{-1}(calC))/(chi_i - s)")
    print("Step 4: N = 1 + 2 epsilon^2 tr[M_a^{-1}(calC)^2]")
    print("Step 5: M^2_{aa,can} = m_a0^2 + 2 epsilon^2 (m_0^2/N_f) c_GA^2 - 2 epsilon^2 [tr(c^A{c^P,chi_0}) - tr(chi_0 c^A c^A)]")
    print("Step 6: ALP inverse propagator with the mesons integrated out (Schur complement of D(s)):")
    print("        D_a(s) = s - M^2_{aa,can} + 2 epsilon^2 tr[calC M_a^{-1}(calC)], det D = -det(M2_PP - s) D_a(s)")
    print("Step 7: N = dD_a/ds")
    print("Step 8c: (c^A, c^P, c_G) -> (c^A + delta, c^P + delta, c_G - tr delta) keeps c^AP, c_GA, calC, N fixed;")
    print("         M^2_{aa,can} -> M^2_{aa,can} - 2 epsilon^2 [tr(delta {c^P, chi_0}) + tr(chi_0 delta^2)], compensated by the a^2 term")
    print("         the transformed spurion acquires (memory/alp/shared/alp-lo-bare-alp-mass-bound.md);")
    print("         this frame has no a^2 spurion term")
    print("Step 9: amplitude numerator V_a + epsilon tr{V [c^A + M_a^{-1}(calC)]}")
    print("\nSymbolic checks (generic Hermitian c^A, c^P, symbolic chi_i, m_0^2, c_G, epsilon, m_a0^2, s):")
    print_results(symbolic_checks())
    print("\nPartially symbolic checks:")
    print_results(partially_symbolic_checks())
    print("\nVerified explicit forms (calC_ii = -chi_i c^AP_ii - (m_0^2/N_f) c_GA):")
    print("  tr[calC M_a^{-1}(calC)] = sum_{i != j} ((chi_i + chi_j)/2)^2 |c^AP_ij|^2/((chi_i + chi_j)/2 - s)")
    print("                            + sum_i calC_ii^2/(chi_i - s)")
    print("                            - (m_0^2/N_f) [sum_i calC_ii/(chi_i - s)]^2/[1 + (m_0^2/N_f) sum_i 1/(chi_i - s)]")
    print("  det D(s) = -det(M2_PP - s) D_a(s),")
    print("  det(M2_PP - s) = prod_{i<j} ((chi_i + chi_j)/2 - s)^2 prod_i (chi_i - s) [1 + (m_0^2/N_f) sum_i 1/(chi_i - s)]")
    for seed in SEEDS:
        print(f"\nRational checks, seed {seed}:")
        print_results(rational_checks(seed))
    print("\nStep 8d: flavor-diagonal sector (P_3, P_8, P_0, a):")
    print_results({"8d quartic": step_flavor_diagonal_quartic()})
    print("\nStep 9:")
    print_results({"9 amplitude": step_one_alp_amplitude()})


if __name__ == "__main__":
    main()
