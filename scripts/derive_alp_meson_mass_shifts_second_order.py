#!/usr/bin/env python3
"""O(epsilon^2) shift of the LO neutral-meson poles under finite ALP mixing, in invariants.

Run:   .venv/bin/python scripts/derive_alp_meson_mass_shifts_second_order.py
Tests: .venv/bin/python -m unittest discover -s scripts -p test_derive_alp_meson_mass_shifts_second_order.py
Record: memory/alp/shared/alp-lo-meson-mass-shift-second-order.md

Conventions (Alda:2026nrd; memory/alp/shared/finite-alp-meson-mixing.md,
memory/alp/shared/alp-lo-weak-mixing-expansion.md): N_f = 3, neutral basis
(P_3, P_8, P_0, P_6, P_7, a), tr(lambda_i lambda_j) = 2 delta_ij,
D(s) = s K - M^2, epsilon = f_pi/f_a. Blocks from `closed_form_blocks` of
scripts/derive_alp_lo_neutral_blocks.py. Neutral space: Hermitian 3x3 matrices commuting
with Q, spanned by lambda_3, lambda_8, lambda_0, lambda_6, lambda_7.

Flavour-space objects (record notation in brackets):
    𝕄_s(X) = {chi_0, X}/2 + (m_0^2/N_f) tr(X) 1 - s X          [\\mathbb M_s]
    𝒞 = -{chi_0, c_q^AP}/2 - (m_0^2/N_f) c_GA 1                [\\mathcal C, Alda:2026nrd tadpole source]
    Sigma(s) = -s c_q^A + {c_q^P, chi_0}/2 - (m_0^2 c_G/N_f) 1
    lambda_n: neutral, 𝕄_0(lambda_n) = m_n0^2 lambda_n, tr(lambda_n lambda_m) = 2 delta_nm.
c_q^AP = c_q^A - c_q^P, c_GA = c_G + tr c_q^A.

Derived and checked here (numbering = steps of the record):
 1. (M^2_PP)_ik = tr(lambda_i 𝕄_0(lambda_k))/2, i.e. M^2_PP is the matrix of 𝕄_0 in the basis
    {lambda_i}, orthonormal under <X, Y> = tr(X Y)/2. 𝕄_0 maps the neutral space into itself and
    tr(Y 𝕄_s(X)) = tr(𝕄_s(Y) X) (symbolic). lambda_n = sum_i v_i lambda_i built from a unit
    eigenvector v of M^2_PP obeys the eigen-equation and tr(lambda_n lambda_m) = 2 delta_nm
    (50 digits).
 2. s K_ia - M^2_ia = epsilon tr(lambda_i Sigma(s)) (symbolic); in the mass basis
    D_na(s) = epsilon tr(lambda_n Sigma(s)) (50 digits, from the blocks).
 3. Sigma(s) = 𝕄_s(c_q^A) + 𝒞 as a matrix identity (symbolic), hence
    D_na(s) = epsilon [(m_n0^2 - s) tr(lambda_n c_q^A) + tr(lambda_n 𝒞)] (50 digits).
 4. D_na(m_n0^2) = epsilon tr(lambda_n 𝒞) exactly (50 digits).
 5. m_n^2 - m_n0^2 = epsilon^2 tr(lambda_n 𝒞)^2 / (m_n0^2 - m_a0^2) + O(epsilon^4); for a degenerate
    eigenvalue (K^0 doublet) the sum of tr(lambda_n 𝒞)^2 over the eigenspace enters for the one coupled
    combination; the others stay at m_n0^2 exactly (solver comparison: residual ~ epsilon^4).
 6. Constant axial redefinition c_q^{A,P} -> c_q^{A,P} + kappa, c_G -> c_G - tr kappa:
    Sigma(s) -> Sigma(s) + 𝕄_s(kappa), 𝒞 unchanged (symbolic).
 7. Y solving 𝕄_{m_a0^2}(Y) = 𝒞 (Alda:2026nrd first-order equation; solved as a matrix equation
    on the neutral space, exact rationals): tr(lambda_n 𝒞) = (m_n0^2 - m_a0^2) tr(lambda_n Y), and
    epsilon tr(lambda_n (c_q^A + Y)) = D_na(m_a0^2)/(m_n0^2 - m_a0^2) (50 digits).

Inputs of the numerical check: illustrative_lo_inputs("neutral") with the bare ALP mass
replaced by m_a0^2 = 3/5 GeV^2 (test value away from the eta and K^0 masses) and
f_a = f_pi/epsilon, epsilon in EPSILON_VALUES. tr S_2 is set to zero; it does not enter
the O(epsilon^2) meson shift (it is O(epsilon^2) in M^2_aa, so it contributes at O(epsilon^4)).

Units: f_pi, f_a in GeV; chi_0, m_0^2, m_a0^2, M^2, s in GeV^2; couplings and epsilon dimensionless.
"""

from __future__ import annotations

import mpmath
import numpy as np
import sympy as sp

from derive_alp_lo_mass_sum_rules import _hermitian_charge_conserving, symbolic_sources
from derive_alp_lo_neutral_blocks import (
    FLAVOR_COUNT,
    NEUTRAL_FIELD_LABELS,
    NEUTRAL_GENERATOR_INDICES,
    closed_form_blocks,
    illustrative_lo_inputs,
)
from solve_alp_meson_spectrum import solve_spectrum_numerical
from symbolic_utils import anticommutator, u3_generators
from verify_alp_weak_mixing_limit import alda_first_order_operator

MESON_COUNT = len(NEUTRAL_GENERATOR_INDICES)
TEST_BARE_ALP_MASS_SQUARED = sp.Rational(3, 5)  # GeV^2, test value
EPSILON_VALUES = tuple(sp.Rational(1, denominator) for denominator in (20, 40, 80, 160, 320))
NEUTRAL_GENERATORS = tuple(u3_generators()[index] for index in NEUTRAL_GENERATOR_INDICES)
CHARGE_MATRIX = sp.diag(sp.Rational(2, 3), sp.Rational(-1, 3), sp.Rational(-1, 3))
HIGH_PRECISION_DIGITS = 50
HIGH_PRECISION_TOLERANCE = mpmath.mpf("1e-40")


def _simplify(expression):
    return sp.simplify(sp.expand(sp.radsimp(sp.expand(expression))))


def _is_zero_vector(vector):
    return all(_simplify(entry) == 0 for entry in vector)


# ------------------------------------------------------------------ building blocks


def neutral_components(matrix):
    """Column (tr(lambda_i X)) over the neutral generators (P_3, P_8, P_0, P_6, P_7)."""
    return sp.Matrix([sp.expand(sp.trace(generator * matrix)) for generator in NEUTRAL_GENERATORS])


def from_neutral_components(components):
    """X = (1/2) sum_i x_i lambda_i, the inverse of `neutral_components` on the neutral space."""
    return sum((component * generator / 2 for component, generator in zip(components, NEUTRAL_GENERATORS)),
               sp.zeros(3))


def meson_mass_operator(argument, quark_spurion, anomaly_mass_squared, squared_momentum=0):
    """𝕄_s(X) = {chi_0, X}/2 + (m_0^2/N_f) tr(X) 1 - s X; s = 0 gives the unmixed meson operator 𝕄_0."""
    return (anticommutator(quark_spurion, argument) / 2
            + anomaly_mass_squared * sp.trace(argument) * sp.eye(FLAVOR_COUNT) / FLAVOR_COUNT
            - squared_momentum * argument)


def mixing_source_matrix(squared_momentum, axial_source, pseudoscalar_source, quark_spurion, gluon_coupling,
                         anomaly_mass_squared):
    """Sigma(s) = -s c_q^A + {c_q^P, chi_0}/2 - (m_0^2 c_G/N_f) 1, with D_ia(s) = epsilon tr(lambda_i Sigma(s))."""
    return (-squared_momentum * axial_source + anticommutator(pseudoscalar_source, quark_spurion) / 2
            - anomaly_mass_squared * gluon_coupling * sp.eye(FLAVOR_COUNT) / FLAVOR_COUNT)


def tadpole_source(axial_minus_pseudoscalar, gluon_axial_combination, quark_spurion, anomaly_mass_squared):
    """𝒞 = -{chi_0, c_q^AP}/2 - (m_0^2/N_f) c_GA 1 (Alda:2026nrd), written in the invariants only."""
    return (-anticommutator(quark_spurion, axial_minus_pseudoscalar) / 2
            - anomaly_mass_squared * gluon_axial_combination * sp.eye(FLAVOR_COUNT) / FLAVOR_COUNT)


def tadpole_source_from_sources(axial_source, pseudoscalar_source, quark_spurion, gluon_coupling,
                                anomaly_mass_squared):
    """𝒞 with c_q^AP = c_q^A - c_q^P and c_GA = c_G + tr c_q^A."""
    return tadpole_source(axial_source - pseudoscalar_source, gluon_coupling + sp.trace(axial_source),
                          quark_spurion, anomaly_mass_squared)


def invariant_mixing_vector(axial_source, pseudoscalar_source, quark_spurion, gluon_coupling, anomaly_mass_squared):
    """Column (tr(lambda_i 𝒞)) over the neutral generators."""
    return neutral_components(tadpole_source_from_sources(axial_source, pseudoscalar_source, quark_spurion,
                                                          gluon_coupling, anomaly_mass_squared))


def solve_neutral_matrix_equation(linear_map, right_side):
    """The neutral matrix Y (Hermitian, commuting with Q) with linear_map(Y) = right_side, solved in flavour space."""
    diagonal = sp.symbols("y_u y_d y_s", real=True)
    real_part, imaginary_part = sp.symbols("y_re y_im", real=True)
    unknown = sp.diag(*diagonal)
    unknown[1, 2] = real_part + sp.I * imaginary_part
    unknown[2, 1] = real_part - sp.I * imaginary_part
    residual = (linear_map(unknown) - right_side).applyfunc(sp.expand)
    equations = [residual[0, 0], residual[1, 1], residual[2, 2], residual[1, 2], residual[2, 1]]
    unknowns = (*diagonal, real_part, imaginary_part)
    (solution,) = sp.linsolve(equations, unknowns)
    return unknown.subs(dict(zip(unknowns, solution))).applyfunc(sp.cancel)


def _block_inputs(sources):
    return {key: sources[key] for key in ("axial_source", "pseudoscalar_source", "quark_spurion",
                                          "gluon_coupling", "anomaly_mass_squared")}


# ------------------------------------------------------------------ high-precision mass basis


def to_high_precision(matrix, digits=HIGH_PRECISION_DIGITS):
    """sympy matrix -> mpmath matrix (complex entries) at the given number of digits."""
    rows, columns = matrix.shape
    result = mpmath.matrix(rows, columns)
    for row in range(rows):
        for column in range(columns):
            value = sp.N(matrix[row, column], digits + 10)
            result[row, column] = mpmath.mpc(mpmath.mpf(str(sp.re(value))), mpmath.mpf(str(sp.im(value))))
    return result


def high_precision_trace(matrix):
    return mpmath.fsum(matrix[index, index] for index in range(matrix.rows))


def high_precision_meson_mass_operator(argument, quark_spurion, anomaly_mass_squared, squared_momentum=0):
    """𝕄_s(X) for mpmath matrices."""
    identity = mpmath.eye(FLAVOR_COUNT)
    return ((quark_spurion * argument + argument * quark_spurion) / 2
            + anomaly_mass_squared * high_precision_trace(argument) * identity / FLAVOR_COUNT
            - squared_momentum * argument)


def mass_eigenstate_generators(sources, digits=HIGH_PRECISION_DIGITS):
    """[(m_n0^2, lambda_n)], lambda_n = sum_i v_i lambda_i with v the unit eigenvectors of M^2_PP.

    Computed with mpmath at `digits` digits; lambda_n are mpmath 3x3 matrices.
    """
    mpmath.mp.dps = digits + 10
    _, mass_squared = closed_form_blocks(**sources)
    eigenvalues, eigenvectors = mpmath.eigsy(mpmath.matrix(mass_squared[:MESON_COUNT, :MESON_COUNT].evalf(digits + 10)))
    generators = [to_high_precision(sp.Matrix(generator), digits) for generator in NEUTRAL_GENERATORS]
    eigenbasis = []
    for column in range(MESON_COUNT):
        generator = mpmath.matrix(FLAVOR_COUNT, FLAVOR_COUNT)
        for row in range(MESON_COUNT):
            generator += eigenvectors[row, column] * generators[row]
        eigenbasis.append((eigenvalues[column], generator))
    return eigenbasis


def _largest_entry(matrix):
    return max(abs(matrix[row, column]) for row in range(matrix.rows) for column in range(matrix.cols))


# ------------------------------------------------------------------ symbolic checks


def symbolic_checks():
    """Steps 1, 2, 3, 6 with symbolic Q-commuting Hermitian sources (chi_0 with a complex d-s entry)."""
    sources = symbolic_sources(diagonal_quark_spurion=False)
    squared_momentum = sp.Symbol("s")
    ratio = sources["pion_decay_constant"] / sources["alp_decay_constant"]
    kinetic, mass_squared = closed_form_blocks(**sources)
    meson_mass_squared = mass_squared[:MESON_COUNT, :MESON_COUNT]
    mixed_entry = squared_momentum * kinetic[:MESON_COUNT, MESON_COUNT] - mass_squared[:MESON_COUNT, MESON_COUNT]
    inputs = _block_inputs(sources)
    quark_spurion, anomaly_mass_squared = sources["quark_spurion"], sources["anomaly_mass_squared"]
    reduced_entry = neutral_components(mixing_source_matrix(squared_momentum, **inputs))
    tadpole = tadpole_source_from_sources(**inputs)

    # step 1: matrix of 𝕄_0 in the basis {lambda_i}
    operator_matrix = sp.Matrix(MESON_COUNT, MESON_COUNT, lambda row, column: sp.trace(
        NEUTRAL_GENERATORS[row] * meson_mass_operator(sp.Matrix(NEUTRAL_GENERATORS[column]), quark_spurion,
                                                      anomaly_mass_squared)) / 2)
    meson_components = sp.Matrix(sp.symbols("x3 x8 x0 x6 x7", real=True))
    operator_half = neutral_components(meson_mass_operator(from_neutral_components(meson_components),
                                                           quark_spurion, anomaly_mass_squared))
    generic_left = sp.Matrix(3, 3, sp.symbols("X0:9"))
    generic_right = sp.Matrix(3, 3, sp.symbols("Y0:9"))
    self_adjoint = sp.expand(
        sp.trace(generic_right * meson_mass_operator(generic_left, quark_spurion, anomaly_mass_squared,
                                                     squared_momentum))
        - sp.trace(meson_mass_operator(generic_right, quark_spurion, anomaly_mass_squared, squared_momentum)
                   * generic_left))
    neutral_argument = _hermitian_charge_conserving("X")
    image = meson_mass_operator(neutral_argument, quark_spurion, anomaly_mass_squared)
    orthonormal = all(sp.simplify(sp.trace(NEUTRAL_GENERATORS[row] * NEUTRAL_GENERATORS[column])
                                  - 2 * (row == column)) == 0
                      for row in range(MESON_COUNT) for column in range(MESON_COUNT))

    # step 3: Sigma(s) = 𝕄_s(c_q^A) + 𝒞; after c_q^P = c_q^A - c_q^AP, c_G = c_GA - tr c_q^A no c_q^A is left in 𝒞
    split_residual = (mixing_source_matrix(squared_momentum, **inputs)
                      - meson_mass_operator(sources["axial_source"], quark_spurion, anomaly_mass_squared,
                                            squared_momentum)
                      - tadpole)
    difference_source = _hermitian_charge_conserving("cAP")
    gluon_axial_symbol = sp.Symbol("c_GA", real=True)
    substituted = dict(inputs, pseudoscalar_source=sources["axial_source"] - difference_source,
                       gluon_coupling=gluon_axial_symbol - sp.trace(sources["axial_source"]))
    substituted_tadpole = tadpole_source_from_sources(**substituted).applyfunc(_simplify)
    axial_symbols = set(sources["axial_source"].free_symbols)

    # step 6: constant axial redefinition
    axial_shift = _hermitian_charge_conserving("kappa")
    shifted = dict(inputs, axial_source=inputs["axial_source"] + axial_shift,
                   pseudoscalar_source=inputs["pseudoscalar_source"] + axial_shift,
                   gluon_coupling=inputs["gluon_coupling"] - sp.trace(axial_shift))

    return {
        "1. tr(lambda_i lambda_k) = 2 delta_ik on the neutral generators": orthonormal,
        "1. (M^2_PP)_ik = tr(lambda_i 𝕄_0(lambda_k))/2": _is_zero_vector(operator_matrix - meson_mass_squared),
        "1. equivalently tr(lambda_i 𝕄_0(X)) = (M^2_PP x)_i for X = (1/2) sum x_k lambda_k": _is_zero_vector(
            operator_half - meson_mass_squared * meson_components),
        "1. tr(Y 𝕄_s(X)) = tr(𝕄_s(Y) X) for generic complex 3x3 X, Y": self_adjoint == 0,
        "1. 𝕄_0 maps the neutral space into itself (Hermitian, commutes with Q)": (
            _is_zero_vector(image - image.H) and _is_zero_vector(image * CHARGE_MATRIX - CHARGE_MATRIX * image)),
        "2. s K_ia - M^2_ia = epsilon tr(lambda_i Sigma(s))": _is_zero_vector(mixed_entry - ratio * reduced_entry),
        "3. Sigma(s) = 𝕄_s(c_q^A) + 𝒞 (matrix identity)": _is_zero_vector(split_residual),
        "3. 𝒞 contains no c_q^A after substituting c_q^AP, c_GA": not (
            set().union(*(entry.free_symbols for entry in substituted_tadpole)) & axial_symbols),
        "3. that 𝒞 equals -{chi_0, c_q^AP}/2 - (m_0^2/N_f) c_GA 1": _is_zero_vector(
            substituted_tadpole - tadpole_source(difference_source, gluon_axial_symbol, quark_spurion,
                                                 anomaly_mass_squared)),
        "6. Sigma(s) -> Sigma(s) + 𝕄_s(kappa) under kappa": _is_zero_vector(
            mixing_source_matrix(squared_momentum, **shifted) - mixing_source_matrix(squared_momentum, **inputs)
            - meson_mass_operator(axial_shift, quark_spurion, anomaly_mass_squared, squared_momentum)),
        "6. 𝒞 unchanged under kappa": _is_zero_vector(tadpole_source_from_sources(**shifted) - tadpole),
    }


def first_order_equation_checks(sources):
    """Step 7 with exact rational sources: Y from 𝕄_{m_a0^2}(Y) = 𝒞 solved in flavour space."""
    bare_alp_mass_squared = sources["bare_alp_mass_squared"]
    quark_spurion, anomaly_mass_squared = sources["quark_spurion"], sources["anomaly_mass_squared"]
    tadpole = tadpole_source_from_sources(**_block_inputs(sources))
    first_order_matrix = solve_neutral_matrix_equation(
        lambda argument: meson_mass_operator(argument, quark_spurion, anomaly_mass_squared, bare_alp_mass_squared),
        tadpole)
    _, mass_squared = closed_form_blocks(**sources)
    meson_mass_squared = mass_squared[:MESON_COUNT, :MESON_COUNT]
    component_solution = from_neutral_components(
        (meson_mass_squared - bare_alp_mass_squared * sp.eye(MESON_COUNT)).LUsolve(
            invariant_mixing_vector(**_block_inputs(sources))))
    residual = (alda_first_order_operator(first_order_matrix, quark_spurion, anomaly_mass_squared,
                                          bare_alp_mass_squared) - tadpole)
    return {
        "7. Y solves the Alda:2026nrd first-order equation 𝕄_{m_a0^2}(Y) = 𝒞": _is_zero_vector(residual),
        "7. same Y as from the generator-basis matrix M^2_PP - m_a0^2": _is_zero_vector(
            first_order_matrix - component_solution),
    }, first_order_matrix


def mass_basis_checks(sources, first_order_matrix, digits=HIGH_PRECISION_DIGITS):
    """Steps 1-4 and 7 in the mass basis with lambda_n at 50 digits (largest differences < 1e-40)."""
    eigenbasis = mass_eigenstate_generators(sources, digits)
    kinetic, mass_squared = closed_form_blocks(**sources)
    epsilon = sources["pion_decay_constant"] / sources["alp_decay_constant"]
    quark_spurion = to_high_precision(sources["quark_spurion"], digits)
    anomaly_mass_squared = mpmath.mpf(str(sp.N(sources["anomaly_mass_squared"], digits + 10)))
    bare_alp_mass_squared = mpmath.mpf(str(sp.N(sources["bare_alp_mass_squared"], digits + 10)))
    epsilon_value = mpmath.mpf(str(sp.N(epsilon, digits + 10)))
    axial = to_high_precision(sources["axial_source"], digits)
    tadpole = to_high_precision(tadpole_source_from_sources(**_block_inputs(sources)), digits)
    first_order = to_high_precision(first_order_matrix, digits)
    squared_momentum = sp.Symbol("s")
    sigma_symbolic = mixing_source_matrix(squared_momentum, **_block_inputs(sources))
    kinetic_column = [mpmath.mpf(str(sp.N(kinetic[row, MESON_COUNT], digits + 10))) for row in range(MESON_COUNT)]
    mass_column = [mpmath.mpf(str(sp.N(mass_squared[row, MESON_COUNT], digits + 10))) for row in range(MESON_COUNT)]
    _, eigenvectors = mpmath.eigsy(mpmath.matrix(mass_squared[:MESON_COUNT, :MESON_COUNT].evalf(digits + 10)))

    def sigma_at(value):
        constant = to_high_precision(sigma_symbolic.subs(squared_momentum, 0), digits)
        return constant - value * axial

    differences = {"eigen": 0, "orthonormal": 0, "mixed entry": 0, "split": 0, "unmixed mass": 0,
                   "first order": 0, "pole vector": 0}
    for row, (eigenvalue_row, generator_row) in enumerate(eigenbasis):
        differences["eigen"] = max(differences["eigen"], _largest_entry(
            high_precision_meson_mass_operator(generator_row, quark_spurion, anomaly_mass_squared)
            - eigenvalue_row * generator_row))
        for column, (_, generator_column) in enumerate(eigenbasis):
            differences["orthonormal"] = max(differences["orthonormal"], abs(
                high_precision_trace(generator_row * generator_column) - 2 * (row == column)))
        trace_axial = high_precision_trace(generator_row * axial)
        trace_tadpole = high_precision_trace(generator_row * tadpole)
        for value in (eigenvalue_row, bare_alp_mass_squared, mpmath.mpf(1) / 3):
            from_blocks = mpmath.fsum(eigenvectors[index, row] * (value * kinetic_column[index] - mass_column[index])
                                      for index in range(MESON_COUNT))
            from_trace = epsilon_value * high_precision_trace(generator_row * sigma_at(value))
            differences["mixed entry"] = max(differences["mixed entry"], abs(from_blocks - from_trace))
            differences["split"] = max(differences["split"], abs(
                from_trace - epsilon_value * ((eigenvalue_row - value) * trace_axial + trace_tadpole)))
        differences["unmixed mass"] = max(differences["unmixed mass"], abs(
            epsilon_value * high_precision_trace(generator_row * sigma_at(eigenvalue_row))
            - epsilon_value * trace_tadpole))
        differences["first order"] = max(differences["first order"], abs(
            trace_tadpole - (eigenvalue_row - bare_alp_mass_squared) * high_precision_trace(generator_row * first_order)))
        differences["pole vector"] = max(differences["pole vector"], abs(
            epsilon_value * high_precision_trace(generator_row * (axial + first_order))
            - epsilon_value * high_precision_trace(generator_row * sigma_at(bare_alp_mass_squared))
            / (eigenvalue_row - bare_alp_mass_squared)))
    labels = {
        "eigen": "1. 𝕄_0(lambda_n) = m_n0^2 lambda_n for lambda_n = sum_i v_i lambda_i (v unit eigenvector of M^2_PP)",
        "orthonormal": "1. tr(lambda_n lambda_m) = 2 delta_nm",
        "mixed entry": "2. v . (s K - M^2)_Pa = epsilon tr(lambda_n Sigma(s)) at s = m_n0^2, m_a0^2, 1/3",
        "split": "3. tr(lambda_n Sigma(s)) = (m_n0^2 - s) tr(lambda_n c_q^A) + tr(lambda_n 𝒞)",
        "unmixed mass": "4. D_na(m_n0^2) = epsilon tr(lambda_n 𝒞)",
        "first order": "7. tr(lambda_n 𝒞) = (m_n0^2 - m_a0^2) tr(lambda_n Y)",
        "pole vector": "7. epsilon tr(lambda_n (c_q^A + Y)) = D_na(m_a0^2)/(m_n0^2 - m_a0^2)",
    }
    return {labels[key] + " (50 digits)": value < HIGH_PRECISION_TOLERANCE for key, value in differences.items()}, \
        max(differences.values())


# ------------------------------------------------------------------ numerical prediction


def mass_basis_groups(meson_mass_squared, relative_tolerance=1e-9):
    """Eigenvalues and eigenvector groups of M^2_PP (float); degenerate eigenvalues grouped."""
    eigenvalues, eigenvectors = np.linalg.eigh(meson_mass_squared)
    groups = []
    for index, value in enumerate(eigenvalues):
        if groups and abs(value - groups[-1][0]) <= relative_tolerance * max(1.0, abs(value)):
            groups[-1][1].append(eigenvectors[:, index])
        else:
            groups.append([value, [eigenvectors[:, index]]])
    return [(value, np.array(vectors).T) for value, vectors in groups]


def second_order_meson_poles(sources, epsilon):
    """Predicted meson poles to O(epsilon^2): list of (m_n0^2, predicted m_n^2) with multiplicity.

    For an eigenspace of dimension g, one pole is m_n0^2 + epsilon^2 sum_n tr(lambda_n 𝒞)^2 / (m_n0^2 - m_a0^2),
    the sum running over an orthonormal basis of the eigenspace (computed as the squared norm of the
    projection of (tr(lambda_i 𝒞)) on the eigenvectors), and g - 1 poles stay at m_n0^2.
    """
    _, mass_squared = closed_form_blocks(**sources)
    meson_mass_squared = np.array(mass_squared[:MESON_COUNT, :MESON_COUNT].evalf(30), dtype=float)
    invariant_vector = np.array(invariant_mixing_vector(**_block_inputs(sources)).evalf(30), dtype=float).ravel()
    bare_alp_mass_squared = float(sources["bare_alp_mass_squared"])
    predictions = []
    for unmixed, vectors in mass_basis_groups(meson_mass_squared):
        projected_norm_squared = float(np.sum((vectors.T @ invariant_vector)**2))
        predictions.append((unmixed, unmixed + float(epsilon)**2 * projected_norm_squared
                            / (unmixed - bare_alp_mass_squared)))
        predictions.extend((unmixed, unmixed) for _ in range(vectors.shape[1] - 1))
    return predictions


def rational_test_sources(epsilon, base_sources=None):
    """illustrative_lo_inputs('neutral') with m_a0^2 = 3/5 GeV^2 and f_a = f_pi/epsilon."""
    sources = dict(base_sources or illustrative_lo_inputs("neutral"))
    sources["bare_alp_mass_squared"] = TEST_BARE_ALP_MASS_SQUARED
    sources["alp_decay_constant"] = sources["pion_decay_constant"] / epsilon
    return sources


def solver_residuals(base_sources=None, epsilon_values=EPSILON_VALUES):
    """Largest |solver pole - prediction| over the meson poles, for the O(eps^2) formula and for m_n0^2."""
    second_order_residuals, zeroth_order_residuals = [], []
    for epsilon in epsilon_values:
        sources = rational_test_sources(epsilon, base_sources)
        kinetic, mass_squared = closed_form_blocks(**sources)
        spectrum = solve_spectrum_numerical(np.array(kinetic.evalf(30), dtype=float),
                                            np.array(mass_squared.evalf(30), dtype=float),
                                            field_labels=NEUTRAL_FIELD_LABELS)
        remaining = list(spectrum.pole_squared_masses)
        second_order, zeroth_order = 0.0, 0.0
        for unmixed, predicted in sorted(second_order_meson_poles(sources, epsilon), key=lambda pair: pair[1]):
            closest = min(range(len(remaining)), key=lambda index: abs(remaining[index] - predicted))
            pole = remaining.pop(closest)
            second_order = max(second_order, abs(pole - predicted))
            zeroth_order = max(zeroth_order, abs(pole - unmixed))
        second_order_residuals.append(second_order)
        zeroth_order_residuals.append(zeroth_order)
    return np.array(second_order_residuals), np.array(zeroth_order_residuals)


def fitted_exponent(epsilon_values, residuals):
    """Slope of log(residual) against log(epsilon)."""
    return float(np.polyfit(np.log([float(value) for value in epsilon_values]), np.log(residuals), 1)[0])


def numerical_checks(base_sources=None, label="neutral", exponent_tolerance=0.15):
    """Step 5: fitted exponents of the solver residuals."""
    second_order, zeroth_order = solver_residuals(base_sources)
    second_exponent = fitted_exponent(EPSILON_VALUES, second_order)
    zeroth_exponent = fitted_exponent(EPSILON_VALUES, zeroth_order)
    return {
        "checks": {
            f"5. {label}: residual of O(eps^2) formula ~ eps^4": abs(second_exponent - 4) < exponent_tolerance,
            f"5. {label}: residual of m_n0^2 ~ eps^2": abs(zeroth_exponent - 2) < exponent_tolerance,
        },
        "second_order_exponent": second_exponent,
        "zeroth_order_exponent": zeroth_exponent,
        "second_order_residuals": second_order,
        "zeroth_order_residuals": zeroth_order,
    }


def derive():
    """Run all checks."""
    neutral = numerical_checks(None, "neutral sources")
    diagonal = numerical_checks(illustrative_lo_inputs("flavor_diagonal"), "flavour-diagonal sources")
    test_sources = rational_test_sources(EPSILON_VALUES[0])
    first_order, first_order_matrix = first_order_equation_checks(test_sources)
    mass_basis, largest = mass_basis_checks(test_sources, first_order_matrix)
    return {
        "symbolic": {"checks": symbolic_checks()},
        "first-order equation": {"checks": first_order},
        "mass basis": {"checks": mass_basis, "largest difference": largest},
        "solver (neutral)": neutral,
        "solver (flavour diagonal)": diagonal,
    }


def main():
    result = derive()
    print("m_n^2 - m_n0^2 = epsilon^2 tr(lambda_n 𝒞)^2 / (m_n0^2 - m_a0^2) + O(epsilon^4),")
    print("  𝒞 = -{chi_0, c_q^AP}/2 - (m_0^2/N_f) c_GA 1,  𝕄_0(lambda_n) = m_n0^2 lambda_n,"
          "  tr(lambda_n lambda_m) = 2 delta_nm.")
    for section, content in result.items():
        for label, passed in content["checks"].items():
            print(f"[{'ok' if passed else 'FAIL'}] {section}: {label}")
    print(f"mass basis: largest difference {mpmath.nstr(result['mass basis']['largest difference'], 3)}")
    for key in ("solver (neutral)", "solver (flavour diagonal)"):
        content = result[key]
        print(f"{key}: epsilon = {[str(value) for value in EPSILON_VALUES]}")
        print(f"  residual O(eps^2) formula: {np.array2string(content['second_order_residuals'], precision=3)}"
              f"  fitted exponent {content['second_order_exponent']:.3f}")
        print(f"  residual m_n0^2:           {np.array2string(content['zeroth_order_residuals'], precision=3)}"
              f"  fitted exponent {content['zeroth_order_exponent']:.3f}")
    passed = all(all(content["checks"].values()) for content in result.values())
    print("All checks pass" if passed else "Some checks FAIL")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
