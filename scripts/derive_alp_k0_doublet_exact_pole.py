#!/usr/bin/env python3
"""K^0 doublet at diagonal chi_0: exact pole at m_K0,0^2, the coupled combination, its contribution to the secular
sum, and the CP-conserving, aligned-phase and CP-violating cases. No expansion in epsilon = f_pi/f_a.

Run:   .venv/bin/python scripts/derive_alp_k0_doublet_exact_pole.py
Tests: .venv/bin/python -m unittest discover -s scripts -p test_derive_alp_k0_doublet_exact_pole.py
Record: memory/alp/shared/alp-lo-k0-doublet-exact-pole.md

Conventions (Alda:2026nrd; memory/alp/shared/alp-lo-exact-pole-equations.md,
memory/alp/shared/alp-lo-neutral-block-reduction.md): chi_0 = diag(x_u, x_d, x_s), m_K0,0^2 = (x_d + x_s)/2,
𝒞 = -{chi_0, c_q^AP}/2 - (m_0^2/N_f) c_GA 1, 𝕄_s(X) = {chi_0, X}/2 + (m_0^2/N_f) tr(X) 1 - s X.
The d-s part of a matrix X is the matrix with only the entries X_ds and X_sd = conj(X_ds).

Derived and checked here:
 1. (symbolic) 𝒞_ds = -m_K0,0^2 (c_q^AP)_ds; 𝕄_s acts on d-s parts as multiplication by m_K0,0^2 - s;
    tr(lambda_6 𝒞)^2 + tr(lambda_7 𝒞)^2 = 4 |𝒞_ds|^2 = 4 m_K0,0^4 |(c_q^AP)_ds|^2; the doublet contribution to
    the secular sum, -2 tr(𝒞_ds-part 𝕄_s^{-1}(𝒞_ds-part)), equals 4 m_K0,0^4 |(c_q^AP)_ds|^2/(s - m_K0,0^2).
 2. (symbolic) Coupled and decoupled generators: the d-s matrices with entry 𝒞_ds/|𝒞_ds| and i 𝒞_ds/|𝒞_ds| at (d, s)
    are orthonormal (tr = 2 delta), eigen-generators of 𝕄_0 with eigenvalue m_K0,0^2, and
    tr(lambda 𝒞) = 2 |𝒞_ds| and 0 respectively. The null vector (-D_7a, D_6a) of D(m_K0,0^2) (neutral-block
    record) is proportional to the decoupled generator.
 3. (exact rationals, epsilon = 1/2 and 2) For CP-conserving, aligned-phase, CP-violating and d-s-free sources:
    (s - m_K0,0^2) divides det(s K - M^2) with multiplicity 1 (2 if (c_q^AP)_ds = 0); the exact null space of
    D(m_K0,0^2) has that dimension and zero ALP components; solve_spectrum_numerical has a pole at m_K0,0^2 with
    ALP component below 1e-10.
 4. CP-conserving (real d-s entries): K_7a = M^2_7a = 0 and tr(lambda_7 𝒞) = 0, so P_7 decouples.
    Aligned phases: the rephasing V = diag(1, 1, w), |w| = 1, makes the d-s entries real and leaves det(s K - M^2)
    unchanged. CP-violating: K_6a M^2_7a - K_7a M^2_6a != 0 (no bare-frame decoupled field), but the canonical
    secular polynomial depends on the d-s sources only through |(c_q^AP)_ds|^2 at fixed M_can,aa: rephasing 𝒞
    leaves det(s - M_can^2) unchanged.

Units: chi_0, m_0^2, masses squared in GeV^2; couplings and epsilon dimensionless. Inputs are test values.
"""

from __future__ import annotations

import numpy as np
import sympy as sp

from derive_alp_lo_neutral_blocks import NEUTRAL_FIELD_LABELS, illustrative_lo_inputs
from derive_alp_meson_mass_shifts_exact import (
    ALP,
    SQUARED_MOMENTUM,
    TEST_CONDENSATE,
    TEST_SECOND_ORDER_SCALAR_TRACE,
    bare_blocks,
    canonical_mass_matrix,
    canonical_mass_matrix_from_sources,
    inverse_meson_operator,
)
from derive_alp_meson_mass_shifts_second_order import (
    MESON_COUNT,
    NEUTRAL_GENERATORS,
    _block_inputs,
    meson_mass_operator,
    mixing_source_matrix,
    neutral_components,
    rational_test_sources,
    tadpole_source,
    tadpole_source_from_sources,
)
from solve_alp_meson_spectrum import inverse_propagator, solve_spectrum_numerical

EPSILON_VALUES = (sp.Rational(1, 2), sp.Integer(2))
ALIGNED_PHASE = (3 + 4 * sp.I) / 5  # |w| = 1, test value
LAMBDA_6, LAMBDA_7 = sp.Matrix(NEUTRAL_GENERATORS[3]), sp.Matrix(NEUTRAL_GENERATORS[4])
GENERATOR_INDEX_6, GENERATOR_INDEX_7 = 3, 4  # positions of P_6, P_7 in (P_3, P_8, P_0, P_6, P_7)


def _is_zero(expression):
    entries = expression if isinstance(expression, sp.MatrixBase) else [expression]
    return all(sp.simplify(sp.expand(entry)) == 0 for entry in entries)


def strange_down_part(matrix):
    """The matrix with only the (d, s) and (s, d) entries of `matrix`."""
    result = sp.zeros(3)
    result[1, 2], result[2, 1] = matrix[1, 2], matrix[2, 1]
    return result


def doublet_generators(tadpole_entry):
    """(coupled, decoupled): d-s matrices with (d, s) entry 𝒞_ds/|𝒞_ds| and i 𝒞_ds/|𝒞_ds|."""
    modulus = sp.sqrt(sp.expand(tadpole_entry * sp.conjugate(tadpole_entry)))
    coupled, decoupled = sp.zeros(3), sp.zeros(3)
    coupled[1, 2], coupled[2, 1] = tadpole_entry / modulus, sp.conjugate(tadpole_entry) / modulus
    decoupled[1, 2], decoupled[2, 1] = sp.I * tadpole_entry / modulus, -sp.I * sp.conjugate(tadpole_entry) / modulus
    return coupled, decoupled


def aligned_phase_sources():
    """CP-conserving test sources with both d-s entries multiplied by the phase w (aligned phases)."""
    sources = dict(illustrative_lo_inputs("cp_conserving"))
    for key in ("axial_source", "pseudoscalar_source"):
        matrix = sp.Matrix(sources[key])
        matrix[1, 2], matrix[2, 1] = matrix[1, 2] * ALIGNED_PHASE, matrix[2, 1] * sp.conjugate(ALIGNED_PHASE)
        sources[key] = matrix
    return sources


def source_cases():
    """Name -> base sources (chi_0 diagonal in all cases)."""
    return {
        "CP-conserving (real d-s entries)": illustrative_lo_inputs("cp_conserving"),
        "aligned phases": aligned_phase_sources(),
        "CP-violating (Im(c^A_ds conj c^P_ds) != 0)": illustrative_lo_inputs("neutral"),
        "no d-s sources": illustrative_lo_inputs("flavor_diagonal"),
    }


# ------------------------------------------------------------------ symbolic


def symbolic_checks():
    """Items 1 and 2 with symbolic diagonal chi_0, complex (c_q^AP)_ds, c_GA, m_0^2, s, epsilon."""
    s = SQUARED_MOMENTUM
    x_u, x_d, x_s, anomaly, epsilon = sp.symbols("x_u x_d x_s m0sq epsilon", positive=True)
    gluon_axial = sp.Symbol("c_GA", real=True)
    real_part, imaginary_part = sp.symbols("Re_cAP_ds Im_cAP_ds", real=True)
    diagonal = sp.symbols("cAP_uu cAP_dd cAP_ss", real=True)
    difference = sp.diag(*diagonal)
    difference[1, 2], difference[2, 1] = real_part + sp.I * imaginary_part, real_part - sp.I * imaginary_part
    quark_spurion = sp.diag(x_u, x_d, x_s)
    kaon = (x_d + x_s) / 2
    tadpole = tadpole_source(difference, gluon_axial, quark_spurion, anomaly)
    doublet_part = strange_down_part(tadpole)
    modulus_squared = real_part**2 + imaginary_part**2
    sum_of_squares = sp.trace(LAMBDA_6 * tadpole)**2 + sp.trace(LAMBDA_7 * tadpole)**2
    contribution = -2 * sp.trace(doublet_part * inverse_meson_operator(doublet_part, quark_spurion, anomaly))
    coupled, decoupled = doublet_generators(tadpole[1, 2])

    # bare null vector at s = m_K0,0^2 (generic sources with this c^AP): D_ja = epsilon tr(lambda_j Sigma(s))
    axial = _symbolic_axial()
    pseudoscalar = axial - difference
    gluon = gluon_axial - sp.trace(axial)
    sigma = mixing_source_matrix(kaon, axial, pseudoscalar, quark_spurion, gluon, anomaly)
    entry_6, entry_7 = epsilon * sp.trace(LAMBDA_6 * sigma), epsilon * sp.trace(LAMBDA_7 * sigma)
    null_generator = -entry_7 * LAMBDA_6 + entry_6 * LAMBDA_7
    ratio = sp.simplify(null_generator[1, 2] / decoupled[1, 2])
    return {
        "1. 𝒞_ds = -m_K0,0^2 (c_q^AP)_ds": _is_zero(tadpole[1, 2] + kaon * difference[1, 2]),
        "1. 𝕄_s(d-s part) = (m_K0,0^2 - s) d-s part": _is_zero(
            meson_mass_operator(doublet_part, quark_spurion, anomaly, s) - (kaon - s) * doublet_part),
        "1. tr(lambda_6 𝒞)^2 + tr(lambda_7 𝒞)^2 = 4|𝒞_ds|^2 = 4 m_K0,0^4 |(c_q^AP)_ds|^2": _is_zero(
            sum_of_squares - 4 * kaon**2 * modulus_squared),
        "1. -2 tr(d-s part 𝕄_s^-1(d-s part)) = 4 m_K0,0^4 |(c_q^AP)_ds|^2/(s - m_K0,0^2)": _is_zero(
            sp.together(contribution - 4 * kaon**2 * modulus_squared / (s - kaon))),
        "2. coupled and decoupled generators: tr = 2 delta": _is_zero(sp.Matrix([
            sp.trace(coupled * coupled) - 2, sp.trace(decoupled * decoupled) - 2, sp.trace(coupled * decoupled)])),
        "2. both are eigen-generators of 𝕄_0 with eigenvalue m_K0,0^2": _is_zero(
            meson_mass_operator(coupled, quark_spurion, anomaly) - kaon * coupled) and _is_zero(
            meson_mass_operator(decoupled, quark_spurion, anomaly) - kaon * decoupled),
        "2. tr(coupled 𝒞)^2 = 4 m_K0,0^4 |(c_q^AP)_ds|^2, tr(decoupled 𝒞) = 0": _is_zero(
            sp.trace(coupled * tadpole)**2 - 4 * kaon**2 * modulus_squared) and _is_zero(
            sp.trace(decoupled * tadpole)),
        "2. null vector (-D_7a, D_6a) of D(m_K0,0^2) is proportional to the decoupled generator": (
            _is_zero(null_generator - ratio * decoupled) and not _is_zero(ratio)),
    }


def _symbolic_axial():
    real_part, imaginary_part = sp.symbols("Re_cA_ds Im_cA_ds", real=True)
    matrix = sp.diag(*sp.symbols("cA_uu cA_dd cA_ss", real=True))
    matrix[1, 2], matrix[2, 1] = real_part + sp.I * imaginary_part, real_part - sp.I * imaginary_part
    return matrix


# ------------------------------------------------------------------ rational


def rational_case_checks(name, base_sources, epsilon):
    """Item 3 (and the case-specific part of item 4) for one source set and one epsilon."""
    s = SQUARED_MOMENTUM
    sources = rational_test_sources(epsilon, base_sources)
    kinetic, mass_squared = bare_blocks(sources, TEST_CONDENSATE, TEST_SECOND_ORDER_SCALAR_TRACE)
    kaon = (sources["quark_spurion"][1, 1] + sources["quark_spurion"][2, 2]) / 2
    difference = sources["axial_source"] - sources["pseudoscalar_source"]
    expected_multiplicity = 1 if sp.expand(difference[1, 2] * difference[2, 1]) != 0 else 2
    polynomial = sp.Poly(sp.expand((s * kinetic - mass_squared).det(method="berkowitz")), s)
    factor, multiplicity, remaining = sp.Poly(s - kaon, s), 0, polynomial
    while remaining.rem(factor).is_zero:
        multiplicity, remaining = multiplicity + 1, remaining.quo(factor)
    null_space = inverse_propagator(kinetic, mass_squared, kaon).nullspace(simplify=True)
    spectrum = solve_spectrum_numerical(np.array(kinetic.evalf(30), dtype=float),
                                        np.array(mass_squared.evalf(30), dtype=float),
                                        field_labels=NEUTRAL_FIELD_LABELS)
    poles = np.array(spectrum.pole_squared_masses)
    at_kaon = [index for index, pole in enumerate(poles) if abs(pole - float(kaon)) < 1e-10]
    alp_components = [abs(spectrum.pole_vector(index)[ALP]) for index in at_kaon]
    checks = {
        f"3. ({name}, epsilon = {epsilon}) (s - m_K0,0^2) divides det(s K - M^2) {expected_multiplicity} time(s)": (
            multiplicity == expected_multiplicity),
        f"3. ({name}, epsilon = {epsilon}) exact null space of D(m_K0,0^2): dimension {expected_multiplicity},"
        " zero ALP component": (len(null_space) == expected_multiplicity
                                and all(sp.simplify(vector[ALP]) == 0 for vector in null_space)),
        f"3. ({name}, epsilon = {epsilon}) solver: pole at m_K0,0^2, ALP component < 1e-10": (
            len(at_kaon) >= 1 and all(value < 1e-10 for value in alp_components)),
    }
    return checks, polynomial


def reduction_checks():
    """Item 4."""
    s = SQUARED_MOMENTUM
    checks = {}
    epsilon = EPSILON_VALUES[0]
    # CP-conserving
    sources = rational_test_sources(epsilon, illustrative_lo_inputs("cp_conserving"))
    kinetic, mass_squared = bare_blocks(sources, TEST_CONDENSATE, TEST_SECOND_ORDER_SCALAR_TRACE)
    tadpole = tadpole_source_from_sources(**_block_inputs(sources))
    checks["4. CP-conserving: K_7a = M^2_7a = 0 and tr(lambda_7 𝒞) = 0 (P_7 decouples)"] = (
        kinetic[GENERATOR_INDEX_7, ALP] == 0 and sp.simplify(mass_squared[GENERATOR_INDEX_7, ALP]) == 0
        and sp.simplify(sp.trace(LAMBDA_7 * tadpole)) == 0)
    # aligned phases: rephasing
    aligned = rational_test_sources(epsilon, aligned_phase_sources())
    rephasing = sp.diag(1, 1, ALIGNED_PHASE)
    rephased = dict(aligned)
    for key in ("axial_source", "pseudoscalar_source", "quark_spurion"):
        rephased[key] = (rephasing * aligned[key] * rephasing.H).applyfunc(sp.simplify)
    aligned_kinetic, aligned_mass = bare_blocks(aligned, TEST_CONDENSATE, TEST_SECOND_ORDER_SCALAR_TRACE)
    aligned_polynomial = sp.Poly(sp.expand((s * aligned_kinetic - aligned_mass).det(method="berkowitz")), s)
    rephased_kinetic, rephased_mass = bare_blocks(rephased, TEST_CONDENSATE, TEST_SECOND_ORDER_SCALAR_TRACE)
    rephased_polynomial = sp.Poly(sp.expand((s * rephased_kinetic - rephased_mass).det(method="berkowitz")), s)
    checks["4. aligned phases: rephasing makes the d-s entries real, chi_0 unchanged"] = (
        sp.im(rephased["axial_source"][1, 2]) == 0 and sp.im(rephased["pseudoscalar_source"][1, 2]) == 0
        and rephased["quark_spurion"] == aligned["quark_spurion"])
    checks["4. aligned phases: det(s K - M^2) unchanged by the rephasing"] = aligned_polynomial == rephased_polynomial
    determinant_aligned = (aligned_kinetic[GENERATOR_INDEX_6, ALP] * aligned_mass[GENERATOR_INDEX_7, ALP]
                           - aligned_kinetic[GENERATOR_INDEX_7, ALP] * aligned_mass[GENERATOR_INDEX_6, ALP])
    checks["4. aligned phases: K_6a M^2_7a - K_7a M^2_6a = 0"] = sp.simplify(determinant_aligned) == 0
    # CP-violating
    violating = rational_test_sources(epsilon, illustrative_lo_inputs("neutral"))
    kinetic, mass_squared = bare_blocks(violating, TEST_CONDENSATE, TEST_SECOND_ORDER_SCALAR_TRACE)
    determinant = (kinetic[GENERATOR_INDEX_6, ALP] * mass_squared[GENERATOR_INDEX_7, ALP]
                   - kinetic[GENERATOR_INDEX_7, ALP] * mass_squared[GENERATOR_INDEX_6, ALP])
    checks["4. CP-violating: K_6a M^2_7a - K_7a M^2_6a != 0 (no decoupled bare-frame field)"] = (
        sp.simplify(determinant) != 0)
    canonical = canonical_mass_matrix_from_sources(violating, TEST_CONDENSATE, TEST_SECOND_ORDER_SCALAR_TRACE)
    tadpole = tadpole_source_from_sources(**_block_inputs(violating))
    entry = tadpole[1, 2]
    phase = sp.conjugate(entry) / sp.sqrt(sp.expand(entry * sp.conjugate(entry)))  # makes the d-s entry real
    rotation = sp.diag(1, 1, sp.conjugate(phase))
    rotated_tadpole = (rotation * tadpole * rotation.H).applyfunc(sp.radsimp)
    rotated = canonical_mass_matrix(canonical[:MESON_COUNT, :MESON_COUNT], neutral_components(rotated_tadpole),
                                    epsilon, canonical[ALP, ALP])
    original_polynomial = sp.expand((s * sp.eye(MESON_COUNT + 1) - canonical).det(method="berkowitz"))
    rotated_polynomial = sp.expand((s * sp.eye(MESON_COUNT + 1) - rotated).det(method="berkowitz"))
    checks["4. CP-violating: rotated 𝒞 has a real d-s entry of modulus |𝒞_ds|"] = (
        sp.simplify(sp.im(rotated_tadpole[1, 2])) == 0
        and sp.simplify(rotated_tadpole[1, 2]**2 - entry * sp.conjugate(entry)) == 0)
    checks["4. CP-violating: det(s - M_can^2) depends on the d-s part only through |(c_q^AP)_ds|^2"] = (
        sp.simplify(original_polynomial - rotated_polynomial) == 0)
    return checks


def derive():
    """Run all checks."""
    checks = dict(symbolic_checks())
    for name, base in source_cases().items():
        for epsilon in EPSILON_VALUES:
            case_checks, _ = rational_case_checks(name, base, epsilon)
            checks.update(case_checks)
    checks.update(reduction_checks())
    return {"checks": checks}


def main():
    result = derive()
    print("K^0 doublet, diagonal chi_0: one pole exactly at m_K0,0^2 (decoupled combination, zero ALP component);")
    print("  the coupled combination enters the secular sum as 4 m_K0,0^4 |(c_q^AP)_ds|^2/(s - m_K0,0^2).")
    for label, passed in result["checks"].items():
        print(f"[{'ok' if passed else 'FAIL'}] {label}")
    passed = all(result["checks"].values())
    print("All checks pass" if passed else "Some checks FAIL")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
