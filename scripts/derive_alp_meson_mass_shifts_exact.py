#!/usr/bin/env python3
"""Exact LO pole equations (no expansion in epsilon = f_pi/f_a) of the neutral mesons and the ALP.

Run:   .venv/bin/python scripts/derive_alp_meson_mass_shifts_exact.py
Tests: .venv/bin/python -m unittest discover -s scripts -p test_derive_alp_meson_mass_shifts_exact.py
Record: memory/alp/shared/alp-lo-exact-pole-equations.md

Conventions (Alda:2026nrd; memory/alp/shared/finite-alp-meson-mixing.md,
memory/alp/shared/alp-lo-mass-sum-rules.md): N_f = 3, neutral basis
(P_3, P_8, P_0, P_6, P_7, a), tr(lambda_i lambda_j) = 2 delta_ij, D(s) = s K - M^2,
epsilon = f_pi/f_a, M^2_aa = m_a0^2 + epsilon^2 (2 m_0^2 c_G^2/N_f - B_0 tr S_2).
Flavour-space objects of the second-order script: 𝕄_s(X) = {chi_0, X}/2 + (m_0^2/N_f) tr(X) 1 - s X,
𝒞 = -{chi_0, c_q^AP}/2 - (m_0^2/N_f) c_GA 1, Sigma(s) with D_ia(s) = epsilon tr(lambda_i Sigma(s)),
lambda_n: neutral, 𝕄_0(lambda_n) = m_n0^2 lambda_n, tr(lambda_n lambda_m) = 2 delta_nm.
𝕄_s^{-1}(X) is the neutral Y with 𝕄_s(Y) = X (solved as a matrix equation); (s - 𝕄_0)^{-1} = -𝕄_s^{-1}.
B_0 tr hat S_2 = B_0 tr S_2 + 2 tr(c_q^A {c_q^P, chi_0}) - 2 tr(c_q^A c_q^A chi_0)  (sum-rule record, (4.1)).
In the code, invariant_vector is the column (tr(lambda_i 𝒞)) over the neutral Gell-Mann generators.

Derived and checked here:
 0. K_ia = -epsilon tr(lambda_i c_q^A), K_aa = 1 + 2 epsilon^2 tr(c_q^A c_q^A), and
    sum_i tr(lambda_i X)^2 = 2 tr(X X) for neutral X (symbolic).
 1. Schur factorization det D(s) = det(s - M^2_PP) D^eff_aa(s),
    D^eff_aa(s) = s K_aa - M^2_aa - D_aP (s - M^2_PP)^{-1} D_Pa (generic symbolic block and exact LO
    blocks with symbolic s); D_aP (s - M^2_PP)^{-1} D_Pa = -2 epsilon^2 tr(Sigma(s) 𝕄_s^{-1}(Sigma(s)))
    (exact rationals, symbolic s).
 2. Canonical frame: with P_i = P_i' + epsilon tr(lambda_i c_q^A) a (det T = 1),
    T^T D(s) T = s 1 - M_can^2 with
        M_can,PP = M^2_PP,  (M_can)_ia = -epsilon tr(lambda_i 𝒞),
        M_can,aa = m_a0^2 + epsilon^2 [ (2 m_0^2/N_f) c_GA^2 - B_0 tr hat S_2 ]
    (symbolic Q-commuting sources, chi_0 with a d-s entry). Every entry is an invariant.
 3. Hence D^eff_aa(s) = s - M_can,aa + 2 epsilon^2 tr(𝒞 𝕄_s^{-1}(𝒞)), equal to the bare Schur form (exact);
    (s - 𝕄_0)^{-1}(𝒞) = -𝕄_s^{-1}(𝒞) (exact, sign of the rewrite);
    sum_ij tr(lambda_i 𝒞) [(s - M^2_PP)^{-1}]_ij tr(lambda_j 𝒞) = -2 tr(𝒞 𝕄_s^{-1}(𝒞)) (exact, symbolic s)
                           = sum_n tr(lambda_n 𝒞)^2/(s - m_n0^2) (50 digits), and
    𝕄_s^{-1}(𝒞) = sum_n lambda_n tr(lambda_n 𝒞)/(2 (m_n0^2 - s)) (50 digits).
 4. Invariance: under c_q^{A,P} -> + kappa, c_G -> c_G - tr kappa, S_2 -> S_2' (from M' = g_R M g_L^dagger)
    the canonical blocks are unchanged (symbolic), and det(s K - M^2) is unchanged (exact rationals).
 5. Meson pole n (non-degenerate eigenvalue m_n0^2):
        (s - m_n0^2) [ s - M_can,aa + 2 epsilon^2 tr(𝒞 𝕄_s^{-1}(𝒞 - tr(lambda_n 𝒞) lambda_n/2)) ]
            = epsilon^2 tr(lambda_n 𝒞)^2,
    the subtracted term equal to sum_{m != n} tr(lambda_m 𝒞)^2/(s - m_m0^2) (50 digits), and equal
    to the bare form (s - m_n0^2) [s K_aa - M^2_aa - sum_{m != n} D_ma^2/(s - m_m0^2)] = D_na(s)^2
    (generic symbolic mass-basis model).
 6. Poles as eigenvalues: det(s - M_can^2) = det(s K - M^2) = det(s - M^2_PP) D^eff_aa(s) (exact rationals);
    arrow form and interlacing in scripts/derive_alp_lo_arrow_matrix_poles.py.
 7. Value at an unmixed mass: det D(m_n0^2) = -epsilon^2 tr(lambda_n 𝒞)^2 prod_{m != n}(m_n0^2 - m_m0^2).
 8. Diagonal chi_0 = diag(x_u, x_d, x_s), mu = m_0^2/N_f, 𝒞_qq = -x_q (c_q^AP)_qq - mu c_GA:
        -2 tr(𝒞 𝕄_s^{-1}(𝒞)) = 2 N(s)/P_eta(s) + 4 m_K0,0^4 |(c_q^AP)_ds|^2 / (s - m_K0,0^2),
        P_eta(s) = prod_q (s - x_q) - mu sum_q prod_{r != q} (s - x_r),
        N(s) = sum_q 𝒞_qq^2 prod_{r != q}(s - x_r) - mu sum_{q<r} (x_q c_qq - x_r c_rr)^2 (s - x_third),
    and the secular polynomial (degree 6, leading coefficient det K = 1)
        det D(s) = (s - m_K0,0^2) { (s - m_K0,0^2) [ P_eta(s) (s - M_can,aa) - 2 epsilon^2 N(s) ]
                                    - 4 epsilon^2 m_K0,0^4 |(c_q^AP)_ds|^2 P_eta(s) }
    (symbolic against det(s - M_can^2); exact rationals against det(s K - M^2) of closed_form_blocks).
 9. Roots: exact real-root isolation of the explicit polynomial (6 real roots; m_K0,0^2 is a root of
    multiplicity 1 if (c_q^AP)_ds != 0 and 2 if it vanishes) compared with solve_spectrum_numerical;
    D^eff_aa vanishes at all other roots.
10. Residues for diagonal chi_0: tr(lambda_n 𝒞)^2 = 2 N(m_n0^2)/P_eta'(m_n0^2) for n = pi^0, eta, eta' and
    the sum of tr(lambda_n 𝒞)^2 over the K^0 doublet = 4 m_K0,0^4 |(c_q^AP)_ds|^2 (50 digits).

Units: f_pi, f_a, B_0, tr S_2 in GeV; chi_0, m_0^2, m_a0^2, M^2, s in GeV^2; couplings dimensionless.
Rational inputs are test values, not fits.
"""

from __future__ import annotations

import mpmath
import numpy as np
import sympy as sp

from derive_alp_lo_mass_sum_rules import (
    _hermitian_charge_conserving,
    add_second_order_scalar_term,
    constant_rotation,
    dressed_second_order_scalar_trace,
    symbolic_second_order_scalar,
    symbolic_sources,
)
from derive_alp_lo_neutral_blocks import FLAVOR_COUNT, NEUTRAL_FIELD_LABELS, closed_form_blocks, illustrative_lo_inputs
from derive_alp_meson_mass_shifts_second_order import (
    HIGH_PRECISION_DIGITS,
    HIGH_PRECISION_TOLERANCE,
    MESON_COUNT,
    NEUTRAL_GENERATORS,
    _block_inputs,
    high_precision_trace,
    invariant_mixing_vector,
    mass_eigenstate_generators,
    meson_mass_operator,
    mixing_source_matrix,
    neutral_components,
    rational_test_sources,
    solve_neutral_matrix_equation,
    tadpole_source_from_sources,
    to_high_precision,
)
from solve_alp_meson_spectrum import solve_spectrum_numerical
from symbolic_utils import anticommutator

TEST_CONDENSATE = sp.Rational(27, 10)  # B_0 in GeV (test value)
TEST_SECOND_ORDER_SCALAR_TRACE = sp.Rational(1, 20)  # tr S_2 in GeV (test value)
TEST_EPSILON = sp.Rational(1, 2)  # f_pi/f_a, finite mixing
SQUARED_MOMENTUM = sp.Symbol("s")
ALP = MESON_COUNT


def _simplify(expression):
    return sp.simplify(sp.expand(sp.radsimp(sp.expand(expression))))


def _is_zero(matrix_or_expression):
    entries = matrix_or_expression if isinstance(matrix_or_expression, sp.MatrixBase) else [matrix_or_expression]
    return all(_simplify(entry) == 0 for entry in entries)


# ------------------------------------------------------------------ building blocks


def bare_blocks(sources, condensate, second_order_scalar_trace):
    """(K, M^2) of closed_form_blocks with -epsilon^2 B_0 tr S_2 added to M^2_aa."""
    kinetic, mass_squared = closed_form_blocks(**sources)
    return kinetic, add_second_order_scalar_term(mass_squared, sources["pion_decay_constant"],
                                                 sources["alp_decay_constant"], condensate,
                                                 second_order_scalar_trace)


def dressed_scalar_trace_times_condensate(sources, condensate, second_order_scalar_trace):
    """B_0 tr hat S_2 = B_0 tr S_2 + 2 tr(c^A {c^P, chi_0}) - 2 tr(c^A c^A chi_0), with m_q = chi_0/(2 B_0)."""
    quark_mass = sources["quark_spurion"] / (2 * condensate)
    return sp.expand(condensate * dressed_second_order_scalar_trace(
        second_order_scalar_trace, sources["axial_source"],
        anticommutator(sources["pseudoscalar_source"], quark_mass), quark_mass))


def canonical_alp_mass_squared(bare_alp_mass_squared, epsilon, anomaly_mass_squared, gluon_axial_combination,
                               condensate_times_dressed_trace):
    """M_can,aa = m_a0^2 + epsilon^2 [(2 m_0^2/N_f) c_GA^2 - B_0 tr hat S_2]."""
    return bare_alp_mass_squared + epsilon**2 * (2 * anomaly_mass_squared * gluon_axial_combination**2 / FLAVOR_COUNT
                                                 - condensate_times_dressed_trace)


def canonical_mass_matrix(meson_mass_squared, invariant_vector, epsilon, alp_mass_squared_canonical):
    """M_can^2 = [[M^2_PP, column -epsilon tr(lambda_i 𝒞)], [row -epsilon tr(lambda_i 𝒞), M_can,aa]]."""
    size = meson_mass_squared.shape[0] + 1
    matrix = sp.zeros(size)
    matrix[:size - 1, :size - 1] = meson_mass_squared
    matrix[:size - 1, size - 1] = -epsilon * invariant_vector
    matrix[size - 1, :size - 1] = -epsilon * invariant_vector.T
    matrix[size - 1, size - 1] = alp_mass_squared_canonical
    return matrix


def canonical_mass_matrix_from_sources(sources, condensate, second_order_scalar_trace):
    """M_can^2 written with c_q^AP, c_GA, tr hat S_2 only (via invariant_mixing_vector and (4.1))."""
    _, mass_squared = closed_form_blocks(**sources)
    epsilon = sources["pion_decay_constant"] / sources["alp_decay_constant"]
    return canonical_mass_matrix(
        mass_squared[:MESON_COUNT, :MESON_COUNT], invariant_mixing_vector(**_block_inputs(sources)), epsilon,
        canonical_alp_mass_squared(sources["bare_alp_mass_squared"], epsilon, sources["anomaly_mass_squared"],
                                   sources["gluon_coupling"] + sp.trace(sources["axial_source"]),
                                   dressed_scalar_trace_times_condensate(sources, condensate,
                                                                         second_order_scalar_trace)))


def kinetic_shift(sources):
    """T = [[1, epsilon t], [0, 1]]: P = P' + epsilon t a, a = a'."""
    epsilon = sources["pion_decay_constant"] / sources["alp_decay_constant"]
    transformation = sp.eye(MESON_COUNT + 1)
    transformation[:MESON_COUNT, ALP] = epsilon * neutral_components(sources["axial_source"])
    return transformation


def schur_complement(inverse_propagator_matrix):
    """D_aa - D_aP D_PP^{-1} D_Pa for the last field."""
    meson_block = inverse_propagator_matrix[:-1, :-1]
    mixed = inverse_propagator_matrix[:-1, -1]
    return inverse_propagator_matrix[-1, -1] - (mixed.T * meson_block.LUsolve(mixed))[0]


def flavour_diagonal_meson_polynomial(quark_spurion_diagonal, anomaly_mass_squared, squared_momentum=SQUARED_MOMENTUM):
    """P_eta(s) = prod_q (s - x_q) - (m_0^2/N_f) sum_q prod_{r != q} (s - x_r)."""
    factors = [squared_momentum - value for value in quark_spurion_diagonal]
    return sp.expand(sp.prod(factors) - anomaly_mass_squared / FLAVOR_COUNT
                     * sum(sp.prod(factors[:q] + factors[q + 1:]) for q in range(3)))


def mixing_numerator_polynomial(quark_spurion_diagonal, difference_diagonal, anomaly_mass_squared,
                                gluon_axial_combination, squared_momentum=SQUARED_MOMENTUM):
    """N(s) = sum_q 𝒞_qq^2 prod_{r != q}(s - x_r) - mu sum_{q<r} (x_q c_qq - x_r c_rr)^2 (s - x_third)."""
    anomaly_term = anomaly_mass_squared / FLAVOR_COUNT
    factors = [squared_momentum - value for value in quark_spurion_diagonal]
    tadpole_diagonal = [-value * difference - anomaly_term * gluon_axial_combination
                        for value, difference in zip(quark_spurion_diagonal, difference_diagonal)]
    first = sum(tadpole_diagonal[q]**2 * sp.prod(factors[:q] + factors[q + 1:]) for q in range(3))
    second = sum((quark_spurion_diagonal[q] * difference_diagonal[q]
                  - quark_spurion_diagonal[r] * difference_diagonal[r])**2 * factors[3 - q - r]
                 for q in range(3) for r in range(q + 1, 3))
    return sp.expand(first - anomaly_term * second)


def explicit_secular_polynomial(quark_spurion_diagonal, difference_source, anomaly_mass_squared,
                                gluon_axial_combination, epsilon, alp_mass_squared_canonical,
                                squared_momentum=SQUARED_MOMENTUM):
    """det D(s) for diagonal chi_0 (item 8)."""
    kaon_mass_squared = (quark_spurion_diagonal[1] + quark_spurion_diagonal[2]) / 2
    difference_diagonal = [difference_source[index, index] for index in range(3)]
    kaon_coupling_squared = sp.expand(difference_source[1, 2] * difference_source[2, 1])
    meson_polynomial = flavour_diagonal_meson_polynomial(quark_spurion_diagonal, anomaly_mass_squared,
                                                         squared_momentum)
    numerator = mixing_numerator_polynomial(quark_spurion_diagonal, difference_diagonal, anomaly_mass_squared,
                                            gluon_axial_combination, squared_momentum)
    kaon_factor = squared_momentum - kaon_mass_squared
    return sp.expand(kaon_factor * (
        kaon_factor * (meson_polynomial * (squared_momentum - alp_mass_squared_canonical) - 2 * epsilon**2 * numerator)
        - 4 * epsilon**2 * kaon_mass_squared**2 * kaon_coupling_squared * meson_polynomial))


def explicit_meson_sum(quark_spurion_diagonal, difference_source, anomaly_mass_squared, gluon_axial_combination,
                       squared_momentum=SQUARED_MOMENTUM):
    """-2 tr(𝒞 𝕄_s^{-1}(𝒞)) = 2 N(s)/P_eta(s) + 4 m_K0,0^4 |c^AP_ds|^2/(s - m_K0,0^2) for diagonal chi_0."""
    kaon_mass_squared = (quark_spurion_diagonal[1] + quark_spurion_diagonal[2]) / 2
    difference_diagonal = [difference_source[index, index] for index in range(3)]
    return (2 * mixing_numerator_polynomial(quark_spurion_diagonal, difference_diagonal, anomaly_mass_squared,
                                            gluon_axial_combination, squared_momentum)
            / flavour_diagonal_meson_polynomial(quark_spurion_diagonal, anomaly_mass_squared, squared_momentum)
            + 4 * kaon_mass_squared**2 * difference_source[1, 2] * difference_source[2, 1]
            / (squared_momentum - kaon_mass_squared))


# ------------------------------------------------------------------ checks


def generic_schur_checks():
    """Item 1 (generic symmetric 4-meson block) and items 3, 5, 7 (generic mass-basis model, 3 mesons)."""
    s = SQUARED_MOMENTUM
    # item 1, generic symmetric block
    size = 4
    entries = sp.symbols(f"M0:{size * (size + 1) // 2}", real=True)
    meson_block = sp.zeros(size)
    position = 0
    for row in range(size):
        for column in range(row, size):
            meson_block[row, column] = meson_block[column, row] = entries[position]
            position += 1
    mixed = sp.Matrix(sp.symbols(f"D0:{size}", real=True))
    alp_entry = sp.Symbol("Daa", real=True)
    full = sp.zeros(size + 1)
    full[:size, :size] = s * sp.eye(size) - meson_block
    full[:size, size] = mixed
    full[size, :size] = mixed.T
    full[size, size] = alp_entry
    factorization = _is_zero(full.det(method="berkowitz")
                             - (s * sp.eye(size) - meson_block).det(method="berkowitz")
                             * (alp_entry - (mixed.T * (s * sp.eye(size) - meson_block).adjugate() * mixed)[0]
                                / (s * sp.eye(size) - meson_block).det(method="berkowitz")))

    # items 3, 5, 7, generic mass basis
    epsilon, alp_canonical_shift, bare_alp_mass_squared = sp.symbols("epsilon delta_aa m_a0sq", real=True)
    count = 3
    unmixed = sp.symbols(f"m0sq_0:{count}", real=True)
    axial = sp.symbols(f"t0:{count}", real=True)
    invariant = sp.symbols(f"c0:{count}", real=True)
    alp_canonical = bare_alp_mass_squared + epsilon**2 * alp_canonical_shift
    mixed_entries = [epsilon * ((unmixed[k] - s) * axial[k] + invariant[k]) for k in range(count)]
    kinetic_alp = 1 + epsilon**2 * sum(value**2 for value in axial)
    mass_alp = (alp_canonical + epsilon**2 * sum(unmixed[k] * axial[k]**2 for k in range(count))
                + 2 * epsilon**2 * sum(axial[k] * invariant[k] for k in range(count)))
    effective_bare = s * kinetic_alp - mass_alp - sum(mixed_entries[k]**2 / (s - unmixed[k]) for k in range(count))
    effective_canonical = s - alp_canonical - epsilon**2 * sum(invariant[k]**2 / (s - unmixed[k]) for k in range(count))
    selected = 0
    bare_equation = ((s - unmixed[selected]) * (s * kinetic_alp - mass_alp - sum(
        mixed_entries[k]**2 / (s - unmixed[k]) for k in range(count) if k != selected)) - mixed_entries[selected]**2)
    canonical_equation = ((s - unmixed[selected]) * (s - alp_canonical - epsilon**2 * sum(
        invariant[k]**2 / (s - unmixed[k]) for k in range(count) if k != selected))
        - epsilon**2 * invariant[selected]**2)

    determinant = sp.prod([s - value for value in unmixed]) * effective_canonical

    return {
        "1. det D = det D_PP (D_aa - D_aP D_PP^-1 D_Pa) (generic 4-meson block)": factorization,
        "3. bare and canonical Schur complements agree (generic mass basis)": _is_zero(
            sp.together(effective_bare - effective_canonical)),
        "5. bare and canonical meson-pole equations agree (mass basis)": _is_zero(sp.together(bare_equation - canonical_equation)),
        "5'. both equal (s - m_n0^2) D^eff_aa(s)": _is_zero(
            sp.together(canonical_equation - (s - unmixed[selected]) * effective_canonical)),
        "7. det D(m_n0^2) = -epsilon^2 tr(lambda_n 𝒞)^2 prod_{m != n}(m_n0^2 - m_m0^2)": _is_zero(
            sp.cancel(determinant).subs(s, unmixed[selected])
            + epsilon**2 * invariant[selected]**2 * sp.prod([unmixed[selected] - unmixed[k]
                                                            for k in range(count) if k != selected])),
    }


def canonical_frame_checks(diagonal_quark_spurion):
    """Items 2 and 4 with symbolic Q-commuting sources."""
    sources = symbolic_sources(diagonal_quark_spurion=diagonal_quark_spurion)
    condensate, second_order_scalar_trace = symbolic_second_order_scalar()
    s = SQUARED_MOMENTUM
    kinetic, mass_squared = bare_blocks(sources, condensate, second_order_scalar_trace)
    transformation = kinetic_shift(sources)
    transformed = (transformation.T * (s * kinetic - mass_squared) * transformation)
    canonical = canonical_mass_matrix_from_sources(sources, condensate, second_order_scalar_trace)

    # constant axial redefinition; S_2 transformed with constant_rotation (M' = g_R M g_L^dagger)
    axial_shift = _hermitian_charge_conserving("kappa")
    quark_mass = sources["quark_spurion"] / (2 * condensate)
    scalar_second = sp.diag(second_order_scalar_trace, 0, 0)  # only tr S_2 enters the blocks
    series = {
        "quark_mass": quark_mass,
        "first_order_scalar": sp.zeros(3),
        "first_order_pseudoscalar": anticommutator(sources["pseudoscalar_source"], quark_mass),
        "second_order_scalar": scalar_second,
        "second_order_pseudoscalar": sp.zeros(3),
        "vector_couplings": (sp.zeros(3), sp.zeros(3)),
        "axial_couplings": (sources["axial_source"], sp.zeros(3)),
    }
    rotated = constant_rotation(series, sp.zeros(3), axial_shift)
    shifted_sources = dict(sources, axial_source=sources["axial_source"] + axial_shift,
                           pseudoscalar_source=sources["pseudoscalar_source"] + axial_shift,
                           gluon_coupling=sources["gluon_coupling"] - sp.trace(axial_shift))
    shifted_trace = sp.expand(rotated["second_order_scalar"].trace())
    shifted_canonical = canonical_mass_matrix_from_sources(shifted_sources, condensate, shifted_trace)
    shifted_kinetic, shifted_mass = bare_blocks(shifted_sources, condensate, shifted_trace)
    shifted_transformation = kinetic_shift(shifted_sources)
    label = "diagonal chi_0" if diagonal_quark_spurion else "chi_0 with d-s entry"
    return {
        f"2. ({label}) T^T D(s) T = s - M_can^2 in invariants": _is_zero(
            transformed - (s * sp.eye(MESON_COUNT + 1) - canonical)),
        f"2'. ({label}) det T = 1": transformation.det() == 1,
        f"4. ({label}) rotated P_1 = {{c^P + kappa, m_q}}": _is_zero(
            rotated["first_order_pseudoscalar"]
            - anticommutator(shifted_sources["pseudoscalar_source"], quark_mass)),
        f"4'. ({label}) tr S_2' = tr S_2 - 4 tr(kappa P_1) - 4 tr(kappa^2 m_q)": _is_zero(
            shifted_trace - second_order_scalar_trace
            + 4 * sp.trace(axial_shift * series["first_order_pseudoscalar"])
            + 4 * sp.trace(axial_shift * axial_shift * quark_mass)),
        f"4''. ({label}) M_can^2 unchanged under kappa": _is_zero(shifted_canonical - canonical),
        f"4'''. ({label}) transformed bare blocks give the same s - M_can^2": _is_zero(
            shifted_transformation.T * (s * shifted_kinetic - shifted_mass) * shifted_transformation
            - (s * sp.eye(MESON_COUNT + 1) - canonical)),
    }


def explicit_form_symbolic_checks():
    """Item 8, symbolic diagonal chi_0 and symbolic c_q^AP, c_GA."""
    s = SQUARED_MOMENTUM
    quark_spurion_diagonal = sp.symbols("x_u x_d x_s", positive=True)
    anomaly_mass_squared, epsilon, alp_canonical = sp.symbols("m0sq epsilon Mcan_aa", positive=True)
    gluon_axial_combination = sp.Symbol("c_GA", real=True)
    difference_source = _hermitian_charge_conserving("cAP")
    quark_spurion = sp.diag(*quark_spurion_diagonal)
    sources = {"axial_source": difference_source, "pseudoscalar_source": sp.zeros(3),
               "quark_spurion": quark_spurion, "gluon_coupling": gluon_axial_combination - sp.trace(difference_source),
               "pion_decay_constant": epsilon, "alp_decay_constant": 1,
               "anomaly_mass_squared": anomaly_mass_squared, "bare_alp_mass_squared": 0}
    # c^A = c^AP, c^P = 0 is one representative; c depends on c^AP, c_GA only (second-order script, item 3)
    _, mass_squared = closed_form_blocks(**sources)
    meson_mass_squared = mass_squared[:MESON_COUNT, :MESON_COUNT]
    invariant_vector = invariant_mixing_vector(**_block_inputs(sources))
    canonical = canonical_mass_matrix(meson_mass_squared, invariant_vector, epsilon, alp_canonical)
    meson_inverse_propagator = s * sp.eye(MESON_COUNT) - meson_mass_squared
    meson_sum = (invariant_vector.T * meson_inverse_propagator.adjugate() * invariant_vector)[0] \
        / meson_inverse_propagator.det(method="berkowitz")
    explicit_sum = explicit_meson_sum(quark_spurion_diagonal, difference_source, anomaly_mass_squared,
                                      gluon_axial_combination)
    explicit_polynomial = explicit_secular_polynomial(quark_spurion_diagonal, difference_source, anomaly_mass_squared,
                                                      gluon_axial_combination, epsilon, alp_canonical)
    kaon_mass_squared = (quark_spurion_diagonal[1] + quark_spurion_diagonal[2]) / 2
    meson_polynomial = flavour_diagonal_meson_polynomial(quark_spurion_diagonal, anomaly_mass_squared)
    return {
        "8. -2 tr(𝒞 𝕄_s^-1(𝒞)) = 2 N/P_eta + 4 m_K0^4 |c^AP_ds|^2/(s - m_K0^2)": _is_zero(
            sp.together(meson_sum - explicit_sum)),
        "8'. det(s - M^2_PP) = (s - m_K0,0^2)^2 P_eta(s)": _is_zero(
            meson_inverse_propagator.det(method="berkowitz") - (s - kaon_mass_squared)**2 * meson_polynomial),
        "8''. det(s - M_can^2) = explicit secular polynomial": _is_zero(
            (s * sp.eye(MESON_COUNT + 1) - canonical).det(method="berkowitz") - explicit_polynomial),
        "8'''. degree 6, leading coefficient 1": (sp.Poly(explicit_polynomial, s).degree() == 6
                                                  and sp.Poly(explicit_polynomial, s).LC() == 1),
    }


def residue_checks(base_sources=None, digits=HIGH_PRECISION_DIGITS):
    """Item 10: for diagonal chi_0, tr(lambda_n 𝒞)^2 = 2 N(m_n0^2) / P_eta'(m_n0^2) for n = pi^0, eta, eta',
    and tr(lambda_n 𝒞)^2 summed over the K^0 doublet = 4 m_K0,0^4 |c^AP_ds|^2 (50 digits)."""
    s = SQUARED_MOMENTUM
    sources = rational_test_sources(TEST_EPSILON, base_sources)
    eigenbasis = mass_eigenstate_generators(sources, digits)
    tadpole = to_high_precision(tadpole_source_from_sources(**_block_inputs(sources)), digits)
    difference_source = sources["axial_source"] - sources["pseudoscalar_source"]
    quark_spurion_diagonal = [sources["quark_spurion"][index, index] for index in range(3)]
    kaon_mass_squared = (quark_spurion_diagonal[1] + quark_spurion_diagonal[2]) / 2
    meson_polynomial = flavour_diagonal_meson_polynomial(quark_spurion_diagonal, sources["anomaly_mass_squared"])
    numerator = mixing_numerator_polynomial(quark_spurion_diagonal,
                                            [difference_source[index, index] for index in range(3)],
                                            sources["anomaly_mass_squared"],
                                            sources["gluon_coupling"] + sp.trace(sources["axial_source"]))
    derivative = sp.diff(meson_polynomial, s)
    largest, kaon_sum = mpmath.mpf(0), mpmath.mpf(0)
    for value, generator in eigenbasis:
        trace_squared = high_precision_trace(generator * tadpole)**2
        if abs(value - mpmath.mpf(str(sp.N(kaon_mass_squared, digits + 10)))) < mpmath.mpf("1e-30"):
            kaon_sum += trace_squared
            continue
        point = sp.Float(str(value), digits + 10)
        closed_form = mpmath.mpf(str(sp.N(2 * numerator.subs(s, point) / derivative.subs(s, point), digits + 10)))
        largest = max(largest, abs(trace_squared - closed_form))
    kaon_closed_form = mpmath.mpf(str(sp.N(4 * kaon_mass_squared**2 * difference_source[1, 2]
                                           * difference_source[2, 1], digits + 10)))
    return {
        "10. tr(lambda_n 𝒞)^2 = 2 N(m_n0^2)/P_eta'(m_n0^2) for pi^0, eta, eta' (50 digits)":
            largest < HIGH_PRECISION_TOLERANCE,
        "10'. K^0 doublet: sum of tr(lambda_n 𝒞)^2 = 4 m_K0,0^4 |c^AP_ds|^2 (50 digits)":
            abs(kaon_sum - kaon_closed_form) < HIGH_PRECISION_TOLERANCE,
    }


def kinetic_trace_checks():
    """Item 0: kinetic blocks in trace form and the completeness relation on the neutral space (symbolic)."""
    sources = symbolic_sources(diagonal_quark_spurion=False)
    kinetic, _ = closed_form_blocks(**sources)
    epsilon = sources["pion_decay_constant"] / sources["alp_decay_constant"]
    axial = sources["axial_source"]
    neutral = _hermitian_charge_conserving("X")
    completeness = sum(sp.trace(generator * neutral)**2 for generator in NEUTRAL_GENERATORS) - 2 * sp.trace(
        neutral * neutral)
    return {
        "0. K_ia = -epsilon tr(lambda_i c_q^A)": _is_zero(kinetic[:MESON_COUNT, ALP]
                                                          + epsilon * neutral_components(axial)),
        "0. K_aa = 1 + 2 epsilon^2 tr(c_q^A c_q^A)": _is_zero(kinetic[ALP, ALP] - 1
                                                              - 2 * epsilon**2 * sp.trace(axial * axial)),
        "0. sum_i tr(lambda_i X)^2 = 2 tr(X X) for neutral X": _is_zero(completeness),
    }


def inverse_meson_operator(argument, quark_spurion, anomaly_mass_squared, squared_momentum=SQUARED_MOMENTUM):
    """𝕄_s^{-1}(X): the neutral Y with 𝕄_s(Y) = X, solved as a matrix equation (s not an eigenvalue of 𝕄_0)."""
    return solve_neutral_matrix_equation(
        lambda unknown: meson_mass_operator(unknown, quark_spurion, anomaly_mass_squared, squared_momentum), argument)


def trace_form_checks(base_sources=None, test_points=(sp.Rational(1, 3), sp.Rational(3, 5)),
                      digits=HIGH_PRECISION_DIGITS):
    """Items 1, 3, 5 in trace form: 𝕄_s^{-1} on the neutral space against the generator-basis matrix
    (exact, symbolic s) and against the sum over lambda_n (50 digits at the test points)."""
    s = SQUARED_MOMENTUM
    sources = rational_test_sources(TEST_EPSILON, base_sources)
    epsilon = TEST_EPSILON
    quark_spurion, anomaly_mass_squared = sources["quark_spurion"], sources["anomaly_mass_squared"]
    _, mass_squared = closed_form_blocks(**sources)
    meson_mass_squared = mass_squared[:MESON_COUNT, :MESON_COUNT]
    meson_inverse = (s * sp.eye(MESON_COUNT) - meson_mass_squared)
    tadpole = tadpole_source_from_sources(**_block_inputs(sources))
    invariant_vector = invariant_mixing_vector(**_block_inputs(sources))
    tadpole_inverse = inverse_meson_operator(tadpole, quark_spurion, anomaly_mass_squared)
    # (s - 𝕄_0)^{-1}(𝒞) solved separately, to fix the sign of the rewrite with 𝕄_s^{-1}
    shifted_inverse = solve_neutral_matrix_equation(
        lambda unknown: s * unknown - meson_mass_operator(unknown, quark_spurion, anomaly_mass_squared), tadpole)
    mixing_source = mixing_source_matrix(s, **_block_inputs(sources))
    mixing_vector = neutral_components(mixing_source)
    trace_form = -2 * sp.trace(tadpole * tadpole_inverse)
    generator_form = (invariant_vector.T * meson_inverse.LUsolve(invariant_vector))[0]
    bare_trace_form = -2 * epsilon**2 * sp.trace(mixing_source * inverse_meson_operator(
        mixing_source, quark_spurion, anomaly_mass_squared))
    bare_generator_form = epsilon**2 * (mixing_vector.T * meson_inverse.LUsolve(mixing_vector))[0]

    eigenbasis = mass_eigenstate_generators(sources, digits)
    tadpole_high = to_high_precision(tadpole, digits)
    spectral_difference, inverse_difference, projected_difference = (mpmath.mpf(0),) * 3
    for point in test_points:
        point_high = mpmath.mpf(str(sp.N(point, digits + 10)))
        inverse_high = to_high_precision(tadpole_inverse.subs(s, point), digits)
        traces = [high_precision_trace(generator * tadpole_high) for _, generator in eigenbasis]
        terms = [trace**2 / (point_high - value) for trace, (value, _) in zip(traces, eigenbasis)]
        spectral = mpmath.fsum(terms)
        spectral_difference = max(spectral_difference,
                                  abs(-2 * high_precision_trace(tadpole_high * inverse_high) - spectral))
        spectral_matrix = mpmath.matrix(3, 3)
        for trace, (value, generator) in zip(traces, eigenbasis):
            spectral_matrix += generator * (trace / (2 * (value - point_high)))
        inverse_difference = max(inverse_difference, max(
            abs(inverse_high[row, column] - spectral_matrix[row, column]) for row in range(3) for column in range(3)))
        # -2 tr(𝒞 𝕄_s^{-1}(𝒞 - tr(lambda_n 𝒞) lambda_n/2)) = sum_{m != n} tr(lambda_m 𝒞)^2/(s - m_m0^2)
        meson_operator_high = mpmath.matrix((meson_mass_squared - point * sp.eye(MESON_COUNT)).evalf(digits + 10))
        generators_high = [to_high_precision(sp.Matrix(generator), digits) for generator in NEUTRAL_GENERATORS]
        for index, (trace, (_, generator)) in enumerate(zip(traces, eigenbasis)):
            subtracted = tadpole_high - generator * (trace / 2)
            subtracted_components = mpmath.matrix([high_precision_trace(basis * subtracted).real
                                                   for basis in generators_high])
            solved = mpmath.lu_solve(meson_operator_high, subtracted_components)
            subtracted_inverse = mpmath.matrix(3, 3)
            for component, basis in zip(solved, generators_high):
                subtracted_inverse += basis * (component / 2)
            projected = -2 * high_precision_trace(tadpole_high * subtracted_inverse)
            projected_difference = max(projected_difference, abs(
                projected - mpmath.fsum(term for other, term in enumerate(terms) if other != index)))
    return {
        "5'. (s - 𝕄_0)^-1(𝒞) = -𝕄_s^-1(𝒞), so 2 tr(𝒞 (s - 𝕄_0)^-1(𝒞)) = -2 tr(𝒞 𝕄_s^-1(𝒞)) (exact)": _is_zero(
            (shifted_inverse + tadpole_inverse).applyfunc(sp.cancel)),
        "1. D_aP (s - M^2_PP)^-1 D_Pa = -2 epsilon^2 tr(Sigma(s) 𝕄_s^-1(Sigma(s))) (exact)": _is_zero(
            sp.cancel(bare_trace_form - bare_generator_form)),
        "3. sum_ij tr(lambda_i 𝒞) [(s - M^2_PP)^-1]_ij tr(lambda_j 𝒞) = -2 tr(𝒞 𝕄_s^-1(𝒞)) (exact)": _is_zero(sp.cancel(trace_form - generator_form)),
        "3. -2 tr(𝒞 𝕄_s^-1(𝒞)) = sum_n tr(lambda_n 𝒞)^2/(s - m_n0^2) (50 digits)":
            spectral_difference < HIGH_PRECISION_TOLERANCE,
        "3. 𝕄_s^-1(𝒞) = sum_n lambda_n tr(lambda_n 𝒞)/(2 (m_n0^2 - s)) (50 digits)":
            inverse_difference < HIGH_PRECISION_TOLERANCE,
        "5. -2 tr(𝒞 𝕄_s^-1(𝒞 - tr(lambda_n 𝒞) lambda_n/2)) = sum_{m != n} (50 digits)":
            projected_difference < HIGH_PRECISION_TOLERANCE,
    }


def rational_checks(epsilon=TEST_EPSILON, condensate=TEST_CONDENSATE,
                    second_order_scalar_trace=TEST_SECOND_ORDER_SCALAR_TRACE, base_sources=None):
    """Items 1, 3, 4, 8, 9 with exact rational sources at finite epsilon."""
    s = SQUARED_MOMENTUM
    sources = rational_test_sources(epsilon, base_sources)
    kinetic, mass_squared = bare_blocks(sources, condensate, second_order_scalar_trace)
    inverse_propagator_matrix = s * kinetic - mass_squared
    bare_polynomial = sp.Poly(sp.expand(inverse_propagator_matrix.det(method="berkowitz")), s)
    canonical = canonical_mass_matrix_from_sources(sources, condensate, second_order_scalar_trace)
    meson_mass_squared = mass_squared[:MESON_COUNT, :MESON_COUNT]
    invariant_vector = invariant_mixing_vector(**_block_inputs(sources))
    effective_bare = schur_complement(inverse_propagator_matrix)
    effective_canonical = (s - canonical[ALP, ALP] - epsilon**2 * (
        invariant_vector.T * (s * sp.eye(MESON_COUNT) - meson_mass_squared).LUsolve(invariant_vector))[0])
    quark_spurion = sources["quark_spurion"]
    tadpole = tadpole_source_from_sources(**_block_inputs(sources))
    effective_operator_form = (s - canonical[ALP, ALP] + 2 * epsilon**2 * sp.trace(
        tadpole * inverse_meson_operator(tadpole, quark_spurion, sources["anomaly_mass_squared"])))
    canonical_polynomial = sp.expand((s * sp.eye(MESON_COUNT + 1) - canonical).det(method="berkowitz"))
    quark_spurion_diagonal = [quark_spurion[index, index] for index in range(3)]
    explicit_polynomial = sp.Poly(explicit_secular_polynomial(
        quark_spurion_diagonal, sources["axial_source"] - sources["pseudoscalar_source"],
        sources["anomaly_mass_squared"], sources["gluon_coupling"] + sp.trace(sources["axial_source"]),
        epsilon, canonical[ALP, ALP]), s)
    kaon_mass_squared = (quark_spurion_diagonal[1] + quark_spurion_diagonal[2]) / 2
    meson_determinant = (s * sp.eye(MESON_COUNT) - meson_mass_squared).det(method="berkowitz")

    # item 4: rational kappa, char polynomial unchanged
    axial_shift = sp.Matrix([[sp.Rational(1, 7), 0, 0],
                             [0, sp.Rational(-2, 5), sp.Rational(1, 3) - sp.I / 6],
                             [0, sp.Rational(1, 3) + sp.I / 6, sp.Rational(3, 4)]])
    quark_mass = quark_spurion / (2 * condensate)
    pseudoscalar_first = anticommutator(sources["pseudoscalar_source"], quark_mass)
    shifted_trace = (second_order_scalar_trace - 4 * sp.trace(axial_shift * pseudoscalar_first)
                     - 4 * sp.trace(axial_shift * axial_shift * quark_mass))
    shifted_sources = dict(sources, axial_source=sources["axial_source"] + axial_shift,
                           pseudoscalar_source=sources["pseudoscalar_source"] + axial_shift,
                           gluon_coupling=sources["gluon_coupling"] - sp.trace(axial_shift))
    shifted_kinetic, shifted_mass = bare_blocks(shifted_sources, condensate, _simplify(shifted_trace))
    shifted_polynomial = sp.Poly(sp.expand((s * shifted_kinetic - shifted_mass).det(method="berkowitz")), s)

    # item 9: exact root isolation and solver
    intervals = explicit_polynomial.intervals(eps=sp.Rational(1, 10**40))
    exact_roots = sorted((lower + upper) / 2 for (lower, upper), multiplicity in intervals
                         for _ in range(multiplicity))
    spectrum = solve_spectrum_numerical(np.array(kinetic.evalf(30), dtype=float),
                                        np.array(mass_squared.evalf(30), dtype=float),
                                        field_labels=NEUTRAL_FIELD_LABELS)
    solver_poles = np.array(spectrum.pole_squared_masses)
    root_differences = np.abs(np.array([float(root) for root in exact_roots]) - np.sort(solver_poles))
    kaon_factor = sp.Poly(s - kaon_mass_squared, s)
    kaon_multiplicity, remaining_polynomial = 0, explicit_polynomial
    while remaining_polynomial.rem(kaon_factor).is_zero:
        kaon_multiplicity, remaining_polynomial = kaon_multiplicity + 1, remaining_polynomial.quo(kaon_factor)
    # roots other than m_K0,0^2 (c^AP_ds != 0: one K^0 root; c^AP_ds = 0: two, and D^eff_aa is regular there)
    other_roots = [root for root in exact_roots if abs(root - kaon_mass_squared) > sp.Rational(1, 10**30)]
    effective_values = [abs(sp.N(effective_canonical.subs(s, root), 50)) for root in other_roots]
    scale = max(abs(float(value)) for value in exact_roots)
    coupled_kaon_norm = sp.expand((sources["axial_source"] - sources["pseudoscalar_source"])[1, 2]
                                  * (sources["axial_source"] - sources["pseudoscalar_source"])[2, 1])
    expected_kaon_multiplicity = 1 if coupled_kaon_norm != 0 else 2
    return {
        "checks": {
            "1. det D(s) = det D_PP(s) D^eff_aa(s) (exact LO blocks)": _is_zero(
                sp.together(bare_polynomial.as_expr() - meson_determinant * effective_bare)),
            "3. bare Schur complement = canonical Schur complement": _is_zero(
                sp.together(effective_bare - effective_canonical)),
            "3. D^eff_aa = s - M_can,aa + 2 epsilon^2 tr(𝒞 𝕄_s^-1(𝒞)) (exact)": _is_zero(
                sp.cancel(effective_bare - effective_operator_form)),
            "6. det(s - M_can^2) = det(s K - M^2) = det(s - M^2_PP) D^eff_aa(s) (exact)": (
                _is_zero(canonical_polynomial - bare_polynomial.as_expr())
                and _is_zero(sp.together(canonical_polynomial - meson_determinant * effective_canonical))),
            "4. det(s K - M^2) unchanged under rational kappa": bare_polynomial == shifted_polynomial,
            "8. explicit secular polynomial = det(s K - M^2) (exact)": explicit_polynomial == bare_polynomial,
            "9. six real roots": len(exact_roots) == 6,
            "9'. m_K0,0^2 is an exact root, multiplicity 1 if c^AP_ds != 0 and 2 if c^AP_ds = 0": (
                kaon_multiplicity == expected_kaon_multiplicity),
            "9''. solver poles = isolated roots (rel. 1e-10)": bool(np.all(root_differences <= 1e-10 * scale)),
            "9'''. D^eff_aa = 0 at the other roots (|value| < 1e-30)": (
                len(other_roots) == 6 - kaon_multiplicity and all(value < 1e-30 for value in effective_values)),
        },
        "exact_roots": exact_roots,
        "solver_poles": solver_poles,
        "largest_root_difference": float(np.max(root_differences)),
        "largest_effective_value": max(effective_values),
        "polynomial": explicit_polynomial,
    }


def derive():
    """Run all checks."""
    return {
        "kinetic blocks": {"checks": kinetic_trace_checks()},
        "generic": {"checks": generic_schur_checks()},
        "trace form (neutral sources)": {"checks": trace_form_checks()},
        "trace form (flavour-diagonal sources)": {"checks": trace_form_checks(illustrative_lo_inputs("flavor_diagonal"))},
        "canonical frame (diagonal chi_0)": {"checks": canonical_frame_checks(True)},
        "canonical frame (chi_0 with d-s entry)": {"checks": canonical_frame_checks(False)},
        "explicit form (symbolic)": {"checks": explicit_form_symbolic_checks()},
        "residues (diagonal chi_0, rational)": {"checks": residue_checks()},
        "rational, epsilon = 1/2": rational_checks(),
        "rational, epsilon = 2 (algebraic test; M^2 not positive, one negative root)": rational_checks(epsilon=sp.Integer(2)),
    }


def main():
    result = derive()
    print("Exact pole equations (canonical frame, no expansion in epsilon):")
    print("  ALP:   s - M_can,aa + 2 epsilon^2 tr(𝒞 𝕄_s^-1(𝒞)) = 0,"
          "  M_can,aa = m_a0^2 + epsilon^2 [(2 m_0^2/N_f) c_GA^2 - B_0 tr hatS_2]")
    print("         -2 tr(𝒞 𝕄_s^-1(𝒞)) = sum_n tr(lambda_n 𝒞)^2/(s - m_n0^2)")
    print("  meson: (s - m_n0^2)[s - M_can,aa - epsilon^2 sum_{m != n} tr(lambda_m 𝒞)^2/(s - m_m0^2)]"
          " = epsilon^2 tr(lambda_n 𝒞)^2")
    for section, content in result.items():
        for label, passed in content["checks"].items():
            print(f"[{'ok' if passed else 'FAIL'}] {section}: {label}")
    for key in ("rational, epsilon = 1/2", "rational, epsilon = 2 (algebraic test; M^2 not positive, one negative root)"):
        content = result[key]
        print(f"{key}: roots {[sp.N(root, 12) for root in content['exact_roots']]}")
        print(f"  largest |solver - root| = {content['largest_root_difference']:.2e},"
              f"  largest |D^eff_aa(root)| = {sp.N(content['largest_effective_value'], 3)}")
    passed = all(all(content["checks"].values()) for content in result.values())
    print("All checks pass" if passed else "Some checks FAIL")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
