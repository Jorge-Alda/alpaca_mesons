#!/usr/bin/env python3
"""Exact LO relations for the meson and ALP poles: charged-meson zero shift,
exact K^0 pole, trace sum rule and product rule.

Run:   .venv/bin/python scripts/derive_alp_lo_mass_sum_rules.py
Tests: .venv/bin/python -m unittest discover -s scripts -p test_derive_alp_lo_mass_sum_rules.py

Conventions (Alda:2026nrd; memory/alp/shared/finite-alp-meson-mixing.md,
memory/alp/shared/alp-lo-charged-meson-decoupling.md,
memory/alp/shared/alp-lo-neutral-block-reduction.md):
N_f = 3, Hermitian U(3) generators lambda_0..lambda_8 with
tr(lambda_i lambda_j) = 2 delta_ij (ordering of scripts/symbolic_utils.py),
fields Phi = (P_i, a), L_2 = 1/2 dPhi^T K dPhi - 1/2 Phi^T M^2 Phi,
D(s) = s K - M^2. With r = f_pi/f_a:

  K_ij = delta_ij,  K_ia = -v_i,  K_aa = 1 + v.v,  v_i = r tr(c_q^A lambda_i),
  M^2_ij = tr({lambda_i, lambda_j} chi_0)/4 + m_0^2 tr(lambda_i) tr(lambda_j)/(2 N_f),
  M^2_ia = -(r/2) tr(lambda_i {c_q^P, chi_0}) + r m_0^2 c_G tr(lambda_i)/N_f,
  M^2_aa = m_a0^2 + 2 r^2 m_0^2 c_G^2/N_f - r^2 B_0 tr S_2,

with chi_0 = 2 B_0 m_q. The blocks without the last term are built by
`closed_form_blocks` of scripts/derive_alp_lo_neutral_blocks.py (not
re-implemented here), which assumes a quark spurion with no second derivative
in a. The term -r^2 B_0 tr S_2 = -(r^2/2) Re tr chi_2 is the only contribution
of the second Taylor coefficient chi_2 = d^2 chi/dx^2 |_{x=0} = 2 B_0 (S_2 + i P_2)
(x = a/f_a, S_2 and P_2 Hermitian) to the quadratic blocks; it is added to
M^2_aa by `add_second_order_scalar_term` and re-derived in section 5 from
`quadratic_blocks` of scripts/derive_alp_lo_joint_vertices.py.

Assumptions: LO chiral Lagrangian plus a canonical bare ALP term, no NLO
terms, no loops. c_q^A, c_q^P, chi_0 are Hermitian and commute with
Q = diag(2/3, -1/3, -1/3): only diagonal entries and a complex d-s entry.
The second-order spurion enters only through the real parameter tr S_2
(tr S_2 = 0 reproduces the spurion linear in a). No expansion in f_pi/f_a is made.

Derived and checked here (symbolic sources unless stated):
 1. Charged mesons (P_1, P_2 = pi^+-, P_4, P_5 = K^+-) have zero K and M^2
    entries to (P_0, P_3, P_6, P_7, P_8, a). Their block is K = 1 and M^2
    independent of c_q^A, c_q^P, c_G, f_a, m_a0^2, m_0^2: zero mass shift.
    For diagonal chi_0 = 2 B_0 diag(m_u, m_d, m_s) the charged block is
    diag(B_0 (m_u+m_d), B_0 (m_u+m_d), B_0 (m_u+m_s), B_0 (m_u+m_s)).
 2. With chi_0 diagonal, s = m_K0,0^2 = [(chi_0)_dd + (chi_0)_ss]/2 is a pole
    of the 6x6 neutral block with null vector (-D_7a, D_6a) on (P_6, P_7) and
    zero ALP component.
 3. Trace sum rule: sum_n m_n^2 = tr(K^{-1} M^2) with
    K^{-1} = [[1 + v v^T, v], [v^T, 1]], and in the sources
      tr(K^{-1} M^2) = (chi_0)_uu + 2 (chi_0)_dd + 2 (chi_0)_ss + m_0^2 + m_a0^2
                       + r^2 [ 2 tr(c_q^A c_q^A chi_0) - 2 tr(c_q^A {c_q^P, chi_0})
                               + (2 m_0^2/N_f) (c_G + tr c_q^A)^2 - B_0 tr S_2 ].
    d tr(K^{-1} M^2) / d tr S_2 = -r^2 B_0 (since (K^{-1})_aa = 1).
 4. Product rule: prod_n m_n^2 = det M^2 / det K, det K = 1, and
      det M^2 = m_K0,0^4 { [det chi_0 + (m_0^2/N_f) e_2(chi_0)]
                            (m_a0^2 - 2 r^2 tr(c_q^P c_q^P chi_0) - r^2 B_0 tr S_2)
                          + 2 r^2 (m_0^2/N_f) det chi_0 (c_G + tr c_q^P)^2 },
    e_2(chi_0) = [(tr chi_0)^2 - tr(chi_0^2)]/2. det M^2 does not depend on c_q^A.
    det M^2 is linear in M^2_aa with cofactor det M^2_PP
    = m_K0,0^4 [det chi_0 + (m_0^2/N_f) e_2(chi_0)], so the tr S_2 term is
    -r^2 B_0 tr S_2 det M^2_PP. In both rules tr S_2 acts as the replacement
    m_a0^2 -> m_a0^2 - r^2 B_0 tr S_2.
    Items 3 and 4 are checked for chi_0 with a complex d-s entry as well.
 5. Taylor check with quadratic_blocks (derive_alp_lo_joint_vertices.py) on
    (P_3, P_8, P_0, P_6, P_7, a), rational chi_0, chi_1, C_R,0..1, C_L,0..1,
    Theta_1..2 and fully symbolic Hermitian S_2, P_2, symbolic f_pi, f_a, B_0:
    K is independent of chi_2; M^2(S_2, P_2) - M^2(0, 0) has only the aa entry
    -r^2 B_0 tr S_2; P_2 does not appear; S_2 = (tr S_2/3) 1, P_2 = 0 gives the
    same blocks as the general S_2, P_2.

 6. Invariant form (memory/alp/shared/alp-dressed-mass-second-order-invariants.md,
    real quark masses: S_0 = m_q, P_0 = 0, S_1 = i[c^S, m_q], P_1 = {c^P, m_q},
    C_V0 = c^V, C_A0 = c^A). hat S_2 is computed from its definition (Hermitian
    part of the second Taylor coefficient of W_R^dagger M W_L) and compared with
    the closed form of the record. With c^AP = c^A - c^P, c_GA = c_G + tr c^A:
      tr hat S_2 = tr S_2 + 4 tr(c^A P_1) - 4 tr(c^A c^A m_q),
      sum_n m_n^2 = chi_uu + 2 chi_dd + 2 chi_ss + m_0^2 + m_a0^2
                    + r^2 [ (2 m_0^2/N_f) c_GA^2 - B_0 tr hat S_2 ],
      det M^2 = m_K0,0^4 { [det chi_0 + (m_0^2/N_f) e_2(chi_0)]
                             (m_a0^2 - r^2 [2 tr(c^AP c^AP chi_0) + B_0 tr hat S_2])
                           + 2 r^2 (m_0^2/N_f) det chi_0 (c_GA - tr c^AP)^2 }.
    Checked against the forms with tr S_2 (symbolic and random rational Hermitian
    inputs) and against the exact Vieta sum and product. Invariance under
    c^V, c^S -> +delta, c^A, c^P -> +kappa, c_G -> c_G - tr kappa (the c_G shift is
    an input, not derived here): the inputs are transformed with
    g_eta = exp(i x theta_eta), theta_R = delta + kappa, theta_L = delta - kappa,
    M' = g_R M g_L^dagger, C_eta' = g_eta C_eta g_eta^dagger + theta_eta.

Comparison with scripts/solve_alp_meson_spectrum.py: the blocks are built
from exact rational sources (illustrative_lo_inputs("neutral") of
derive_alp_lo_neutral_blocks.py plus a variant with a d-s entry in chi_0);
the poles are obtained with `solve_spectrum_numerical` (float64, because the
exact solver cannot find all poles of a generic 6x6 system in closed form)
and compared with the exact closed forms at relative tolerance 1e-10. The
sum and product of the poles are also checked exactly by Vieta's formulas on
the exact characteristic polynomial det(s K - M^2). This is done for tr S_2 = 0
and for the rational test values B_0 = 27/10 GeV, tr S_2 = 1/20 GeV; for the
latter the exact differences of the Vieta sum and product with respect to
tr S_2 = 0 are compared with -r^2 B_0 tr S_2 and -r^2 B_0 tr S_2 det M^2_PP.

Units: f_pi, f_a, B_0, tr S_2 in GeV; chi_0, m_0^2, m_a0^2, M^2, s in GeV^2;
K and the couplings are dimensionless. The rational inputs are test values, not fits.
"""

from __future__ import annotations

import random

import numpy as np
import sympy as sp

from derive_alp_lo_neutral_blocks import (
    ALL_GENERATOR_INDICES,
    CHARGED_GENERATOR_INDICES,
    FLAVOR_COUNT,
    NEUTRAL_GENERATOR_INDICES,
    closed_form_blocks,
    illustrative_lo_inputs,
)
from derive_alp_lo_joint_vertices import quadratic_blocks
from solve_alp_meson_spectrum import inverse_propagator, solve_spectrum_numerical
from symbolic_utils import u3_generators

KAON_POSITIONS_IN_NEUTRAL_BASIS = (3, 4)  # (P_6, P_7) in (P_3, P_8, P_0, P_6, P_7, a)
ALP_POSITION_IN_NEUTRAL_BASIS = 5


def _hermitian_charge_conserving(name, diagonal_only=False):
    """Symbolic Hermitian 3x3 matrix with real diagonal entries and a complex d-s entry."""
    up, down, strange = sp.symbols(f"{name}_uu {name}_dd {name}_ss", real=True)
    matrix = sp.diag(up, down, strange)
    if not diagonal_only:
        real_part, imaginary_part = sp.symbols(f"Re{name}_ds Im{name}_ds", real=True)
        matrix[1, 2] = real_part + sp.I * imaginary_part
        matrix[2, 1] = real_part - sp.I * imaginary_part
    return matrix


def symbolic_sources(diagonal_quark_spurion=True):
    """Symbolic charge-conserving LO inputs (keyword arguments of closed_form_blocks)."""
    return {
        "axial_source": _hermitian_charge_conserving("cA"),
        "pseudoscalar_source": _hermitian_charge_conserving("cP"),
        "quark_spurion": _hermitian_charge_conserving("chi0", diagonal_only=diagonal_quark_spurion),
        "gluon_coupling": sp.Symbol("c_G", real=True),
        "pion_decay_constant": sp.Symbol("f_pi", positive=True),
        "alp_decay_constant": sp.Symbol("f_a", positive=True),
        "anomaly_mass_squared": sp.Symbol("m0sq", positive=True),
        "bare_alp_mass_squared": sp.Symbol("ma0sq", positive=True),
    }


def symbolic_second_order_scalar():
    """Symbols (B_0, tr S_2) for the second-order scalar spurion term."""
    return sp.Symbol("B_0", positive=True), sp.Symbol("trS_2", real=True)


RATIONAL_CONDENSATE = sp.Rational(27, 10)  # B_0 in GeV (test value)
RATIONAL_SECOND_ORDER_SCALAR_TRACE = sp.Rational(1, 20)  # tr S_2 in GeV (test value)


def add_second_order_scalar_term(mass_squared, pion_decay_constant, alp_decay_constant,
                                 condensate, second_order_scalar_trace):
    """Copy of M^2 with -(f_pi/f_a)^2 B_0 tr S_2 added to the aa entry (last row and column)."""
    shifted = sp.Matrix(mass_squared)
    shifted[-1, -1] += -(pion_decay_constant / alp_decay_constant)**2 * condensate * second_order_scalar_trace
    return shifted


def blocks_with_second_order_scalar(sources, condensate, second_order_scalar_trace,
                                    generator_indices=NEUTRAL_GENERATOR_INDICES):
    """(K, M^2) from closed_form_blocks plus the tr S_2 term in M^2_aa."""
    kinetic, mass_squared = closed_form_blocks(**sources, generator_indices=generator_indices)
    return kinetic, add_second_order_scalar_term(mass_squared, sources["pion_decay_constant"],
                                                 sources["alp_decay_constant"], condensate,
                                                 second_order_scalar_trace)


def rational_sources(quark_spurion_ds_entry=0):
    """Exact rational test inputs: illustrative_lo_inputs("neutral"), optionally with a d-s entry in chi_0."""
    inputs = dict(illustrative_lo_inputs("neutral"))
    quark_spurion = sp.Matrix(inputs["quark_spurion"])
    quark_spurion[1, 2] = quark_spurion_ds_entry
    quark_spurion[2, 1] = sp.conjugate(quark_spurion_ds_entry)
    inputs["quark_spurion"] = quark_spurion
    return inputs


def _simplify(expression):
    return sp.simplify(sp.expand(sp.radsimp(sp.expand(expression))))


def _is_zero_matrix(matrix):
    return all(_simplify(entry) == 0 for entry in matrix)


# ------------------------------------------------------------ 1. charged mesons


def charged_meson_checks():
    """Charged-neutral entries vanish; charged block K = 1 and M^2 depends only on chi_0."""
    sources = symbolic_sources(diagonal_quark_spurion=False)
    kinetic, mass_squared = closed_form_blocks(**sources, generator_indices=ALL_GENERATOR_INDICES)
    charged = list(CHARGED_GENERATOR_INDICES)  # row index = generator index in the 10x10 basis
    neutral = [index for index in range(10) if index not in charged]
    charged_block_mass = mass_squared.extract(charged, charged).applyfunc(_simplify)
    alp_parameters = {sources[key] for key in ("gluon_coupling", "alp_decay_constant",
                                               "anomaly_mass_squared", "bare_alp_mass_squared")}
    alp_parameters |= sources["axial_source"].free_symbols | sources["pseudoscalar_source"].free_symbols
    checks = {
        "K charged-neutral entries vanish": _is_zero_matrix(kinetic.extract(charged, neutral)),
        "M^2 charged-neutral entries vanish": _is_zero_matrix(mass_squared.extract(charged, neutral)),
        "K charged block is the unit matrix": _is_zero_matrix(kinetic.extract(charged, charged) - sp.eye(4)),
        "M^2 charged block independent of ALP inputs and m_0^2":
            not (charged_block_mass.free_symbols & alp_parameters),
    }
    up_mass, down_mass, strange_mass, condensate = sp.symbols("m_u m_d m_s B_0", positive=True)
    quark_spurion = 2 * condensate * sp.diag(up_mass, down_mass, strange_mass)
    unmixed_charged_block = charged_block_mass.subs({
        sources["quark_spurion"][0, 0]: quark_spurion[0, 0],
        sources["quark_spurion"][1, 1]: quark_spurion[1, 1],
        sources["quark_spurion"][2, 2]: quark_spurion[2, 2],
        sp.Symbol("Rechi0_ds", real=True): 0, sp.Symbol("Imchi0_ds", real=True): 0,
    }).applyfunc(_simplify)
    charged_pion_mass_squared = condensate * (up_mass + down_mass)
    charged_kaon_mass_squared = condensate * (up_mass + strange_mass)
    checks["charged block = diag(B0(mu+md), B0(mu+md), B0(mu+ms), B0(mu+ms))"] = _is_zero_matrix(
        unmixed_charged_block - sp.diag(charged_pion_mass_squared, charged_pion_mass_squared,
                                        charged_kaon_mass_squared, charged_kaon_mass_squared))
    return {
        "checks": checks,
        "charged_block_mass_squared": charged_block_mass,
        "charged_pion_mass_squared": charged_pion_mass_squared,
        "charged_kaon_mass_squared": charged_kaon_mass_squared,
    }


# --------------------------------------------------------------- 2. K^0 pole


def neutral_kaon_mass_squared(quark_spurion):
    """m_K0,0^2 = [(chi_0)_dd + (chi_0)_ss]/2."""
    return (quark_spurion[1, 1] + quark_spurion[2, 2]) / 2


def neutral_kaon_null_vector(kinetic, mass_squared, quark_spurion):
    """Vector (-D_7a, D_6a) on (P_6, P_7), zero elsewhere, at s = m_K0,0^2."""
    propagator = inverse_propagator(kinetic, mass_squared, neutral_kaon_mass_squared(quark_spurion))
    position_6, position_7 = KAON_POSITIONS_IN_NEUTRAL_BASIS
    null_vector = sp.zeros(6, 1)
    null_vector[position_6] = -propagator[position_7, ALP_POSITION_IN_NEUTRAL_BASIS]
    null_vector[position_7] = propagator[position_6, ALP_POSITION_IN_NEUTRAL_BASIS]
    return propagator, null_vector


def neutral_kaon_pole_checks():
    """Exact K^0 pole for diagonal chi_0 and symbolic complex d-s entries of c_q^A, c_q^P."""
    sources = symbolic_sources(diagonal_quark_spurion=True)
    kinetic, mass_squared = closed_form_blocks(**sources)
    propagator, null_vector = neutral_kaon_null_vector(kinetic, mass_squared, sources["quark_spurion"])
    return {
        "checks": {
            "D(m_K0,0^2) (-D_7a, D_6a, 0...) = 0": _is_zero_matrix(propagator * null_vector),
            "null vector is nonzero for generic sources": not _is_zero_matrix(null_vector),
            "null vector has zero ALP component": null_vector[ALP_POSITION_IN_NEUTRAL_BASIS] == 0,
        },
        "null_vector": null_vector.applyfunc(_simplify),
    }


# ------------------------------------------------------ 3. trace sum rule


def alp_kinetic_vector(kinetic):
    """v_i = -K_ia (column, neutral mesons only)."""
    return -kinetic[:5, ALP_POSITION_IN_NEUTRAL_BASIS]


def kinetic_inverse_closed_form(kinetic_vector):
    """K^{-1} = [[1 + v v^T, v], [v^T, 1]]."""
    upper_left = sp.eye(5) + kinetic_vector * kinetic_vector.T
    return sp.Matrix(sp.BlockMatrix([[upper_left, kinetic_vector], [kinetic_vector.T, sp.Matrix([[1]])]]))


def trace_sum_rule_block_form(kinetic, mass_squared):
    """tr M^2_PP + v^T M^2_PP v + 2 v.M^2_Pa + M^2_aa."""
    kinetic_vector = alp_kinetic_vector(kinetic)
    meson_block = mass_squared[:5, :5]
    mixed_column = mass_squared[:5, ALP_POSITION_IN_NEUTRAL_BASIS]
    return (meson_block.trace() + (kinetic_vector.T * meson_block * kinetic_vector)[0]
            + 2 * (kinetic_vector.T * mixed_column)[0] + mass_squared[5, 5])


def trace_sum_rule_closed_form(axial_source, pseudoscalar_source, quark_spurion, gluon_coupling,
                               pion_decay_constant, alp_decay_constant, anomaly_mass_squared,
                               bare_alp_mass_squared, condensate=0, second_order_scalar_trace=0):
    """tr(K^{-1} M^2) in the sources (neutral 6x6 block), including -r^2 B_0 tr S_2."""
    ratio = pion_decay_constant / alp_decay_constant
    unmixed_meson_sum = quark_spurion[0, 0] + 2 * quark_spurion[1, 1] + 2 * quark_spurion[2, 2] + anomaly_mass_squared
    mixing_term = ratio**2 * (
        2 * sp.trace(axial_source * axial_source * quark_spurion)
        - 2 * sp.trace(axial_source * (pseudoscalar_source * quark_spurion + quark_spurion * pseudoscalar_source))
        + 2 * anomaly_mass_squared / FLAVOR_COUNT * (gluon_coupling + sp.trace(axial_source))**2
        - condensate * second_order_scalar_trace)
    return unmixed_meson_sum + bare_alp_mass_squared + mixing_term


def trace_sum_rule_checks(diagonal_quark_spurion):
    sources = symbolic_sources(diagonal_quark_spurion)
    condensate, second_order_scalar_trace = symbolic_second_order_scalar()
    kinetic, mass_squared = blocks_with_second_order_scalar(sources, condensate, second_order_scalar_trace)
    kinetic_inverse = kinetic_inverse_closed_form(alp_kinetic_vector(kinetic))
    trace_exact = (kinetic_inverse * mass_squared).trace()
    ratio = sources["pion_decay_constant"] / sources["alp_decay_constant"]
    return {
        "K K^{-1} = 1 with K^{-1} = [[1+vv^T, v],[v^T, 1]]": _is_zero_matrix(kinetic * kinetic_inverse - sp.eye(6)),
        "tr(K^{-1} M^2) = block form": _simplify(trace_exact - trace_sum_rule_block_form(kinetic, mass_squared)) == 0,
        "tr(K^{-1} M^2) = closed form in sources (with tr S_2)": _simplify(
            trace_exact - trace_sum_rule_closed_form(**sources, condensate=condensate,
                                                     second_order_scalar_trace=second_order_scalar_trace)) == 0,
        "d tr(K^{-1} M^2) / d trS_2 = -r^2 B_0": _simplify(
            sp.diff(trace_exact, second_order_scalar_trace) + ratio**2 * condensate) == 0,
    }


# --------------------------------------------------------- 4. product rule


def second_symmetric_polynomial(matrix):
    """e_2 = [(tr X)^2 - tr(X^2)]/2."""
    return (sp.trace(matrix)**2 - sp.trace(matrix * matrix)) / 2


def product_rule_closed_form(axial_source, pseudoscalar_source, quark_spurion, gluon_coupling,
                             pion_decay_constant, alp_decay_constant, anomaly_mass_squared,
                             bare_alp_mass_squared, condensate=0, second_order_scalar_trace=0):
    """det M^2 of the neutral 6x6 block in the sources (independent of c_q^A), including tr S_2."""
    del axial_source  # det M^2 does not involve K
    ratio = pion_decay_constant / alp_decay_constant
    anomaly_term = anomaly_mass_squared / FLAVOR_COUNT
    flavor_diagonal_determinant = quark_spurion.det() + anomaly_term * second_symmetric_polynomial(quark_spurion)
    return neutral_kaon_mass_squared(quark_spurion)**2 * (
        flavor_diagonal_determinant * (bare_alp_mass_squared
                                       - 2 * ratio**2 * sp.trace(pseudoscalar_source**2 * quark_spurion)
                                       - ratio**2 * condensate * second_order_scalar_trace)
        + 2 * ratio**2 * anomaly_term * quark_spurion.det() * (gluon_coupling + sp.trace(pseudoscalar_source))**2)


def product_rule_checks(diagonal_quark_spurion):
    sources = symbolic_sources(diagonal_quark_spurion)
    condensate, second_order_scalar_trace = symbolic_second_order_scalar()
    kinetic, mass_squared = blocks_with_second_order_scalar(sources, condensate, second_order_scalar_trace)
    ratio = sources["pion_decay_constant"] / sources["alp_decay_constant"]
    mass_determinant = sp.expand(mass_squared.det(method="berkowitz"))
    meson_determinant = mass_squared[:5, :5].det(method="berkowitz")
    alp_entry = sp.Symbol("Maa_generic")
    generic_alp_entry = sp.Matrix(mass_squared)
    generic_alp_entry[5, 5] = alp_entry
    generic_determinant = sp.expand(generic_alp_entry.det(method="berkowitz"))
    flavor_diagonal_determinant = (sources["quark_spurion"].det()
                                   + sources["anomaly_mass_squared"] / FLAVOR_COUNT
                                   * second_symmetric_polynomial(sources["quark_spurion"]))
    return {
        "det K = 1": _simplify(kinetic.det(method="berkowitz") - 1) == 0,
        "det M^2_PP = m_K0,0^4 [det chi_0 + (m_0^2/N_f) e_2(chi_0)]": _simplify(
            meson_determinant - neutral_kaon_mass_squared(sources["quark_spurion"])**2 * flavor_diagonal_determinant) == 0,
        "det M^2 = closed form in sources (with tr S_2)": _simplify(
            mass_determinant - product_rule_closed_form(**sources, condensate=condensate,
                                                        second_order_scalar_trace=second_order_scalar_trace)) == 0,
        "det M^2 independent of c_q^A": not (mass_determinant.free_symbols
                                             & sources["axial_source"].free_symbols),
        "det M^2 linear in M^2_aa with coefficient det M^2_PP": (
            sp.diff(generic_determinant, alp_entry, 2) == 0
            and _simplify(sp.diff(generic_determinant, alp_entry) - meson_determinant) == 0),
        "d det M^2 / d trS_2 = -r^2 B_0 det M^2_PP, second derivative 0": (
            sp.diff(mass_determinant, second_order_scalar_trace, 2) == 0
            and _simplify(sp.diff(mass_determinant, second_order_scalar_trace)
                          + ratio**2 * condensate * meson_determinant) == 0),
    }


# ------------------------------------------------- comparison with the solver


def _exact_pole_sum_and_product(kinetic, mass_squared):
    """Vieta sum and product of the roots of det(s K - M^2) (degree 6: (-1)^6 = +1)."""
    squared_momentum = sp.Symbol("s")
    characteristic = sp.Poly(sp.expand(
        inverse_propagator(kinetic, mass_squared, squared_momentum).det(method="berkowitz")), squared_momentum)
    coefficients = characteristic.all_coeffs()
    return (_simplify(-coefficients[1] / coefficients[0]), _simplify(coefficients[-1] / coefficients[0]),
            characteristic)


def solver_comparison(sources, condensate=0, second_order_scalar_trace=0, relative_tolerance=1e-10):
    """Sum and product of poles: exact Vieta and float64 solver vs closed forms."""
    kinetic, mass_squared = blocks_with_second_order_scalar(sources, condensate, second_order_scalar_trace)
    exact_sum, exact_product, characteristic = _exact_pole_sum_and_product(kinetic, mass_squared)
    second_order = {"condensate": condensate, "second_order_scalar_trace": second_order_scalar_trace}
    expected_sum = _simplify(trace_sum_rule_closed_form(**sources, **second_order))
    expected_product = _simplify(product_rule_closed_form(**sources, **second_order))
    spectrum = solve_spectrum_numerical(np.array(kinetic.evalf(30), dtype=float),
                                        np.array(mass_squared.evalf(30), dtype=float))
    poles = np.array(spectrum.pole_squared_masses)
    kaon_pole = float(neutral_kaon_mass_squared(sources["quark_spurion"]))
    diagonal_quark_spurion = sources["quark_spurion"][1, 2] == 0
    checks = {
        "Vieta sum = trace closed form (exact)": _simplify(exact_sum - expected_sum) == 0,
        "Vieta product = det closed form (exact)": _simplify(exact_product - expected_product) == 0,
        "solver sum = trace closed form": abs(poles.sum() - float(expected_sum)) <= relative_tolerance * abs(float(expected_sum)),
        "solver product = det closed form": abs(np.prod(poles) - float(expected_product)) <= relative_tolerance * abs(float(expected_product)),
    }
    if second_order_scalar_trace != 0:
        reference_sum, reference_product, _ = _exact_pole_sum_and_product(*closed_form_blocks(**sources))
        ratio = sources["pion_decay_constant"] / sources["alp_decay_constant"]
        shift = -ratio**2 * condensate * second_order_scalar_trace
        checks["Vieta sum(trS_2) - sum(0) = -r^2 B_0 trS_2 (exact)"] = _simplify(exact_sum - reference_sum - shift) == 0
        checks["Vieta product(trS_2) - product(0) = -r^2 B_0 trS_2 det M^2_PP (exact)"] = _simplify(
            exact_product - reference_product - shift * mass_squared[:5, :5].det(method="berkowitz")) == 0
    if diagonal_quark_spurion:
        closest = int(np.argmin(np.abs(poles - kaon_pole)))
        checks["det D(m_K0,0^2) = 0 (exact)"] = characteristic.eval(neutral_kaon_mass_squared(sources["quark_spurion"])) == 0
        checks["solver has a pole at m_K0,0^2"] = abs(poles[closest] - kaon_pole) <= relative_tolerance * kaon_pole
        checks["that pole vector has zero ALP component"] = abs(spectrum.pole_vectors[ALP_POSITION_IN_NEUTRAL_BASIS, closest]) <= 1e-10
    return {"checks": checks, "poles": poles, "exact_sum": exact_sum, "exact_product": exact_product}


def charged_solver_comparison(sources, relative_tolerance=1e-10):
    """10x10 system: the charged poles are (chi_uu+chi_dd)/2 and (chi_uu+chi_ss)/2 twice, with no ALP component."""
    kinetic, mass_squared = closed_form_blocks(**sources, generator_indices=ALL_GENERATOR_INDICES)
    spectrum = solve_spectrum_numerical(np.array(kinetic.evalf(30), dtype=float),
                                        np.array(mass_squared.evalf(30), dtype=float))
    quark_spurion = sources["quark_spurion"]
    expected = [float((quark_spurion[0, 0] + quark_spurion[1, 1]) / 2), float((quark_spurion[0, 0] + quark_spurion[2, 2]) / 2)]
    poles = np.array(spectrum.pole_squared_masses)
    checks = {}
    for label, value in zip(("pi^+-", "K^+-"), expected):
        matching = [index for index, pole in enumerate(poles) if abs(pole - value) <= relative_tolerance * value]
        checks[f"{label} pole at unmixed value, twice"] = len(matching) == 2
        checks[f"{label} pole vectors have zero ALP and neutral components"] = all(
            np.max(np.abs(spectrum.pole_vectors[[0, 3, 6, 7, 8, 9], index])) <= 1e-10 for index in matching)
    return {"checks": checks}


# ------------------------------- 5. second-order spurion from quadratic_blocks


def _symbolic_hermitian(name):
    """Fully symbolic Hermitian 3x3 matrix."""
    matrix = sp.zeros(3)
    for row in range(3):
        matrix[row, row] = sp.Symbol(f"{name}_{row}{row}", real=True)
        for column in range(row + 1, 3):
            real_part = sp.Symbol(f"Re{name}_{row}{column}", real=True)
            imaginary_part = sp.Symbol(f"Im{name}_{row}{column}", real=True)
            matrix[row, column] = real_part + sp.I * imaginary_part
            matrix[column, row] = real_part - sp.I * imaginary_part
    return matrix


def second_order_spurion_taylor_checks():
    """Dependence of the exact quadratic blocks on chi_2 = 2 B_0 (S_2 + i P_2), from quadratic_blocks."""
    rational = sp.Rational
    pion_decay_constant, alp_decay_constant, condensate = sp.symbols("f_pi f_a B_0", positive=True)
    ratio = pion_decay_constant / alp_decay_constant
    generators = u3_generators()
    neutral_generators = [sp.Matrix(generators[index]) for index in NEUTRAL_GENERATOR_INDICES]
    quark_spurion = rational_sources(rational(1, 50) + sp.I * rational(1, 70))["quark_spurion"]
    first_derivative = sp.Matrix([[rational(1, 3), rational(1, 7) + sp.I / 5, 0],
                                  [sp.I / 9, rational(-1, 4) + sp.I / 2, rational(2, 5)],
                                  [rational(1, 6), -sp.I / 3, rational(3, 8)]])
    right_derivatives = [sp.diag(1, rational(1, 2), rational(-1, 3)), sp.diag(rational(1, 5), 0, rational(1, 4)), sp.zeros(3)]
    left_derivatives = [sp.diag(rational(1, 4), rational(-1, 2), rational(2, 3)), sp.diag(0, rational(1, 3), 0), sp.zeros(3)]
    theta_derivatives = [0, rational(1, 2), rational(3, 7)]
    scalar_second, pseudoscalar_second = _symbolic_hermitian("S2"), _symbolic_hermitian("P2")
    trace_only_scalar = sp.trace(scalar_second) / 3 * sp.eye(3)

    def blocks(second_derivative):
        return quadratic_blocks(neutral_generators, [quark_spurion, first_derivative, second_derivative],
                                right_derivatives, left_derivatives, theta_derivatives,
                                pion_decay_constant, alp_decay_constant, rational(4, 5), rational(3, 10))

    kinetic_general, mass_general = blocks(2 * condensate * (scalar_second + sp.I * pseudoscalar_second))
    kinetic_zero, mass_zero = blocks(sp.zeros(3))
    kinetic_trace_only, mass_trace_only = blocks(2 * condensate * trace_only_scalar)
    expected_difference = sp.zeros(6)
    expected_difference[5, 5] = -ratio**2 * condensate * sp.trace(scalar_second)
    mass_symbols = set().union(*(entry.free_symbols for entry in mass_general))
    return {
        "checks": {
            "K independent of chi_2": _is_zero_matrix(kinetic_general - kinetic_zero),
            "M^2(S_2, P_2) - M^2(0, 0) = -r^2 B_0 tr S_2 in aa only": _is_zero_matrix(
                mass_general - mass_zero - expected_difference),
            "P_2 does not appear in K or M^2": not (mass_symbols & pseudoscalar_second.free_symbols),
            "S_2 = (tr S_2/3) 1, P_2 = 0 gives the same K and M^2": (
                _is_zero_matrix(kinetic_general - kinetic_trace_only)
                and _is_zero_matrix(mass_general - mass_trace_only)),
        },
    }


# ----------------------- 6. invariant form with the dressed source hat S_2


def _commutator(first, second):
    return first * second - second * first


def _anticommutator(first, second):
    return first * second + second * first


def dressed_second_order_scalar_from_definition(quark_mass, first_order_scalar, first_order_pseudoscalar,
                                                second_order_scalar, second_order_pseudoscalar,
                                                vector_couplings, axial_couplings):
    """hat S_2 = Hermitian part of the second Taylor coefficient of W_R^dagger M W_L.

    M = S + i P with S_0 = m_q, P_0 = 0; vector_couplings = (C_V,0, C_V,1),
    axial_couplings = (C_A,0, C_A,1); C_R = C_V + C_A, C_L = C_V - C_A;
    W_1 = i C_0, W_2 = i C_1 - C_0^2 (memory/alp/shared/alp-dressed-mass-second-order-invariants.md).
    Leibniz rule for the product of three series, coefficients = derivatives at x = 0.
    """
    right = [vector + axial for vector, axial in zip(vector_couplings, axial_couplings)]
    left = [vector - axial for vector, axial in zip(vector_couplings, axial_couplings)]
    right_first_adjoint = -sp.I * right[0]
    right_second_adjoint = -sp.I * right[1] - right[0] * right[0]
    left_first = sp.I * left[0]
    left_second = sp.I * left[1] - left[0] * left[0]
    mass_first = first_order_scalar + sp.I * first_order_pseudoscalar
    mass_second = second_order_scalar + sp.I * second_order_pseudoscalar
    dressed_second = (right_second_adjoint * quark_mass + 2 * right_first_adjoint * mass_first + mass_second
                      + 2 * right_first_adjoint * quark_mass * left_first + quark_mass * left_second
                      + 2 * mass_first * left_first)
    return ((dressed_second + dressed_second.H) / 2).applyfunc(sp.expand)


def dressed_second_order_scalar_closed_form(quark_mass, first_order_scalar, first_order_pseudoscalar,
                                            second_order_scalar, vector_couplings, axial_couplings):
    """Specialization to real quark masses of memory/alp/shared/alp-dressed-mass-second-order-invariants.md:
    hat S_2 = S_2 - 2i[C_V0, S_1] + 2{C_A0, P_1} - i[C_V1, m_q] - [C_V0,[C_V0,m_q]] - {C_A0,{C_A0,m_q}}."""
    vector_0, vector_1 = vector_couplings
    axial_0 = axial_couplings[0]
    return (second_order_scalar - 2 * sp.I * _commutator(vector_0, first_order_scalar)
            + 2 * _anticommutator(axial_0, first_order_pseudoscalar) - sp.I * _commutator(vector_1, quark_mass)
            - _commutator(vector_0, _commutator(vector_0, quark_mass))
            - _anticommutator(axial_0, _anticommutator(axial_0, quark_mass))).applyfunc(sp.expand)


def dressed_second_order_scalar_trace(second_order_scalar_trace, axial_coupling, first_order_pseudoscalar, quark_mass):
    """tr hat S_2 = tr S_2 + 4 tr(C_A0 P_1) - 4 tr(C_A0 C_A0 m_q) (claim 1)."""
    return (second_order_scalar_trace + 4 * sp.trace(axial_coupling * first_order_pseudoscalar)
            - 4 * sp.trace(axial_coupling * axial_coupling * quark_mass))


def trace_sum_rule_invariant_form(axial_source, pseudoscalar_source, quark_spurion, gluon_coupling,
                                  pion_decay_constant, alp_decay_constant, anomaly_mass_squared,
                                  bare_alp_mass_squared, condensate, dressed_scalar_trace):
    """sum_n m_n^2 = chi_uu + 2 chi_dd + 2 chi_ss + m_0^2 + m_a0^2
                     + r^2 [ (2 m_0^2/N_f) c_GA^2 - B_0 tr hat S_2 ],  c_GA = c_G + tr c_q^A."""
    del pseudoscalar_source  # enters only through tr hat S_2
    ratio = pion_decay_constant / alp_decay_constant
    gluon_axial_combination = gluon_coupling + sp.trace(axial_source)
    return (quark_spurion[0, 0] + 2 * quark_spurion[1, 1] + 2 * quark_spurion[2, 2] + anomaly_mass_squared
            + bare_alp_mass_squared
            + ratio**2 * (2 * anomaly_mass_squared / FLAVOR_COUNT * gluon_axial_combination**2
                          - condensate * dressed_scalar_trace))


def product_rule_invariant_form(axial_source, pseudoscalar_source, quark_spurion, gluon_coupling,
                                pion_decay_constant, alp_decay_constant, anomaly_mass_squared,
                                bare_alp_mass_squared, condensate, dressed_scalar_trace):
    """det M^2 = m_K0,0^4 { [det chi_0 + (m_0^2/N_f) e_2(chi_0)]
                    (m_a0^2 - r^2 [2 tr(c^AP c^AP chi_0) + B_0 tr hat S_2])
                  + 2 r^2 (m_0^2/N_f) det chi_0 (c_GA - tr c^AP)^2 },
    c^AP = c_q^A - c_q^P, c_GA = c_G + tr c_q^A (so c_GA - tr c^AP = c_G + tr c_q^P)."""
    ratio = pion_decay_constant / alp_decay_constant
    anomaly_term = anomaly_mass_squared / FLAVOR_COUNT
    axial_minus_pseudoscalar = axial_source - pseudoscalar_source
    gluon_axial_combination = gluon_coupling + sp.trace(axial_source)
    flavor_diagonal_determinant = quark_spurion.det() + anomaly_term * second_symmetric_polynomial(quark_spurion)
    return neutral_kaon_mass_squared(quark_spurion)**2 * (
        flavor_diagonal_determinant * (bare_alp_mass_squared - ratio**2 * (
            2 * sp.trace(axial_minus_pseudoscalar * axial_minus_pseudoscalar * quark_spurion)
            + condensate * dressed_scalar_trace))
        + 2 * ratio**2 * anomaly_term * quark_spurion.det()
        * (gluon_axial_combination - sp.trace(axial_minus_pseudoscalar))**2)


def random_rational_hermitian(generator, charge_conserving=False):
    """Exact rational Hermitian 3x3 matrix (entries p/q, |p| <= 9, 1 <= q <= 9).
    charge_conserving: only diagonal entries and a complex d-s entry."""
    def number():
        return sp.Rational(generator.randint(-9, 9), generator.randint(1, 9))
    matrix = sp.zeros(3)
    for row in range(3):
        matrix[row, row] = number()
        for column in range(row + 1, 3):
            if charge_conserving and row == 0:
                continue
            matrix[row, column] = number() + sp.I * number()
            matrix[column, row] = sp.conjugate(matrix[row, column])
    return matrix


def first_and_second_order_inputs(axial_source, pseudoscalar_source, quark_mass, scalar_source, vector_source,
                                  vector_coupling_first, axial_coupling_first, second_order_scalar,
                                  second_order_pseudoscalar):
    """Inputs of hat S_2 in the Alda:2026nrd first-order parametrization:
    S_1 = i[c^S, m_q], P_1 = {c^P, m_q}, C_V0 = c^V, C_A0 = c^A; C_V1, C_A1, S_2, P_2 free."""
    return {
        "quark_mass": quark_mass,
        "first_order_scalar": sp.I * _commutator(scalar_source, quark_mass),
        "first_order_pseudoscalar": _anticommutator(pseudoscalar_source, quark_mass),
        "second_order_scalar": second_order_scalar,
        "second_order_pseudoscalar": second_order_pseudoscalar,
        "vector_couplings": (vector_source, vector_coupling_first),
        "axial_couplings": (axial_source, axial_coupling_first),
    }


def constant_rotation(series, vector_shift, axial_shift):
    """Inputs after g_eta = exp(i x theta_eta), theta_R = delta + kappa, theta_L = delta - kappa:
    M' = g_R M g_L^dagger and C_eta' = g_eta C_eta g_eta^dagger + theta_eta, to second (M) and first (C) order."""
    right_generator, left_generator = vector_shift + axial_shift, vector_shift - axial_shift
    quark_mass = series["quark_mass"]
    mass_first = series["first_order_scalar"] + sp.I * series["first_order_pseudoscalar"]
    mass_second = series["second_order_scalar"] + sp.I * series["second_order_pseudoscalar"]
    # derivatives of g_R = 1 + i x theta_R - x^2 theta_R^2/2 + ... and of g_L^dagger
    right_1, right_2 = sp.I * right_generator, -right_generator**2
    left_adjoint_1, left_adjoint_2 = -sp.I * left_generator, -left_generator**2
    rotated_first = right_1 * quark_mass + mass_first + quark_mass * left_adjoint_1
    rotated_second = (right_2 * quark_mass + 2 * right_1 * mass_first + mass_second
                      + 2 * right_1 * quark_mass * left_adjoint_1 + quark_mass * left_adjoint_2
                      + 2 * mass_first * left_adjoint_1)
    rotated_first, rotated_second = rotated_first.applyfunc(sp.expand), rotated_second.applyfunc(sp.expand)
    # C_eta,0' = C_eta,0 + theta_eta, C_eta,1' = C_eta,1 + i[theta_eta, C_eta,0]; C_V,A = (C_R +- C_L)/2
    (vector_0, vector_1), (axial_0, axial_1) = series["vector_couplings"], series["axial_couplings"]
    right = (vector_0 + axial_0 + right_generator,
             vector_1 + axial_1 + sp.I * _commutator(right_generator, vector_0 + axial_0))
    left = (vector_0 - axial_0 + left_generator,
            vector_1 - axial_1 + sp.I * _commutator(left_generator, vector_0 - axial_0))
    rotated_couplings = {
        "vector_couplings": tuple(((right_term + left_term) / 2).applyfunc(sp.expand)
                                  for right_term, left_term in zip(right, left)),
        "axial_couplings": tuple(((right_term - left_term) / 2).applyfunc(sp.expand)
                                 for right_term, left_term in zip(right, left)),
    }
    return {
        "quark_mass": quark_mass,
        "first_order_scalar": ((rotated_first + rotated_first.H) / 2).applyfunc(sp.expand),
        "first_order_pseudoscalar": ((rotated_first - rotated_first.H) / (2 * sp.I)).applyfunc(sp.expand),
        "second_order_scalar": ((rotated_second + rotated_second.H) / 2).applyfunc(sp.expand),
        "second_order_pseudoscalar": ((rotated_second - rotated_second.H) / (2 * sp.I)).applyfunc(sp.expand),
        **rotated_couplings,
    }


RANDOM_SEEDS = (11, 23, 37)


def _random_series(generator, charge_conserving=False):
    """Random exact rational c^A, c^P, c^S, c^V, m_q and second-order inputs: (hat S_2 inputs, couplings)."""
    axial_source, pseudoscalar_source, scalar_source, vector_source, quark_mass = (
        random_rational_hermitian(generator, charge_conserving) for _ in range(5))
    series = first_and_second_order_inputs(
        axial_source, pseudoscalar_source, quark_mass, scalar_source, vector_source,
        *(random_rational_hermitian(generator, charge_conserving) for _ in range(4)))
    return series, {"axial_source": axial_source, "pseudoscalar_source": pseudoscalar_source,
                    "scalar_source": scalar_source, "vector_source": vector_source}


def _random_rule_sources(generator, couplings, quark_mass):
    """Sources of the closed forms with chi_0 = 2 B_0 m_q and rational B_0, c_G (test values)."""
    condensate = sp.Rational(generator.randint(1, 9), generator.randint(1, 9))
    return condensate, {
        "axial_source": couplings["axial_source"], "pseudoscalar_source": couplings["pseudoscalar_source"],
        "quark_spurion": 2 * condensate * quark_mass,
        "gluon_coupling": sp.Rational(generator.randint(-9, 9), generator.randint(1, 9)),
        "pion_decay_constant": sp.Rational(13, 100), "alp_decay_constant": sp.Rational(7, 3),
        "anomaly_mass_squared": sp.Rational(4, 5), "bare_alp_mass_squared": sp.Rational(3, 10),
    }


def dressed_scalar_trace_checks():
    """Claim 1: definition of hat S_2 = record closed form, and its trace identity."""
    symbolic_inputs = {
        "quark_mass": _symbolic_hermitian("mq"),
        "first_order_scalar": _symbolic_hermitian("S1"),
        "first_order_pseudoscalar": _symbolic_hermitian("P1"),
        "second_order_scalar": _symbolic_hermitian("S2"),
        "second_order_pseudoscalar": _symbolic_hermitian("P2"),
        "vector_couplings": (_symbolic_hermitian("CV0"), _symbolic_hermitian("CV1")),
        "axial_couplings": (_symbolic_hermitian("CA0"), _symbolic_hermitian("CA1")),
    }
    checks = {}

    def compare(label, inputs):
        from_definition = dressed_second_order_scalar_from_definition(**inputs)
        closed_inputs = {key: value for key, value in inputs.items() if key != "second_order_pseudoscalar"}
        checks[f"{label}: hat S_2 definition = record closed form"] = _is_zero_matrix(
            from_definition - dressed_second_order_scalar_closed_form(**closed_inputs))
        checks[f"{label}: tr hat S_2 = tr S_2 + 4 tr(C_A0 P_1) - 4 tr(C_A0 C_A0 m_q)"] = sp.expand(
            from_definition.trace() - dressed_second_order_scalar_trace(
                inputs["second_order_scalar"].trace(), inputs["axial_couplings"][0],
                inputs["first_order_pseudoscalar"], inputs["quark_mass"])) == 0

    compare("symbolic generic Hermitian inputs", symbolic_inputs)
    for seed in RANDOM_SEEDS:
        series, _ = _random_series(random.Random(seed))
        compare(f"rational seed {seed}, P_1 = {{c^P, m_q}}", series)
    return checks


def invariant_form_symbolic_checks(diagonal_quark_spurion):
    """Claims 2-4 with symbolic charge-conserving sources; tr hat S_2 from claim 1 with P_1 = {c^P, m_q}."""
    sources = symbolic_sources(diagonal_quark_spurion)
    condensate, second_order_scalar_trace = symbolic_second_order_scalar()
    axial_source, pseudoscalar_source, quark_spurion = (
        sources["axial_source"], sources["pseudoscalar_source"], sources["quark_spurion"])
    quark_mass = quark_spurion / (2 * condensate)
    dressed_trace = dressed_second_order_scalar_trace(
        second_order_scalar_trace, axial_source, _anticommutator(pseudoscalar_source, quark_mass), quark_mass)
    ratio_squared = (sources["pion_decay_constant"] / sources["alp_decay_constant"])**2
    old = {"condensate": condensate, "second_order_scalar_trace": second_order_scalar_trace}
    new = {"condensate": condensate, "dressed_scalar_trace": dressed_trace}
    return {
        "r^2 B0 trS2 = r^2 [B0 tr hatS2 - 2 tr(cA{cP,chi0}) + 2 tr(cA cA chi0)]": _simplify(
            ratio_squared * condensate * second_order_scalar_trace - ratio_squared * (
                condensate * dressed_trace
                - 2 * sp.trace(axial_source * _anticommutator(pseudoscalar_source, quark_spurion))
                + 2 * sp.trace(axial_source * axial_source * quark_spurion))) == 0,
        "trace rule: invariant form = form with tr S_2": _simplify(
            trace_sum_rule_invariant_form(**sources, **new) - trace_sum_rule_closed_form(**sources, **old)) == 0,
        "product rule: invariant form = form with tr S_2": _simplify(
            product_rule_invariant_form(**sources, **new) - product_rule_closed_form(**sources, **old)) == 0,
    }


def invariant_form_random_checks():
    """Claims 3-4 for general (not charge-conserving) rational Hermitian inputs, tr hat S_2 from the definition."""
    checks = {}
    for seed in RANDOM_SEEDS:
        generator = random.Random(seed)
        series, couplings = _random_series(generator)
        condensate, sources = _random_rule_sources(generator, couplings, series["quark_mass"])
        dressed_trace = dressed_second_order_scalar_from_definition(**series).trace()
        old = {"condensate": condensate, "second_order_scalar_trace": series["second_order_scalar"].trace()}
        new = {"condensate": condensate, "dressed_scalar_trace": dressed_trace}
        checks[f"seed {seed}: trace rule invariant form = form with tr S_2"] = sp.expand(
            trace_sum_rule_invariant_form(**sources, **new) - trace_sum_rule_closed_form(**sources, **old)) == 0
        checks[f"seed {seed}: product rule invariant form = form with tr S_2"] = sp.expand(
            product_rule_invariant_form(**sources, **new) - product_rule_closed_form(**sources, **old)) == 0
    return checks


def solver_test_series(sources, condensate, second_order_scalar_trace):
    """Fixed rational charge-conserving inputs of hat S_2 with m_q = chi_0/(2 B_0), c^A, c^P from sources
    and tr S_2 = second_order_scalar_trace (test values, not fits)."""
    rational = sp.Rational
    scalar_source = sp.Matrix([[rational(1, 5), 0, 0], [0, rational(-2, 3), rational(1, 4) + sp.I / 6],
                               [0, rational(1, 4) - sp.I / 6, rational(1, 7)]])
    vector_source = sp.Matrix([[rational(-1, 2), 0, 0], [0, rational(1, 3), rational(-2, 5) + sp.I / 3],
                               [0, rational(-2, 5) - sp.I / 3, rational(3, 4)]])
    vector_first = sp.Matrix([[rational(2, 7), 0, 0], [0, rational(1, 9), sp.I / 5], [0, -sp.I / 5, rational(-1, 6)]])
    axial_first = sp.diag(rational(1, 8), rational(-3, 5), rational(1, 2))
    second_scalar = sp.Matrix([[rational(1, 3), 0, 0], [0, rational(-1, 7), rational(1, 9) - sp.I / 4],
                               [0, rational(1, 9) + sp.I / 4, 0]])
    second_scalar[2, 2] = second_order_scalar_trace - second_scalar[0, 0] - second_scalar[1, 1]
    second_pseudoscalar = sp.diag(rational(1, 6), rational(2, 5), rational(-1, 3))
    return first_and_second_order_inputs(
        sources["axial_source"], sources["pseudoscalar_source"], sources["quark_spurion"] / (2 * condensate),
        scalar_source, vector_source, vector_first, axial_first, second_scalar, second_pseudoscalar)


def invariant_form_solver_comparison(sources, condensate, second_order_scalar_trace):
    """Exact Vieta sum and product of det(s K - M^2) (blocks with tr S_2) vs the invariant forms."""
    series = solver_test_series(sources, condensate, second_order_scalar_trace)
    dressed_trace = dressed_second_order_scalar_from_definition(**series).trace()
    kinetic, mass_squared = blocks_with_second_order_scalar(sources, condensate, second_order_scalar_trace)
    exact_sum, exact_product, _ = _exact_pole_sum_and_product(kinetic, mass_squared)
    new = {"condensate": condensate, "dressed_scalar_trace": dressed_trace}
    return {
        "checks": {
            "inputs: tr S_2 != 0, c^P != 0, tr hat S_2 != tr S_2": (
                series["second_order_scalar"].trace() != 0 and sources["pseudoscalar_source"] != sp.zeros(3)
                and _simplify(dressed_trace - second_order_scalar_trace) != 0),
            "Vieta sum = trace rule invariant form (exact)": _simplify(
                exact_sum - trace_sum_rule_invariant_form(**sources, **new)) == 0,
            "Vieta product = product rule invariant form (exact)": _simplify(
                exact_product - product_rule_invariant_form(**sources, **new)) == 0,
        },
        "dressed_scalar_trace": _simplify(dressed_trace),
    }


def constant_shift_invariance_checks():
    """Claim 5: c^V, c^S -> +delta; c^A, c^P -> +kappa; c_G -> c_G - tr kappa. The second-order inputs are
    transformed with constant_rotation (transformation laws), not assumed."""
    shift_sets = []
    for seed in RANDOM_SEEDS:
        generator = random.Random(seed + 1000)
        shift_sets.append((f"seed {seed}", seed, random_rational_hermitian(generator), random_rational_hermitian(generator)))
    shift_sets.append(("symbolic charge-conserving delta, kappa (seed 11 inputs)", RANDOM_SEEDS[0],
                       _hermitian_charge_conserving("delta"), _hermitian_charge_conserving("kappa")))
    checks = {}
    for label, seed, vector_shift, axial_shift in shift_sets:
        generator = random.Random(seed)
        series, couplings = _random_series(generator)
        condensate, original_sources = _random_rule_sources(generator, couplings, series["quark_mass"])
        rotated = constant_rotation(series, vector_shift, axial_shift)
        shifted = {key: couplings[key] + (axial_shift if key in ("axial_source", "pseudoscalar_source") else vector_shift)
                   for key in couplings}
        shifted_sources = dict(original_sources, axial_source=shifted["axial_source"],
                               pseudoscalar_source=shifted["pseudoscalar_source"],
                               gluon_coupling=original_sources["gluon_coupling"] - sp.trace(axial_shift))
        quark_mass = series["quark_mass"]
        dressed = dressed_second_order_scalar_from_definition(**series)
        dressed_rotated = dressed_second_order_scalar_from_definition(**rotated)
        original_new = {"condensate": condensate, "dressed_scalar_trace": dressed.trace()}
        shifted_new = {"condensate": condensate, "dressed_scalar_trace": dressed_rotated.trace()}
        original_old = {"condensate": condensate, "second_order_scalar_trace": series["second_order_scalar"].trace()}
        shifted_old = {"condensate": condensate, "second_order_scalar_trace": rotated["second_order_scalar"].trace()}

        def unchanged(function, shifted_arguments, original_arguments):
            return _simplify(function(**shifted_sources, **shifted_arguments)
                             - function(**original_sources, **original_arguments)) == 0

        checks.update({
            f"{label}: rotated S_1 = i[c^S+delta, m_q], P_1 = {{c^P+kappa, m_q}}, C_V0 = c^V+delta, C_A0 = c^A+kappa": (
                _is_zero_matrix(rotated["first_order_scalar"] - sp.I * _commutator(shifted["scalar_source"], quark_mass))
                and _is_zero_matrix(rotated["first_order_pseudoscalar"]
                                    - _anticommutator(shifted["pseudoscalar_source"], quark_mass))
                and _is_zero_matrix(rotated["vector_couplings"][0] - shifted["vector_source"])
                and _is_zero_matrix(rotated["axial_couplings"][0] - shifted["axial_source"])),
            f"{label}: c^AP unchanged": _is_zero_matrix(
                (shifted["axial_source"] - shifted["pseudoscalar_source"])
                - (couplings["axial_source"] - couplings["pseudoscalar_source"])),
            f"{label}: c_GA = c_G + tr c^A unchanged": _simplify(
                shifted_sources["gluon_coupling"] + sp.trace(shifted["axial_source"])
                - original_sources["gluon_coupling"] - sp.trace(couplings["axial_source"])) == 0,
            f"{label}: raw tr S_2 changes": _simplify(
                rotated["second_order_scalar"].trace() - series["second_order_scalar"].trace()) != 0,
            f"{label}: hat S_2 (matrix) unchanged": _is_zero_matrix(dressed_rotated - dressed),
            f"{label}: tr hat S_2 unchanged": _simplify(dressed_rotated.trace() - dressed.trace()) == 0,
            f"{label}: trace rule invariant form unchanged": unchanged(
                trace_sum_rule_invariant_form, shifted_new, original_new),
            f"{label}: product rule invariant form unchanged": unchanged(
                product_rule_invariant_form, shifted_new, original_new),
            f"{label}: trace rule with transformed raw tr S_2 unchanged": unchanged(
                trace_sum_rule_closed_form, shifted_old, original_old),
            f"{label}: product rule with transformed raw tr S_2 unchanged": unchanged(
                product_rule_closed_form, shifted_old, original_old),
        })
    return checks


def derive():
    """Run all checks; return a dict of sections with their checks and results."""
    with_second_order = {"condensate": RATIONAL_CONDENSATE,
                         "second_order_scalar_trace": RATIONAL_SECOND_ORDER_SCALAR_TRACE}
    return {
        "charged mesons": charged_meson_checks(),
        "K0 pole": neutral_kaon_pole_checks(),
        "trace rule (diagonal chi_0)": {"checks": trace_sum_rule_checks(True)},
        "trace rule (chi_0 with d-s entry)": {"checks": trace_sum_rule_checks(False)},
        "product rule (diagonal chi_0)": {"checks": product_rule_checks(True)},
        "product rule (chi_0 with d-s entry)": {"checks": product_rule_checks(False)},
        "solver (diagonal chi_0)": solver_comparison(rational_sources()),
        "solver (chi_0 with d-s entry)": solver_comparison(
            rational_sources(sp.Rational(1, 50) + sp.I * sp.Rational(1, 70))),
        "solver (diagonal chi_0, trS_2 != 0)": solver_comparison(rational_sources(), **with_second_order),
        "solver (chi_0 with d-s entry, trS_2 != 0)": solver_comparison(
            rational_sources(sp.Rational(1, 50) + sp.I * sp.Rational(1, 70)), **with_second_order),
        "solver charged (10x10)": charged_solver_comparison(rational_sources()),
        "second-order spurion (quadratic_blocks)": second_order_spurion_taylor_checks(),
        "claim 1: tr hat S_2": {"checks": dressed_scalar_trace_checks()},
        "claims 2-4: invariant forms (symbolic, diagonal chi_0)": {"checks": invariant_form_symbolic_checks(True)},
        "claims 2-4: invariant forms (symbolic, chi_0 with d-s entry)": {"checks": invariant_form_symbolic_checks(False)},
        "claims 3-4: invariant forms (random rational Hermitian)": {"checks": invariant_form_random_checks()},
        "claim 4: invariant forms vs Vieta (diagonal chi_0)": invariant_form_solver_comparison(
            rational_sources(), **with_second_order),
        "claim 4: invariant forms vs Vieta (chi_0 with d-s entry)": invariant_form_solver_comparison(
            rational_sources(sp.Rational(1, 50) + sp.I * sp.Rational(1, 70)), **with_second_order),
        "claim 5: constant shift invariance": {"checks": constant_shift_invariance_checks()},
    }


def main():
    result = derive()
    print("Neutral basis (P_3, P_8, P_0, P_6, P_7, a); r = f_pi/f_a; N_f = 3.")
    print(f"m_pi+^2 = {result['charged mesons']['charged_pion_mass_squared']},"
          f"  m_K+^2 = {result['charged mesons']['charged_kaon_mass_squared']}")
    print("K0 null vector at s = m_K0,0^2 (symbolic sources):")
    sp.pprint(result["K0 pole"]["null_vector"].T)
    print("Trace sum rule: sum m_n^2 = tr(K^-1 M^2) =")
    print("  chi_uu + 2 chi_dd + 2 chi_ss + m0^2 + ma0^2 + r^2 [2 tr(cA cA chi0) - 2 tr(cA {cP, chi0})"
          " + (2 m0^2/N_f)(c_G + tr cA)^2 - B0 trS2]")
    print("Product rule: prod m_n^2 = det M^2 (det K = 1) =")
    print("  m_K0,0^4 { [det chi0 + (m0^2/N_f) e_2(chi0)] (ma0^2 - 2 r^2 tr(cP cP chi0) - r^2 B0 trS2)"
          " + 2 r^2 (m0^2/N_f) det chi0 (c_G + tr cP)^2 }")
    print("  i.e. tr S_2 contributes -r^2 B0 trS2 det M^2_PP, det M^2_PP = m_K0,0^4 [det chi0 + (m0^2/N_f) e_2(chi0)]")
    print("Invariant forms (c^AP = c^A - c^P, c_GA = c_G + tr c^A, hat S_2 dressed second-order source):")
    print("  tr hatS2 = trS2 + 4 tr(cA P1) - 4 tr(cA cA m_q),  P1 = {cP, m_q}")
    print("  r^2 B0 trS2 = r^2 [B0 tr hatS2 - 2 tr(cA {cP, chi0}) + 2 tr(cA cA chi0)]")
    print("  sum m_n^2 = chi_uu + 2 chi_dd + 2 chi_ss + m0^2 + ma0^2 + r^2 [(2 m0^2/N_f) c_GA^2 - B0 tr hatS2]")
    print("  det M^2 = m_K0,0^4 { [det chi0 + (m0^2/N_f) e_2(chi0)] (ma0^2 - r^2 [2 tr(cAP cAP chi0) + B0 tr hatS2])"
          " + 2 r^2 (m0^2/N_f) det chi0 (c_GA - tr cAP)^2 }")
    for section in ("claim 4: invariant forms vs Vieta (diagonal chi_0)",
                    "claim 4: invariant forms vs Vieta (chi_0 with d-s entry)"):
        print(f"{section}: tr hatS2 = {result[section]['dressed_scalar_trace']} GeV")
    print(f"Rational test values for the second-order spurion: B0 = {RATIONAL_CONDENSATE} GeV,"
          f" trS2 = {RATIONAL_SECOND_ORDER_SCALAR_TRACE} GeV")
    for section in ("solver (diagonal chi_0)", "solver (chi_0 with d-s entry)",
                    "solver (diagonal chi_0, trS_2 != 0)", "solver (chi_0 with d-s entry, trS_2 != 0)"):
        print(f"{section}: exact sum = {result[section]['exact_sum']}, exact product = {result[section]['exact_product']}")
    all_passed = True
    for section, content in result.items():
        for label, passed in content["checks"].items():
            all_passed &= bool(passed)
            print(f"[{'PASS' if passed else 'FAIL'}] {section}: {label}")
    print("ALL CHECKS PASS" if all_passed else "SOME CHECKS FAIL")


if __name__ == "__main__":
    main()
