#!/usr/bin/env python3
"""Analytical lower bound on the bare ALP mass from M^2 >= 0 at LO.

Run:   .venv/bin/python scripts/derive_alp_lo_stability_bound.py
Tests: .venv/bin/python -m unittest discover -s scripts -p test_derive_alp_lo_stability_bound.py

Conventions (Alda:2026nrd; memory/alp/shared/finite-alp-meson-mixing.md and
scripts/derive_alp_lo_neutral_blocks.py): N_f = 3, Hermitian U(3) generators
with tr(lambda_i lambda_j) = 2 delta_ij, t_i = tr(lambda_i), chi_0 = 2 B_0 m_q
diagonal in the quark basis (chi_u, chi_d, chi_s > 0), epsilon = f_pi/f_a,
L_2 = 1/2 dPhi^T K dPhi - 1/2 Phi^T M^2 Phi on Phi = (P_i, a), with

  M^2_ij = tr({lambda_i, lambda_j} chi_0)/4 + m_0^2 t_i t_j/(2 N_f),
  M^2_ia = -(epsilon/2) tr(lambda_i {c_q^P, chi_0}) + epsilon m_0^2 c_G t_i/N_f,
  M^2_aa = m_a0^2 + 2 epsilon^2 m_0^2 c_G^2/N_f.

The quark spurion is linear in a (finite-mixing record); no a^2 term from the
quark spurion is included. Units: chi_q, m_0^2, m_a0^2, M^2 in GeV^2; f_pi, f_a
in GeV; couplings dimensionless.

Derivation (each step is computed below, not assumed):
 1. M^2 depends on c_q^P, c_G, chi_0, m_0^2, epsilon only, not on c_q^A.
 2. The topological part of M^2 is (m_0^2/(2 N_f)) u u^T, u = (t_i, 2 epsilon c_G),
    a rank-one positive semidefinite matrix.
 3. With M^2_PP > 0 (chi_q > 0, m_0^2 >= 0), M^2 >= 0 <=> Schur complement
    M^2_aa - b^T (M^2_PP)^{-1} b >= 0, b_i = M^2_ia, i.e.
        m_a0^2 >= b^T (M^2_PP)^{-1} b - 2 epsilon^2 m_0^2 c_G^2/N_f.
 4. Flavour-diagonal sources: the neutral diagonal mesons P_3, P_8, P_0 are
    traded for quark-basis coordinates phi_q (diag(phi_u, phi_d, phi_s) =
    sum P_i lambda_i). Then
        P^T M^2_PP P = phi^T A phi, A = diag(chi_q)/2 + (m_0^2/(2 N_f)) 1 1^T,
        2 a b.P = 2 a beta.phi, beta_q = epsilon (-c_q chi_q + m_0^2 c_G/N_f),
    and A^{-1} follows from the Sherman-Morrison formula.
 5. Result (exact):
        m_a0^2 >= epsilon^2 B,
        B = 2 sum_q (c_q^P)^2 chi_q
            - 2 (m_0^2/N_f) (tr c_q^P + c_G)^2 / (1 + (m_0^2/N_f) sum_q 1/chi_q).
    For generic Hermitian c_q^P (d-s entries and also charged entries) the
    same form holds with sum_q (c_q^P)^2 chi_q -> tr(c_q^P chi_0 c_q^P).
 6. Sufficient bound: dropping the topological square gives
    B <= B_no_top = 2 tr(c_q^P chi_0 c_q^P).
 7. Limits: m_0^2 -> infinity; c_G = 0; vanishing mass mixing b = 0; B <= 0.

Verification: symbolic identities; independent constrained minimisation for
m_0^2 -> infinity; exact Schur complement of the generator-basis blocks; the
illustrative values of scripts/diagnose_alp_meson_spectrum.py; eigenvalues of
M^2 just above and just below the bound at random points.

Section 4 of the output evaluates B for explicit inputs (stated assumptions)
only as an order-of-magnitude illustration.
"""

from __future__ import annotations

import numpy as np
import sympy as sp

from derive_alp_lo_neutral_blocks import (
    ALL_GENERATOR_INDICES,
    FLAVOR_COUNT,
    NEUTRAL_GENERATOR_INDICES,
    closed_form_blocks,
    illustrative_lo_inputs,
)
from derive_alp_charged_meson_decoupling import symbolic_hermitian
from symbolic_utils import u3_generators

DIAGONAL_GENERATOR_INDICES = (3, 8, 0)

chi_u, chi_d, chi_s = sp.symbols("chi_u chi_d chi_s", positive=True)
c_u, c_d, c_s = sp.symbols("c_u c_d c_s", real=True)
gluon_coupling = sp.Symbol("c_G", real=True)
anomaly_mass_squared = sp.Symbol("m_0_sq", nonnegative=True)
bare_alp_mass_squared = sp.Symbol("m_a0_sq", real=True)
pion_decay_constant, alp_decay_constant = sp.symbols("f_pi f_a", positive=True)
ratio_epsilon = sp.Symbol("epsilon", positive=True)

QUARK_SPURION = (chi_u, chi_d, chi_s)
PSEUDOSCALAR_DIAGONAL = (c_u, c_d, c_s)


# --------------------------------------------------------------- closed forms


def bound_coefficient(pseudoscalar_source, quark_spurion_diagonal, anomaly_mass, gluon):
    """B such that M^2 >= 0 <=> m_a0^2 >= (f_pi/f_a)^2 B (chi_0 diagonal, any Hermitian c_q^P).

    B = 2 tr(c chi c) - 2 (m_0^2/N_f)(tr c + c_G)^2 / (1 + (m_0^2/N_f) sum_q 1/chi_q).
    """
    chi_matrix = sp.diag(*quark_spurion_diagonal)
    singlet = anomaly_mass / FLAVOR_COUNT
    inverse_sum = sum(1 / value for value in quark_spurion_diagonal)
    return (2 * (pseudoscalar_source * chi_matrix * pseudoscalar_source).trace()
            - 2 * singlet * (pseudoscalar_source.trace() + gluon)**2 / (1 + singlet * inverse_sum))


def bound_coefficient_without_topological_term(pseudoscalar_source, quark_spurion_diagonal):
    """B_no_top = 2 tr(c_q^P chi_0 c_q^P): the bound with the topological square dropped."""
    return 2 * (pseudoscalar_source * sp.diag(*quark_spurion_diagonal) * pseudoscalar_source).trace()


def bound_coefficient_large_singlet_mass(pseudoscalar_source, quark_spurion_diagonal, gluon):
    """m_0^2 -> infinity limit of B: 2 tr(c chi c) - 2 (tr c + c_G)^2 / sum_q 1/chi_q."""
    inverse_sum = sum(1 / value for value in quark_spurion_diagonal)
    return (bound_coefficient_without_topological_term(pseudoscalar_source, quark_spurion_diagonal)
            - 2 * (pseudoscalar_source.trace() + gluon)**2 / inverse_sum)


# --------------------------------------------------------- generator-basis Schur


def lo_inputs(axial_source, pseudoscalar_source, quark_spurion_diagonal, gluon, f_pi, f_a, anomaly_mass,
              bare_mass):
    return {
        "axial_source": axial_source, "pseudoscalar_source": pseudoscalar_source,
        "quark_spurion": sp.diag(*quark_spurion_diagonal), "gluon_coupling": gluon,
        "pion_decay_constant": f_pi, "alp_decay_constant": f_a,
        "anomaly_mass_squared": anomaly_mass, "bare_alp_mass_squared": bare_mass,
    }


def schur_bound_from_blocks(inputs, generator_indices=NEUTRAL_GENERATOR_INDICES):
    """b^T (M^2_PP)^{-1} b - (M^2_aa - m_a0^2), computed directly from the LO blocks."""
    _, mass_squared = closed_form_blocks(**inputs, generator_indices=generator_indices)
    meson = mass_squared[:-1, :-1]
    mixing = mass_squared[:-1, -1]
    solution = meson.LUsolve(mixing)
    induced = mass_squared[-1, -1] - inputs["bare_alp_mass_squared"]
    return (mixing.T * solution)[0] - induced


# ------------------------------------------------------------------ step 1, 2


def axial_source_independence():
    """dM^2/d(c_q^A entries) = 0 for symbolic charge-conserving Hermitian c_q^A, c_q^P."""
    axial = symbolic_hermitian("cA")
    inputs = lo_inputs(axial, symbolic_hermitian("cP"), QUARK_SPURION, gluon_coupling,
                       pion_decay_constant, alp_decay_constant, anomaly_mass_squared, bare_alp_mass_squared)
    _, mass_squared = closed_form_blocks(**inputs)
    return all(sp.simplify(mass_squared.diff(symbol)) == sp.zeros(*mass_squared.shape)
               for symbol in axial.free_symbols)


def topological_part_is_rank_one_square():
    """M^2(m_0^2) - M^2(m_0^2 = 0) = (m_0^2/(2 N_f)) u u^T with u = (t_i, 2 epsilon c_G)."""
    inputs = lo_inputs(sp.zeros(3), symbolic_hermitian("cP"), QUARK_SPURION, gluon_coupling,
                       pion_decay_constant, alp_decay_constant, anomaly_mass_squared, bare_alp_mass_squared)
    _, with_topological = closed_form_blocks(**inputs)
    _, without_topological = closed_form_blocks(**dict(inputs, anomaly_mass_squared=0))
    generators = u3_generators()
    direction = sp.Matrix([generators[index].trace() for index in NEUTRAL_GENERATOR_INDICES]
                          + [2 * pion_decay_constant / alp_decay_constant * gluon_coupling])
    expected = anomaly_mass_squared / (2 * FLAVOR_COUNT) * direction * direction.T
    return sp.simplify(with_topological - without_topological - expected) == sp.zeros(6, 6)


# ------------------------------------------------------- step 3, 4, 5 (diagonal)


def quark_basis_map():
    """Matrix T with P_i = T_iq phi_q for i in (3, 8, 0): P_i = tr(lambda_i diag(phi))/2."""
    generators = u3_generators()
    return sp.Matrix([[generators[index][flavour, flavour] / 2 for flavour in range(3)]
                      for index in DIAGONAL_GENERATOR_INDICES])


def flavour_diagonal_derivation():
    """Steps 3-5 for c_q^P = diag(c_u, c_d, c_s); returns intermediate and final expressions."""
    inputs = lo_inputs(sp.zeros(3), sp.diag(*PSEUDOSCALAR_DIAGONAL), QUARK_SPURION, gluon_coupling,
                       pion_decay_constant, alp_decay_constant, anomaly_mass_squared, bare_alp_mass_squared)
    _, mass_squared = closed_form_blocks(**inputs, generator_indices=DIAGONAL_GENERATOR_INDICES)
    meson, mixing = mass_squared[:3, :3], mass_squared[:3, 3]
    induced = sp.simplify(mass_squared[3, 3] - bare_alp_mass_squared)
    epsilon = pion_decay_constant / alp_decay_constant
    singlet_half = anomaly_mass_squared / (2 * FLAVOR_COUNT)
    ones = sp.ones(3, 1)

    transform = quark_basis_map()
    quark_meson = sp.simplify(transform.T * meson * transform)
    quark_mixing = sp.simplify(transform.T * mixing)
    diagonal_part = sp.diag(*QUARK_SPURION) / 2
    expected_quark_meson = diagonal_part + singlet_half * ones * ones.T
    expected_quark_mixing = sp.Matrix([epsilon * (-coupling * chi + anomaly_mass_squared * gluon_coupling / FLAVOR_COUNT)
                                       for coupling, chi in zip(PSEUDOSCALAR_DIAGONAL, QUARK_SPURION)])

    diagonal_inverse = diagonal_part.inv()
    sherman_morrison_inverse = (diagonal_inverse - singlet_half * diagonal_inverse * ones * ones.T * diagonal_inverse
                                / (1 + singlet_half * (ones.T * diagonal_inverse * ones)[0]))

    schur_quark_basis = sp.factor(sp.simplify(
        (quark_mixing.T * sherman_morrison_inverse * quark_mixing)[0] - induced))
    schur_generator_basis = sp.simplify((mixing.T * meson.inv() * mixing)[0] - induced)
    closed_form = epsilon**2 * bound_coefficient(sp.diag(*PSEUDOSCALAR_DIAGONAL), QUARK_SPURION,
                                                 anomaly_mass_squared, gluon_coupling)
    scaled = sp.simplify(schur_generator_basis * alp_decay_constant**2 / pion_decay_constant**2)
    return {
        "quark_meson": quark_meson, "quark_mixing": quark_mixing, "induced": induced,
        "sherman_morrison_inverse": sherman_morrison_inverse, "schur_quark_basis": schur_quark_basis,
        "closed_form": closed_form, "coefficient": scaled,
        "checks": {
            "M^2_PP > 0 (leading minors positive for chi_q > 0, m_0^2 >= 0)": all(
                sp.factor(meson[:order, :order].det()).is_positive for order in (1, 2, 3)),
            "P^T M^2_PP P = phi^T A phi with A = diag(chi)/2 + (m_0^2/(2 N_f)) 1 1^T":
                sp.simplify(quark_meson - expected_quark_meson) == sp.zeros(3, 3),
            "b.P = beta.phi with beta_q = epsilon(-c_q chi_q + m_0^2 c_G/N_f)":
                sp.simplify(quark_mixing - expected_quark_mixing) == sp.zeros(3, 1),
            "M^2_aa - m_a0^2 = 2 epsilon^2 m_0^2 c_G^2/N_f":
                sp.simplify(induced - 2 * epsilon**2 * anomaly_mass_squared * gluon_coupling**2 / FLAVOR_COUNT) == 0,
            "Sherman-Morrison inverse: A A^{-1} = 1":
                sp.simplify(expected_quark_meson * sherman_morrison_inverse) == sp.eye(3),
            "Schur complement equal in quark and generator bases":
                sp.simplify(schur_quark_basis - schur_generator_basis) == 0,
            "bound = epsilon^2 B (closed form)": sp.simplify(schur_generator_basis - closed_form) == 0,
            "bound (f_a/f_pi)^2 independent of f_pi and f_a":
                not ({pion_decay_constant, alp_decay_constant} & scaled.free_symbols),
        },
    }


# --------------------------------------------------------- step 5, d-s and charged


def flavour_violating_derivation():
    """Generic neutral c_q^P (complex d-s entry), 5 neutral mesons, symbolic; returns checks."""
    pseudoscalar = symbolic_hermitian("cP")
    inputs = lo_inputs(sp.zeros(3), pseudoscalar, QUARK_SPURION, gluon_coupling, pion_decay_constant,
                       alp_decay_constant, anomaly_mass_squared, bare_alp_mass_squared)
    schur = schur_bound_from_blocks(inputs)
    epsilon = pion_decay_constant / alp_decay_constant
    closed_form = epsilon**2 * bound_coefficient(pseudoscalar, QUARK_SPURION, anomaly_mass_squared, gluon_coupling)
    diagonal_part = pseudoscalar.copy()
    diagonal_part[1, 2] = diagonal_part[2, 1] = 0
    off_diagonal_contribution = sp.simplify(
        bound_coefficient(pseudoscalar, QUARK_SPURION, anomaly_mass_squared, gluon_coupling)
        - bound_coefficient(diagonal_part, QUARK_SPURION, anomaly_mass_squared, gluon_coupling))
    absolute_square = sp.Symbol("cP_ds_re", real=True)**2 + sp.Symbol("cP_ds_im", real=True)**2
    return {"checks": {
        "generic neutral c_q^P: Schur bound = epsilon^2 B with tr(c chi c)": sp.simplify(schur - closed_form) == 0,
        "d-s entry adds 2 (chi_d + chi_s) |c_ds^P|^2 to B":
            sp.simplify(off_diagonal_contribution - 2 * (chi_d + chi_s) * absolute_square) == 0,
    }}


def random_hermitian(generator, charge_conserving):
    matrix = sp.zeros(3)
    for row in range(3):
        matrix[row, row] = sp.Rational(int(generator.integers(-9, 10)), int(generator.integers(1, 7)))
        for column in range(row + 1, 3):
            if charge_conserving and row == 0:
                continue
            entry = (sp.Rational(int(generator.integers(-9, 10)), int(generator.integers(1, 7)))
                     + sp.I * sp.Rational(int(generator.integers(-9, 10)), int(generator.integers(1, 7))))
            matrix[row, column], matrix[column, row] = entry, sp.conjugate(entry)
    return matrix


def random_exact_point(generator, charge_conserving=True):
    """Random exact rational inputs (GeV units) with chi_q > 0 and m_0^2 > 0."""
    spurion = tuple(sp.Rational(int(generator.integers(1, 60)), 100) for _ in range(3))
    return lo_inputs(random_hermitian(generator, True), random_hermitian(generator, charge_conserving), spurion,
                     sp.Rational(int(generator.integers(-6, 7)), 4), sp.Rational(23, 250),
                     sp.Rational(int(generator.integers(1, 40)), 40), sp.Rational(int(generator.integers(1, 20)), 10),
                     0)


def charged_entries_check(seed=11, sample_count=3):
    """All nine generators with charge-violating c_q^P at random exact points: Schur = epsilon^2 B."""
    generator = np.random.default_rng(seed)
    agreements = []
    for _ in range(sample_count):
        inputs = random_exact_point(generator, charge_conserving=False)
        epsilon = inputs["pion_decay_constant"] / inputs["alp_decay_constant"]
        schur = schur_bound_from_blocks(inputs, generator_indices=ALL_GENERATOR_INDICES)
        closed = epsilon**2 * bound_coefficient(inputs["pseudoscalar_source"], tuple(inputs["quark_spurion"].diagonal()),
                                                inputs["anomaly_mass_squared"], inputs["gluon_coupling"])
        agreements.append(sp.simplify(schur - closed) == 0)
    return all(agreements)


# ---------------------------------------------------------------- step 6, 7


def sufficient_bound_and_limits():
    """Step 6 and the limits of step 7 for diagonal c_q^P; returns expressions and checks."""
    pseudoscalar = sp.diag(*PSEUDOSCALAR_DIAGONAL)
    full = bound_coefficient(pseudoscalar, QUARK_SPURION, anomaly_mass_squared, gluon_coupling)
    without_topological = bound_coefficient_without_topological_term(pseudoscalar, QUARK_SPURION)
    large_singlet = bound_coefficient_large_singlet_mass(pseudoscalar, QUARK_SPURION, gluon_coupling)
    singlet = anomaly_mass_squared / FLAVOR_COUNT
    inverse_sum = sum(1 / value for value in QUARK_SPURION)
    trace = sum(PSEUDOSCALAR_DIAGONAL)
    gap = sp.simplify(without_topological - full)
    expected_gap = 2 * singlet * (trace + gluon_coupling)**2 / (1 + singlet * inverse_sum)

    # Independent route for m_0^2 -> infinity: minimise the quadratic form without the
    # topological square on the surface u = sum_q phi_q + 2 epsilon c_G a = 0.
    phi = sp.symbols("phi_u phi_d phi_s", real=True)
    alp, multiplier, epsilon = sp.symbols("a lambda_L epsilon", real=True)
    form = (sum(chi * value**2 for chi, value in zip(QUARK_SPURION, phi)) / 2
            - 2 * epsilon * alp * sum(c * chi * value for c, chi, value in zip(PSEUDOSCALAR_DIAGONAL, QUARK_SPURION, phi)))
    constraint = sum(phi) + 2 * epsilon * gluon_coupling * alp
    lagrangian_function = form - multiplier * constraint
    solution = sp.solve([lagrangian_function.diff(value) for value in phi] + [constraint], list(phi) + [multiplier],
                        dict=True)[0]
    minimum = sp.simplify(form.subs(solution))
    constrained_bound = sp.simplify(-minimum / alp**2)

    # Lagrange identity: (sum 1/chi)(sum c^2 chi) - (sum c)^2 = sum_{q<q'} (c_q chi_q - c_q' chi_q')^2/(chi_q chi_q')
    pairs = ((0, 1), (0, 2), (1, 2))
    lagrange = sum((PSEUDOSCALAR_DIAGONAL[i] * QUARK_SPURION[i] - PSEUDOSCALAR_DIAGONAL[j] * QUARK_SPURION[j])**2
                   / (QUARK_SPURION[i] * QUARK_SPURION[j]) for i, j in pairs)
    numerator_without_gluon = sp.simplify(full.subs(gluon_coupling, 0) * (1 + singlet * inverse_sum) / 2)
    aligned = {coupling: singlet * gluon_coupling / chi for coupling, chi in zip(PSEUDOSCALAR_DIAGONAL, QUARK_SPURION)}
    proportional = {coupling: sp.Symbol("alpha", real=True) / chi
                    for coupling, chi in zip(PSEUDOSCALAR_DIAGONAL, QUARK_SPURION)}
    gradient = [sp.simplify(full.diff(coupling)) for coupling in PSEUDOSCALAR_DIAGONAL]
    return {
        "full": full, "without_topological": without_topological, "large_singlet": large_singlet,
        "checks": {
            "B_no_top equals the Schur bound at m_0^2 = 0": sp.simplify(full.subs(anomaly_mass_squared, 0)
                                                                       - without_topological) == 0,
            "B_no_top - B = 2 (m_0^2/N_f)(tr c + c_G)^2/(1 + (m_0^2/N_f) sum 1/chi) >= 0":
                sp.simplify(gap - expected_gap) == 0 and expected_gap.is_nonnegative,
            "(i) lim_{m_0^2 -> oo} B = 2 sum c^2 chi - 2 (tr c + c_G)^2/sum 1/chi":
                sp.simplify(sp.limit(full, anomaly_mass_squared, sp.oo) - large_singlet) == 0,
            "(i) independent check: minimum on the surface u = 0 gives epsilon^2 B_oo":
                sp.simplify(constrained_bound - epsilon**2 * large_singlet) == 0,
            "(i) B decreases monotonically in m_0^2 (dB/dm_0^2 <= 0)":
                sp.simplify(full.diff(anomaly_mass_squared)).is_nonpositive,
            "(i) Lagrange identity for B_oo at c_G = 0":
                sp.simplify(large_singlet.subs(gluon_coupling, 0) * inverse_sum / 2 - lagrange) == 0,
            "(i) c_G = 0, c_q = alpha/chi_q: B_oo = 0": sp.simplify(
                large_singlet.subs(gluon_coupling, 0).subs(proportional)) == 0,
            "(ii) c_G = 0: (1 + m_0^2/N_f sum 1/chi) B/2 = sum c^2 chi + (m_0^2/N_f) * Lagrange sum":
                sp.simplify(numerator_without_gluon
                            - sum(c**2 * chi for c, chi in zip(PSEUDOSCALAR_DIAGONAL, QUARK_SPURION))
                            - singlet * lagrange) == 0,
            "(ii) stationary point of B in c_q^P is c_q = m_0^2 c_G/(N_f chi_q)":
                all(sp.simplify(entry.subs(aligned)) == 0 for entry in gradient),
            "(ii) there the mass mixing b vanishes and B = -2 m_0^2 c_G^2/N_f":
                sp.simplify(full.subs(aligned) + 2 * singlet * gluon_coupling**2) == 0,
            "(ii) c_q^P = 0: B = -2 (m_0^2/N_f) c_G^2/(1 + (m_0^2/N_f) sum 1/chi)":
                sp.simplify(full.subs({c: 0 for c in PSEUDOSCALAR_DIAGONAL})
                            + 2 * singlet * gluon_coupling**2 / (1 + singlet * inverse_sum)) == 0,
        },
    }


# ------------------------------------------------------ numerical eigenvalue tests


def smallest_mass_eigenvalue(inputs, generator_indices=NEUTRAL_GENERATOR_INDICES):
    _, mass_squared = closed_form_blocks(**inputs, generator_indices=generator_indices)
    return float(np.linalg.eigvalsh(np.array(mass_squared.evalf(30), dtype=float))[0])


def eigenvalue_sign_checks(seed=3, sample_count=6, relative_offset=1e-6):
    """At random points: smallest eigenvalue of M^2 >0 just above and <0 just below epsilon^2 B."""
    generator = np.random.default_rng(seed)
    results = []
    for index in range(sample_count):
        charge_conserving = index % 2 == 0
        inputs = random_exact_point(generator, charge_conserving=True)
        if not charge_conserving:
            inputs["pseudoscalar_source"] = random_hermitian(generator, False)
        indices = NEUTRAL_GENERATOR_INDICES if charge_conserving else ALL_GENERATOR_INDICES
        epsilon = inputs["pion_decay_constant"] / inputs["alp_decay_constant"]
        threshold = epsilon**2 * bound_coefficient(inputs["pseudoscalar_source"], tuple(inputs["quark_spurion"].diagonal()),
                                                   inputs["anomaly_mass_squared"], inputs["gluon_coupling"])
        offset = relative_offset * (abs(threshold) + 1)
        above = smallest_mass_eigenvalue(dict(inputs, bare_alp_mass_squared=threshold + offset), indices)
        below = smallest_mass_eigenvalue(dict(inputs, bare_alp_mass_squared=threshold - offset), indices)
        results.append((float(threshold), above, below))
    return results


def light_eigenvalue_check(ratio_values=(sp.Rational(1, 10), sp.Rational(1, 100), sp.Rational(1, 1000))):
    """c_q^A = c_q^P = 0, c_G = 1/2, m_a0^2 = 0: smallest eigenvalue of M^2 (= lightest pole, K = 1)
    divided by -epsilon^2 B tends to 1 as epsilon -> 0 (illustrative rational chi_q, m_0^2)."""
    ratios = []
    for epsilon in ratio_values:
        inputs = lo_inputs(sp.zeros(3), sp.zeros(3), (sp.Rational(3, 200), sp.Rational(1, 40), sp.Rational(47, 100)),
                           sp.Rational(1, 2), epsilon, 1, sp.Rational(4, 5), 0)
        coefficient = bound_coefficient(sp.zeros(3), tuple(inputs["quark_spurion"].diagonal()),
                                        inputs["anomaly_mass_squared"], inputs["gluon_coupling"])
        ratios.append(smallest_mass_eigenvalue(inputs) / float(-epsilon**2 * coefficient))
    return ratios


def completed_square_remark():
    """Remark on the source convention (diagonal c_q^P): adding 2 epsilon^2 tr(c chi c) a^2/2 to the
    potential, i.e. epsilon^2 B_no_top to M^2_aa, turns the chi_0 part into the square
    (1/2) sum_q chi_q (phi_q - 2 epsilon c_q a)^2 (in the normalisation 2 V_2 used here), so the bound
    becomes epsilon^2 (B - B_no_top) <= 0. Such an a^2 term is NOT part of the linear spurion convention."""
    phi = sp.symbols("phi_u phi_d phi_s", real=True)
    alp, epsilon = sp.symbols("a epsilon", real=True)
    linear_spurion = (sum(chi * value**2 for chi, value in zip(QUARK_SPURION, phi)) / 2
                      - 2 * epsilon * alp * sum(c * chi * value for c, chi, value
                                                in zip(PSEUDOSCALAR_DIAGONAL, QUARK_SPURION, phi)))
    added = epsilon**2 * alp**2 * bound_coefficient_without_topological_term(sp.diag(*PSEUDOSCALAR_DIAGONAL),
                                                                              QUARK_SPURION)
    square = sum(chi * (value - 2 * epsilon * c * alp)**2
                 for chi, value, c in zip(QUARK_SPURION, phi, PSEUDOSCALAR_DIAGONAL)) / 2
    return sp.expand(linear_spurion + added - square) == 0


# ------------------------------------------------------- illustrative evaluation

# Assumed inputs (order-of-magnitude illustration only; not fits):
# masses in GeV from ParticleDataGroup:2024cfk (rounded), f_pi as requested by the plan.
ILLUSTRATION_INPUTS = {
    "f_pi [GeV]": 0.0924,
    "m_pi0 [GeV]": 0.1349768,
    "m_K (isospin average) [GeV]": 0.4956,
    "m_eta [GeV]": 0.547862,
    "m_eta' [GeV]": 0.95778,
    "m_u/m_d": 0.474,
    "m_s/m_ud, m_ud = (m_u + m_d)/2": 27.33,
}


def illustration_spurion(values=ILLUSTRATION_INPUTS):
    """chi_q = 2 B_0 m_q from 2 B_0 m_ud = m_pi0^2 (LO, isospin-averaged) and the quark-mass ratios;
    m_0^2 from the LO relation m_eta^2 + m_eta'^2 = 2 m_K^2 + m_0^2 (this convention, N_f = 3)."""
    light = values["m_pi0 [GeV]"]**2
    ratio = values["m_u/m_d"]
    chi = (2 * light * ratio / (1 + ratio), 2 * light / (1 + ratio), values["m_s/m_ud, m_ud = (m_u + m_d)/2"] * light)
    singlet = values["m_eta [GeV]"]**2 + values["m_eta' [GeV]"]**2 - 2 * values["m_K (isospin average) [GeV]"]**2
    return chi, singlet


ILLUSTRATION_BENCHMARKS = {
    "c^P = diag(1, 1, 1), c_G = 0": (sp.diag(1, 1, 1), 0),
    "c^P = diag(1, 1, 1), c_G = 1": (sp.diag(1, 1, 1), 1),
    "c^P = diag(1, -1/2, -1/2), c_G = 1": (sp.diag(1, sp.Rational(-1, 2), sp.Rational(-1, 2)), 1),
    "c^P = 0, c_G = 1": (sp.zeros(3), 1),
}
ILLUSTRATION_DECAY_CONSTANTS = (1.0, 10.0, 100.0)


def illustration():
    chi, singlet = illustration_spurion()
    f_pi = ILLUSTRATION_INPUTS["f_pi [GeV]"]
    pion_mass_squared = ILLUSTRATION_INPUTS["m_pi0 [GeV]"]**2
    rows = []
    for label, (pseudoscalar, gluon) in ILLUSTRATION_BENCHMARKS.items():
        full = float(bound_coefficient(pseudoscalar, chi, singlet, gluon))
        without_topological = float(bound_coefficient_without_topological_term(pseudoscalar, chi))
        crossing = np.sqrt(full / pion_mass_squared) if full > 0 else None
        bounds = [(f_a, (f_pi / f_a)**2 * full) for f_a in ILLUSTRATION_DECAY_CONSTANTS]
        rows.append({"label": label, "B": full, "B_no_top": without_topological,
                     "f_a/f_pi at bound = m_pi^2": crossing, "bounds": bounds})
    return chi, singlet, rows


# --------------------------------------------------------------------- run


def illustrative_values_check():
    """Reproduce the exact bounds of scripts/diagnose_alp_meson_spectrum.py from epsilon^2 B."""
    expected = {"flavor_diagonal": sp.Rational(1613899121, 18280800000),
                "neutral": sp.Rational(14983933993, 91404000000)}
    checks = {}
    for structure, value in expected.items():
        inputs = illustrative_lo_inputs(structure)
        epsilon = inputs["pion_decay_constant"] / inputs["alp_decay_constant"]
        closed = epsilon**2 * bound_coefficient(inputs["pseudoscalar_source"], tuple(inputs["quark_spurion"].diagonal()),
                                                inputs["anomaly_mass_squared"], inputs["gluon_coupling"])
        checks[f"epsilon^2 B reproduces the diagnose-script value ({structure}) = {value}"] = (
            sp.simplify(closed - value) == 0)
    return checks


def main():
    print("Step 1-2: structure of M^2")
    print(f"  [{'pass' if axial_source_independence() else 'FAIL'}] M^2 independent of c_q^A")
    print(f"  [{'pass' if topological_part_is_rank_one_square() else 'FAIL'}] "
          "topological part = (m_0^2/(2 N_f)) u u^T, u = (t_i, 2 epsilon c_G)")
    diagonal = flavour_diagonal_derivation()
    print("Step 3-5: flavour-diagonal sources, quark basis")
    print(f"  A^-1 (Sherman-Morrison) = {diagonal['sherman_morrison_inverse']}")
    print("  m_a0^2 >= (f_pi/f_a)^2 B,  B = 2 sum_q (c_q^P)^2 chi_q"
          " - 2 (m_0^2/N_f)(tr c^P + c_G)^2 / (1 + (m_0^2/N_f) sum_q 1/chi_q)")
    print(f"  sympy form of B: {sp.factor(diagonal['coefficient'])}")
    for label, passed in diagonal["checks"].items():
        print(f"  [{'pass' if passed else 'FAIL'}] {label}")
    print("Step 5: flavour-violating entries")
    for label, passed in flavour_violating_derivation()["checks"].items():
        print(f"  [{'pass' if passed else 'FAIL'}] {label}")
    print(f"  [{'pass' if charged_entries_check() else 'FAIL'}] charge-violating c_q^P, 9 generators, "
          "random exact points: Schur = epsilon^2 B")
    print("Step 6-7: sufficient bound and limits")
    for label, passed in sufficient_bound_and_limits()["checks"].items():
        print(f"  [{'pass' if passed else 'FAIL'}] {label}")
    for label, passed in illustrative_values_check().items():
        print(f"  [{'pass' if passed else 'FAIL'}] {label}")
    print("Eigenvalue signs of M^2 at m_a0^2 = epsilon^2 B (1 +- 1e-6), random points:")
    for threshold, above, below in eigenvalue_sign_checks():
        print(f"  [{'pass' if above > 0 > below else 'FAIL'}] epsilon^2 B = {threshold:.6g} GeV^2: "
              f"smallest eigenvalue above {above:.3g}, below {below:.3g}")
    ratios = light_eigenvalue_check()
    print(f"  c^A = c^P = 0, m_a0^2 = 0: lightest pole / (-epsilon^2 B) at epsilon = 0.1, 0.01, 0.001: "
          f"{', '.join(f'{value:.8f}' for value in ratios)}")
    print(f"  [{'pass' if completed_square_remark() else 'FAIL'}] remark: an extra a^2 term epsilon^2 B_no_top "
          "(not in the linear spurion convention) would complete the chi_0 square and give B - B_no_top <= 0")

    chi, singlet, rows = illustration()
    print("Illustration (order of magnitude only; assumed inputs):")
    for name, value in ILLUSTRATION_INPUTS.items():
        print(f"  {name} = {value}")
    print(f"  derived: chi_u, chi_d, chi_s = {chi[0]:.5g}, {chi[1]:.5g}, {chi[2]:.5g} GeV^2; "
          f"m_0^2 = {singlet:.4g} GeV^2")
    for row in rows:
        print(f"  {row['label']}: B = {row['B']:.4g} GeV^2, B_no_top = {row['B_no_top']:.4g} GeV^2, "
              f"f_a/f_pi at epsilon^2 B = m_pi0^2: {row['f_a/f_pi at bound = m_pi^2']}")
        for f_a, bound in row["bounds"]:
            root = f"{1e3 * np.sqrt(bound):.3g} MeV" if bound >= 0 else "none (bound negative)"
            print(f"    f_a = {f_a:g} GeV: m_a0^2 >= {bound:.3g} GeV^2, sqrt = {root}")


if __name__ == "__main__":
    main()
