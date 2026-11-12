#!/usr/bin/env python3
"""LO neutral poles as eigenvalues of the canonical mass matrix: arrow form, interlacing, and closed-form
approximations for an isolated meson pole and for an isolated ALP pole (no expansion in epsilon = f_pi/f_a).

Run:   .venv/bin/python scripts/derive_alp_lo_arrow_matrix_poles.py
Tests: .venv/bin/python -m unittest discover -s scripts -p test_derive_alp_lo_arrow_matrix_poles.py
Records: memory/alp/shared/alp-lo-arrow-matrix-spectrum.md (items 1, 2)
         memory/alp/shared/alp-lo-isolated-meson-pole.md (item 3)
         memory/alp/shared/alp-lo-isolated-alp-pole.md (item 4)

Conventions (Alda:2026nrd; memory/alp/shared/alp-lo-exact-pole-equations.md): N_f = 3, neutral basis
(P_3, P_8, P_0, P_6, P_7, a), canonical frame K = 1, M_can^2 from `canonical_mass_matrix_from_sources`.
lambda_n: generator of the meson mass eigenstate n, 𝕄_0(lambda_n) = m_n0^2 lambda_n, tr(lambda_n lambda_m) = 2 delta_nm;
𝒞 = -{chi_0, c_q^AP}/2 - (m_0^2/N_f) c_GA 1; M = M_can,aa (record: M^2_{\\hat a\\hat a,can}).
In the code, tr_lambda_C is the list of tr(lambda_n 𝒞) over the mass eigenstates n (GeV^2).

Derived and checked here:
 1. Rotating the meson block of M_can^2 to the mass basis gives the arrow matrix
        A = [[diag(m_n0^2), column -epsilon tr(lambda_n 𝒞)], [row -epsilon tr(lambda_n 𝒞), M]]
    (only the diagonal, the last row and the last column are nonzero); its eigenvalues are the eigenvalues of
    M_can^2 and the generalized eigenvalues of (K, M^2) of the bare frame, and
        det(s - A) = (s - M) prod_n (s - m_n0^2) - epsilon^2 sum_n tr(lambda_n 𝒞)^2 prod_{n' != n} (s - m_n'0^2)
    (50 digits, LO blocks at epsilon = 1/2, 1, 2; the identity also at s = m_n0^2).
 2. Cauchy interlacing, unmixed masses sorted increasingly m_1 <= ... <= m_5, poles mu_1 <= ... <= mu_6:
        mu_k <= m_k <= mu_{k+1}   (k = 1..5),
    and the opposite direction m_k <= mu_k fails (LO blocks at several epsilon; arrow models).
 3. Isolated meson n: exact (s - m_n0^2)(s - Mt_n(s)) = epsilon^2 tr(lambda_n 𝒞)^2, Mt_n(s) = M + epsilon^2 sum_{n' != n}
    tr(lambda_n' 𝒞)^2/(s - m_n'0^2). Freezing Mt_n at s = m_n0^2 gives
        s_pm = (m_n0^2 + Mt_n0)/2 pm sqrt((Mt_n0 - m_n0^2)^2/4 + epsilon^2 tr(lambda_n 𝒞)^2).
    Exact bound: every pole mu lies within |Mt_n(mu) - Mt_n0| of s_+ or s_-. Error scan against the
    exact eigenvalues (arrow model, several spacings Delta and epsilon up to 2): error ~ Delta^-2 at fixed
    epsilon. Applied to the LO blocks for each meson.
 4. Isolated ALP: exact s = g(s), g(s) = M + epsilon^2 sum_n tr(lambda_n 𝒞)^2/(s - m_n0^2) = M - 2 epsilon^2 tr(𝒞 𝕄_s^{-1}(𝒞)).
    Frozen value s^(1) = g(M). For M above (below) all m_n0^2 and q = epsilon^2 sum_n tr(lambda_n 𝒞)^2/(M - m_n0^2)^2 < 1:
    the iteration converges, the pole lies between M and s^(1), and |s^(1) - s_a| <= q/(1+q) |s^(1) - M|.
    For |g'(s_a)| > 1 the iteration does not converge to s_a. Error scan against the exact eigenvalue.

Units: all masses squared, tr(lambda_n 𝒞) and M in GeV^2; epsilon dimensionless. Model inputs are test values, not fits.
"""

from __future__ import annotations

import mpmath
import numpy as np
import sympy as sp

from derive_alp_lo_neutral_blocks import illustrative_lo_inputs
from derive_alp_meson_mass_shifts_exact import (
    TEST_CONDENSATE,
    TEST_SECOND_ORDER_SCALAR_TRACE,
    bare_blocks,
    canonical_mass_matrix_from_sources,
)
from derive_alp_meson_mass_shifts_second_order import (
    HIGH_PRECISION_DIGITS,
    HIGH_PRECISION_TOLERANCE,
    MESON_COUNT,
    NEUTRAL_GENERATORS,
    _block_inputs,
    high_precision_trace,
    mass_eigenstate_generators,
    rational_test_sources,
    tadpole_source_from_sources,
    to_high_precision,
)

LO_EPSILON_VALUES = (sp.Rational(1, 2), sp.Integer(1), sp.Integer(2))
MODEL_EPSILON_VALUES = ("0.1", "0.5", "1", "2")
MODEL_SPACINGS = (2, 4, 8, 16, 32, 64)
DEGENERACY_TOLERANCE = mpmath.mpf("1e-30")


def _mp(value):
    return mpmath.mpf(str(sp.N(value, HIGH_PRECISION_DIGITS + 10)))


# ------------------------------------------------------------------ arrow matrix


def arrow_matrix(unmixed_masses_squared, tr_lambda_C, epsilon, alp_entry):
    """A = [[diag(m_n0^2), column -epsilon tr(lambda_n 𝒞)], [row -epsilon tr(lambda_n 𝒞), M]] (mpmath)."""
    size = len(unmixed_masses_squared) + 1
    matrix = mpmath.matrix(size, size)
    for index, (mass_squared, tr_lambda_n_C) in enumerate(zip(unmixed_masses_squared, tr_lambda_C)):
        matrix[index, index] = mass_squared
        matrix[index, size - 1] = matrix[size - 1, index] = -epsilon * tr_lambda_n_C
    matrix[size - 1, size - 1] = alp_entry
    return matrix


def sorted_eigenvalues(matrix):
    eigenvalues, _ = mpmath.eigsy(matrix)
    return sorted(eigenvalues[index] for index in range(matrix.rows))


def arrow_determinant(squared_momentum, unmixed_masses_squared, tr_lambda_C, epsilon, alp_entry):
    """(s - M) prod_n (s - m_n0^2) - epsilon^2 sum_n tr(lambda_n 𝒞)^2 prod_{n' != n} (s - m_n'0^2)."""
    factors = [squared_momentum - value for value in unmixed_masses_squared]
    return ((squared_momentum - alp_entry) * mpmath.fprod(factors)
            - epsilon**2 * mpmath.fsum(tr_lambda_n_C**2 * mpmath.fprod(factors[:index] + factors[index + 1:])
                                       for index, tr_lambda_n_C in enumerate(tr_lambda_C)))


def lo_arrow_inputs(sources, condensate=TEST_CONDENSATE, second_order_scalar_trace=TEST_SECOND_ORDER_SCALAR_TRACE,
                    digits=HIGH_PRECISION_DIGITS):
    """LO arrow-matrix inputs: unmixed m_n0^2, tr(lambda_n 𝒞), epsilon, M = M_can,aa, and M_can^2 (sympy).

    Degenerate eigenvalues (K^0 doublet) are rotated so that at most one generator of the eigenspace has
    tr(lambda_n 𝒞) != 0: that generator gets sqrt(sum of tr(lambda 𝒞)^2 over the eigenspace) and the others 0.
    """
    canonical = canonical_mass_matrix_from_sources(sources, condensate, second_order_scalar_trace)
    tadpole = to_high_precision(tadpole_source_from_sources(**_block_inputs(sources)), digits)
    eigenbasis = mass_eigenstate_generators(sources, digits)
    unmixed, tr_lambda_C = [], []
    for value, generator in eigenbasis:
        tr_lambda_n_C = high_precision_trace(generator * tadpole).real
        if unmixed and abs(value - unmixed[-1]) < DEGENERACY_TOLERANCE:
            tr_lambda_C[-1] = mpmath.sqrt(tr_lambda_C[-1]**2 + tr_lambda_n_C**2)
            unmixed.append(value)
            tr_lambda_C.append(mpmath.mpf(0))
        else:
            unmixed.append(value)
            tr_lambda_C.append(tr_lambda_n_C)
    epsilon = _mp(sources["pion_decay_constant"] / sources["alp_decay_constant"])
    return {"unmixed": unmixed, "tr_lambda_C": tr_lambda_C, "epsilon": epsilon,
            "alp_entry": _mp(canonical[MESON_COUNT, MESON_COUNT]), "canonical": canonical}


def arrow_checks(base_sources=None, label="neutral sources", epsilon_values=LO_EPSILON_VALUES):
    """Item 1 on the LO blocks."""
    rotation_difference, eigenvalue_difference, bare_difference, determinant_difference = (mpmath.mpf(0),) * 4
    for epsilon in epsilon_values:
        sources = rational_test_sources(epsilon, base_sources)
        inputs = lo_arrow_inputs(sources)
        canonical = to_high_precision(inputs["canonical"]).apply(mpmath.re)
        # rotation of the meson block: columns = coefficients of lambda_n on the Gell-Mann generators
        generators = [to_high_precision(sp.Matrix(generator)) for generator in NEUTRAL_GENERATORS]
        rotation = mpmath.eye(MESON_COUNT + 1)
        eigenbasis = mass_eigenstate_generators(sources)
        for column, (_, generator) in enumerate(eigenbasis):
            for row, basis in enumerate(generators):
                rotation[row, column] = (high_precision_trace(basis * generator) / 2).real
        rotated = rotation.T * canonical * rotation
        # compare with the arrow matrix built with the unrotated eigenbasis tr_lambda_C
        tadpole = to_high_precision(tadpole_source_from_sources(**_block_inputs(sources)))
        raw_tr_lambda_C = [high_precision_trace(generator * tadpole).real for _, generator in eigenbasis]
        arrow = arrow_matrix([value for value, _ in eigenbasis], raw_tr_lambda_C, inputs["epsilon"],
                             inputs["alp_entry"])
        rotation_difference = max(rotation_difference, max(abs(rotated[row, column] - arrow[row, column])
                                                           for row in range(6) for column in range(6)))
        poles = sorted_eigenvalues(canonical)
        arrow_poles = sorted_eigenvalues(arrow_matrix(inputs["unmixed"], inputs["tr_lambda_C"], inputs["epsilon"],
                                                      inputs["alp_entry"]))
        eigenvalue_difference = max(eigenvalue_difference, max(abs(a - b) for a, b in zip(poles, arrow_poles)))
        kinetic, mass_squared = bare_blocks(sources, TEST_CONDENSATE, TEST_SECOND_ORDER_SCALAR_TRACE)
        cholesky = mpmath.cholesky(to_high_precision(kinetic).apply(mpmath.re))
        inverse = cholesky**-1
        bare_poles = sorted_eigenvalues(inverse * to_high_precision(mass_squared).apply(mpmath.re) * inverse.T)
        bare_difference = max(bare_difference, max(abs(a - b) for a, b in zip(poles, bare_poles)))
        for point in list(inputs["unmixed"]) + [mpmath.mpf(1) / 3]:
            shifted = point * mpmath.eye(6) - arrow
            direct = mpmath.mpf(str(sp.Matrix(6, 6, lambda row, column: sp.Float(
                str(shifted[row, column]), HIGH_PRECISION_DIGITS + 10)).det(method="berkowitz")))
            formula = arrow_determinant(point, [value for value, _ in eigenbasis], raw_tr_lambda_C,
                                        inputs["epsilon"], inputs["alp_entry"])
            determinant_difference = max(determinant_difference, abs(direct - formula))
    return {
        f"1. ({label}) rotated M_can^2 = arrow matrix (50 digits)": rotation_difference < HIGH_PRECISION_TOLERANCE,
        f"1. ({label}) eigenvalues of M_can^2 = eigenvalues of the arrow matrix": (
            eigenvalue_difference < HIGH_PRECISION_TOLERANCE),
        f"1. ({label}) = generalized eigenvalues of the bare (K, M^2)": bare_difference < HIGH_PRECISION_TOLERANCE,
        f"1. ({label}) det(s - A) formula, also at s = m_n0^2": determinant_difference < HIGH_PRECISION_TOLERANCE,
    }


# ------------------------------------------------------------------ interlacing


def interlacing_holds(unmixed_sorted, poles_sorted, tolerance=HIGH_PRECISION_TOLERANCE):
    """mu_k <= m_k <= mu_{k+1} for k = 1..len(unmixed)."""
    return all(poles_sorted[index] <= value + tolerance and value <= poles_sorted[index + 1] + tolerance
               for index, value in enumerate(unmixed_sorted))


def opposite_direction_holds(unmixed_sorted, poles_sorted, tolerance=HIGH_PRECISION_TOLERANCE):
    """m_k <= mu_k <= m_{k+1}: the wrong direction (should fail)."""
    return all(value <= poles_sorted[index] + tolerance for index, value in enumerate(unmixed_sorted))


def interlacing_checks():
    """Item 2 on LO blocks (neutral and flavour-diagonal sources, epsilon = 1/2, 1, 2) and on arrow models."""
    correct, opposite = True, False
    for base in (None, illustrative_lo_inputs("flavor_diagonal")):
        for epsilon in LO_EPSILON_VALUES:
            inputs = lo_arrow_inputs(rational_test_sources(epsilon, base))
            poles = sorted_eigenvalues(to_high_precision(inputs["canonical"]).apply(mpmath.re))
            correct = correct and interlacing_holds(sorted(inputs["unmixed"]), poles)
            opposite = opposite or opposite_direction_holds(sorted(inputs["unmixed"]), poles)
    for alp_entry in ("-3", "0.5", "2.5", "9"):
        for epsilon in MODEL_EPSILON_VALUES:
            unmixed = [mpmath.mpf(value) for value in ("0", "1", "2", "3", "4")]
            tr_lambda_C = [mpmath.mpf(value) for value in ("1", "-0.5", "0.7", "0.3", "-1.2")]
            poles = sorted_eigenvalues(arrow_matrix(unmixed, tr_lambda_C, mpmath.mpf(epsilon), mpmath.mpf(alp_entry)))
            correct = correct and interlacing_holds(unmixed, poles)
            opposite = opposite or opposite_direction_holds(unmixed, poles)
    return {"2. mu_k <= m_n'0^2 <= mu_{k+1} in all cases": correct,
            "2. the opposite direction m_n'0^2 <= mu_k fails in every case": not opposite}


# ------------------------------------------------------------------ isolated meson pole


def energy_dependent_alp_mass(squared_momentum, meson, unmixed, tr_lambda_C, epsilon, alp_entry):
    """Mt_n(s) = M + epsilon^2 sum_{n' != n, tr(lambda_n' 𝒞) != 0} tr(lambda_n' 𝒞)^2/(s - m_n'0^2)."""
    return alp_entry + epsilon**2 * mpmath.fsum(
        tr_lambda_C[index]**2 / (squared_momentum - unmixed[index])
        for index in range(len(unmixed)) if index != meson and tr_lambda_C[index] != 0)


def isolated_meson_poles(meson, unmixed, tr_lambda_C, epsilon, alp_entry):
    """(s_-, s_+, Mt_n0, meson root): two-level roots with Mt_n frozen at s = m_n0^2."""
    frozen = energy_dependent_alp_mass(unmixed[meson], meson, unmixed, tr_lambda_C, epsilon, alp_entry)
    centre = (unmixed[meson] + frozen) / 2
    half_splitting = mpmath.sqrt((frozen - unmixed[meson])**2 / 4 + epsilon**2 * tr_lambda_C[meson]**2)
    lower, upper = centre - half_splitting, centre + half_splitting
    return lower, upper, frozen, (lower if frozen >= unmixed[meson] else upper)


def isolated_meson_model(spacing, epsilon):
    """Arrow model: meson n = 0 at 1 GeV^2, M = 11/10 GeV^2, other mesons at 1 + spacing (1, 3/2, 2, 3),
    all tr(lambda_n 𝒞) = 1 GeV^2."""
    unmixed = [mpmath.mpf(1)] + [1 + spacing * mpmath.mpf(ratio) for ratio in ("1", "1.5", "2", "3")]
    return unmixed, [mpmath.mpf(1)] * 5, mpmath.mpf(epsilon), mpmath.mpf("1.1")


def isolated_meson_scan():
    """Item 3 on the arrow model: errors of s_-, s_+ against mu_1, mu_2 and the exact bound."""
    rows, bound_holds = [], True
    for epsilon in MODEL_EPSILON_VALUES:
        for spacing in MODEL_SPACINGS:
            unmixed, tr_lambda_C, epsilon_value, alp_entry = isolated_meson_model(spacing, epsilon)
            poles = sorted_eigenvalues(arrow_matrix(unmixed, tr_lambda_C, epsilon_value, alp_entry))
            lower, upper, frozen, _ = isolated_meson_poles(0, unmixed, tr_lambda_C, epsilon_value, alp_entry)
            errors = []
            for approximation, pole in ((lower, poles[0]), (upper, poles[1])):
                error = abs(approximation - pole)
                bound = abs(energy_dependent_alp_mass(pole, 0, unmixed, tr_lambda_C, epsilon_value, alp_entry) - frozen)
                bound_holds = bound_holds and error <= bound + HIGH_PRECISION_TOLERANCE
                errors.append(error)
            rows.append((epsilon, spacing, poles[0], poles[1], errors[0], errors[1]))
    exponents = {}
    for epsilon in MODEL_EPSILON_VALUES:
        selected = [row for row in rows if row[0] == epsilon and row[1] >= 16]
        exponents[epsilon] = float(np.polyfit(np.log([row[1] for row in selected]),
                                              np.log([float(max(row[4], row[5])) for row in selected]), 1)[0])
    return rows, bound_holds, exponents


def isolated_meson_lo(base_sources=None):
    """Item 3 on the LO blocks: for each meson with tr(lambda_n 𝒞) != 0, error of the meson root and the exact bound."""
    rows, bound_holds = [], True
    for epsilon in LO_EPSILON_VALUES:
        inputs = lo_arrow_inputs(rational_test_sources(epsilon, base_sources))
        unmixed, tr_lambda_C = inputs["unmixed"], inputs["tr_lambda_C"]
        poles = sorted_eigenvalues(arrow_matrix(unmixed, tr_lambda_C, inputs["epsilon"], inputs["alp_entry"]))
        # remove the exact poles of decoupled mesons (tr(lambda_n 𝒞) = 0, e.g. one K^0 combination) before matching
        for value, tr_lambda_n_C in zip(unmixed, tr_lambda_C):
            if tr_lambda_n_C == 0:
                poles.remove(min(poles, key=lambda pole: abs(pole - value)))
        for meson in range(len(unmixed)):
            if tr_lambda_C[meson] == 0:
                continue
            lower, upper, frozen, meson_root = isolated_meson_poles(meson, unmixed, tr_lambda_C, inputs["epsilon"],
                                                                    inputs["alp_entry"])
            pole = min(poles, key=lambda value: abs(value - meson_root))
            bound = abs(energy_dependent_alp_mass(pole, meson, unmixed, tr_lambda_C, inputs["epsilon"],
                                                  inputs["alp_entry"]) - frozen)
            bound_holds = bound_holds and min(abs(lower - pole), abs(upper - pole)) <= bound + HIGH_PRECISION_TOLERANCE
            rows.append((epsilon, unmixed[meson], pole, abs(meson_root - pole), bound))
    return rows, bound_holds


# ------------------------------------------------------------------ isolated ALP pole


def alp_fixed_point_map(squared_momentum, unmixed, tr_lambda_C, epsilon, alp_entry):
    """g(s) = M + epsilon^2 sum_n tr(lambda_n 𝒞)^2/(s - m_n0^2)."""
    return alp_entry + epsilon**2 * mpmath.fsum(tr_lambda_n_C**2 / (squared_momentum - value)
                                               for value, tr_lambda_n_C in zip(unmixed, tr_lambda_C) if tr_lambda_n_C != 0)


def alp_map_slope(squared_momentum, unmixed, tr_lambda_C, epsilon):
    """|g'(s)| = epsilon^2 sum_n tr(lambda_n 𝒞)^2/(s - m_n0^2)^2."""
    return epsilon**2 * mpmath.fsum(tr_lambda_n_C**2 / (squared_momentum - value)**2
                                    for value, tr_lambda_n_C in zip(unmixed, tr_lambda_C) if tr_lambda_n_C != 0)


def alp_pole_iteration(unmixed, tr_lambda_C, epsilon, alp_entry, tolerance=mpmath.mpf("1e-45"), maximum=2000):
    """s^(k+1) = g(s^(k)) from s^(0) = M; returns (last iterate, number of iterations, converged)."""
    current = alp_entry
    for iteration in range(1, maximum + 1):
        following = alp_fixed_point_map(current, unmixed, tr_lambda_C, epsilon, alp_entry)
        if abs(following - current) < tolerance:
            return following, iteration, True
        current = following
    return current, maximum, False


def isolated_alp_model(distance, epsilon, side="above"):
    """Arrow model: mesons at 40..44 GeV^2, tr(lambda_n 𝒞) = (1, -1/2, 7/10, 3/10, -6/5) GeV^2, M = 44 + d (above),
    40 - d (below) or 42 + d (between, d < 1)."""
    unmixed = [mpmath.mpf(value) for value in (40, 41, 42, 43, 44)]
    tr_lambda_C = [mpmath.mpf(value) for value in ("1", "-0.5", "0.7", "0.3", "-1.2")]
    alp_entry = {"above": 44, "below": 40, "between": 42}[side] + (-1 if side == "below" else 1) * mpmath.mpf(distance)
    return unmixed, tr_lambda_C, mpmath.mpf(epsilon), alp_entry


def iteration_consistent(iterate, converged, poles, unmixed, tr_lambda_C, epsilon):
    """If the iteration converges, it converges to a pole with |g'| <= 1."""
    if not converged:
        return True
    pole = min(poles, key=lambda value: abs(value - iterate))
    return abs(iterate - pole) < mpmath.mpf("1e-40") and alp_map_slope(pole, unmixed, tr_lambda_C, epsilon) <= 1


def isolated_alp_scan():
    """Item 4 on the arrow model (ALP above, below and between the mesons)."""
    rows, bound_holds, convergence_correct, local_bound_holds = [], True, True, True
    cases = [(side, distance) for side in ("above", "below") for distance in (1, 2, 4, 8, 16, 32)]
    cases += [("between", distance) for distance in ("0.5",)]
    for side, distance in cases:
        for epsilon in MODEL_EPSILON_VALUES:
            unmixed, tr_lambda_C, epsilon_value, alp_entry = isolated_alp_model(distance, epsilon, side)
            poles = sorted_eigenvalues(arrow_matrix(unmixed, tr_lambda_C, epsilon_value, alp_entry))
            frozen = alp_fixed_point_map(alp_entry, unmixed, tr_lambda_C, epsilon_value, alp_entry)
            iterate, iterations, converged = alp_pole_iteration(unmixed, tr_lambda_C, epsilon_value, alp_entry)
            convergence_correct = convergence_correct and iteration_consistent(iterate, converged, poles, unmixed,
                                                                               tr_lambda_C, epsilon_value)
            ratio = alp_map_slope(alp_entry, unmixed, tr_lambda_C, epsilon_value)
            if side == "between":
                pole = min(poles, key=lambda value: abs(value - iterate)) if converged else poles[3]
            else:
                pole = poles[-1] if side == "above" else poles[0]
                # |g'(s_a)| <= |s_a - M| / |s_a - m_nearest| < 1 for M outside the meson masses
                nearest = unmixed[-1] if side == "above" else unmixed[0]
                local_bound_holds = local_bound_holds and (
                    alp_map_slope(pole, unmixed, tr_lambda_C, epsilon_value)
                    <= abs(pole - alp_entry) / abs(pole - nearest) + HIGH_PRECISION_TOLERANCE < 1)
                convergence_correct = convergence_correct and converged and abs(iterate - pole) < mpmath.mpf("1e-40")
                if ratio < 1:
                    bound_holds = bound_holds and (
                        abs(frozen - pole) <= ratio / (1 + ratio) * abs(frozen - alp_entry) + HIGH_PRECISION_TOLERANCE
                        and min(alp_entry, frozen) <= pole <= max(alp_entry, frozen))
            rows.append((side, epsilon, distance, pole, abs(frozen - pole), ratio,
                         alp_map_slope(pole, unmixed, tr_lambda_C, epsilon_value), iterations, converged))
    return rows, bound_holds, convergence_correct, local_bound_holds


def isolated_alp_lo(base_sources=None):
    """Item 4 on the LO blocks: ALP corner M, the pole closest to the frozen value, |g'| at that pole, iteration."""
    rows = []
    for epsilon in LO_EPSILON_VALUES:
        inputs = lo_arrow_inputs(rational_test_sources(epsilon, base_sources))
        arguments = (inputs["unmixed"], inputs["tr_lambda_C"], inputs["epsilon"])
        frozen = alp_fixed_point_map(inputs["alp_entry"], *arguments, inputs["alp_entry"])
        poles = sorted_eigenvalues(arrow_matrix(*arguments, inputs["alp_entry"]))
        pole = min(poles, key=lambda value: abs(value - frozen))
        iterate, iterations, converged = alp_pole_iteration(*arguments, inputs["alp_entry"])
        rows.append((epsilon, inputs["alp_entry"], frozen, pole, alp_map_slope(pole, *arguments), iterate, converged))
    return rows


def isolated_alp_lo_consistent(rows):
    """The iteration converges, and to the pole nearest to s^(1), which has |g'| < 1 (M above all mesons here)."""
    return all(converged and abs(iterate - pole) < mpmath.mpf("1e-40") and slope < 1
               for _, _, _, pole, slope, iterate, converged in rows)


# ------------------------------------------------------------------ driver


def derive():
    """Run all checks; return checks and tables."""
    mpmath.mp.dps = HIGH_PRECISION_DIGITS + 10
    meson_rows, meson_bound, exponents = isolated_meson_scan()
    meson_lo_rows, meson_lo_bound = isolated_meson_lo()
    alp_rows, alp_bound, alp_convergence, alp_local = isolated_alp_scan()
    alp_lo_rows = isolated_alp_lo()
    checks = {}
    checks.update(arrow_checks(None, "neutral sources"))
    checks.update(arrow_checks(illustrative_lo_inputs("flavor_diagonal"), "flavour-diagonal sources"))
    checks.update(interlacing_checks())
    checks.update({
        "3. model: each pole within |Mt_n(mu) - Mt_n0| of s_+ or s_-": meson_bound,
        "3. model: error ~ spacing^-2 at fixed epsilon (fitted exponent within 0.1 of -2, spacing >= 16)": all(
            abs(value + 2) < 0.1 for value in exponents.values()),
        "3. LO blocks: bound holds for every meson, epsilon = 1/2, 1, 2": meson_lo_bound,
        "4. model: pole between M and s^(1), |s^(1) - s_a| <= q/(1+q)|s^(1) - M| when q < 1": alp_bound,
        "4. model: M outside the meson masses: |g'(s_a)| <= |s_a - M|/|s_a - m_nearest| < 1": alp_local,
        "4. model: iteration from M converges to s_a (M outside), and only ever to poles with |g'| <= 1": alp_convergence,
        "4. LO blocks: iteration converges to the pole nearest s^(1), with |g'| < 1": isolated_alp_lo_consistent(
            alp_lo_rows),
    })
    return {"checks": checks, "meson_rows": meson_rows, "exponents": exponents, "meson_lo_rows": meson_lo_rows,
            "alp_rows": alp_rows, "alp_lo_rows": alp_lo_rows}


def main():
    result = derive()
    short = lambda value: mpmath.nstr(value, 4)
    print("Item 3, arrow model (m_n0^2 = 1, M = 1.1, others at 1 + spacing (1, 1.5, 2, 3), tr(lambda_n 𝒞) = 1; GeV^2):")
    print("  epsilon  spacing   mu_1        mu_2        |s_- - mu_1|  |s_+ - mu_2|")
    for epsilon, spacing, first, second, lower_error, upper_error in result["meson_rows"]:
        print(f"  {epsilon:>7}  {spacing:>7}   {short(first):<10}  {short(second):<10}  "
              f"{short(lower_error):<12}  {short(upper_error)}")
    print("  fitted exponent of the error in the spacing (spacing >= 16): "
          + ", ".join(f"epsilon = {key}: {value:.3f}" for key, value in result["exponents"].items()))
    print("Item 3, LO blocks (neutral sources, m_a0^2 = 3/5 GeV^2):")
    print("  epsilon  m_n0^2      pole        |error|     bound")
    for epsilon, unmixed, pole, error, bound in result["meson_lo_rows"]:
        print(f"  {str(epsilon):>7}  {short(unmixed):<10}  {short(pole):<10}  {short(error):<10}  {short(bound)}")
    print("Item 4, arrow model (mesons at 40..44 GeV^2; M = 44 + d above, 40 - d below, 42 + d between):")
    print("  side     epsilon  d    s_a          |s^(1) - s_a|  q        |g'(s_a)|  iterations")
    for side, epsilon, distance, pole, error, ratio, slope, iterations, converged in result["alp_rows"]:
        print(f"  {side:<8} {epsilon:>7}  {distance:<3}  {short(pole):<11}  {short(error):<13}  "
              f"{short(ratio):<7}  {short(slope):<9}  {iterations if converged else 'no convergence'}")
    print("Item 4, LO blocks (neutral sources):")
    print("  epsilon  M           s^(1)       nearest pole  |g'|      iteration")
    for epsilon, alp_entry, frozen, pole, slope, iterate, converged in result["alp_lo_rows"]:
        print(f"  {str(epsilon):>7}  {short(alp_entry):<10}  {short(frozen):<10}  {short(pole):<12}  {short(slope):<8}  "
              f"{short(iterate) if converged else 'no convergence'}")
    for label, passed in result["checks"].items():
        print(f"[{'ok' if passed else 'FAIL'}] {label}")
    passed = all(result["checks"].values())
    print("All checks pass" if passed else "Some checks FAIL")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
