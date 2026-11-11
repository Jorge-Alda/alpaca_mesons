#!/usr/bin/env python3
"""Checks of the general pole solver on exact, random and LO (P_i, a) data.

Run:   .venv/bin/python scripts/verify_alp_lo_spectrum.py
Tests: .venv/bin/python -m unittest discover -s scripts -p test_verify_alp_lo_spectrum.py

Conventions as in scripts/solve_alp_meson_spectrum.py: D(s) = s K - M^2,
pole vectors v^(n) with v^(n)T K v^(n) = 1, residues R_n = v^(n) v^(n)T, ALP
as the last field. The single-pole quantities of
memory/alp/shared/finite-alp-meson-mixing.md are
  w^(a)_i(s) = [(M^2_PP - s K_PP)^{-1}]_ij (s K_ja - M^2_ja),
  N_a(s)     = K_aa + 2 K_ia w_i + K_ij w_i w_j   (the pole normalization),
and the LSZ factor is N_a(m_a^2)^{-1/2}. With v^(n) normalized by K one
expects v^(n) = +-(w^(n), 1)/sqrt(N(m_n^2)), i.e.
  w^(n)_i = v^(n)_i / v^(n)_a,   v^(n)_a^2 N(m_n^2) = 1,   (R_n)_aa = 1/N(m_n^2).

Checks, grouped by the tasks of openspec change solve-lo-finite-mixing-spectrum:
  1.2  null vectors, K-orthonormality, sum_n v v^T = K^{-1}, and
       D^{-1}(s) = sum_n R_n/(s - m_n^2) against a direct inverse, for
       (a) exact rational data with known rational poles (random rational S,
       K = S^T S, M^2 = S^T diag(m^2) S, and the LO-form data of
       derive_alp_inverse_propagator_congruence.py), solved exactly and in
       floating point; (b) random real symmetric (K, M^2), seeded;
       (c) the LO blocks of derive_alp_lo_neutral_blocks.py on the neutral
       6x6 basis (P_3, P_8, P_0, P_6, P_7, a), floating point, for
       flavour-diagonal, CP-conserving (real d-s entries) and generic neutral
       (complex d-s entries) sources.
  6x6 vs 10x10: for charge-conserving sources the 10x10 spectrum on
       (P_0..P_8, a) is the union of the 6x6 spectrum and of the charged
       (P_1, P_2, P_4, P_5) block, and the neutral restriction of the 10x10
       D^{-1}(s) equals the 6x6 D^{-1}(s). The exact pole at
       m_K0^2 = (chi_d + chi_s)/2 with zero ALP and zero (P_3, P_8, P_0)
       components (scripts/derive_alp_neutral_block_reduction.py) is found.
  1.3  agreement with w^(a)(m_a^2) and N_a of the finite-mixing record, both
       from the explicit formulas and from the symbolic expressions of
       derive_alp_finite_mixing.py (including d/ds of the Schur complement).
  1.4  invariance under Phi = T Phi' with a random rational invertible T:
       same poles, v' = +-T^{-1} v, R' = T^{-1} R T^{-T}.

Units: K dimensionless, M^2 and s in a common squared-mass unit (GeV^2 for
the LO blocks). All test numbers are illustrative.
"""

from __future__ import annotations

import random

import numpy as np
import sympy as sp

from derive_alp_finite_mixing import derive as derive_finite_mixing
from derive_alp_inverse_propagator_congruence import generic_invertible_matrix, neutral_test_quadratic_form
from derive_alp_lo_neutral_blocks import (
    ALL_GENERATOR_INDICES,
    CHARGED_GENERATOR_INDICES,
    NEUTRAL_FIELD_LABELS,
    NEUTRAL_GENERATOR_INDICES,
    closed_form_blocks,
    illustrative_lo_inputs,
)
from solve_alp_meson_spectrum import inverse_propagator, solve_spectrum_exact, solve_spectrum_numerical

NUMERICAL_TOLERANCE = 1e-9
SQUARED_MOMENTUM = sp.Symbol("s")


# ------------------------------------------------------------------- inputs


def random_rational_quadratic_form(size=4, seed=20261111):
    """K = S^T S and M^2 = S^T diag(m^2) S with random rational S and distinct rational poles.

    The expected K-normalized pole vectors are the columns of S^{-1} (up to sign).
    """
    generator = random.Random(seed)
    conversion = generic_invertible_matrix(size=size, seed=seed)
    poles = []
    while len(poles) < size:
        candidate = sp.Rational(generator.randint(1, 99), generator.randint(1, 20))
        if candidate not in poles:
            poles.append(candidate)
    return {
        "kinetic": conversion.T * conversion,
        "mass_squared": conversion.T * sp.diag(*poles) * conversion,
        "pole_squared_masses": sorted(poles),
        "expected_vectors": conversion.inv(),
        "expected_poles_unsorted": poles,
    }


def random_symmetric_quadratic_form(size, seed):
    """Random real symmetric (K, M^2) with K positive definite; NumPy floats."""
    generator = np.random.default_rng(seed)
    square_root = generator.normal(size=(size, size))
    kinetic = square_root @ square_root.T + 0.5 * np.eye(size)
    mass = generator.normal(size=(size, size))
    return kinetic, (mass + mass.T) / 2


LO_SOURCE_STRUCTURES = ("flavor_diagonal", "cp_conserving", "neutral")


def lo_blocks_numerical(source_structure="neutral", generator_indices=NEUTRAL_GENERATOR_INDICES):
    """LO K and M^2 of derive_alp_lo_neutral_blocks.py as floats; default neutral 6x6."""
    kinetic, mass_squared = closed_form_blocks(**illustrative_lo_inputs(source_structure),
                                               generator_indices=generator_indices)
    return np.array(kinetic.evalf(30), dtype=float), np.array(mass_squared.evalf(30), dtype=float)


# --------------------------------------------------------- task 1.2 checks


def same_up_to_sign_exact(left, right):
    """True if the exact column vectors satisfy left = +right or left = -right."""
    zero = sp.zeros(left.rows, 1)
    return ((left - right).applyfunc(sp.simplify) == zero
            or (left + right).applyfunc(sp.simplify) == zero)


def exact_spectrum_checks(spectrum, expected_poles=None, symbolic_propagator=True):
    """Exact checks of a PoleSpectrum returned by solve_spectrum_exact."""
    kinetic, mass_squared = spectrum.kinetic, spectrum.mass_squared
    size = spectrum.size
    vectors = spectrum.pole_vectors
    zero = sp.zeros(size)
    checks = {
        "D(m_n^2) v^(n) = 0": all(
            (inverse_propagator(kinetic, mass_squared, pole) * spectrum.pole_vector(index)).applyfunc(sp.simplify)
            == sp.zeros(size, 1)
            for index, pole in enumerate(spectrum.pole_squared_masses)),
        "v^(m)T K v^(n) = delta_mn": (vectors.T * kinetic * vectors).applyfunc(sp.simplify) == sp.eye(size),
        "sum_n v v^T = K^-1": (spectrum.completeness_sum() - kinetic.inv()).applyfunc(sp.simplify) == zero,
    }
    if expected_poles is not None:
        checks["poles equal the constructed values"] = list(spectrum.pole_squared_masses) == sorted(expected_poles)
    test_points = [sp.Rational(7, 3), sp.Rational(-1, 2)]
    checks["D^-1(s) = sum R/(s-m^2) at rational s"] = all(
        (spectrum.propagator_from_poles(point) - inverse_propagator(kinetic, mass_squared, point).inv())
        .applyfunc(sp.simplify) == zero for point in test_points)
    if symbolic_propagator:
        difference = (inverse_propagator(kinetic, mass_squared, SQUARED_MOMENTUM)
                      * spectrum.propagator_from_poles(SQUARED_MOMENTUM)).applyfunc(sp.cancel)
        checks["D(s) sum R/(s-m^2) = 1 as rational functions of s"] = (
            difference - sp.eye(size)).applyfunc(sp.simplify) == zero
    return checks


def numerical_spectrum_checks(spectrum, test_points=(0.37, -1.3, 2.0 + 0.5j)):
    """Floating-point checks of a PoleSpectrum; relative tolerance NUMERICAL_TOLERANCE."""
    kinetic, mass_squared = spectrum.kinetic, spectrum.mass_squared
    vectors = spectrum.pole_vectors
    scale = max(1.0, np.max(np.abs(mass_squared)), np.max(np.abs(spectrum.pole_squared_masses)))
    inverse_kinetic = np.linalg.inv(kinetic)

    def close(left, right, reference=1.0):
        return bool(np.max(np.abs(left - right)) <= NUMERICAL_TOLERANCE * max(1.0, reference))

    return {
        "M^2 v = m^2 K v": close(mass_squared @ vectors, kinetic @ vectors * np.array(spectrum.pole_squared_masses), scale),
        "v^(m)T K v^(n) = delta_mn": close(vectors.T @ kinetic @ vectors, np.eye(spectrum.size)),
        "sum_n v v^T = K^-1": close(spectrum.completeness_sum(), inverse_kinetic, np.max(np.abs(inverse_kinetic))),
        "D^-1(s) = sum R/(s-m^2) vs direct inverse": all(
            close(spectrum.propagator_from_poles(point),
                  np.linalg.inv(inverse_propagator(kinetic, mass_squared, point)),
                  np.max(np.abs(np.linalg.inv(inverse_propagator(kinetic, mass_squared, point)))))
            for point in test_points),
    }


def exact_and_numerical_agree(exact_spectrum):
    """Compare the exact solution with the floating-point solver on the same data."""
    numerical = solve_spectrum_numerical(np.array(exact_spectrum.kinetic.evalf(30), dtype=float),
                                         np.array(exact_spectrum.mass_squared.evalf(30), dtype=float))
    exact_poles = np.array([float(pole) for pole in exact_spectrum.pole_squared_masses])
    exact_residues = [np.array(residue.evalf(30), dtype=float) for residue in exact_spectrum.residues]
    return {
        "poles agree": bool(np.allclose(exact_poles, numerical.pole_squared_masses, rtol=1e-10, atol=1e-12)),
        "residues agree": all(np.allclose(left, right, rtol=1e-9, atol=1e-11)
                              for left, right in zip(exact_residues, numerical.residues)),
        "pole vectors agree (same sign convention)": bool(np.allclose(
            np.array(exact_spectrum.pole_vectors.evalf(30), dtype=float), numerical.pole_vectors, rtol=1e-9, atol=1e-11)),
    }


# --------------------------------------------------------- task 1.3 checks


def single_pole_quantities(kinetic, mass_squared, squared_momentum):
    """w^(a)(s) and N_a(s) of the finite-mixing record (ALP is the last field)."""
    exact = isinstance(kinetic, sp.MatrixBase)
    meson = slice(0, kinetic.shape[0] - 1)
    alp = kinetic.shape[0] - 1
    meson_block = mass_squared[meson, meson] - squared_momentum * kinetic[meson, meson]
    mixed = squared_momentum * kinetic[meson, alp] - mass_squared[meson, alp]
    if exact:
        mixed = sp.Matrix(mixed)
        overlap = meson_block.inv() * mixed
        normalization = (kinetic[alp, alp] + 2 * (kinetic[meson, alp].T * overlap)[0]
                         + (overlap.T * kinetic[meson, meson] * overlap)[0])
        return overlap.applyfunc(sp.simplify), sp.simplify(normalization)
    overlap = np.linalg.solve(meson_block, mixed)
    normalization = kinetic[alp, alp] + 2 * kinetic[meson, alp] @ overlap + overlap @ kinetic[meson, meson] @ overlap
    return overlap, float(normalization)


def single_pole_record_quantities(kinetic, mass_squared, squared_momentum):
    """w and N from the symbolic expressions of derive_alp_finite_mixing.derive().

    Also returns d/ds of its Schur complement at s, which equals N at a pole.
    """
    result = derive_finite_mixing()
    meson = slice(0, kinetic.shape[0] - 1)
    alp = kinetic.shape[0] - 1
    substitution = {
        result["kinetic_meson"]: sp.Matrix(kinetic[meson, meson]),
        result["mass_meson"]: sp.Matrix(mass_squared[meson, meson]),
        result["kinetic_mixing"]: sp.Matrix(kinetic[meson, alp]),
        result["mass_mixing"]: sp.Matrix(mass_squared[meson, alp]),
        result["kinetic_alp"]: kinetic[alp, alp],
        result["mass_alp"]: mass_squared[alp, alp],
    }

    def evaluate(expression):
        value = expression.subs(substitution).subs(result["s"], squared_momentum).doit()
        return sp.Matrix(value.as_explicit() if hasattr(value, "as_explicit") else value).applyfunc(sp.simplify)

    # There the overlap is written -(s K_P - M_P^2)^{-1}(s K_Pa - M_Pa^2), equal to w(s).
    overlap = evaluate(result["unnormalized_meson_overlap"])
    return overlap, evaluate(result["pole_norm"])[0], evaluate(result["derivative_schur"])[0]


def single_pole_checks_exact(spectrum):
    """Task 1.3 on an exact spectrum: every pole with a nonzero ALP component."""
    entries = []
    alp = spectrum.size - 1
    for index, pole in enumerate(spectrum.pole_squared_masses):
        vector = spectrum.pole_vector(index)
        if vector[alp] == 0:
            continue
        meson_block = spectrum.mass_squared[:alp, :alp] - pole * spectrum.kinetic[:alp, :alp]
        if meson_block.det() == 0:
            continue
        overlap, normalization = single_pole_quantities(spectrum.kinetic, spectrum.mass_squared, pole)
        record_overlap, record_normalization, schur_derivative = single_pole_record_quantities(
            spectrum.kinetic, spectrum.mass_squared, pole)
        entries.append({
            "pole": pole,
            "v_P / v_a = w(m^2)": (vector[:alp, 0] / vector[alp] - overlap).applyfunc(sp.simplify) == sp.zeros(alp, 1),
            "v_a^2 N(m^2) = 1": sp.simplify(vector[alp]**2 * normalization - 1) == 0,
            "(R_n)_aa = 1/N(m^2)": sp.simplify(spectrum.residues[index][alp, alp] * normalization - 1) == 0,
            "w agrees with derive_alp_finite_mixing": (record_overlap - overlap).applyfunc(sp.simplify) == sp.zeros(alp, 1),
            "N agrees with derive_alp_finite_mixing": sp.simplify(record_normalization - normalization) == 0,
            "d/ds Schur complement = N at the pole": sp.simplify(schur_derivative - normalization) == 0,
        })
    return entries


def alp_like_pole(spectrum):
    """Index of the pole with the largest |v^(n)_a| (ALP is the last field)."""
    components = [abs(complex(spectrum.pole_vector(index)[spectrum.size - 1])) for index in range(spectrum.size)]
    return int(np.argmax(components))


def single_pole_checks_numerical(spectrum, relative_tolerance=1e-9):
    """Task 1.3 in floating point for the ALP-like pole of a NumPy spectrum."""
    alp = spectrum.size - 1
    index = alp_like_pole(spectrum)
    pole = spectrum.pole_squared_masses[index]
    vector = spectrum.pole_vector(index)
    overlap, normalization = single_pole_quantities(spectrum.kinetic, spectrum.mass_squared, pole)
    return {
        "pole": pole,
        "v_P / v_a = w(m^2)": bool(np.allclose(vector[:alp] / vector[alp], overlap, rtol=relative_tolerance, atol=1e-12)),
        "v_a^2 N(m^2) = 1": bool(abs(vector[alp]**2 * normalization - 1) < relative_tolerance),
        "(R_n)_aa = 1/N(m^2)": bool(abs(spectrum.residues[index][alp, alp] * normalization - 1) < relative_tolerance),
        "LSZ factor N^(-1/2)": normalization**-0.5,
        "ALP component v_a": vector[alp],
    }


# --------------------------------------------------------- task 1.4 checks


def pole_groups(poles, tolerance=0):
    """Index lists of poles that coincide within `tolerance` (exact equality for 0).

    The solver returns sorted poles, so neighbouring poles are compared. The
    individual vectors of a degenerate pole are basis dependent; only the sum
    of their residues is compared there.
    """
    groups = [[0]]
    for index in range(1, len(poles)):
        if abs(poles[index] - poles[groups[-1][-1]]) <= tolerance:
            groups[-1].append(index)
        else:
            groups.append([index])
    return groups


def field_change_checks_exact(kinetic, mass_squared, seed):
    """Exact poles, vectors and residues before and after Phi = T Phi'."""
    kinetic, mass_squared = sp.Matrix(kinetic), sp.Matrix(mass_squared)
    size = kinetic.rows
    conversion = generic_invertible_matrix(size=size, seed=seed)
    inverse_conversion = conversion.inv()
    before = solve_spectrum_exact(kinetic, mass_squared)
    after = solve_spectrum_exact(conversion.T * kinetic * conversion, conversion.T * mass_squared * conversion)
    groups = pole_groups(before.pole_squared_masses)
    return {
        "poles unchanged": before.pole_squared_masses == after.pole_squared_masses,
        "v' = +-T^-1 v (simple poles)": all(
            same_up_to_sign_exact(after.pole_vector(group[0]), inverse_conversion * before.pole_vector(group[0]))
            for group in groups if len(group) == 1),
        "R' = T^-1 R T^-T (summed over each pole)": all(
            (sum((after.residues[index] for index in group), sp.zeros(size))
             - inverse_conversion * sum((before.residues[index] for index in group), sp.zeros(size))
             * inverse_conversion.T).applyfunc(sp.simplify) == sp.zeros(size)
            for group in groups),
    }


def field_change_checks_numerical(kinetic, mass_squared, seed, relative_tolerance=1e-8):
    """Floating-point version (e.g. LO blocks), with an exact rational T converted to floats.

    Poles closer than 1e-8 times the largest |m^2| (or 1e-8) are grouped as
    degenerate, e.g. the K0 pair (P_6, P_7) for flavour-diagonal sources.
    """
    size = kinetic.shape[0]
    conversion = np.array(generic_invertible_matrix(size=size, seed=seed).evalf(), dtype=float)
    inverse_conversion = np.linalg.inv(conversion)
    before = solve_spectrum_numerical(kinetic, mass_squared)
    after = solve_spectrum_numerical(conversion.T @ kinetic @ conversion, conversion.T @ mass_squared @ conversion)
    scale = max(1.0, np.max(np.abs(before.pole_squared_masses)))
    groups = pole_groups(before.pole_squared_masses, tolerance=1e-8 * scale)

    def close(left, right):
        return bool(np.max(np.abs(left - right)) <= relative_tolerance * max(1.0, np.max(np.abs(right))))

    def summed(residues, group):
        return np.sum([residues[index] for index in group], axis=0)

    simple = [group[0] for group in groups if len(group) == 1]
    mapped_vectors = inverse_conversion @ before.pole_vectors
    return {
        "number of degenerate poles": sum(1 for group in groups if len(group) > 1),
        "poles unchanged": bool(np.allclose(before.pole_squared_masses, after.pole_squared_masses,
                                            rtol=relative_tolerance, atol=relative_tolerance * scale)),
        "v' = +-T^-1 v (simple poles)": all(
            close(after.pole_vectors[:, index], mapped_vectors[:, index])
            or close(after.pole_vectors[:, index], -mapped_vectors[:, index]) for index in simple),
        "R' = T^-1 R T^-T (summed over each pole)": all(
            close(summed(after.residues, group),
                  inverse_conversion @ summed(before.residues, group) @ inverse_conversion.T)
            for group in groups),
    }


# ------------------------------------------------- neutral 6x6 vs full 10x10


def neutral_versus_full_checks(source_structure, test_points=(0.37, -1.3, 2.0 + 0.5j), relative_tolerance=1e-9):
    """Compare the 6x6 neutral system with the full 10x10 system (floating point)."""
    neutral = solve_spectrum_numerical(*lo_blocks_numerical(source_structure), field_labels=NEUTRAL_FIELD_LABELS)
    full_kinetic, full_mass = lo_blocks_numerical(source_structure, ALL_GENERATOR_INDICES)
    full = solve_spectrum_numerical(full_kinetic, full_mass)
    charged = list(CHARGED_GENERATOR_INDICES)
    charged_poles = np.linalg.eigvalsh(full_mass[np.ix_(charged, charged)])
    neutral_positions = list(NEUTRAL_GENERATOR_INDICES) + [len(ALL_GENERATOR_INDICES)]
    scale = max(1.0, np.max(np.abs(full.pole_squared_masses)))

    def restricted_propagator_agrees(point):
        full_inverse = np.linalg.inv(inverse_propagator(full_kinetic, full_mass, point))
        restricted = full_inverse[np.ix_(neutral_positions, neutral_positions)]
        neutral_inverse = np.linalg.inv(inverse_propagator(neutral.kinetic, neutral.mass_squared, point))
        return bool(np.max(np.abs(restricted - neutral_inverse)) <= relative_tolerance * max(1.0, np.max(np.abs(neutral_inverse))))

    return {
        "K_PP of the charged block is the unit matrix":
            bool(np.allclose(full_kinetic[np.ix_(charged, charged)], np.eye(4), atol=1e-14)),
        "10x10 poles = 6x6 poles and charged-block poles": bool(np.allclose(
            np.sort(full.pole_squared_masses), np.sort(np.concatenate([neutral.pole_squared_masses, charged_poles])),
            rtol=relative_tolerance, atol=relative_tolerance * scale)),
        "neutral restriction of the 10x10 D^-1(s) = 6x6 D^-1(s)": all(restricted_propagator_agrees(point) for point in test_points),
    }


def neutral_kaon_pole_checks(spectrum, inputs, relative_tolerance=1e-9):
    """Pole at m_K0^2 = (chi_d + chi_s)/2 with pole vector only on (P_6, P_7) (6x6 numerical spectrum)."""
    spurion = inputs["quark_spurion"]
    kaon_mass_squared = float((spurion[1, 1] + spurion[2, 2]) / 2)
    poles = np.array(spectrum.pole_squared_masses)
    matching = [index for index in range(spectrum.size)
                if abs(poles[index] - kaon_mass_squared) <= relative_tolerance * max(1.0, kaon_mass_squared)]
    return {
        "number of poles at m_K0^2": len(matching),
        "a pole at m_K0^2 exists": bool(matching),
        "its pole vectors have zero ALP and (P_3, P_8, P_0) components": bool(matching) and all(
            np.max(np.abs(spectrum.pole_vector(index)[[0, 1, 2, 5]])) <= 1e-9 for index in matching),
    }


# --------------------------------------------------------------------- run


def derive():
    """Run every check; returns a dict of named groups of boolean checks."""
    groups = {}
    for seed in (20261111, 7):
        data = random_rational_quadratic_form(size=4, seed=seed)
        spectrum = solve_spectrum_exact(data["kinetic"], data["mass_squared"])
        checks = exact_spectrum_checks(spectrum, data["pole_squared_masses"])
        order = sorted(range(4), key=lambda index: data["expected_poles_unsorted"][index])
        checks["v^(n) = +-S^-1 e_n"] = all(
            same_up_to_sign_exact(spectrum.pole_vector(position), data["expected_vectors"][:, original])
            for position, original in enumerate(order))
        checks.update(exact_and_numerical_agree(spectrum))
        groups[f"1.2 exact random rational S, seed {seed}"] = checks
    lo_form = neutral_test_quadratic_form()
    lo_form_spectrum = solve_spectrum_exact(lo_form["kinetic"], lo_form["mass_squared"])
    checks = exact_spectrum_checks(lo_form_spectrum, lo_form["pole_squared_masses"])
    checks.update(exact_and_numerical_agree(lo_form_spectrum))
    groups["1.2 exact LO-form data of the congruence script"] = checks
    for size, seed in ((4, 1), (10, 2), (10, 3)):
        kinetic, mass_squared = random_symmetric_quadratic_form(size, seed)
        groups[f"1.2 random symmetric, size {size}, seed {seed}"] = numerical_spectrum_checks(
            solve_spectrum_numerical(kinetic, mass_squared))
    for source_structure in LO_SOURCE_STRUCTURES:
        label = f"LO 6x6 blocks, {source_structure} sources"
        kinetic, mass_squared = lo_blocks_numerical(source_structure)
        spectrum = solve_spectrum_numerical(kinetic, mass_squared, field_labels=NEUTRAL_FIELD_LABELS)
        groups[f"1.2 {label}"] = numerical_spectrum_checks(spectrum)
        groups[f"1.3 {label}, ALP-like pole"] = single_pole_checks_numerical(spectrum)
        groups[f"1.4 {label}"] = field_change_checks_numerical(kinetic, mass_squared, seed=20261111)
        groups[f"6x6 vs 10x10, {source_structure} sources"] = neutral_versus_full_checks(source_structure)
        groups[f"{label}, exact pole at m_K0^2"] = neutral_kaon_pole_checks(
            spectrum, illustrative_lo_inputs(source_structure))
    for entry in single_pole_checks_exact(lo_form_spectrum):
        groups[f"1.3 exact LO-form data, pole {entry['pole']}"] = entry
    for seed in (20261111, 7, 123):
        groups[f"1.4 exact LO-form data, T seed {seed}"] = field_change_checks_exact(
            lo_form["kinetic"], lo_form["mass_squared"], seed)
    return groups


def main():
    for group, checks in derive().items():
        print(group)
        for label, value in checks.items():
            if isinstance(value, bool):
                print(f"  [{'pass' if value else 'FAIL'}] {label}")
            else:
                print(f"  {label} = {value}")


if __name__ == "__main__":
    main()
