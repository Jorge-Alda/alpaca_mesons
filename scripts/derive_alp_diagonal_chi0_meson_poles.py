#!/usr/bin/env python3
"""Explicit exact LO pole equations of pi^0, eta, eta' for diagonal chi_0 = diag(x_u, x_d, x_s), general m_q,
no expansion in epsilon = f_pi/f_a and no isospin-limit rotation.

Run:   .venv/bin/python scripts/derive_alp_diagonal_chi0_meson_poles.py
Tests: .venv/bin/python -m unittest discover -s scripts -p test_derive_alp_diagonal_chi0_meson_poles.py
Record: memory/alp/shared/alp-lo-diagonal-chi0-meson-poles.md

Conventions (Alda:2026nrd; memory/alp/shared/alp-lo-exact-pole-equations.md): N_f = 3, neutral basis
(P_3, P_8, P_0, P_6, P_7, a), D(s) = s K - M^2, epsilon = f_pi/f_a, 𝒞 = -{chi_0, c_q^AP}/2 - (m_0^2/N_f) c_GA 1,
𝕄_s(X) = {chi_0, X}/2 + (m_0^2/N_f) tr(X) 1 - s X, lambda_n the generator of the mass eigenstate n,
M = M_can,aa (record: M^2_{\\hat a\\hat a,can}). The building blocks are reused from
scripts/derive_alp_meson_mass_shifts_exact.py:
    P(s) = (s - m_pi0,0^2)(s - m_eta,0^2)(s - m_eta',0^2) = prod_q (s - x_q) - (m_0^2/N_f) sum_q prod_{r != q}(s - x_r)
        (`flavour_diagonal_meson_polynomial`),
    N(s) (`mixing_numerator_polynomial`), det D(s) (`explicit_secular_polynomial`).

Derived and checked here:
 1. N(s) has degree <= 2 in s (symbolic x_q, c_q^AP, c_GA, m_0^2), so 2 N(s)/P(s) has the exact partial
    fractions sum_n R_n/(s - m_n0^2) with R_n = 2 N(m_n0^2)/prod_{n' != n}(m_n0^2 - m_n'0^2) (symbolic roots).
 2. R_n = tr(lambda_n 𝒞)^2 for n = pi^0, eta, eta' (50 digits, lambda_n from the eigenvectors of M^2_PP).
 3. Explicit per-meson equation, n = pi^0, eta, eta':
        (s - m_n0^2) [ s - M - epsilon^2 sum_{n' != n} R_n'/(s - m_n'0^2)
                       - 4 epsilon^2 m_K0,0^4 |(c_q^AP)_ds|^2/(s - m_K0,0^2) ] = epsilon^2 R_n,
    multiplied by (s - m_K0,0^2) prod_{n' != n}(s - m_n'0^2), equals the secular polynomial divided by
    (s - m_K0,0^2) (symbolic roots, generic residues).
 4. Exact rationals (neutral, CP-conserving and flavour-diagonal test sources; epsilon = 1/2 and 2): explicit
    polynomial = det(s K - M^2); exact real-root isolation (to 1e-50) against solve_spectrum_numerical; the
    per-meson equation of item 3 vanishes at every isolated root other than m_K0,0^2 (|residual| < 1e-30, 50 digits).

Units: chi_0, m_0^2, M, s, masses squared in GeV^2; couplings and epsilon dimensionless. Inputs are test values.
"""

from __future__ import annotations

import mpmath
import numpy as np
import sympy as sp

from derive_alp_lo_mass_sum_rules import _hermitian_charge_conserving
from derive_alp_lo_neutral_blocks import NEUTRAL_FIELD_LABELS, illustrative_lo_inputs
from derive_alp_meson_mass_shifts_exact import (
    ALP,
    SQUARED_MOMENTUM,
    TEST_CONDENSATE,
    TEST_SECOND_ORDER_SCALAR_TRACE,
    bare_blocks,
    canonical_mass_matrix_from_sources,
    explicit_secular_polynomial,
    flavour_diagonal_meson_polynomial,
    mixing_numerator_polynomial,
)
from derive_alp_meson_mass_shifts_second_order import (
    HIGH_PRECISION_DIGITS,
    _block_inputs,
    high_precision_trace,
    mass_eigenstate_generators,
    rational_test_sources,
    tadpole_source_from_sources,
    to_high_precision,
)
from solve_alp_meson_spectrum import solve_spectrum_numerical

EPSILON_VALUES = (sp.Rational(1, 2), sp.Integer(2))
SOURCE_STRUCTURES = ("neutral", "cp_conserving", "flavor_diagonal")
RESIDUAL_TOLERANCE = mpmath.mpf("1e-30")


def _is_zero(expression):
    return sp.simplify(sp.together(sp.expand(expression))) == 0


# ------------------------------------------------------------------ explicit forms


def explicit_inputs(sources, condensate=TEST_CONDENSATE, second_order_scalar_trace=TEST_SECOND_ORDER_SCALAR_TRACE):
    """x_q, (c_q^AP)_qq, (c_q^AP)_ds, c_GA, m_0^2, epsilon and M = M_can,aa for diagonal chi_0."""
    difference = sources["axial_source"] - sources["pseudoscalar_source"]
    return {
        "quark_spurion_diagonal": [sources["quark_spurion"][index, index] for index in range(3)],
        "difference_source": difference,
        "anomaly_mass_squared": sources["anomaly_mass_squared"],
        "gluon_axial_combination": sources["gluon_coupling"] + sp.trace(sources["axial_source"]),
        "epsilon": sources["pion_decay_constant"] / sources["alp_decay_constant"],
        "alp_entry": canonical_mass_matrix_from_sources(sources, condensate, second_order_scalar_trace)[ALP, ALP],
    }


def kaon_mass_squared(inputs):
    """m_K0,0^2 = (x_d + x_s)/2."""
    return (inputs["quark_spurion_diagonal"][1] + inputs["quark_spurion_diagonal"][2]) / 2


def kaon_weight(inputs):
    """4 m_K0,0^4 |(c_q^AP)_ds|^2."""
    difference = inputs["difference_source"]
    return sp.expand(4 * kaon_mass_squared(inputs)**2 * difference[1, 2] * difference[2, 1])


def unmixed_flavour_diagonal_masses(inputs, digits=HIGH_PRECISION_DIGITS):
    """m_pi0,0^2 < m_eta,0^2 < m_eta',0^2: the roots of P(s) (mpmath, `digits` digits)."""
    polynomial = sp.Poly(flavour_diagonal_meson_polynomial(inputs["quark_spurion_diagonal"],
                                                           inputs["anomaly_mass_squared"]), SQUARED_MOMENTUM)
    return sorted(mpmath.mpf(str(sp.re(root))) for root in polynomial.nroots(n=digits + 10))


def explicit_residues(inputs, unmixed, digits=HIGH_PRECISION_DIGITS):
    """R_n = 2 N(m_n0^2)/prod_{n' != n}(m_n0^2 - m_n'0^2) for n = pi^0, eta, eta' (mpmath)."""
    numerator = mixing_numerator_polynomial(
        inputs["quark_spurion_diagonal"], [inputs["difference_source"][index, index] for index in range(3)],
        inputs["anomaly_mass_squared"], inputs["gluon_axial_combination"])
    residues = []
    for index, value in enumerate(unmixed):
        at_root = mpmath.mpf(str(sp.N(numerator.subs(SQUARED_MOMENTUM, sp.Float(str(value), digits + 10)),
                                      digits + 10)))
        residues.append(2 * at_root / mpmath.fprod(value - other for other_index, other in enumerate(unmixed)
                                                   if other_index != index))
    return residues


def meson_equation_residual(squared_momentum, meson, unmixed, residues, inputs):
    """Left side minus right side of the explicit per-meson equation (item 3), mpmath."""
    epsilon = mpmath.mpf(str(sp.N(inputs["epsilon"], HIGH_PRECISION_DIGITS + 10)))
    alp_entry = mpmath.mpf(str(sp.N(inputs["alp_entry"], HIGH_PRECISION_DIGITS + 10)))
    kaon = mpmath.mpf(str(sp.N(kaon_mass_squared(inputs), HIGH_PRECISION_DIGITS + 10)))
    weight = mpmath.mpf(str(sp.N(kaon_weight(inputs), HIGH_PRECISION_DIGITS + 10)))
    others = mpmath.fsum(residues[other] / (squared_momentum - unmixed[other])
                         for other in range(3) if other != meson)
    return ((squared_momentum - unmixed[meson])
            * (squared_momentum - alp_entry - epsilon**2 * others - epsilon**2 * weight / (squared_momentum - kaon))
            - epsilon**2 * residues[meson])


# ------------------------------------------------------------------ symbolic checks


def symbolic_checks():
    """Items 1 and 3 with symbolic inputs."""
    s = SQUARED_MOMENTUM
    quark_spurion_diagonal = sp.symbols("x_u x_d x_s", positive=True)
    anomaly_mass_squared = sp.Symbol("m0sq", positive=True)
    gluon_axial_combination = sp.Symbol("c_GA", real=True)
    difference = _hermitian_charge_conserving("cAP")
    numerator = mixing_numerator_polynomial(quark_spurion_diagonal, [difference[index, index] for index in range(3)],
                                            anomaly_mass_squared, gluon_axial_combination)
    degree = sp.Poly(numerator, s).degree()

    roots = sp.symbols("r_1 r_2 r_3", real=True)
    coefficients = sp.symbols("a_0 a_1 a_2", real=True)
    generic_numerator = coefficients[0] + coefficients[1] * s + coefficients[2] * s**2
    product = sp.prod([s - root for root in roots])
    residues = [2 * generic_numerator.subs(s, root) / sp.prod([root - other for other in roots if other != root])
                for root in roots]
    partial_fractions = _is_zero(2 * generic_numerator / product
                                 - sum(residue / (s - root) for residue, root in zip(residues, roots)))

    epsilon, alp_entry, kaon, weight = sp.symbols("epsilon M m_K0sq W", real=True)
    polynomial = sp.expand((s - kaon) * ((s - kaon) * ((s - alp_entry) * product - 2 * epsilon**2 * generic_numerator)
                                         - epsilon**2 * weight * product))
    equations = []
    for meson in range(3):
        others = sum(residues[other] / (s - roots[other]) for other in range(3) if other != meson)
        equation = ((s - roots[meson]) * (s - alp_entry - epsilon**2 * others - epsilon**2 * weight / (s - kaon))
                    - epsilon**2 * residues[meson])
        cofactor = (s - kaon) * sp.prod([s - roots[other] for other in range(3) if other != meson])
        equations.append(_is_zero(equation * cofactor - polynomial / (s - kaon)))
    return {
        "1. N(s) has degree <= 2 in s": degree <= 2,
        "1. 2 N(s)/P(s) = sum_n R_n/(s - m_n0^2), R_n = 2 N(m_n0^2)/prod_{n' != n}(m_n0^2 - m_n'0^2)": partial_fractions,
        "3. per-meson equation x (s - m_K0,0^2) prod_{n' != n}(s - m_n'0^2) = det D(s)/(s - m_K0,0^2), n = 1, 2, 3": all(
            equations),
    }


# ------------------------------------------------------------------ rational checks


def residue_trace_checks(structure, digits=HIGH_PRECISION_DIGITS):
    """Item 2: R_n = tr(lambda_n 𝒞)^2 for pi^0, eta, eta' (50 digits)."""
    mpmath.mp.dps = digits + 10
    sources = rational_test_sources(EPSILON_VALUES[0], illustrative_lo_inputs(structure))
    inputs = explicit_inputs(sources)
    unmixed = unmixed_flavour_diagonal_masses(inputs, digits)
    residues = explicit_residues(inputs, unmixed, digits)
    tadpole = to_high_precision(tadpole_source_from_sources(**_block_inputs(sources)), digits)
    largest = mpmath.mpf(0)
    eigenbasis = mass_eigenstate_generators(sources, digits)
    for value, residue in zip(unmixed, residues):
        value_generator = [generator for eigenvalue, generator in eigenbasis
                           if abs(eigenvalue - value) < mpmath.mpf("1e-40")]
        trace_squared = high_precision_trace(value_generator[0] * tadpole).real**2
        largest = max(largest, abs(trace_squared - residue) + (len(value_generator) != 1))
    return largest


def rational_checks(structure, epsilon, digits=HIGH_PRECISION_DIGITS):
    """Item 4 for one source structure and one epsilon."""
    mpmath.mp.dps = digits + 10
    s = SQUARED_MOMENTUM
    sources = rational_test_sources(epsilon, illustrative_lo_inputs(structure))
    inputs = explicit_inputs(sources)
    kinetic, mass_squared = bare_blocks(sources, TEST_CONDENSATE, TEST_SECOND_ORDER_SCALAR_TRACE)
    bare_polynomial = sp.Poly(sp.expand((s * kinetic - mass_squared).det(method="berkowitz")), s)
    explicit_polynomial = sp.Poly(explicit_secular_polynomial(
        inputs["quark_spurion_diagonal"], inputs["difference_source"], inputs["anomaly_mass_squared"],
        inputs["gluon_axial_combination"], inputs["epsilon"], inputs["alp_entry"]), s)
    intervals = explicit_polynomial.intervals(eps=sp.Rational(1, 10**50))
    exact_roots = sorted((lower + upper) / 2 for (lower, upper), multiplicity in intervals for _ in range(multiplicity))
    spectrum = solve_spectrum_numerical(np.array(kinetic.evalf(30), dtype=float),
                                        np.array(mass_squared.evalf(30), dtype=float),
                                        field_labels=NEUTRAL_FIELD_LABELS)
    solver_difference = float(np.max(np.abs(np.array([float(root) for root in exact_roots])
                                            - np.sort(spectrum.pole_squared_masses))))
    unmixed = unmixed_flavour_diagonal_masses(inputs, digits)
    residues = explicit_residues(inputs, unmixed, digits)
    kaon = kaon_mass_squared(inputs)
    largest_residual = mpmath.mpf(0)
    for root in exact_roots:
        if abs(root - kaon) < sp.Rational(1, 10**40):
            continue
        point = mpmath.mpf(str(sp.N(root, digits + 10)))
        for meson in range(3):
            largest_residual = max(largest_residual,
                                   abs(meson_equation_residual(point, meson, unmixed, residues, inputs)))
    scale = max(abs(float(root)) for root in exact_roots)
    return {
        "polynomial": explicit_polynomial == bare_polynomial,
        "roots": len(exact_roots) == 6,
        "solver_difference": solver_difference,
        "solver": solver_difference <= 1e-10 * max(1.0, scale),
        "largest_residual": largest_residual,
        "residual": largest_residual < RESIDUAL_TOLERANCE,
        "exact_roots": exact_roots,
    }


def derive():
    """Run all checks; return checks and per-case numbers."""
    checks = dict(symbolic_checks())
    numbers = {}
    for structure in SOURCE_STRUCTURES:
        largest = residue_trace_checks(structure)
        checks[f"2. ({structure}) R_n = tr(lambda_n 𝒞)^2 for pi^0, eta, eta' (50 digits)"] = largest < mpmath.mpf(
            "1e-40")
        numbers[(structure, "residue")] = largest
        for epsilon in EPSILON_VALUES:
            result = rational_checks(structure, epsilon)
            label = f"({structure}, epsilon = {epsilon})"
            checks[f"4. {label} explicit polynomial = det(s K - M^2) (exact)"] = result["polynomial"]
            checks[f"4. {label} six real roots, solver agrees (rel. 1e-10)"] = result["roots"] and result["solver"]
            checks[f"4. {label} per-meson equations vanish at the roots (< 1e-30)"] = result["residual"]
            numbers[(structure, epsilon)] = result
    return {"checks": checks, "numbers": numbers}


def main():
    result = derive()
    for key, value in result["numbers"].items():
        if key[1] == "residue":
            print(f"{key[0]}: largest |R_n - tr(lambda_n 𝒞)^2| = {mpmath.nstr(value, 3)}")
        else:
            print(f"{key[0]}, epsilon = {key[1]}: roots {[sp.N(root, 10) for root in value['exact_roots']]}")
            print(f"  largest |solver - root| = {value['solver_difference']:.2e} GeV^2,"
                  f"  largest per-meson residual = {mpmath.nstr(value['largest_residual'], 3)}")
    for label, passed in result["checks"].items():
        print(f"[{'ok' if passed else 'FAIL'}] {label}")
    passed = all(result["checks"].values())
    print("All checks pass" if passed else "Some checks FAIL")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
