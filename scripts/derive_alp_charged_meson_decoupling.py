#!/usr/bin/env python3
"""Decoupling of the charged mesons from the neutral mesons and the ALP at LO.

Run:   .venv/bin/python scripts/derive_alp_charged_meson_decoupling.py
Tests: .venv/bin/python -m unittest discover -s scripts -p test_derive_alp_charged_meson_decoupling.py

Conventions as in scripts/derive_alp_lo_neutral_blocks.py (Alda:2026nrd):
U(3) generators lambda_0..lambda_8 with tr(lambda_i lambda_j) = 2 delta_ij,
L_2 = 1/2 dPhi^T K dPhi - 1/2 Phi^T M^2 Phi on (P_0, ..., P_8, a), quark
charge matrix Q = diag(2/3, -1/3, -1/3) in units of e.

Statement. Let c_q^A, c_q^P and chi_0 be Hermitian 3x3 matrices that commute
with Q, i.e. their u-d and u-s entries vanish; diagonal entries and the d-s
off-diagonal entry are arbitrary (chi_0 is not required to be diagonal).
Then every entry of K and M^2 between a charged meson (P_1, P_2 = pi^+-,
P_4, P_5 = K^+-) and a neutral field (P_0, P_3, P_6, P_7, P_8 or a) vanishes,
including the m_0^2 and c_G terms. Without Q-commutation this fails.

Algebraic reason (verified, not assumed, below). The neutral generators are
those commuting with Q: span(lambda_0, lambda_3, lambda_6, lambda_7, lambda_8).
The charged generators lie in the image of ad_Q = [Q, .]. For Y commuting
with Q, tr([Q, Z] Y) = tr(Z [Y, Q]) = 0, so tr(lambda_c Y) = 0. The K_ia,
M^2_ia entries are traces of lambda_c with c_q^A or {c_q^P, chi_0}, which
commute with Q. Since ad_Q is a derivation, {[Q, Z], lambda_n} = [Q, {Z, lambda_n}]
for neutral lambda_n, so {lambda_c, lambda_n} also lies in the image of ad_Q
and M^2_cn = tr({lambda_c, lambda_n} chi_0)/4 vanishes. The m_0^2 and c_G
terms contain tr(lambda_c) = 0.

Method. (1) Classify the generators by ad_Q. (2) Build the 10x10 Hessian of
the explicit LO Lagrangian (blocks_from_explicit_lagrangian) with fully
symbolic Q-commuting sources and symbolic f_pi, f_a, m_0^2, m_a0^2, c_G;
check that it equals the closed forms and that the charged-neutral entries
vanish identically. (3) Repeat with fully general symbolic Hermitian sources
and list the charged-neutral entries: they are linear in the u-d and u-s
entries of c_q^A, {c_q^P, chi_0} and chi_0, so they vanish if and only if
those entries vanish. (4) Explicit counterexamples with a single u-d or u-s
entry.

Units: f_pi, f_a in GeV; chi_0, m_0^2, m_a0^2 in GeV^2; couplings
dimensionless. Assumptions: those of the LO blocks script (LO chiral
Lagrangian, canonical bare ALP term, CP-conserving stationary vacuum,
Hermitian sources, quark spurion linear in a, no NLO, no loops).
"""

from __future__ import annotations

import sympy as sp

from derive_alp_lo_neutral_blocks import (
    ALL_GENERATOR_INDICES,
    CHARGED_GENERATOR_INDICES,
    NEUTRAL_GENERATOR_INDICES,
    blocks_from_explicit_lagrangian,
    charged_neutral_entries,
    closed_form_blocks,
)
from symbolic_utils import anticommutator, commutator, u3_generators

QUARK_CHARGE = sp.diag(sp.Rational(2, 3), sp.Rational(-1, 3), sp.Rational(-1, 3))
FLAVORS = ("u", "d", "s")


def symbolic_hermitian(name, charge_conserving=True):
    """Hermitian 3x3 matrix of real symbols; u-d and u-s entries zero if charge_conserving.

    Diagonal entries {name}_uu etc.; off-diagonal entry (j, k), j < k, is
    {name}_jk_re + i {name}_jk_im.
    """
    matrix = sp.zeros(3)
    for row in range(3):
        matrix[row, row] = sp.Symbol(f"{name}_{FLAVORS[row]}{FLAVORS[row]}", real=True)
        for column in range(row + 1, 3):
            if charge_conserving and row == 0:
                continue
            pair = FLAVORS[row] + FLAVORS[column]
            entry = sp.Symbol(f"{name}_{pair}_re", real=True) + sp.I * sp.Symbol(f"{name}_{pair}_im", real=True)
            matrix[row, column], matrix[column, row] = entry, sp.conjugate(entry)
    return matrix


def symbolic_scalars():
    """Symbolic f_pi, f_a > 0 and real m_0^2, m_a0^2, c_G (keyword names of the blocks functions)."""
    return {
        "pion_decay_constant": sp.Symbol("f_pi", positive=True),
        "alp_decay_constant": sp.Symbol("f_a", positive=True),
        "anomaly_mass_squared": sp.Symbol("m_0_sq", real=True),
        "bare_alp_mass_squared": sp.Symbol("m_a0_sq", real=True),
        "gluon_coupling": sp.Symbol("c_G", real=True),
    }


def symbolic_inputs(charge_conserving=True):
    """Symbolic sources (c_q^A, c_q^P, chi_0) and scalars for the blocks functions."""
    return {
        "axial_source": symbolic_hermitian("cA", charge_conserving),
        "pseudoscalar_source": symbolic_hermitian("cP", charge_conserving),
        "quark_spurion": symbolic_hermitian("chi", charge_conserving),
        **symbolic_scalars(),
    }


def classify_generators():
    """Indices of generators commuting with Q, and check that the others lie in the image of ad_Q."""
    generators = u3_generators()
    commuting = tuple(index for index, generator in enumerate(generators)
                      if commutator(QUARK_CHARGE, generator) == sp.zeros(3))
    coefficients = sp.symbols("z_0:9")
    trial = sum((coefficient * generator for coefficient, generator in zip(coefficients, generators)), sp.zeros(3))
    image = commutator(QUARK_CHARGE, trial)
    in_image = {}
    for index in ALL_GENERATOR_INDICES:
        if index in commuting:
            continue
        equations = list(image - generators[index])
        in_image[index] = bool(sp.solve(equations, coefficients, dict=True))
    return commuting, in_image


def anticommutator_lemma():
    """{lambda_c, lambda_n} is orthogonal (trace) to every Q-commuting matrix, for charged c, neutral n."""
    generators = u3_generators()
    test_matrix = symbolic_hermitian("Y", charge_conserving=True)
    return all(sp.expand(sp.trace(anticommutator(generators[charged], generators[neutral]) * test_matrix)) == 0
               for charged in CHARGED_GENERATOR_INDICES for neutral in NEUTRAL_GENERATOR_INDICES)


def _difference_vanishes(left, right):
    return (left - right).applyfunc(lambda entry: sp.simplify(sp.expand(entry))) == sp.zeros(*left.shape)


def symbolic_decoupling(charge_conserving=True):
    """10x10 Lagrangian Hessian with symbolic sources; closed-form agreement and charged-neutral entries."""
    inputs = symbolic_inputs(charge_conserving)
    kinetic, mass_squared = blocks_from_explicit_lagrangian(**inputs)
    closed_kinetic, closed_mass = closed_form_blocks(**inputs, generator_indices=ALL_GENERATOR_INDICES)
    return {
        "inputs": inputs,
        "kinetic": kinetic,
        "mass_squared": mass_squared,
        "closed forms agree": _difference_vanishes(kinetic, closed_kinetic) and _difference_vanishes(mass_squared, closed_mass),
        "charged-neutral kinetic entries": charged_neutral_entries(kinetic).applyfunc(sp.expand),
        "charged-neutral mass entries": charged_neutral_entries(mass_squared).applyfunc(sp.expand),
    }


def charge_violating_entries_vanish_only_without_u_entries(general_result):
    """For general sources: set the u-d, u-s entries of c_q^A, c_q^P, chi_0 to zero and recover decoupling."""
    u_entry_symbols = [sp.Symbol(f"{matrix}_{pair}_{part}", real=True)
                       for matrix in ("cA", "cP", "chi") for pair in ("ud", "us") for part in ("re", "im")]
    entries = list(general_result["charged-neutral kinetic entries"]) + list(general_result["charged-neutral mass entries"])
    vanish_after = all(sp.expand(entry.subs({symbol: 0 for symbol in u_entry_symbols})) == 0 for entry in entries)
    return {"vanish when u-d, u-s entries vanish": vanish_after,
            "number of nonzero charged-neutral entries for general sources": sum(1 for entry in entries if entry != 0)}


def counterexamples():
    """Single charge-violating entries (rational) that couple a charged meson to a neutral field."""
    rational = sp.Rational
    base = {
        "axial_source": sp.diag(1, rational(1, 2), rational(-1, 3)),
        "pseudoscalar_source": sp.diag(rational(1, 4), rational(-1, 2), rational(2, 3)),
        "quark_spurion": sp.diag(rational(3, 200), rational(1, 40), rational(47, 100)),
        "gluon_coupling": rational(1, 2),
        "pion_decay_constant": rational(23, 250),
        "alp_decay_constant": rational(1, 5),
        "anomaly_mass_squared": rational(4, 5),
        "bare_alp_mass_squared": rational(3, 10),
    }
    up_down = sp.zeros(3)
    up_down[0, 1] = up_down[1, 0] = rational(1, 3)
    up_strange = sp.zeros(3)
    up_strange[0, 2] = up_strange[2, 0] = rational(1, 4)
    cases = {
        "u-d entry 1/3 in c_q^A: K_{P_1 a}": ({"axial_source": base["axial_source"] + up_down}, "kinetic", 1, 9),
        "u-s entry 1/4 in c_q^P: M^2_{P_4 a}": ({"pseudoscalar_source": base["pseudoscalar_source"] + up_strange}, "mass", 4, 9),
        # M^2_{P_1 P_3} stays zero here; the u-d entry of chi_0 couples P_1 to P_8, P_0 and a.
        "u-d entry 1/100 in chi_0: M^2_{P_1 P_8}": ({"quark_spurion": base["quark_spurion"] + up_down * 3 / 100}, "mass", 1, 8),
    }
    values = {}
    for label, (change, which, row, column) in cases.items():
        kinetic, mass_squared = blocks_from_explicit_lagrangian(**{**base, **change})
        values[label] = (kinetic if which == "kinetic" else mass_squared)[row, column]
    unchanged_kinetic, unchanged_mass = blocks_from_explicit_lagrangian(**base)
    values["reference without the entry: all charged-neutral entries"] = (
        charged_neutral_entries(unchanged_kinetic) == sp.zeros(4, 6) and charged_neutral_entries(unchanged_mass) == sp.zeros(4, 6))
    return values


def derive():
    """Run the classification, the symbolic proof and the counterexamples."""
    commuting, in_image = classify_generators()
    conserving = symbolic_decoupling(charge_conserving=True)
    general = symbolic_decoupling(charge_conserving=False)
    checks = {
        "generators commuting with Q are lambda_0, lambda_3, lambda_6, lambda_7, lambda_8":
            commuting == tuple(sorted(NEUTRAL_GENERATOR_INDICES)),
        "lambda_1, lambda_2, lambda_4, lambda_5 lie in the image of ad_Q":
            tuple(in_image) == CHARGED_GENERATOR_INDICES and all(in_image.values()),
        "tr({lambda_c, lambda_n} Y) = 0 for every Q-commuting Y": anticommutator_lemma(),
        "symbolic Q-commuting sources: 10x10 Lagrangian Hessian equals the closed forms": conserving["closed forms agree"],
        "symbolic Q-commuting sources: charged-neutral K entries vanish":
            conserving["charged-neutral kinetic entries"] == sp.zeros(4, 6),
        "symbolic Q-commuting sources: charged-neutral M^2 entries vanish (incl. m_0^2, c_G)":
            conserving["charged-neutral mass entries"] == sp.zeros(4, 6),
        "general symbolic Hermitian sources: 10x10 Lagrangian Hessian equals the closed forms": general["closed forms agree"],
    }
    general_summary = charge_violating_entries_vanish_only_without_u_entries(general)
    checks["general sources: charged-neutral entries vanish once u-d, u-s entries vanish"] = \
        general_summary["vanish when u-d, u-s entries vanish"]
    examples = counterexamples()
    checks["counterexample entries are nonzero"] = all(
        value != 0 for label, value in examples.items() if not label.startswith("reference"))
    checks["reference (diagonal sources) is decoupled"] = examples["reference without the entry: all charged-neutral entries"]
    return {"checks": checks, "conserving": conserving, "general": general,
            "general_summary": general_summary, "counterexamples": examples}


def main():
    result = derive()
    print("Charged-meson decoupling at LO; sources commuting with Q = diag(2/3, -1/3, -1/3).")
    for label, passed in result["checks"].items():
        print(f"[{'pass' if passed else 'FAIL'}] {label}")
    print("General Hermitian sources, charged-neutral entries (rows P_1, P_2, P_4, P_5; columns P_3, P_8, P_0, P_6, P_7, a):")
    print("K:", result["general"]["charged-neutral kinetic entries"])
    print("M^2:", result["general"]["charged-neutral mass entries"])
    print("Counterexamples:")
    for label, value in result["counterexamples"].items():
        print(f"  {label}: {value}")


if __name__ == "__main__":
    main()
