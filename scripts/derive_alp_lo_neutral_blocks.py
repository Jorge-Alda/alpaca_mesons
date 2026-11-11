#!/usr/bin/env python3
"""LO kinetic matrix K and mass matrix M^2 of the neutral mesons and the ALP.

Run:   .venv/bin/python scripts/derive_alp_lo_neutral_blocks.py
Tests: .venv/bin/python -m unittest discover -s scripts -p test_derive_alp_lo_neutral_blocks.py

Conventions (Alda:2026nrd; memory/alp/shared/finite-alp-meson-mixing.md):
N_f = 3, meson matrix sum_i P_i lambda_i with Hermitian U(3) generators
normalized by tr(lambda_i lambda_j) = 2 delta_ij, indexed (lambda_0, ..., lambda_8)
as in scripts/symbolic_utils.py. Quadratic Lagrangian
L_2 = 1/2 dPhi^T K dPhi - 1/2 Phi^T M^2 Phi. The closed forms are

  K_ij = delta_ij,  K_ia = -(f_pi/f_a) tr(c_q^A lambda_i),  K_aa = 1 + sum_i K_ia K_ia,
  M^2_ij = tr({lambda_i, lambda_j} chi_0)/4 + m_0^2 tr(lambda_i) tr(lambda_j)/(2 N_f),
  M^2_ia = -(f_pi/(2 f_a)) tr(lambda_i {c_q^P, chi_0}) + f_pi m_0^2 c_G tr(lambda_i)/(N_f f_a),
  M^2_aa = m_a0^2 + 2 f_pi^2 m_0^2 c_G^2/(N_f f_a^2),

with chi_0 = 2 B_0 m_q. In K_aa the sum runs over all nine generators.

Default basis. The ALP is electrically neutral. When c_q^A, c_q^P and chi_0
commute with the quark charge matrix Q = diag(2/3, -1/3, -1/3) (diagonal
entries plus d-s off-diagonal entries only), the charged generators
lambda_1, lambda_2 (pi^+-) and lambda_4, lambda_5 (K^+-) have no K or M^2
entries with lambda_0, lambda_3, lambda_6, lambda_7, lambda_8 or with a; this
is proven for symbolic sources in scripts/derive_alp_charged_meson_decoupling.py.
The default LO input is therefore the neutral 6x6 system on
    (P_3, P_8, P_0, P_6, P_7, a).
K_aa still contains tr(c_q^A c_q^A) over all generators; for charge-conserving
c_q^A the charged K_ia vanish, so the 6x6 K_aa equals 1 + sum over the five
neutral K_ia (checked below).

The full 10x10 system (P_0, ..., P_8, a) is kept only as a cross-check: it is
compared with (i) the generic Hessian entries of
scripts/derive_alp_lo_finite_mixing.py evaluated with explicit generators and
(ii) the Hessian of the LO Lagrangian expanded explicitly in all ten fields.
The 6x6 blocks must equal the restriction of that 10x10 Hessian, and the
charged-neutral entries of the 10x10 Hessian must vanish.

Assumptions as in the finite-mixing record: LO chiral Lagrangian plus a
canonical bare ALP term, CP-conserving stationary vacuum, Hermitian c_q^A,
c_q^P, quark spurion linear in a, no NLO terms and no loops. No expansion in
f_pi/f_a is made.

Units: f_pi, f_a in GeV; chi_0, m_0^2, m_a0^2 and M^2 in GeV^2; K and the
couplings c_q^A, c_q^P, c_G are dimensionless. The input values below are
exact rationals chosen for testing only; they are not fitted values.
"""

from __future__ import annotations

import sympy as sp

from derive_alp_lo_finite_mixing import derive as derive_lo_hessian_entries
from symbolic_utils import anticommutator, u3_generators

FLAVOR_COUNT = 3
ALL_GENERATOR_INDICES = tuple(range(9))
NEUTRAL_GENERATOR_INDICES = (3, 8, 0, 6, 7)
CHARGED_GENERATOR_INDICES = (1, 2, 4, 5)
FULL_FIELD_LABELS = tuple([f"P_{index}" for index in ALL_GENERATOR_INDICES] + ["a"])
NEUTRAL_FIELD_LABELS = tuple([f"P_{index}" for index in NEUTRAL_GENERATOR_INDICES] + ["a"])
SOURCE_STRUCTURES = ("flavor_diagonal", "cp_conserving", "neutral", "charge_violating")


def field_labels(generator_indices=NEUTRAL_GENERATOR_INDICES):
    """Labels (P_i for the chosen generators, then a) of the rows of K and M^2."""
    return tuple([f"P_{index}" for index in generator_indices] + ["a"])


def illustrative_lo_inputs(source_structure="neutral"):
    """Exact rational test inputs (GeV units); not physical fits. chi_0 is diagonal.

    source_structure:
      "flavor_diagonal"  c_q^A, c_q^P diagonal: only P_3, P_8, P_0 mix with a.
      "cp_conserving"    diagonal plus real d-s entries: P_7 decouples.
      "neutral"          diagonal plus complex d-s entries with
                         Im(c_A,ds conj(c_P,ds)) != 0: generic neutral sources.
      "charge_violating" general Hermitian c_q^A, c_q^P with u-d and u-s
                         entries; used only to show that the charged mesons
                         then do not decouple (10x10 system needed).
    """
    rational = sp.Rational
    if source_structure == "flavor_diagonal":
        axial_source = sp.diag(1, rational(1, 2), rational(-1, 3))
        pseudoscalar_source = sp.diag(rational(1, 4), rational(-1, 2), rational(2, 3))
    elif source_structure in ("cp_conserving", "neutral"):
        axial_entry, pseudoscalar_entry = rational(1, 3), rational(1, 2)
        if source_structure == "neutral":
            axial_entry, pseudoscalar_entry = rational(1, 3) + sp.I / 4, rational(1, 2) - sp.I / 3
        axial_source = sp.Matrix([
            [1, 0, 0],
            [0, rational(1, 2), axial_entry],
            [0, sp.conjugate(axial_entry), rational(-1, 3)],
        ])
        pseudoscalar_source = sp.Matrix([
            [rational(1, 4), 0, 0],
            [0, rational(-1, 2), pseudoscalar_entry],
            [0, sp.conjugate(pseudoscalar_entry), rational(2, 3)],
        ])
    elif source_structure == "charge_violating":
        axial_source = sp.Matrix([
            [1, rational(1, 3) + sp.I / 4, rational(-1, 5)],
            [rational(1, 3) - sp.I / 4, rational(1, 2), sp.I / 7],
            [rational(-1, 5), -sp.I / 7, rational(-1, 3)],
        ])
        pseudoscalar_source = sp.Matrix([
            [rational(1, 4), rational(-1, 6), rational(1, 2) - sp.I / 3],
            [rational(-1, 6), rational(-1, 2), rational(2, 9)],
            [rational(1, 2) + sp.I / 3, rational(2, 9), rational(2, 3)],
        ])
    else:
        raise ValueError(f"unknown source_structure {source_structure!r}; choose from {SOURCE_STRUCTURES}")
    return {
        "axial_source": axial_source,
        "pseudoscalar_source": pseudoscalar_source,
        "quark_spurion": sp.diag(rational(3, 200), rational(1, 40), rational(47, 100)),
        "gluon_coupling": rational(1, 2),
        "pion_decay_constant": rational(23, 250),
        "alp_decay_constant": rational(1, 5),
        "anomaly_mass_squared": rational(4, 5),
        "bare_alp_mass_squared": rational(3, 10),
    }


def _clean(expression):
    return sp.radsimp(sp.expand(expression))


def closed_form_blocks(axial_source, pseudoscalar_source, quark_spurion, gluon_coupling,
                       pion_decay_constant, alp_decay_constant, anomaly_mass_squared,
                       bare_alp_mass_squared, generator_indices=NEUTRAL_GENERATOR_INDICES):
    """K and M^2 from the closed forms on (P_i for i in generator_indices, a).

    Default: the neutral 6x6 basis (P_3, P_8, P_0, P_6, P_7, a). Pass
    ALL_GENERATOR_INDICES for the full 10x10 system. K_aa is the closed form
    1 + 2 (f_pi/f_a)^2 tr(c_q^A c_q^A), i.e. the sum over all nine generators.
    """
    generators = u3_generators()
    selected = [generators[index] for index in generator_indices]
    size = len(selected) + 1
    alp = size - 1
    ratio = pion_decay_constant / alp_decay_constant
    kinetic = sp.eye(size)
    mass_squared = sp.zeros(size)
    for row, generator_row in enumerate(selected):
        kinetic[row, alp] = kinetic[alp, row] = _clean(-ratio * sp.trace(axial_source * generator_row))
        for column, generator_column in enumerate(selected):
            mass_squared[row, column] = _clean(
                sp.trace(anticommutator(generator_row, generator_column) * quark_spurion) / 4
                + anomaly_mass_squared * sp.trace(generator_row) * sp.trace(generator_column) / (2 * FLAVOR_COUNT))
        mass_squared[row, alp] = mass_squared[alp, row] = _clean(
            -ratio / 2 * sp.trace(generator_row * anticommutator(pseudoscalar_source, quark_spurion))
            + ratio * anomaly_mass_squared * gluon_coupling * sp.trace(generator_row) / FLAVOR_COUNT)
    kinetic[alp, alp] = _clean(1 + 2 * ratio**2 * sp.trace(axial_source * axial_source))
    mass_squared[alp, alp] = _clean(
        bare_alp_mass_squared
        + 2 * ratio**2 * anomaly_mass_squared * gluon_coupling**2 / FLAVOR_COUNT)
    return kinetic, mass_squared


def blocks_from_generic_hessian_entries(axial_source, pseudoscalar_source, quark_spurion, gluon_coupling,
                                        pion_decay_constant, alp_decay_constant, anomaly_mass_squared,
                                        bare_alp_mass_squared, generator_indices=ALL_GENERATOR_INDICES):
    """Evaluate the Hessian entries of derive_alp_lo_finite_mixing.derive() for all pairs.

    That script differentiates the LO Lagrangian for two meson coordinates
    P_i, P_j with abstract generators; its bilinear entries (i, j), (i, a),
    (a, a) are evaluated here with the explicit generators, including i = j.
    """
    result = derive_lo_hessian_entries()
    generators = u3_generators()
    selected = [generators[index] for index in generator_indices]
    size = len(selected) + 1
    alp = size - 1
    scalars = {
        result["f_pi"]: pion_decay_constant,
        result["f_a"]: alp_decay_constant,
        result["anomaly_mass"]: anomaly_mass_squared,
        result["bare_alp_mass"]: bare_alp_mass_squared,
        result["gluon_coupling"]: gluon_coupling,
    }
    sources = {
        result["axial_source"]: axial_source,
        result["pseudoscalar_source"]: pseudoscalar_source,
        result["quark_spurion"]: quark_spurion,
    }

    def evaluate(entry, generator_i, generator_j):
        replaced = entry.subs({**sources, result["generator_i"]: generator_i, result["generator_j"]: generator_j})
        replaced = replaced.subs(result["flavor_count"], FLAVOR_COUNT).subs(scalars).doit()
        return _clean(replaced)

    kinetic_meson, kinetic_mixed, kinetic_alp = result["kinetic_entries"]
    mass_meson, mass_mixed, mass_alp = result["mass_entries"]
    kinetic = sp.zeros(size)
    mass_squared = sp.zeros(size)
    for row, generator_row in enumerate(selected):
        for column, generator_column in enumerate(selected):
            kinetic[row, column] = evaluate(kinetic_meson, generator_row, generator_column)
            mass_squared[row, column] = evaluate(mass_meson, generator_row, generator_column)
        kinetic[row, alp] = kinetic[alp, row] = evaluate(kinetic_mixed, generator_row, generator_row)
        mass_squared[row, alp] = mass_squared[alp, row] = evaluate(mass_mixed, generator_row, generator_row)
    kinetic[alp, alp] = evaluate(kinetic_alp, generators[0], generators[0])
    mass_squared[alp, alp] = evaluate(mass_alp, generators[0], generators[0])
    return kinetic, mass_squared


def blocks_from_explicit_lagrangian(axial_source, pseudoscalar_source, quark_spurion, gluon_coupling,
                                    pion_decay_constant, alp_decay_constant, anomaly_mass_squared,
                                    bare_alp_mass_squared, generator_indices=ALL_GENERATOR_INDICES):
    """Hessians of the LO Lagrangian written out in the fields (P_i for i in generator_indices, a).

    Same operators as derive_alp_lo_finite_mixing.py:
      L_kin  = (da)^2/2 + f_pi^2/4 tr(u u),  u = -dPhi/f_pi + 2 da c_q^A/f_a,
      L_mass = -tr(Phi Phi chi_0)/4 + f_pi/2 tr(Phi a {c_q^P, chi_0})/f_a
               - f_pi^2 m_0^2 theta^2/(4 N_f) - m_a0^2 a^2/2,
      theta  = -tr(Phi)/f_pi - 2 c_G a/f_a,  Phi = sum_i P_i lambda_i.
    K = Hessian of L_kin in the field derivatives, M^2 = -Hessian of L_mass.
    Works for symbolic or rational sources. With the default (all nine
    generators) this is the full 10x10 Hessian.
    """
    generators = u3_generators()
    selected = [generators[index] for index in generator_indices]
    count = len(selected)
    fields = sp.symbols(f"P_0:{count}", real=True) + (sp.Symbol("a", real=True),)
    field_derivatives = sp.symbols(f"dP_0:{count}", real=True) + (sp.Symbol("da", real=True),)
    meson_matrix = sum((field * generator for field, generator in zip(fields[:count], selected)), sp.zeros(3))
    derivative_matrix = sum((field * generator for field, generator in zip(field_derivatives[:count], selected)), sp.zeros(3))
    alp, derivative_alp = fields[count], field_derivatives[count]
    vielbein = -derivative_matrix / pion_decay_constant + 2 * derivative_alp * axial_source / alp_decay_constant
    kinetic_lagrangian = derivative_alp**2 / 2 + pion_decay_constant**2 * sp.trace(vielbein * vielbein) / 4
    topological_invariant = -sp.trace(meson_matrix) / pion_decay_constant - 2 * gluon_coupling * alp / alp_decay_constant
    mass_lagrangian = (
        -sp.trace(meson_matrix * meson_matrix * quark_spurion) / 4
        + pion_decay_constant / 2 * sp.trace(meson_matrix * alp * anticommutator(pseudoscalar_source, quark_spurion)) / alp_decay_constant
        - pion_decay_constant**2 * anomaly_mass_squared * topological_invariant**2 / (4 * FLAVOR_COUNT)
        - bare_alp_mass_squared * alp**2 / 2
    )
    kinetic = sp.hessian(sp.expand(kinetic_lagrangian), field_derivatives).applyfunc(_clean)
    mass_squared = (-sp.hessian(sp.expand(mass_lagrangian), fields)).applyfunc(_clean)
    return kinetic, mass_squared


def restrict_to_generators(full_matrix, generator_indices=NEUTRAL_GENERATOR_INDICES):
    """Rows and columns (P_i for i in generator_indices, a) of a 10x10 matrix on (P_0..P_8, a)."""
    positions = list(generator_indices) + [len(ALL_GENERATOR_INDICES)]
    return full_matrix.extract(positions, positions)


def charged_neutral_entries(full_matrix):
    """Entries of a 10x10 matrix between charged (P_1, P_2, P_4, P_5) and (neutral P_i, a)."""
    neutral_positions = list(NEUTRAL_GENERATOR_INDICES) + [len(ALL_GENERATOR_INDICES)]
    return full_matrix.extract(list(CHARGED_GENERATOR_INDICES), neutral_positions)


def derive(source_structure="neutral"):
    """Build the 6x6 neutral blocks; cross-check with the full 10x10 Lagrangian Hessians.

    For charge-conserving structures the charged-neutral entries of the
    10x10 Hessian must vanish; for "charge_violating" they must not.
    """
    inputs = illustrative_lo_inputs(source_structure)
    kinetic, mass_squared = closed_form_blocks(**inputs)
    full_kinetic, full_mass = closed_form_blocks(**inputs, generator_indices=ALL_GENERATOR_INDICES)
    kinetic_generic, mass_generic = blocks_from_generic_hessian_entries(**inputs)
    kinetic_explicit, mass_explicit = blocks_from_explicit_lagrangian(**inputs)
    zero_full, zero_neutral = sp.zeros(10), sp.zeros(6)
    neutral_kinetic_mixing = kinetic[:5, 5]
    charge_conserving = source_structure != "charge_violating"
    decoupled = (charged_neutral_entries(kinetic_explicit) == sp.zeros(4, 6)
                 and charged_neutral_entries(mass_explicit) == sp.zeros(4, 6))
    checks = {
        "6x6 K and M^2 are real": all(entry.is_real for entry in list(kinetic) + list(mass_squared)),
        "6x6 K and M^2 are symmetric": kinetic == kinetic.T and mass_squared == mass_squared.T,
        "10x10 closed-form K matches the generic Hessian entries": (full_kinetic - kinetic_generic).applyfunc(_clean) == zero_full,
        "10x10 closed-form M^2 matches the generic Hessian entries": (full_mass - mass_generic).applyfunc(_clean) == zero_full,
        "10x10 closed-form K matches the explicit Lagrangian Hessian": (full_kinetic - kinetic_explicit).applyfunc(_clean) == zero_full,
        "10x10 closed-form M^2 matches the explicit Lagrangian Hessian": (full_mass - mass_explicit).applyfunc(_clean) == zero_full,
        "6x6 K equals the restriction of the 10x10 Lagrangian Hessian":
            (kinetic - restrict_to_generators(kinetic_explicit)).applyfunc(_clean) == zero_neutral,
        "6x6 M^2 equals the restriction of the 10x10 Lagrangian Hessian":
            (mass_squared - restrict_to_generators(mass_explicit)).applyfunc(_clean) == zero_neutral,
        ("charged-neutral entries of the 10x10 K and M^2 vanish" if charge_conserving
         else "charged-neutral entries of the 10x10 K and M^2 do not all vanish"): decoupled == charge_conserving,
    }
    if charge_conserving:
        checks["6x6 K_aa = 1 + sum over the five neutral K_ia^2"] = _clean(
            kinetic[5, 5] - 1 - (neutral_kinetic_mixing.T * neutral_kinetic_mixing)[0]) == 0
    return {"inputs": inputs, "kinetic": kinetic, "mass_squared": mass_squared,
            "full_kinetic": full_kinetic, "full_mass_squared": full_mass, "checks": checks}


def main():
    for source_structure in SOURCE_STRUCTURES:
        result = derive(source_structure)
        print(f"LO blocks on {NEUTRAL_FIELD_LABELS}, N_f = 3, {source_structure} sources; "
              "illustrative exact rationals (GeV units).")
        print(f"K_ia = {list(result['kinetic'][:5, 5])}")
        print(f"K_aa = {result['kinetic'][5, 5]}")
        print(f"M^2_ia = {list(result['mass_squared'][:5, 5])}")
        print(f"M^2_aa = {result['mass_squared'][5, 5]}")
        for check, passed in result["checks"].items():
            print(f"[{'pass' if passed else 'FAIL'}] {check}")
        print()


if __name__ == "__main__":
    main()
