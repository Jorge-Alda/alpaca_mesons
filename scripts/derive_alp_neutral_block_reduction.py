#!/usr/bin/env python3
"""How far the LO neutral system (P_3, P_8, P_0, P_6, P_7, a) can be reduced.

Run:   .venv/bin/python scripts/derive_alp_neutral_block_reduction.py
Tests: .venv/bin/python -m unittest discover -s scripts -p test_derive_alp_neutral_block_reduction.py

Conventions as in scripts/derive_alp_lo_neutral_blocks.py (Alda:2026nrd):
K and M^2 on (P_3, P_8, P_0, P_6, P_7, a) from the closed forms, which
scripts/derive_alp_charged_meson_decoupling.py shows to be equal to the LO
Lagrangian Hessian for symbolic sources commuting with Q. The sources
c_q^A, c_q^P, chi_0 are Hermitian with zero u-d and u-s entries; their d-s
entries are written c_ds = Re c_ds + i Im c_ds. Positions in the 6x6:
0 = P_3, 1 = P_8, 2 = P_0, 3 = P_6, 4 = P_7, 5 = a. A flavour rotation V
(vector-like, so it does not change the anomaly term) acts by X -> V X V^dagger
on c_q^A, c_q^P and chi_0; V must commute with Q, so it acts on the d-s plane.

Results checked here (all symbolic, or exact rationals for examples):
 (a) With diagonal chi_0 = diag(chi_u, chi_d, chi_s): K_PP = 1, the M^2 entries
     between (P_3, P_8, P_0) and (P_6, P_7) vanish, and the (P_6, P_7) block of
     M^2 is m_K0^2 times the unit matrix with m_K0^2 = (chi_d + chi_s)/2 = B_0 (m_d + m_s).
     The flavour-diagonal mesons and the K0 doublet mix only through a. With a
     d-s entry in chi_0 the (P_3, P_8, P_0)-(P_6, P_7) entries are nonzero; the
     (P_6, P_7) block stays proportional to the unit matrix (no
     K0-K0bar-type entry is generated at LO).
 (b) The ALP couples to the K0 doublet through the kinetic vector
     k = (K_{P_6 a}, K_{P_7 a}) and the mass vector mu = (M^2_{P_6 a}, M^2_{P_7 a}),
     with k_6 mu_7 - k_7 mu_6 = 2 (f_pi/f_a)^2 (chi_d + chi_s) Im(c^A_ds conj(c^P_ds)).
     A K0 combination n with n.k = n.mu = 0 (decoupled field) exists if and only
     if this vanishes. Real d-s entries (CP conservation in the convention
     Phi -> -Phi^T, under which P_7 is CP-even and P_3, P_8, P_0, P_6, a are
     CP-odd) give k_7 = mu_7 = 0: P_7 decouples. For aligned phases the
     s-quark rephasing V = diag(1, 1, e^{i phi}), which rotates (P_6, P_7) by a
     constant angle and keeps chi_0 diagonal, makes both entries real.
     Independently of phases, s = m_K0^2 is always an exact pole: the vector
     with components (-w_7, w_6) on (P_6, P_7), w_j = m_K0^2 K_ja - M^2_ja, and
     zero elsewhere is a null vector of D(m_K0^2). Its ALP component is zero,
     so it contributes nothing to the a-a propagator.
 (c) A d-s rotation removes all entries between (P_3, P_8, P_0) and
     (P_6, P_7, coupled through a) if and only if, in the rotated basis, the
     d-s entries of chi_0, c_q^A, c_q^P all vanish, i.e. if and only if the d-s
     blocks of the three matrices commute pairwise. For m_d != m_s the only
     rotations keeping chi_0 diagonal are phases. A rotation diagonalizing
     the d-s block of c_q^A (when it does not commute with chi_0) produces a
     d-s entry in chi_0 and hence direct M^2 entries between (P_3, P_8, P_0)
     and (P_6, P_7).

Units: f_pi, f_a in GeV; chi, m_0^2, m_a0^2, M^2, s in GeV^2. Assumptions as
in the LO blocks script (LO, canonical bare ALP term, CP-conserving stationary
vacuum, Hermitian sources commuting with Q, no NLO, no loops).
"""

from __future__ import annotations

import sympy as sp

from derive_alp_charged_meson_decoupling import symbolic_hermitian, symbolic_scalars
from derive_alp_lo_neutral_blocks import closed_form_blocks, illustrative_lo_inputs
from symbolic_utils import commutator

POSITION = {"P_3": 0, "P_8": 1, "P_0": 2, "P_6": 3, "P_7": 4, "a": 5}
FLAVOR_DIAGONAL_POSITIONS = [0, 1, 2]
NEUTRAL_KAON_POSITIONS = [3, 4]
ALP_POSITION = 5
SQUARED_MOMENTUM = sp.Symbol("s")
# CP in the convention Phi -> -Phi^T, a -> -a: sign of each field in the 6x6 basis.
CP_SIGNS = sp.diag(-1, -1, -1, -1, 1, -1)


def _vanishes(matrix):
    return matrix.applyfunc(lambda entry: sp.simplify(sp.expand(entry))) == sp.zeros(*matrix.shape)


def diagonal_quark_spurion():
    return sp.diag(*sp.symbols("chi_uu chi_dd chi_ss", real=True))


def symbolic_neutral_inputs(diagonal_spurion=True):
    """Symbolic Q-commuting c_q^A, c_q^P; chi_0 diagonal or with a d-s entry."""
    return {
        "axial_source": symbolic_hermitian("cA"),
        "pseudoscalar_source": symbolic_hermitian("cP"),
        "quark_spurion": diagonal_quark_spurion() if diagonal_spurion else symbolic_hermitian("chi"),
        **symbolic_scalars(),
    }


def ds_entry(matrix):
    return matrix[1, 2]


def rotate_sources(inputs, rotation):
    """Apply X -> V X V^dagger to c_q^A, c_q^P and chi_0 with V = rotation (3x3)."""
    rotated = dict(inputs)
    for name in ("axial_source", "pseudoscalar_source", "quark_spurion"):
        rotated[name] = (rotation * inputs[name] * rotation.H).applyfunc(sp.expand)
    return rotated


def ds_rotation(cosine, sine):
    """Real rotation in the d-s plane (commutes with Q)."""
    return sp.Matrix([[1, 0, 0], [0, cosine, sine], [0, -sine, cosine]])


# ------------------------------------------------------------------- (a)


def part_a():
    """Structure of the 6x6 for diagonal and for d-s off-diagonal chi_0."""
    kinetic, mass_squared = closed_form_blocks(**symbolic_neutral_inputs(diagonal_spurion=True))
    chi = diagonal_quark_spurion()
    kaon_mass_squared = (chi[1, 1] + chi[2, 2]) / 2
    off_kinetic, off_mass = closed_form_blocks(**symbolic_neutral_inputs(diagonal_spurion=False))
    off_chi = symbolic_hermitian("chi")
    return {
        "checks": {
            "K_PP = 1 (diagonal chi_0)": kinetic[:5, :5] == sp.eye(5),
            "M^2 between (P_3, P_8, P_0) and (P_6, P_7) vanishes (diagonal chi_0)":
                _vanishes(mass_squared.extract(FLAVOR_DIAGONAL_POSITIONS, NEUTRAL_KAON_POSITIONS)),
            "M^2 on (P_6, P_7) = (chi_d + chi_s)/2 times unit matrix (diagonal chi_0)":
                _vanishes(mass_squared.extract(NEUTRAL_KAON_POSITIONS, NEUTRAL_KAON_POSITIONS) - kaon_mass_squared * sp.eye(2)),
            "M^2 between (P_3, P_8, P_0) and (P_6, P_7) is nonzero with a d-s entry in chi_0":
                not _vanishes(off_mass.extract(FLAVOR_DIAGONAL_POSITIONS, NEUTRAL_KAON_POSITIONS)),
            "with a d-s entry in chi_0, M^2 on (P_6, P_7) is still (chi_dd + chi_ss)/2 times unit matrix":
                _vanishes(off_mass.extract(NEUTRAL_KAON_POSITIONS, NEUTRAL_KAON_POSITIONS)
                          - (off_chi[1, 1] + off_chi[2, 2]) / 2 * sp.eye(2)),
        },
        "kinetic": kinetic,
        "mass_squared": mass_squared,
        "kaon_mass_squared": kaon_mass_squared,
        "off_diagonal_spurion_coupling": off_mass.extract(FLAVOR_DIAGONAL_POSITIONS, NEUTRAL_KAON_POSITIONS).applyfunc(sp.simplify),
    }


# ------------------------------------------------------------------- (b)


def part_b():
    """ALP couplings to the K0 doublet, decoupling condition, CP, rephasing, exact pole."""
    inputs = symbolic_neutral_inputs(diagonal_spurion=True)
    kinetic, mass_squared = closed_form_blocks(**inputs)
    kinetic_vector = kinetic.extract(NEUTRAL_KAON_POSITIONS, [ALP_POSITION])
    mass_vector = mass_squared.extract(NEUTRAL_KAON_POSITIONS, [ALP_POSITION])
    axial_ds, pseudoscalar_ds = ds_entry(inputs["axial_source"]), ds_entry(inputs["pseudoscalar_source"])
    ratio = inputs["pion_decay_constant"] / inputs["alp_decay_constant"]
    chi = inputs["quark_spurion"]
    cross = sp.expand(kinetic_vector[0] * mass_vector[1] - kinetic_vector[1] * mass_vector[0])
    expected_cross = 2 * ratio**2 * (chi[1, 1] + chi[2, 2]) * sp.im(sp.expand(axial_ds * sp.conjugate(pseudoscalar_ds)))
    # A K0 combination n = (n_6, n_7) decoupled from a needs n.k = 0 and n.mu = 0.
    combination = sp.Matrix(sp.symbols("n_6 n_7"))
    decoupled_solutions = sp.linsolve([(combination.T * kinetic_vector)[0], (combination.T * mass_vector)[0]], list(combination))

    # CP: real sources are invariant; complex d-s entries are not.
    real_entries = {sp.Symbol(f"{name}_ds_im", real=True): 0 for name in ("cA", "cP")}
    real_kinetic, real_mass = kinetic.subs(real_entries), mass_squared.subs(real_entries)
    p7_offdiagonal = [index for index in range(6) if index != POSITION["P_7"]]

    # Rephasing V = diag(1, 1, w), |w| = 1, w = (3 + 4i)/5, acts as an O(2) rotation on (P_6, P_7).
    phase = (3 + 4 * sp.I) / 5
    rephasing = sp.diag(1, 1, phase)
    rephased_kinetic, rephased_mass = closed_form_blocks(**rotate_sources(inputs, rephasing))
    field_rotation = None
    for sine in (sp.Rational(4, 5), sp.Rational(-4, 5)):
        candidate = sp.eye(6)
        candidate[3:5, 3:5] = sp.Matrix([[sp.Rational(3, 5), -sine], [sine, sp.Rational(3, 5)]])
        if _vanishes(rephased_kinetic - candidate.T * kinetic * candidate) and _vanishes(rephased_mass - candidate.T * mass_squared * candidate):
            field_rotation = candidate

    # Aligned phases: c^A_ds and c^P_ds both proportional to w = (3+4i)/5; rephasing by w makes them real.
    aligned = illustrative_lo_inputs("cp_conserving")
    for name in ("axial_source", "pseudoscalar_source"):
        aligned[name] = aligned[name].copy()
        aligned[name][1, 2] = aligned[name][1, 2] * phase
        aligned[name][2, 1] = sp.conjugate(aligned[name][1, 2])
    aligned_kinetic, aligned_mass = closed_form_blocks(**aligned)
    realigned_kinetic, realigned_mass = closed_form_blocks(**rotate_sources(aligned, rephasing))

    # Exact pole at s = m_K0^2 with zero ALP component, for all symbolic sources.
    kaon_mass_squared = (chi[1, 1] + chi[2, 2]) / 2
    alp_row = kaon_mass_squared * kinetic_vector - mass_vector
    null_vector = sp.zeros(6, 1)
    null_vector[3], null_vector[4] = -alp_row[1], alp_row[0]
    residual = (kaon_mass_squared * kinetic - mass_squared) * null_vector

    def p7_decoupled(kinetic_matrix, mass_matrix):
        return all(sp.simplify(matrix[POSITION["P_7"], index]) == 0
                   for matrix in (kinetic_matrix, mass_matrix) for index in p7_offdiagonal)

    return {
        "checks": {
            "k_6 mu_7 - k_7 mu_6 = 2 (f_pi/f_a)^2 (chi_d + chi_s) Im(c^A_ds conj c^P_ds)": sp.simplify(cross - expected_cross) == 0,
            "a decoupled K0 combination n != 0 exists only if k, mu are parallel":
                decoupled_solutions == sp.FiniteSet((0, 0)),
            "real d-s entries: P_7 decouples (5x5 plus P_7)": p7_decoupled(real_kinetic, real_mass),
            "real d-s entries: K, M^2 invariant under CP (Phi -> -Phi^T, a -> -a)":
                _vanishes(CP_SIGNS * real_kinetic * CP_SIGNS - real_kinetic) and _vanishes(CP_SIGNS * real_mass * CP_SIGNS - real_mass),
            "complex d-s entries: K not CP invariant": not _vanishes(CP_SIGNS * kinetic * CP_SIGNS - kinetic),
            "s-quark rephasing = constant O(2) rotation of (P_6, P_7), other fields unchanged": field_rotation is not None,
            "aligned complex phases: P_7 does not decouple before rephasing": not p7_decoupled(aligned_kinetic, aligned_mass),
            "aligned complex phases: P_7 decouples after rephasing": p7_decoupled(realigned_kinetic, realigned_mass),
            "D(m_K0^2) annihilates (0, 0, 0, -w_7, w_6, 0) (exact pole, zero ALP component)": _vanishes(residual),
        },
        "kinetic_vector": kinetic_vector.applyfunc(sp.simplify),
        "mass_vector": mass_vector.applyfunc(sp.simplify),
        "cross": sp.factor(cross),
        "null_vector": null_vector.applyfunc(sp.simplify),
    }


# ------------------------------------------------------------------- (c)


def ds_blocks_commute(*matrices):
    """Pairwise commutation of the d-s blocks (rows/columns 1, 2) of 3x3 matrices."""
    blocks = [matrix[1:3, 1:3] for matrix in matrices]
    return all(commutator(left, right).applyfunc(sp.expand) == sp.zeros(2)
               for position, left in enumerate(blocks) for right in blocks[position + 1:])


def flavor_diagonal_kaon_entries(kinetic, mass_squared):
    """All entries that couple (P_3, P_8, P_0) to (P_6, P_7) directly or through a."""
    return list(mass_squared.extract(FLAVOR_DIAGONAL_POSITIONS, NEUTRAL_KAON_POSITIONS)) + \
        list(kinetic.extract(NEUTRAL_KAON_POSITIONS, [ALP_POSITION])) + \
        list(mass_squared.extract(NEUTRAL_KAON_POSITIONS, [ALP_POSITION]))


def part_c():
    """Condition for separating the flavour-diagonal part from the K0 part by a d-s rotation."""
    inputs = symbolic_neutral_inputs(diagonal_spurion=False)
    kinetic, mass_squared = closed_form_blocks(**inputs)
    ds_symbols = [sp.Symbol(f"{name}_ds_{part}", real=True) for name in ("chi", "cA", "cP") for part in ("re", "im")]
    conditions = [sp.expand(entry) for entry in flavor_diagonal_kaon_entries(kinetic, mass_squared)]
    # sympy.solve treats chi_dd + chi_ss as generic (nonzero); m_K0^2 > 0 excludes the special case.
    solutions = sp.solve(conditions, ds_symbols, dict=True)
    zero_entries = {symbol: 0 for symbol in ds_symbols}

    # Rotations keeping diag(chi_d, chi_s) diagonal, chi_d != chi_s.
    rotation_entries = sp.symbols("v_11 v_12 v_21 v_22")
    general_rotation = sp.Matrix(2, 2, rotation_entries)
    chi_d, chi_s = sp.symbols("chi_dd chi_ss", real=True)
    keeping = sp.solve(list(commutator(general_rotation, sp.diag(chi_d, chi_s))), [rotation_entries[1], rotation_entries[2]], dict=True)

    rational = sp.Rational
    rotation = ds_rotation(rational(3, 5), rational(4, 5))
    # Commuting example: all sources diagonal, then rotated into a basis with d-s entries.
    diagonal = illustrative_lo_inputs("flavor_diagonal")
    commuting = rotate_sources(diagonal, rotation)
    commuting_kinetic, commuting_mass = closed_form_blocks(**commuting)
    back_kinetic, back_mass = closed_form_blocks(**rotate_sources(commuting, rotation.T))
    # Non-commuting example: c_q^A d-s block = O^T diag(1/2, -1/3) O, chi_0 diagonal with m_d != m_s.
    noncommuting = illustrative_lo_inputs("flavor_diagonal")
    noncommuting["axial_source"] = rotation.T * noncommuting["axial_source"] * rotation
    noncommuting_kinetic, noncommuting_mass = closed_form_blocks(**noncommuting)
    diagonalized = rotate_sources(noncommuting, rotation)
    diagonalized_kinetic, diagonalized_mass = closed_form_blocks(**diagonalized)

    def separated(kinetic_matrix, mass_matrix):
        return all(sp.simplify(entry) == 0 for entry in flavor_diagonal_kaon_entries(kinetic_matrix, mass_matrix))

    return {
        "checks": {
            "separation conditions hold when the d-s entries of chi_0, c_q^A, c_q^P vanish":
                all(sp.simplify(condition.subs(zero_entries)) == 0 for condition in conditions),
            "separation conditions force all d-s entries to vanish (for chi_d + chi_s != 0)": solutions == [zero_entries],
            "for chi_d != chi_s, rotations keeping chi_0 diagonal are diagonal (phases)":
                keeping == [{rotation_entries[1]: 0, rotation_entries[2]: 0}],
            "commuting example: d-s blocks commute": ds_blocks_commute(
                commuting["axial_source"], commuting["pseudoscalar_source"], commuting["quark_spurion"]),
            "commuting example: entries present in the rotated basis": not separated(commuting_kinetic, commuting_mass),
            "commuting example: separated after rotating back (4x4 plus K0 pair)": separated(back_kinetic, back_mass),
            "non-commuting example: d-s blocks of c_q^A and chi_0 do not commute": not ds_blocks_commute(
                noncommuting["axial_source"], noncommuting["quark_spurion"]),
            "non-commuting example: K_{P_6 a} or K_{P_7 a} nonzero in the mass basis": not _vanishes(
                noncommuting_kinetic.extract(NEUTRAL_KAON_POSITIONS, [ALP_POSITION])),
            "non-commuting example: diagonalizing c_q^A removes K_{P_6 a}, K_{P_7 a}": _vanishes(
                diagonalized_kinetic.extract(NEUTRAL_KAON_POSITIONS, [ALP_POSITION])),
            "non-commuting example: ... but produces a d-s entry in chi_0": ds_entry(diagonalized["quark_spurion"]) != 0,
            "non-commuting example: ... and direct M^2 entries (P_3, P_8, P_0)-(P_6, P_7)": not _vanishes(
                diagonalized_mass.extract(FLAVOR_DIAGONAL_POSITIONS, NEUTRAL_KAON_POSITIONS)),
            "non-commuting example: ... while the (P_6, P_7) block of M^2 stays proportional to the unit matrix":
                diagonalized_mass[3, 3] == diagonalized_mass[4, 4] and diagonalized_mass[3, 4] == 0,
        },
        "conditions": conditions,
        "solutions": solutions,
        "rotated_spurion": diagonalized["quark_spurion"],
        "rotated_coupling": diagonalized_mass.extract(FLAVOR_DIAGONAL_POSITIONS, NEUTRAL_KAON_POSITIONS),
    }


CONCLUSION = """Minimal LO block sizes in the quark-mass basis (chi_0 diagonal; any Hermitian
chi_0 commuting with Q is brought there by a d-s rotation), m_d != m_s:
  * generic neutral sources, Im(c^A_ds conj c^P_ds) != 0: 6x6 on (P_3, P_8, P_0, P_6, P_7, a).
    (P_3, P_8, P_0) and (P_6, P_7) have no direct K or M^2 entries; they mix only through a.
    The (P_6, P_7) block is m_K0^2 1, m_K0^2 = (chi_d + chi_s)/2. One pole is always exactly
    at m_K0^2 with pole vector in the (P_6, P_7) plane and zero ALP component; the other five
    poles are the zeros of det D(s)/(s - m_K0^2).
  * Im(c^A_ds conj c^P_ds) = 0 (includes real d-s entries, CP conservation, and c^A_ds = 0 or
    c^P_ds = 0): after the s-quark rephasing, 5x5 on (P_3, P_8, P_0, P_6', a) plus a decoupled
    K0 combination P_7' with m^2 = m_K0^2.
  * c^A_ds = c^P_ds = 0 (flavour-diagonal sources): 4x4 on (P_3, P_8, P_0, a) plus the
    degenerate K0 pair (P_6, P_7) at m_K0^2.
  A d-s rotation cannot do better unless the d-s blocks of c_q^A, c_q^P, chi_0 commute; then the
  rotated sources are flavour diagonal and the last case applies."""


def derive():
    return {"a": part_a(), "b": part_b(), "c": part_c(), "conclusion": CONCLUSION}


def main():
    result = derive()
    for part in ("a", "b", "c"):
        print(f"({part})")
        for label, passed in result[part]["checks"].items():
            print(f"  [{'pass' if passed else 'FAIL'}] {label}")
    print("k = (K_{P_6 a}, K_{P_7 a}) =", list(result["b"]["kinetic_vector"]))
    print("mu = (M^2_{P_6 a}, M^2_{P_7 a}) =", list(result["b"]["mass_vector"]))
    print("k_6 mu_7 - k_7 mu_6 =", result["b"]["cross"])
    print("M^2 between (P_3, P_8, P_0) and (P_6, P_7) with a d-s entry in chi_0:", result["a"]["off_diagonal_spurion_coupling"])
    print("non-commuting example, chi_0 after diagonalizing c_q^A:", result["c"]["rotated_spurion"])
    print(result["conclusion"])


if __name__ == "__main__":
    main()
