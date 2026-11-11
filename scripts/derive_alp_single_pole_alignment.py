#!/usr/bin/env python3
"""Coordinates in which the ALP decouples from the mesons at one selected pole.

Run:   .venv/bin/python scripts/derive_alp_single_pole_alignment.py
Tests: .venv/bin/python -m unittest discover -s scripts -p test_derive_alp_single_pole_alignment.py

Conventions as in derive_alp_inverse_propagator_congruence.py and Alda:2026nrd:
Phi = (P_i, a), D(s) = s K - M^2, pole vector v^(k) = (w^{(k)}, 1). The shift
P_i = P_i' + w_i^{(k)}(m_k^2) a', a = a' is the field change Phi = T Phi' with
T = [[1, w], [0, 1]]. Squared masses are in an arbitrary common unit.

Checks:
1. Regular pole (meson block m_k^2 K_PP - M^2_PP invertible): the shift gives
   m_k^2 K'_{ia} - M^2'_{ia} = 0, K'_{ia} and M^2'_{ia} do not vanish
   separately, D'_{ia}(s) does not vanish away from the pole, and the result
   agrees with the formulas of scripts/derive_alp_finite_mixing.py.
2. A one-dimensional null space of the meson block allows a zero of det D
   only if the meson null vector u obeys D_{aP} u = 0 (symbolic check).
3. Exceptional poles: (a) the pole vector has no ALP component (no shift
   exists); (b) a degenerate pole, where w is not unique and every choice
   aligns. The generalized eigenproblem, i.e. the null space of D(m_k^2),
   handles both.
"""

import sympy as sp

from derive_alp_finite_mixing import derive as derive_finite_mixing_formulas
from derive_alp_inverse_propagator_congruence import (
    SQUARED_MASS,
    inverse_propagator,
    neutral_test_quadratic_form,
    pole_normalization,
    rational_orthogonal_matrix,
)


def shift_matrix(meson_shift):
    """T for P_i = P_i' + w_i a', a = a'."""
    size = meson_shift.rows + 1
    conversion = sp.eye(size)
    conversion[0:size - 1, size - 1] = meson_shift
    return conversion


def pole_meson_overlap(kinetic, mass_squared, pole):
    """w^{(k)}(s) = [(M^2 - sK)^{-1}]_{PP} (s K_Pa - M^2_Pa) of the finite-mixing record."""
    meson = slice(0, kinetic.rows - 1)
    alp = kinetic.rows - 1
    meson_block = pole * kinetic[meson, meson] - mass_squared[meson, meson]
    mixed_block = pole * kinetic[meson, alp] - mass_squared[meson, alp]
    return -meson_block.inv() * mixed_block


def regular_pole_checks():
    data = neutral_test_quadratic_form()
    kinetic, mass_squared = data["kinetic"], data["mass_squared"]
    inverse = inverse_propagator(kinetic, mass_squared)
    finite_mixing = derive_finite_mixing_formulas()
    results = []
    for pole, vector in zip(data["pole_squared_masses"], data["pole_vectors"]):
        meson_shift = pole_meson_overlap(kinetic, mass_squared, pole)
        conversion = shift_matrix(meson_shift)
        kinetic_new = conversion.T * kinetic * conversion
        mass_squared_new = conversion.T * mass_squared * conversion
        inverse_new = conversion.T * inverse * conversion
        mixed_on_pole = pole * kinetic_new[0:3, 3] - mass_squared_new[0:3, 3]
        mixed_off_pole = sp.expand(inverse_new[0:3, 3])
        substitution = {
            finite_mixing["kinetic_meson"]: sp.ImmutableMatrix(kinetic[0:3, 0:3]),
            finite_mixing["mass_meson"]: sp.ImmutableMatrix(mass_squared[0:3, 0:3]),
            finite_mixing["kinetic_mixing"]: sp.ImmutableMatrix(kinetic[0:3, 3]),
            finite_mixing["mass_mixing"]: sp.ImmutableMatrix(mass_squared[0:3, 3]),
            finite_mixing["kinetic_alp"]: kinetic[3, 3],
            finite_mixing["mass_alp"]: mass_squared[3, 3],
            finite_mixing["s"]: pole,
        }
        overlap_from_record = sp.Matrix(finite_mixing["unnormalized_meson_overlap"].subs(substitution).doit())
        # derive_alp_finite_mixing.py uses P = P' - rho a, so rho = -w.
        alp_block_from_record = sp.Matrix(
            finite_mixing["converted_alp"].subs({**substitution, finite_mixing["meson_shift"]: sp.ImmutableMatrix(-meson_shift)}).doit()
        )[0]
        null_vectors = inverse.subs(SQUARED_MASS, pole).nullspace()
        generalized_vector = null_vectors[0] / null_vectors[0][3]
        results.append({
            "pole": pole,
            "w": meson_shift,
            "meson block invertible": (pole * kinetic[0:3, 0:3] - mass_squared[0:3, 0:3]).det() != 0,
            "w equals generalized-eigenproblem vector": generalized_vector[0:3, 0] == meson_shift and len(null_vectors) == 1,
            "m_k^2 K'_ia - M^2'_ia = 0": mixed_on_pole == sp.zeros(3, 1),
            "K'_ia nonzero": kinetic_new[0:3, 3] != sp.zeros(3, 1),
            "M^2'_ia nonzero": mass_squared_new[0:3, 3] != sp.zeros(3, 1),
            "D'_ia(s) nonzero away from the pole": mixed_off_pole != sp.zeros(3, 1),
            "D'_aa(m_k^2) = 0": inverse_new[3, 3].subs(SQUARED_MASS, pole) == 0,
            "primed pole vector is (0,0,0,1)": inverse_new.subs(SQUARED_MASS, pole) * sp.Matrix([0, 0, 0, 1]) == sp.zeros(4, 1),
            "K'_aa equals pole normalization": kinetic_new[3, 3] == pole_normalization(vector, kinetic),
            "agrees with derive_alp_finite_mixing overlap": overlap_from_record == meson_shift,
            "agrees with derive_alp_finite_mixing ALP block": sp.expand(alp_block_from_record - inverse_new[3, 3].subs(SQUARED_MASS, pole)) == 0,
        })
    return results


def singular_meson_block_determinant():
    """det D at s_k when the meson block has a one-dimensional null space.

    The meson block can be diagonalized by an orthogonal meson-only rotation,
    a special case of Phi = T Phi'; in that basis it is diag(0, d_2, d_3).
    """
    mixed_first, mixed_second, mixed_third, alp_entry, second, third = sp.symbols("b_1 b_2 b_3 c d_2 d_3")
    inverse_at_pole = sp.Matrix([
        [0, 0, 0, mixed_first],
        [0, second, 0, mixed_second],
        [0, 0, third, mixed_third],
        [mixed_first, mixed_second, mixed_third, alp_entry],
    ])
    determinant = sp.factor(inverse_at_pole.det())
    return determinant, sp.expand(determinant + mixed_first**2 * second * third) == 0


def decoupled_meson_example(degenerate):
    """K = 1, P_3 decoupled from (P_1, P_2, a), which mix through a rational rotation.

    degenerate=False: m_3^2 is not a pole of the mixed system; the pole s = m_3^2
    has pole vector (0, 0, 1, 0) only.
    degenerate=True: m_3^2 equals a pole of the mixed system; the null space of
    D(m_3^2) is two-dimensional.
    """
    rotation = rational_orthogonal_matrix(sp.Matrix([
        [0, sp.Rational(1, 2), sp.Rational(-1, 3)],
        [sp.Rational(-1, 2), 0, sp.Rational(1, 4)],
        [sp.Rational(1, 3), sp.Rational(-1, 4), 0],
    ]))
    mixed_poles = [sp.Rational(1, 10), sp.Rational(2, 5), sp.Rational(4, 5)]
    mixed_mass = rotation.T * sp.diag(*mixed_poles) * rotation
    decoupled_pole = mixed_poles[1] if degenerate else sp.Rational(3, 5)
    mass_squared = sp.zeros(4)
    order = [0, 1, 3]
    for row, full_row in enumerate(order):
        for column, full_column in enumerate(order):
            mass_squared[full_row, full_column] = mixed_mass[row, column]
    mass_squared[2, 2] = decoupled_pole
    kinetic = sp.eye(4)
    inverse = inverse_propagator(kinetic, mass_squared)
    inverse_at_pole = inverse.subs(SQUARED_MASS, decoupled_pole)
    meson_block = inverse_at_pole[0:3, 0:3]
    null_vectors = inverse_at_pole.nullspace()
    alp_components = [vector[3] for vector in null_vectors]
    result = {
        "pole": decoupled_pole,
        "meson block singular": meson_block.det() == 0,
        "null space dimension": len(null_vectors),
        "ALP components of null vectors": alp_components,
    }
    free_component = sp.Symbol("w_3")
    if degenerate:
        aligned = [vector / vector[3] for vector in null_vectors if vector[3] != 0][0]
        meson_shift = sp.Matrix([aligned[0], aligned[1], free_component])
        conversion = shift_matrix(meson_shift)
        inverse_new = conversion.T * inverse_at_pole * conversion
        result["shift aligns for every w_3"] = sp.expand(inverse_new[0:3, 3]) == sp.zeros(3, 1)
        # det D has a double zero here, so use the limit of (s - s_k) D^{-1}(s).
        residue = ((SQUARED_MASS - decoupled_pole) * inverse.adjugate() / inverse.det()).applyfunc(
            lambda entry: sp.cancel(entry).subs(SQUARED_MASS, decoupled_pole)
        )
        orthonormal = sp.GramSchmidt(null_vectors, orthonormal=True)  # K = 1 here
        result["residue rank"] = residue.rank()
        result["residue is sum over orthonormal pole vectors"] = sp.simplify(
            residue - sum((vector * vector.T for vector in orthonormal), sp.zeros(4))
        ) == sp.zeros(4)
    else:
        # The pole equation for w can be solved in the meson rows (w_3 free) but
        # the ALP row then fails: (w, 1) is not a null vector of D(m_k^2).
        meson_shift = sp.Matrix(sp.symbols("w_1 w_2 w_3"))
        meson_rows = sp.solve(list(meson_block * meson_shift + inverse_at_pole[0:3, 3]), list(meson_shift), dict=True)
        candidate = sp.Matrix([*meson_shift.subs(meson_rows[0]), 1])
        alp_row = sp.simplify((inverse_at_pole * candidate)[3])
        result["meson rows solvable with free w_3"] = len(meson_rows) == 1 and free_component not in meson_rows[0]
        result["ALP row of D(w,1) at the pole"] = alp_row
        result["no shift aligns this pole"] = alp_row != 0 and all(component == 0 for component in alp_components)
    return result


def main():
    print("Single-pole alignment P_i = P_i' + w_i^(k)(m_k^2) a'; exact rationals; N_f = 3 neutral mesons.")
    for entry in regular_pole_checks():
        print(f"pole s = {entry['pole']}: w = {list(entry['w'])}")
        for label, passed in entry.items():
            if isinstance(passed, bool):
                print(f"  [{'pass' if passed else 'FAIL'}] {label}")
    determinant, passed = singular_meson_block_determinant()
    print(f"det D(s_k) with meson block diag(0, d_2, d_3): {determinant}  [{'pass' if passed else 'FAIL'}]")
    for degenerate in (False, True):
        print("Exceptional pole, " + ("degenerate" if degenerate else "pole vector without ALP component") + ":")
        for label, value in decoupled_meson_example(degenerate).items():
            print(f"  {label}: {value}")


if __name__ == "__main__":
    main()
