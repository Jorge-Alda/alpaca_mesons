#!/usr/bin/env python3
"""Alignment of several poles of the ALP-meson inverse propagator.

Run:   .venv/bin/python scripts/derive_alp_multi_pole_alignment.py
Tests: .venv/bin/python -m unittest discover -s scripts -p test_derive_alp_multi_pole_alignment.py

Conventions as in derive_alp_inverse_propagator_congruence.py and Alda:2026nrd:
Phi = (P_i, a), D(s) = s K - M^2, Phi = T Phi', D' = T^T D T. Pole k is
aligned with coordinate J of Phi' when the column J of D'(s_k) vanishes,
i.e. column J of T is a pole vector v^{(k)}. For the ALP coordinate
(J = a) and a shift T = [[1, w], [0, 1]] this is D'_{ia}(s_k) = 0.

Checks:
1. Shift-type maps: two poles are aligned by one w only if
   w^{(1)}(m_1^2) = w^{(2)}(m_2^2); with K positive definite and m_1^2 != m_2^2
   this never happens (exact counterexample, rank test, symbolic identity).
   A shift also cannot diagonalize D(s) at all s.
2. General T: pole vectors of distinct poles are K-orthogonal; T built from
   all of them gives T^T K T = diag and T^T M^2 T = diag, so D'(s) is diagonal
   at every s. Aligning a subset decouples those coordinates from each other
   at every s, but not from the rest unless the other columns are
   K-orthogonal to them.
3. The claim needs K positive definite (or at least a pencil without
   coinciding Jordan structure): a counterexample with indefinite K has a
   double pole in D^{-1}, which no T can diagonalize.
4. Floating-point check (scipy.linalg.eigh) of the pencil diagonalization for
   generic symmetric M^2 and positive definite K with irrational poles.
"""

import numpy as np
import scipy.linalg
import sympy as sp

from derive_alp_inverse_propagator_congruence import (
    SQUARED_MASS,
    generic_invertible_matrix,
    inverse_propagator,
    neutral_test_quadratic_form,
    pole_normalization,
)
from derive_alp_single_pole_alignment import pole_meson_overlap, shift_matrix


def shift_type_two_pole_counterexample(first=0, second=3):
    data = neutral_test_quadratic_form()
    kinetic, mass_squared = data["kinetic"], data["mass_squared"]
    poles = data["pole_squared_masses"]
    first_overlap = pole_meson_overlap(kinetic, mass_squared, poles[first])
    second_overlap = pole_meson_overlap(kinetic, mass_squared, poles[second])
    # Stack the two alignment conditions D_PP(s_k) w = -D_Pa(s_k), k = 1, 2.
    coefficient_rows, right_hand_side = [], []
    for pole in (poles[first], poles[second]):
        inverse_at_pole = inverse_propagator(kinetic, mass_squared).subs(SQUARED_MASS, pole)
        coefficient_rows.append(inverse_at_pole[0:3, 0:3])
        right_hand_side.append(-inverse_at_pole[0:3, 3])
    coefficients = sp.Matrix.vstack(*coefficient_rows)
    augmented = coefficients.row_join(sp.Matrix.vstack(*right_hand_side))
    # A common (w, 1) would be a null vector of D(s_1) - D(s_2) = (s_1 - s_2) K.
    vector = sp.Matrix([*first_overlap, 1])
    difference = sp.expand(
        inverse_propagator(kinetic, mass_squared).subs(SQUARED_MASS, poles[first])
        - inverse_propagator(kinetic, mass_squared).subs(SQUARED_MASS, poles[second])
        - (poles[first] - poles[second]) * kinetic
    )
    return {
        "w first": first_overlap,
        "w second": second_overlap,
        "overlaps differ": first_overlap != second_overlap,
        "rank of conditions": coefficients.rank(),
        "rank with right-hand side": augmented.rank(),
        "no common shift": augmented.rank() > coefficients.rank(),
        "D(s_1) - D(s_2) = (s_1 - s_2) K": difference == sp.zeros(4),
        "v^T K v > 0 forbids a common pole vector": bool(pole_normalization(vector, kinetic) > 0),
        "shift aligning first pole misaligns second": (
            shift_matrix(first_overlap).T * inverse_propagator(kinetic, mass_squared).subs(SQUARED_MASS, poles[second]) * shift_matrix(first_overlap)
        )[0:3, 3] != sp.zeros(3, 1),
    }


def shift_cannot_diagonalize():
    """A shift keeps K_PP, M^2_PP; removing both mixed blocks needs two different w."""
    data = neutral_test_quadratic_form()
    kinetic, mass_squared = data["kinetic"], data["mass_squared"]
    kinetic_removing_shift = -kinetic[0:3, 0:3].inv() * kinetic[0:3, 3]
    mass_removing_shift = -mass_squared[0:3, 0:3].inv() * mass_squared[0:3, 3]
    meson_mass_block = mass_squared[0:3, 0:3]
    return {
        "shift removing K_ia (canonical shift w = v)": kinetic_removing_shift,
        "shift removing M^2_ia": mass_removing_shift,
        "the two shifts differ": kinetic_removing_shift != mass_removing_shift,
        "canonical shift equals axial overlap v": kinetic_removing_shift == data["axial_overlap"],
        "M^2_PP is not diagonal": meson_mass_block != sp.diag(*[meson_mass_block[i, i] for i in range(3)]),
    }


def general_conversion_checks():
    data = neutral_test_quadratic_form()
    kinetic, mass_squared = data["kinetic"], data["mass_squared"]
    poles, vectors = data["pole_squared_masses"], data["pole_vectors"]
    inverse = inverse_propagator(kinetic, mass_squared)
    gram = sp.Matrix(4, 4, lambda k, l: (vectors[k].T * kinetic * vectors[l])[0])
    normalized = [vector / sp.sqrt(pole_normalization(vector, kinetic)) for vector in vectors]
    full_alignment = sp.Matrix.hstack(*normalized)
    kinetic_new = sp.simplify(full_alignment.T * kinetic * full_alignment)
    mass_squared_new = sp.simplify(full_alignment.T * mass_squared * full_alignment)
    inverse_new = sp.simplify(full_alignment.T * inverse * full_alignment)

    # Align only poles 0 and 3 (columns 3, 4); columns 1, 2 generic.
    generic = generic_invertible_matrix(seed=99)
    partial = sp.Matrix.hstack(generic[:, 0], generic[:, 1], vectors[0], vectors[3])
    partial_inverse = sp.expand(partial.T * inverse * partial)
    # K-orthogonal completion: project the generic columns away from v_0, v_3.
    completion = []
    for column in (generic[:, 0], generic[:, 1]):
        projected = column
        for vector in (vectors[0], vectors[3]):
            projected = projected - vector * (vector.T * kinetic * column)[0] / pole_normalization(vector, kinetic)
        completion.append(projected)
    orthogonal_partial = sp.Matrix.hstack(*completion, vectors[0], vectors[3])
    orthogonal_inverse = sp.expand(orthogonal_partial.T * inverse * orthogonal_partial)
    return {
        "pole vectors K-orthogonal": all(gram[k, l] == 0 for k in range(4) for l in range(4) if k != l),
        "T^T K T = 1": kinetic_new == sp.eye(4),
        "T^T M^2 T = diag(m_k^2)": mass_squared_new == sp.diag(*poles),
        "D'(s) diagonal at every s": all(inverse_new[k, l] == 0 for k in range(4) for l in range(4) if k != l),
        "partial: T invertible": partial.det() != 0,
        "partial: aligned columns vanish at their poles": (
            partial_inverse.subs(SQUARED_MASS, poles[0])[:, 2] == sp.zeros(4, 1)
            and partial_inverse.subs(SQUARED_MASS, poles[3])[:, 3] == sp.zeros(4, 1)
        ),
        "partial: aligned coordinates decoupled at every s": partial_inverse[2, 3] == 0,
        "partial: aligned and generic coordinates still mix off pole": partial_inverse[0:2, 2:4] != sp.zeros(2, 2),
        "partial, K-orthogonal completion: block diagonal at every s": orthogonal_inverse[0:2, 2:4] == sp.zeros(2, 2),
    }


def indefinite_kinetic_counterexample():
    """K = [[0,1],[1,0]] (not positive definite), M^2 = diag(1, 0)."""
    kinetic = sp.Matrix([[0, 1], [1, 0]])
    mass_squared = sp.diag(1, 0)
    inverse = inverse_propagator(kinetic, mass_squared)
    propagator = inverse.inv().applyfunc(sp.cancel)
    double_pole = propagator.applyfunc(lambda entry: sp.cancel(entry * SQUARED_MASS**2).subs(SQUARED_MASS, 0))
    null_vector = inverse.subs(SQUARED_MASS, 0).nullspace()[0]
    return {
        "det D": sp.factor(inverse.det()),
        "pole vector": list(null_vector),
        "v^T K v": pole_normalization(null_vector, kinetic),
        "coefficient of 1/s^2 in D^{-1}": double_pole,
        "double pole present, so no T makes D'(s) diagonal": double_pole != sp.zeros(2),
    }


def floating_point_pencil_check(seed=4):
    """Generic K > 0 and symmetric M^2 with irrational poles (scipy.linalg.eigh)."""
    generator = np.random.default_rng(seed)
    factor = generator.normal(size=(4, 4))
    kinetic = factor.T @ factor + np.eye(4)
    symmetric = generator.normal(size=(4, 4))
    mass_squared = symmetric @ symmetric.T + 0.1 * np.eye(4)
    poles, vectors = scipy.linalg.eigh(mass_squared, kinetic)
    kinetic_new = vectors.T @ kinetic @ vectors
    mass_squared_new = vectors.T @ mass_squared @ vectors
    return {
        "poles": poles,
        "max |T^T K T - 1|": float(np.max(np.abs(kinetic_new - np.eye(4)))),
        "max |T^T M^2 T - diag|": float(np.max(np.abs(mass_squared_new - np.diag(poles)))),
    }


def main():
    print("Several-pole alignment; exact rationals for N_f = 3 neutral mesons unless stated.")
    sections = [
        ("Shift-type map, two poles", shift_type_two_pole_counterexample()),
        ("Shift-type map and diagonalization at all s", shift_cannot_diagonalize()),
        ("General T built from pole vectors", general_conversion_checks()),
        ("Indefinite K (outside the physical case)", indefinite_kinetic_counterexample()),
        ("Floating-point pencil check, generic data", floating_point_pencil_check()),
    ]
    for title, entries in sections:
        print(title + ":")
        for label, value in entries.items():
            if isinstance(value, bool):
                print(f"  [{'pass' if value else 'FAIL'}] {label}")
            else:
                print(f"  {label}: {value if not isinstance(value, sp.MatrixBase) else value.tolist()}")


if __name__ == "__main__":
    main()
