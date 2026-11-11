#!/usr/bin/env python3
"""Inverse propagator of the ALP and neutral mesons under a linear field change.

Run:   .venv/bin/python scripts/derive_alp_inverse_propagator_congruence.py
Tests: .venv/bin/python -m unittest discover -s scripts -p test_derive_alp_inverse_propagator_congruence.py

Conventions (Alda:2026nrd, memory/alp/shared/finite-alp-meson-mixing.md):
fields Phi_I = (P_1, P_2, P_3, a) for N_f = 3 neutral mesons and the ALP,
L_2 = 1/2 dPhi^T K dPhi - 1/2 Phi^T M^2 Phi, inverse propagator
D(s) = s K - M^2. The field change is Phi = T Phi' with constant invertible T.
A pole vector at s_k is a null vector v of D(s_k); in Alda:2026nrd's
normalization v = (w^{(k)}, 1). The amputated vertex linear in the external
leg is V_I Phi_I, so V' = T^T V.

Units: K is dimensionless, M^2 and s are squared masses in an arbitrary common
unit; the rational test values have no physical meaning.

Assumptions: K symmetric positive definite, M^2 symmetric, simple real poles.
The test data keep the leading-order kinetic form K_ij = delta_ij,
K_ia = -v_i, K_aa = 1 + v.v and use an arbitrary symmetric M^2 with rational
poles, built from a rational orthogonal matrix (Cayley transform).
All arithmetic is exact (sympy rationals).
"""

import random

import sympy as sp

SQUARED_MASS = sp.Symbol("s")


def rational_orthogonal_matrix(antisymmetric):
    """Cayley transform (1 - A)(1 + A)^{-1}: orthogonal with rational entries."""
    identity = sp.eye(antisymmetric.rows)
    return (identity - antisymmetric) * (identity + antisymmetric).inv()


def neutral_test_quadratic_form():
    """Exact rational K and M^2 on (P_1, P_2, P_3, a) with known poles.

    K = S^T S with S the leading-order kinetic shift, M^2 = S^T O^T diag(m^2) O S.
    The pole vector of pole k is S^{-1} O^T e_k, rescaled to unit ALP component.
    """
    axial_overlap = sp.Matrix([sp.Rational(1, 3), sp.Rational(-1, 5), sp.Rational(1, 2)])
    kinetic_shift = sp.eye(4)
    kinetic_shift[0:3, 3] = -axial_overlap
    antisymmetric = sp.Matrix([
        [0, sp.Rational(1, 2), sp.Rational(-1, 3), sp.Rational(1, 4)],
        [sp.Rational(-1, 2), 0, sp.Rational(1, 5), sp.Rational(-2, 3)],
        [sp.Rational(1, 3), sp.Rational(-1, 5), 0, sp.Rational(1, 2)],
        [sp.Rational(-1, 4), sp.Rational(2, 3), sp.Rational(-1, 2), 0],
    ])
    rotation = rational_orthogonal_matrix(antisymmetric)
    pole_squared_masses = [sp.Rational(1, 50), sp.Rational(3, 10), sp.Rational(9, 10), sp.Rational(1, 5)]
    kinetic = kinetic_shift.T * kinetic_shift
    mass_squared = kinetic_shift.T * rotation.T * sp.diag(*pole_squared_masses) * rotation * kinetic_shift
    pole_vectors = []
    for index in range(4):
        vector = kinetic_shift.inv() * rotation.T[:, index]
        pole_vectors.append(vector / vector[3])
    return {
        "axial_overlap": axial_overlap,
        "kinetic": kinetic,
        "mass_squared": mass_squared,
        "pole_squared_masses": pole_squared_masses,
        "pole_vectors": pole_vectors,
    }


def generic_invertible_matrix(size=4, seed=20261111):
    """Deterministic rational matrix with no structure; checked to be invertible."""
    generator = random.Random(seed)
    while True:
        matrix = sp.Matrix(size, size, lambda row, column: sp.Rational(generator.randint(-9, 9), generator.randint(1, 7)))
        if matrix.det() != 0:
            return matrix


def inverse_propagator(kinetic, mass_squared):
    return SQUARED_MASS * kinetic - mass_squared


def propagator_residue(inverse, pole):
    """Residue of D(s)^{-1} at a simple zero of det D: adj D(s_k) / det'D(s_k)."""
    determinant = inverse.det()
    return inverse.subs(SQUARED_MASS, pole).adjugate() / sp.diff(determinant, SQUARED_MASS).subs(SQUARED_MASS, pole)


def pole_normalization(vector, kinetic):
    """v^T K v; for v = (w, 1) this is K_aa + 2 K_ia w_i + K_ij w_i w_j."""
    return (vector.T * kinetic * vector)[0]


def one_leg_amplitude(vertex, vector, kinetic):
    """LSZ amplitude V.v / sqrt(v^T K v) of the finite-mixing record."""
    return (vertex.T * vector)[0] / sp.sqrt(pole_normalization(vector, kinetic))


def derive(seed=20261111):
    """Compare poles, pole vectors, residues and amplitudes before and after Phi = T Phi'."""
    data = neutral_test_quadratic_form()
    kinetic, mass_squared = data["kinetic"], data["mass_squared"]
    conversion = generic_invertible_matrix(seed=seed)
    vertex = sp.Matrix([sp.Rational(2, 3), sp.Rational(-5, 4), sp.Rational(1, 7), sp.Rational(3, 2)])

    inverse = inverse_propagator(kinetic, mass_squared)
    kinetic_new = conversion.T * kinetic * conversion
    mass_squared_new = conversion.T * mass_squared * conversion
    inverse_new = inverse_propagator(kinetic_new, mass_squared_new)
    vertex_new = conversion.T * vertex

    checks = {}
    checks["congruence"] = sp.expand(inverse_new - conversion.T * inverse * conversion) == sp.zeros(4)
    determinant = sp.expand(inverse.det())
    determinant_new = sp.expand(inverse_new.det())
    checks["determinant ratio is det(T)^2"] = sp.expand(determinant_new - conversion.det() ** 2 * determinant) == 0
    expected_determinant = kinetic.det() * sp.prod([SQUARED_MASS - pole for pole in data["pole_squared_masses"]])
    checks["zeros of det D are the chosen poles"] = sp.expand(determinant - expected_determinant) == 0
    checks["kinetic and mass blocks change individually"] = kinetic_new != kinetic and mass_squared_new != mass_squared

    per_pole = []
    for pole, vector in zip(data["pole_squared_masses"], data["pole_vectors"]):
        vector_new = conversion.inv() * vector
        residue = propagator_residue(inverse, pole)
        residue_new = propagator_residue(inverse_new, pole)
        amplitude = one_leg_amplitude(vertex, vector, kinetic)
        amplitude_new = one_leg_amplitude(vertex_new, vector_new, kinetic_new)
        rescaled_new = vector_new / vector_new[3]
        amplitude_rescaled = one_leg_amplitude(vertex_new, rescaled_new, kinetic_new)
        per_pole.append({
            "pole": pole,
            "null vector before": inverse.subs(SQUARED_MASS, pole) * vector == sp.zeros(4, 1),
            "null vector after (T^-1 v)": inverse_new.subs(SQUARED_MASS, pole) * vector_new == sp.zeros(4, 1),
            "residue is v v^T / v^T K v": residue == vector * vector.T / pole_normalization(vector, kinetic),
            "residue maps as T^-1 R T^-T": residue_new == conversion.inv() * residue * conversion.inv().T,
            "pole normalization invariant": pole_normalization(vector_new, kinetic_new) == pole_normalization(vector, kinetic),
            "V^T R V invariant": sp.simplify((vertex_new.T * residue_new * vertex_new)[0] - (vertex.T * residue * vertex)[0]) == 0,
            "LSZ amplitude invariant": sp.simplify(amplitude_new - amplitude) == 0,
            "amplitude squared invariant after a-normalization": sp.simplify(amplitude_rescaled**2 - amplitude**2) == 0,
            "amplitude": amplitude,
        })
    return {"data": data, "conversion": conversion, "checks": checks, "per_pole": per_pole}


def main():
    result = derive()
    print("Field change Phi = T Phi' on (P_1, P_2, P_3, a); exact rationals; generic T.")
    print(f"T = {result['conversion'].tolist()}")
    print(f"poles s_k = {result['data']['pole_squared_masses']}")
    for label, passed in result["checks"].items():
        print(f"[{'pass' if passed else 'FAIL'}] {label}")
    for entry in result["per_pole"]:
        print(f"pole s = {entry['pole']}: amplitude = {entry['amplitude']}")
        for label, passed in entry.items():
            if isinstance(passed, bool):
                print(f"  [{'pass' if passed else 'FAIL'}] {label}")


if __name__ == "__main__":
    main()
