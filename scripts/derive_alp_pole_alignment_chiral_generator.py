#!/usr/bin/env python3
"""The pole-alignment shift as a constant-generator chiral redefinition.

Run:   .venv/bin/python scripts/derive_alp_pole_alignment_chiral_generator.py
Tests: .venv/bin/python -m unittest discover -s scripts -p test_derive_alp_pole_alignment_chiral_generator.py

Conventions (Alda:2026nrd; memory/alp/shared/alp-chiral-source-transformations.md,
alp-general-chiral-redefinition-invariants.md, finite-alp-meson-mixing.md):
x = a/f_a, U = exp(i Phi/f_pi), Phi = P_i lambda_i with the nine U(3)
generators, tr(lambda_i lambda_j) = 2 delta_ij, lambda_0 = sqrt(2/3) 1.
g_R = exp(i x G_R), g_L = exp(i x G_L), U' = g_R U g_L^dagger,
chi' = g_R chi g_L^dagger, r_mu = (d_mu x) C_R, l_mu = (d_mu x) C_L,
C_A = (C_R - C_L)/2, Theta' = Theta - (1/2i) log(det g_R / det g_L).
Leading-order Lagrangian, as in scripts/derive_alp_lo_finite_mixing.py:
  1/2 (da)^2 + f_pi^2/4 <u u> + f_pi^2/4 <chi U^dag + U chi^dag>
  - f_pi^2 m_0^2/(4 N_f) X_TC^2 - 1/2 m_a0^2 a^2,
with u = -dPhi/f_pi + 2 (dx) C_A at linear order, chi(x) = chi_0 + i x {c_P, chi_0}
and X_TC = -<Phi>/f_pi - 2 Theta(x), Theta = c_G x.

Generator choice: G_R = -G_L = -(f_a/(2 f_pi)) w_i lambda_i, constant.
Then g_R = exp[-i a w_i lambda_i/(2 f_pi)] is the element called L_a in the
finite-mixing record, and g_L = L_a^dagger = R_a.

Checks (exact, symbolic in all nine w_i; rational source data):
1. Hessian of the original Lagrangian reproduces K and M^2 of the
   finite-mixing record.
2. U' = g_R U g_L^dagger has linear part P_i' = P_i - w_i a, a' = a.
3. The Hessian computed from the transformed sources C_A', chi'(x) (to x^2)
   and Theta' equals T^T K T and T^T M^2 T with T the shift P = P' + w a'.
4. c_q^A - c_q^P is unchanged, and c_G' = c_G + (f_a/(2 f_pi)) w_i <lambda_i>.
Units: f_pi, f_a in an arbitrary mass unit; chi_0, m_0^2, m_a0^2 in its square.
"""

import sympy as sp

from derive_alp_single_pole_alignment import shift_matrix

FLAVOR_COUNT = 3


def u3_generators():
    """Gell-Mann matrices lambda_1..lambda_8 and lambda_0 = sqrt(2/3) 1, ordered (0, 1, ..., 8)."""
    i = sp.I
    gell_mann = [
        sp.sqrt(sp.Rational(2, 3)) * sp.eye(3),
        sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]]),
        sp.Matrix([[0, -i, 0], [i, 0, 0], [0, 0, 0]]),
        sp.Matrix([[1, 0, 0], [0, -1, 0], [0, 0, 0]]),
        sp.Matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]]),
        sp.Matrix([[0, 0, -i], [0, 0, 0], [i, 0, 0]]),
        sp.Matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]]),
        sp.Matrix([[0, 0, 0], [0, 0, -i], [0, i, 0]]),
        sp.Matrix([[1, 0, 0], [0, 1, 0], [0, 0, -2]]) / sp.sqrt(3),
    ]
    return gell_mann


def source_data():
    i = sp.I
    return {
        "f_pi": sp.Rational(23, 250),
        "f_a": sp.Rational(3, 25),
        "chi_0": sp.diag(sp.Rational(1, 50), sp.Rational(1, 40), sp.Rational(1, 2)),
        "c_A": sp.Matrix([
            [sp.Rational(1, 2), sp.Rational(1, 3) + i / 5, 0],
            [sp.Rational(1, 3) - i / 5, sp.Rational(-1, 4), sp.Rational(1, 6) * i],
            [0, -sp.Rational(1, 6) * i, sp.Rational(2, 3)],
        ]),
        "c_P": sp.Matrix([
            [sp.Rational(-1, 3), 0, sp.Rational(1, 5) - i / 7],
            [0, sp.Rational(1, 2), sp.Rational(1, 4)],
            [sp.Rational(1, 5) + i / 7, sp.Rational(1, 4), sp.Rational(1, 7)],
        ]),
        "c_G": sp.Rational(1, 3),
        "m_0_sq": sp.Rational(4, 5),
        "m_a0_sq": sp.Rational(1, 4),
    }


def truncated_exponential(matrix, order):
    return sum((matrix**n / sp.factorial(n) for n in range(order + 1)), sp.zeros(matrix.rows))


def quadratic_hessian(data, axial_source, mass_source_coefficients, topological_slope):
    """K and M^2 of the LO Lagrangian on (P_0..P_8, a).

    mass_source_coefficients = [chi_0, chi_1, chi_2] with chi(x) = chi_0 + x chi_1 + x^2 chi_2;
    topological_slope is dTheta/dx (Theta linear in x for constant generators).
    """
    generators = u3_generators()
    f_pi, f_a = data["f_pi"], data["f_a"]
    mesons = sp.symbols("P_0:9", real=True)
    meson_derivatives = sp.symbols("dP_0:9", real=True)
    alp, alp_derivative = sp.symbols("a da", real=True)
    fields = [*mesons, alp]
    derivatives = [*meson_derivatives, alp_derivative]

    derivative_meson_matrix = sum((d * g for d, g in zip(meson_derivatives, generators)), sp.zeros(3))
    vielbein = -derivative_meson_matrix / f_pi + 2 * alp_derivative / f_a * axial_source
    kinetic_lagrangian = alp_derivative**2 / 2 + f_pi**2 / 4 * (vielbein * vielbein).trace()

    meson_matrix = sum((p * g for p, g in zip(mesons, generators)), sp.zeros(3))
    goldstone = truncated_exponential(sp.I * meson_matrix / f_pi, 2)
    x = alp / f_a
    chi = mass_source_coefficients[0] + x * mass_source_coefficients[1] + x**2 * mass_source_coefficients[2]
    topological = -meson_matrix.trace() / f_pi - 2 * topological_slope * x
    mass_lagrangian = (
        f_pi**2 / 4 * (chi * goldstone.H + goldstone * chi.H).trace()
        - f_pi**2 * data["m_0_sq"] / (4 * FLAVOR_COUNT) * topological**2
        - data["m_a0_sq"] * alp**2 / 2
    )
    zero_fields = {field: 0 for field in fields}
    kinetic = sp.Matrix(10, 10, lambda r, c: sp.diff(kinetic_lagrangian, derivatives[r], derivatives[c]))
    mass_squared = sp.Matrix(10, 10, lambda r, c: -sp.diff(mass_lagrangian, fields[r], fields[c]).subs(zero_fields))
    return kinetic.applyfunc(sp.expand), mass_squared.applyfunc(sp.expand)


def record_blocks(data):
    """K and M^2 blocks as written in memory/alp/shared/finite-alp-meson-mixing.md."""
    generators = u3_generators()
    f_pi, f_a, chi_0 = data["f_pi"], data["f_a"], data["chi_0"]
    axial_overlap = [f_pi / f_a * (data["c_A"] * g).trace() for g in generators]
    kinetic = sp.eye(10)
    mass_squared = sp.zeros(10)
    for i, gi in enumerate(generators):
        kinetic[i, 9] = kinetic[9, i] = -axial_overlap[i]
        for j, gj in enumerate(generators):
            mass_squared[i, j] = ((gi * gj + gj * gi) * chi_0).trace() / 4 + data["m_0_sq"] / (2 * FLAVOR_COUNT) * gi.trace() * gj.trace()
        mass_squared[i, 9] = mass_squared[9, i] = (
            -f_pi / (2 * f_a) * (gi * (data["c_P"] * chi_0 + chi_0 * data["c_P"])).trace()
            + f_pi * data["m_0_sq"] * data["c_G"] / (FLAVOR_COUNT * f_a) * gi.trace()
        )
    kinetic[9, 9] = 1 + sum(v * v for v in axial_overlap)
    mass_squared[9, 9] = data["m_a0_sq"] + 2 * f_pi**2 * data["m_0_sq"] * data["c_G"]**2 / (FLAVOR_COUNT * f_a**2)
    return kinetic.applyfunc(sp.expand), mass_squared.applyfunc(sp.expand)


def derive():
    data = source_data()
    generators = u3_generators()
    f_pi, f_a, chi_0 = data["f_pi"], data["f_a"], data["chi_0"]
    meson_shift = sp.Matrix(sp.symbols("w_0:9", real=True))
    generator = -f_a / (2 * f_pi) * sum((w * g for w, g in zip(meson_shift, generators)), sp.zeros(3))
    right_generator, left_generator = generator, -generator
    x, alp = sp.symbols("x a", real=True)

    checks = {}
    # 1. Original Hessian versus the finite-mixing record.
    chi_1 = sp.I * (data["c_P"] * chi_0 + chi_0 * data["c_P"])
    kinetic, mass_squared = quadratic_hessian(data, data["c_A"], [chi_0, chi_1, sp.zeros(3)], data["c_G"])
    record_kinetic, record_mass = record_blocks(data)
    checks["original K equals finite-mixing record"] = (kinetic - record_kinetic).applyfunc(sp.simplify) == sp.zeros(10)
    checks["original M^2 equals finite-mixing record"] = (mass_squared - record_mass).applyfunc(sp.simplify) == sp.zeros(10)

    # 2. g_R is L_a of the finite-mixing record; linear part of U' = g_R U g_L^dagger.
    exponent_from_record = -sp.I * alp / (2 * f_pi) * sum((w * g for w, g in zip(meson_shift, generators)), sp.zeros(3))
    checks["i x G_R = -i a w.lambda/(2 f_pi) (g_R = L_a)"] = (
        (sp.I * (alp / f_a) * right_generator - exponent_from_record).applyfunc(sp.expand) == sp.zeros(3)
    )
    checks["G_L = -G_R (g_L = L_a^dagger = R_a)"] = left_generator == -right_generator
    order_parameter = sp.Symbol("epsilon")
    mesons = sp.Matrix(sp.symbols("P_0:9", real=True))
    meson_matrix = sum((p * g for p, g in zip(mesons, generators)), sp.zeros(3))
    goldstone = sp.eye(3) + sp.I * order_parameter * meson_matrix / f_pi
    right_element = sp.eye(3) + sp.I * order_parameter * (alp / f_a) * right_generator
    left_element_dagger = sp.eye(3) - sp.I * order_parameter * (alp / f_a) * left_generator
    transformed = (right_element * goldstone * left_element_dagger).applyfunc(
        lambda entry: sp.expand(entry).coeff(order_parameter, 1)
    )
    transformed_mesons = [sp.expand(f_pi * (transformed * g).trace() / (2 * sp.I)) for g in generators]
    checks["linear map P_i' = P_i - w_i a"] = all(
        sp.expand(new - old + w * alp) == 0 for new, old, w in zip(transformed_mesons, mesons, meson_shift)
    )

    # 3. Transformed sources and their Hessian.
    right_element = truncated_exponential(sp.I * x * right_generator, 2)
    left_element = truncated_exponential(sp.I * x * left_generator, 2)
    transformed_chi = (right_element * (chi_0 + x * chi_1) * left_element.H).applyfunc(sp.expand)
    transformed_chi_coefficients = [transformed_chi.applyfunc(lambda entry: entry.coeff(x, n)) for n in range(3)]
    transformed_axial = data["c_A"] + (right_generator - left_generator) / 2
    trace_difference = (right_generator - left_generator).trace()
    transformed_topological_slope = data["c_G"] - trace_difference / 2
    kinetic_new, mass_squared_new = quadratic_hessian(
        data, transformed_axial, transformed_chi_coefficients, transformed_topological_slope
    )
    conversion = shift_matrix(meson_shift)
    checks["transformed K equals T^T K T"] = (kinetic_new - conversion.T * kinetic * conversion).applyfunc(sp.simplify) == sp.zeros(10)
    checks["transformed M^2 equals T^T M^2 T"] = (mass_squared_new - conversion.T * mass_squared * conversion).applyfunc(sp.simplify) == sp.zeros(10)
    checks["transformed theory has an a^2 quark-source term"] = transformed_chi_coefficients[2] != sp.zeros(3)

    # 4. Source shifts.
    transformed_pseudoscalar_anticommutator = -sp.I * transformed_chi_coefficients[1]
    shifted_pseudoscalar = data["c_P"] + right_generator
    checks["c_P' = c_P + G_R (from chi' at order x)"] = (
        transformed_pseudoscalar_anticommutator - (shifted_pseudoscalar * chi_0 + chi_0 * shifted_pseudoscalar)
    ).applyfunc(sp.expand) == sp.zeros(3)
    checks["c_A' - c_P' = c_A - c_P"] = (
        (transformed_axial - shifted_pseudoscalar) - (data["c_A"] - data["c_P"])
    ).applyfunc(sp.expand) == sp.zeros(3)
    expected_gluon_shift = f_a / (2 * f_pi) * sum(w * g.trace() for w, g in zip(meson_shift, generators))
    checks["c_G' = c_G + (f_a/2 f_pi) w_i <lambda_i>"] = sp.expand(transformed_topological_slope - data["c_G"] - expected_gluon_shift) == 0
    return {
        "checks": checks,
        "trace of G_R - G_L": sp.simplify(trace_difference),
        "gluon coefficient shift": sp.simplify(expected_gluon_shift),
    }


def general_generator_linear_map():
    """Linear part of U' = g_R U g_L^dagger for arbitrary g_eta(x) with g_eta(0) = 1.

    Only the first Taylor coefficients G_eta1 = -i dg_eta/dx(0) enter. The
    result is P_i' = P_i - w_i a, a' = a with
    w_i = -(f_pi/(2 f_a)) <lambda_i (G_R1 - G_L1)>: always shift type, fixed by
    the axial part; the vector part and higher Taylor terms start at aP or a^2.
    """
    generators = u3_generators()
    f_pi, f_a, alp, order_parameter = sp.symbols("f_pi f_a a epsilon", positive=True)
    right_first = sum((c * g for c, g in zip(sp.symbols("r_0:9", real=True), generators)), sp.zeros(3))
    left_first = sum((c * g for c, g in zip(sp.symbols("l_0:9", real=True), generators)), sp.zeros(3))
    mesons = sp.symbols("P_0:9", real=True)
    meson_matrix = sum((p * g for p, g in zip(mesons, generators)), sp.zeros(3))
    goldstone = sp.eye(3) + sp.I * order_parameter * meson_matrix / f_pi
    right_element = sp.eye(3) + sp.I * order_parameter * alp / f_a * right_first
    left_element = sp.eye(3) + sp.I * order_parameter * alp / f_a * left_first
    transformed = (right_element * goldstone * left_element.H).applyfunc(
        lambda entry: sp.expand(entry).coeff(order_parameter, 1)
    )
    transformed_mesons = [sp.expand(f_pi * (transformed * g).trace() / (2 * sp.I)) for g in generators]
    shifts = [-f_pi / (2 * f_a) * ((right_first - left_first) * g).trace() for g in generators]
    return all(
        sp.expand(new - old + shift * alp) == 0
        for new, old, shift in zip(transformed_mesons, mesons, shifts)
    )


def main():
    print(f"[{'pass' if general_generator_linear_map() else 'FAIL'}] any g(x) with g(0)=1: linear map is a shift, w_i = -(f_pi/2f_a)<lambda_i (G_R1 - G_L1)>")
    result = derive()
    print("Constant generators G_R = -G_L = -(f_a/2f_pi) w_i lambda_i; U(3), all nine w_i symbolic.")
    for label, passed in result["checks"].items():
        print(f"[{'pass' if passed else 'FAIL'}] {label}")
    print(f"tr(G_R - G_L) = {result['trace of G_R - G_L']}")
    print(f"c_G' - c_G = {result['gluon coefficient shift']}")


if __name__ == "__main__":
    main()
