#!/usr/bin/env python
"""Derive tree-level neutral-pseudoscalar resonance masses and mixing.

Invocation: .venv/bin/python scripts/derive_dynamical_pseudoscalar_masses.py
Tests:      .venv/bin/python scripts/derive_dynamical_pseudoscalar_masses.py --test

Natural units; chi_0 and squared masses have units GeV^2, f_pi and d_m
have units GeV. Real canonical resonance components use R=lambda_i R_i/sqrt(2).
Inputs: memory/nlo/quadratic-lagrangian.md and neutral-meson-matrices.md.
The p^4 matching omits resonance-mass spurion corrections and loops.
"""

import sys
import unittest

import sympy as sp


def show(label, expression):
    print(f"{label}:\n{sp.latex(expression)}\n")


def derive():
    """Calculate Hessians, poles and eigenvectors from quadratic Lagrangians."""
    chi_u, chi_d, chi_s = sp.symbols("chi_u chi_d chi_s", real=True)
    f_pi, d_m, resonance_mass_squared = sp.symbols(
        "f_pi d_m M_R^2", nonzero=True, real=True
    )
    meson_fields = sp.Matrix(sp.symbols("pi_3 eta_8 eta_0", real=True))
    resonance_fields = sp.Matrix(sp.symbols("R_3 R_8 R_0", real=True))
    generators = [sp.diag(1, -1, 0), sp.diag(1, 1, -2) / sp.sqrt(3),
                  sp.sqrt(sp.Rational(2, 3)) * sp.eye(3)]
    chi_0 = sp.diag(chi_u, chi_d, chi_s)
    phi = sum((g * p for g, p in zip(generators, meson_fields)), sp.zeros(3))
    resonance = sum((g * r / sp.sqrt(2) for g, r in
                     zip(generators, resonance_fields)), sp.zeros(3))
    chi_minus = -sp.I * (phi * chi_0 + chi_0 * phi) / f_pi
    interaction = sp.expand(sp.I * d_m * sp.trace(resonance * chi_minus))
    bilinear = sp.simplify(sp.Matrix([
        [sp.diff(interaction, r, p) for p in meson_fields]
        for r in resonance_fields
    ]))
    mass_correction = sp.simplify(-bilinear.T * bilinear / resonance_mass_squared)
    kinetic_correction = sp.simplify(bilinear.T * bilinear / resonance_mass_squared**2)
    local_lagrangian = sp.simplify(
        -d_m**2 * sp.trace(chi_minus * chi_minus) / (2 * resonance_mass_squared)
    )
    local_mass_hessian = sp.simplify(-sp.hessian(local_lagrangian, meson_fields))

    # Derive a two-light-state, one-resonance example in the LO mass eigenbasis.
    s, expansion = sp.symbols("s epsilon", real=True)
    mass_1, mass_2, mass_r = sp.symbols("m_1^2 m_2^2 M^2", real=True)
    mixing_1, mixing_2 = sp.symbols("b_1 b_2", real=True)
    mass_matrix = sp.Matrix([[mass_1, 0, -expansion * mixing_1],
                            [0, mass_2, -expansion * mixing_2],
                            [-expansion * mixing_1, -expansion * mixing_2, mass_r]])
    inverse_propagator = s * sp.eye(3) - mass_matrix
    effective = sp.simplify(inverse_propagator[:2, :2] -
                           inverse_propagator[:2, 2:] *
                           inverse_propagator[2:, 2:].inv() *
                           inverse_propagator[2:, :2])
    pole_shift = sp.symbols("delta_m_1^2", real=True)
    determinant = inverse_propagator.det()
    pole_equation = sp.expand(determinant.subs(s, mass_1 + expansion**2 * pole_shift))
    pole_shift_value = sp.factor(sp.solve(pole_equation.coeff(expansion, 2), pole_shift)[0])
    light_mixing, heavy_mixing = sp.symbols("v_2 v_R", real=True)
    eigenvector = sp.Matrix([1, expansion**2 * light_mixing, expansion * heavy_mixing])
    residual = (mass_matrix - (mass_1 + expansion**2 * pole_shift_value) * sp.eye(3)) * eigenvector
    heavy_value = sp.factor(sp.solve(sp.expand(residual[2]).coeff(expansion, 1), heavy_mixing)[0])
    light_value = sp.factor(sp.solve(
        sp.expand(residual[1].subs(heavy_mixing, heavy_value)).coeff(expansion, 2), light_mixing
    )[0])
    norm = sp.sqrt((eigenvector.T * eigenvector)[0])
    normalized = eigenvector.applyfunc(lambda entry: sp.factor(sp.series(
        (entry / norm).subs({light_mixing: light_value, heavy_mixing: heavy_value}),
        expansion, 0, 3).removeO()))

    # Derive the generalized light-field eigenproblem with residual local LECs.
    local_order = sp.symbols("t", real=True)
    kinetic_11, kinetic_12, kinetic_22 = sp.symbols("delta_K_11 delta_K_12 delta_K_22", real=True)
    mass_11, mass_12, mass_22 = sp.symbols("delta_M_11^2 delta_M_12^2 delta_M_22^2", real=True)
    delta_kinetic = sp.Matrix([[kinetic_11, kinetic_12], [kinetic_12, kinetic_22]])
    delta_mass = sp.Matrix([[mass_11, mass_12], [mass_12, mass_22]])
    local_mass_shift, local_mixing, local_normalization = sp.symbols("delta_m_local^2 v_local n_local", real=True)
    local_vector = sp.Matrix([1 + local_order * local_normalization, local_order * local_mixing])
    local_residual = ((sp.diag(mass_1, mass_2) + local_order * delta_mass) -
                      (mass_1 + local_order * local_mass_shift) *
                      (sp.eye(2) + local_order * delta_kinetic)) * local_vector
    local_normalization_equation = (local_vector.T * (sp.eye(2) + local_order * delta_kinetic) * local_vector)[0] - 1
    local_solution = sp.solve([sp.expand(entry).coeff(local_order, 1) for entry in local_residual] +
                              [sp.expand(local_normalization_equation).coeff(local_order, 1)],
                              [local_mass_shift, local_mixing, local_normalization])

    # One meson plus one resonance: derive the exact roots before expanding.
    two_mass_matrix = mass_matrix.extract([0, 2], [0, 2])
    exact_poles = sp.solve((s * sp.eye(2) - two_mass_matrix).det(), s)
    # Low-energy Taylor coefficients of the full resonance propagator.
    inverse_propagator_series = sp.series(1 / (mass_r - s), s, 0, 2).removeO()
    return locals()


class DynamicalMassTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.calculation = derive()

    def test_hessian_matching_to_full_nonet(self):
        c = self.calculation
        self.assertEqual(sp.simplify(c["mass_correction"] - c["local_mass_hessian"]), sp.zeros(3))

    def test_pole_and_eigenvector_orders(self):
        c = self.calculation
        residual = c["residual"].subs({c["light_mixing"]: c["light_value"],
                                       c["heavy_mixing"]: c["heavy_value"]})
        for entry in residual:
            self.assertEqual(sp.series(entry, c["expansion"], 0, 3).removeO().simplify(), 0)

    def test_exact_two_state_poles(self):
        c = self.calculation
        characteristic = (c["s"] * sp.eye(2) - c["two_mass_matrix"]).det()
        for pole in c["exact_poles"]:
            self.assertEqual(sp.simplify(characteristic.subs(c["s"], pole)), 0)

    def test_isospin_limit(self):
        c = self.calculation
        bilinear = c["bilinear"].subs(c["chi_d"], c["chi_u"])
        self.assertEqual(bilinear[0, 1], 0)
        self.assertEqual(bilinear[0, 2], 0)

    def test_schur_complement_determinant(self):
        c = self.calculation
        self.assertEqual(sp.factor(c["inverse_propagator"].det() -
                                   (c["s"] - c["mass_r"]) * c["effective"].det()), 0)

    def test_generalized_eigenproblem_and_normalization(self):
        c = self.calculation
        for entry in c["local_residual"]:
            self.assertEqual(sp.series(entry.subs(c["local_solution"]), c["local_order"],
                                       0, 2).removeO().simplify(), 0)
        self.assertEqual(sp.series(c["local_normalization_equation"].subs(c["local_solution"]),
                                   c["local_order"], 0, 2).removeO().simplify(), 0)


def main():
    c = derive()
    print("Tree level, real fields, nondegenerate LO poles; epsilon counts bilinear mixing.")
    print("Mass shifts first occur at epsilon^2; heavy admixtures at epsilon.\n")
    for label, key in [
        ("chi_minus at first field order", "chi_minus"),
        ("i d_m trace(R chi_minus)", "interaction"),
        ("B_ij from the mixed Hessian; M_PR^2=-B^T", "bilinear"),
        ("Delta M^2=-B^T B/M_R^2", "mass_correction"),
        ("Delta K=B^T B/M_R^4 (higher chiral order for B=O(p^2))", "kinetic_correction"),
        ("Local exchange Lagrangian", "local_lagrangian"),
        ("Explicit three-state mass matrix", "mass_matrix"),
        ("Schur complement inverse propagator", "effective"),
        ("Light-pole correction coefficient delta_m_1^2", "pole_shift_value"),
        ("Resonance admixture coefficient v_R", "heavy_value"),
        ("Other-light-state admixture coefficient v_2", "light_value"),
        ("Normalized light eigenvector through epsilon^2", "normalized"),
        ("Exact two-state squared masses", "exact_poles"),
        ("Low-energy resonance propagator", "inverse_propagator_series"),
        ("Generalized light-field local shift, mixing and normalization", "local_solution"),
    ]:
        show(label, c[key])


if __name__ == "__main__":
    if "--test" in sys.argv:
        unittest.main(argv=[sys.argv[0]])
    else:
        main()
