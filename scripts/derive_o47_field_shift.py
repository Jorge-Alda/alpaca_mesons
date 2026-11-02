#!/usr/bin/env python3
"""Derive the quadratic meson-field shift for the O47 kinetic contribution.

Run: .venv/bin/python scripts/derive_o47_field_shift.py

Inputs: memory/nlo/singlet-source-derivative-basis.md,
memory/nlo/neutral-meson-matrices.md, and memory/nlo/lsz-invariant-split.md.
Conventions: Phi=sum_i lambda_i P_i, tr(lambda_i lambda_j)=2 delta_ij,
lambda_0=sqrt(2/3) I_3. All fields, f_pi, and f_a have mass dimension one.
K and rho are dimensionless; M^2 has dimension two. The coefficient g_47
multiplies O47 and is dimensionless. Retain quadratic fields and first order
in 1/f_a. No further chiral truncation is applied to the abstract K formula.
"""

import sys
import unittest

import sympy as sp


def show(label, expression):
    print(f"{label}:\n{sp.latex(expression)}\n")


def derive():
    f_pi, f_a = sp.symbols("f_pi f_a", real=True, nonzero=True)
    c_GA, g_47 = sp.symbols("c_GA g_47", real=True)
    derivative_a = sp.Symbol("partial_a", real=True)
    derivative_pi3, derivative_eta8, derivative_eta0 = sp.symbols(
        "partial_pi_3 partial_eta_8 partial_eta_0", real=True
    )
    mass_spurion = sp.MatrixSymbol("chi_0", 3, 3)
    generators = (
        sp.diag(1, -1, 0),
        sp.diag(1, 1, -2) / sp.sqrt(3),
        sp.sqrt(sp.Rational(2, 3)) * sp.eye(3),
    )
    derivative_fields = (
        derivative_pi3,
        derivative_eta8,
        derivative_eta0,
    )
    derivative_phi = sum(
        (generator * derivative for generator, derivative in zip(
            generators, derivative_fields
        )),
        sp.zeros(3),
    )
    show("Trace of the meson derivative matrix", sp.trace(derivative_phi))
    mixed_lagrangian = (
        -4 * g_47 * c_GA / (f_a * f_pi)
        * derivative_a * sp.trace(derivative_phi) * sp.Trace(mass_spurion)
    )
    show("O47 quadratic mixed Lagrangian", mixed_lagrangian)
    kinetic_contribution = sp.Matrix([
        sp.diff(mixed_lagrangian, derivative_a, derivative)
        for derivative in derivative_fields
    ])
    show("O47 mixed kinetic block in the (3,8,0) basis", kinetic_contribution)
    show("Singlet shift with the leading meson kinetic block I_3",
         sp.Eq(sp.Symbol("rho_0"), kinetic_contribution[2]))

    dimension = sp.Symbol("n", integer=True, positive=True)
    kinetic = sp.MatrixSymbol("K", dimension, dimension)
    mass_squared = sp.MatrixSymbol("M^2", dimension, dimension)
    mixed_kinetic = sp.MatrixSymbol("K_mixed", dimension, 1)
    mixed_mass = sp.MatrixSymbol("M^2_mixed", dimension, 1)
    rho = sp.MatrixSymbol("rho", dimension, 1)
    alp_kinetic = sp.MatrixSymbol("K_aa", 1, 1)
    alp_mass_squared = sp.MatrixSymbol("M^2_aa", 1, 1)
    transformation = sp.BlockMatrix([
        [sp.Identity(dimension), -rho],
        [sp.ZeroMatrix(1, dimension), sp.Identity(1)],
    ])
    full_kinetic = sp.BlockMatrix([
        [kinetic, mixed_kinetic],
        [mixed_kinetic.T, alp_kinetic],
    ])
    full_mass_squared = sp.BlockMatrix([
        [mass_squared, mixed_mass],
        [mixed_mass.T, alp_mass_squared],
    ])
    transformed_kinetic = sp.block_collapse(
        transformation.T * full_kinetic * transformation
    )
    transformed_mass = sp.block_collapse(
        transformation.T * full_mass_squared * transformation
    )
    derived_mixed_kinetic = transformed_kinetic.blocks[0, 1]
    derived_mixed_mass = transformed_mass.blocks[0, 1]
    show("Full field-coordinate transformation matrix", transformation)
    show("Shifted mixed kinetic block", derived_mixed_kinetic)
    show("Shifted mixed mass block", derived_mixed_mass)

    contribution_47 = sp.MatrixSymbol("Delta_K_47", dimension, 1)
    aligned_mixed = sp.MatrixSymbol("K_aligned_mixed", dimension, 1)
    selected_rho = kinetic.I * contribution_47
    show("Shift that removes the O47 kinetic contribution", selected_rho)
    remaining_kinetic = derived_mixed_kinetic.subs(
        mixed_kinetic, aligned_mixed + contribution_47
    ).subs(rho, selected_rho).doit()
    show("Remaining mixed kinetic block", remaining_kinetic)
    show("Induced change of the mixed mass block",
         (derived_mixed_mass - mixed_mass).subs(rho, selected_rho).doit())

    source_row = sp.MatrixSymbol("V", 1, dimension)
    direct_alp_source = sp.MatrixSymbol("V_a", 1, 1)
    source_coefficients = sp.BlockMatrix([[source_row, direct_alp_source]])
    transformed_source = sp.block_collapse(source_coefficients * transformation)
    shifted_direct_source = transformed_source.blocks[0, 1]
    show("Direct ALP source after the meson-field shift", shifted_direct_source)
    show("Direct ALP source change for the O47 shift",
         (shifted_direct_source - direct_alp_source).subs(
             rho, selected_rho
         ).doit())


class FieldShiftTests(unittest.TestCase):
    def test_direct_substitution_with_nondiagonal_kinetic_matrix(self):
        """Compare the block result with scalar-polynomial Hessians."""
        kinetic = sp.Matrix([[2, 1], [1, 3]])
        mass_squared = sp.Matrix([[7, 2], [2, 11]])
        extra_mixing = sp.Matrix([0, sp.Symbol("d")])
        selected_rho = kinetic.inv() * extra_mixing
        meson1, meson2, alp = sp.symbols("P1 P2 a")
        meson_fields = sp.Matrix([meson1, meson2])
        shifted_fields = meson_fields - selected_rho * alp
        kinetic_polynomial = (
            (shifted_fields.T * kinetic * shifted_fields)[0] / 2
            + (shifted_fields.T * extra_mixing)[0] * alp
        )
        mass_polynomial = -(shifted_fields.T * mass_squared * shifted_fields)[0] / 2
        derived_kinetic_mixing = sp.Matrix([
            sp.diff(kinetic_polynomial, field, alp)
            for field in meson_fields
        ])
        derived_mass_mixing = sp.Matrix([
            -sp.diff(mass_polynomial, field, alp)
            for field in meson_fields
        ])
        self.assertEqual(sp.simplify(derived_kinetic_mixing), sp.zeros(2, 1))
        self.assertEqual(
            sp.simplify(derived_mass_mixing + mass_squared * selected_rho),
            sp.zeros(2, 1),
        )


if __name__ == "__main__":
    if "--test" in sys.argv:
        unittest.main(argv=[sys.argv[0]])
    else:
        derive()
