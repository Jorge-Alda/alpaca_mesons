#!/usr/bin/env python3
"""Derive CCWZ restrictions on extra ALP kinetic mixing and LSZ covariance.

Run: .venv/bin/python scripts/derive_ccwz_kinetic_covariance.py
Tests: .venv/bin/python scripts/derive_ccwz_kinetic_covariance.py --test

Inputs: memory/nlo/kinetic-block-alignment.md,
memory/nlo/singlet-source-derivative-basis.md, and
memory/nlo/lsz-invariant-split.md. Use Phi=lambda_i P_i with
tr(lambda_i lambda_j)=2 delta_ij and the (3,8,0) neutral basis.
Fields and decay constants have mass dimension one; chi_0 and M^2 have
dimension two. K, v, Delta, and LSZ overlaps are dimensionless.
Retain quadratic fields and first order in 1/f_a, about a real-mass vacuum.
The abstract K and M^2 are symmetric and invertible where inverted.
"""

import sys
import unittest

import sympy as sp


def show(label, expression):
    print(f"{label}:\n{sp.latex(expression)}\n")


def source_kinetic_terms():
    """Expand V4, O46 and O47 using the CCWZ derivative building blocks."""
    f_pi, f_a = sp.symbols("f_pi f_a", nonzero=True, real=True)
    c_GA, g_46, g_47, V_4 = sp.symbols("c_GA g_46 g_47 V_4", real=True)
    derivative_a = sp.Symbol("partial_a", real=True)
    derivative_fields = sp.symbols("partial_pi_3 partial_eta_8 partial_eta_0")
    axial_coefficients = sp.symbols("c_A_u c_A_d c_A_s", real=True)
    mass_spurion = sp.diag(*sp.symbols("chi_u chi_d chi_s", real=True))
    generators = (
        sp.diag(1, -1, 0),
        sp.diag(1, 1, -2) / sp.sqrt(3),
        sp.sqrt(sp.Rational(2, 3)) * sp.eye(3),
    )
    derivative_phi = sum((generator * derivative for generator, derivative
                          in zip(generators, derivative_fields)), sp.zeros(3))
    vielbein = -derivative_phi / f_pi + 2 * derivative_a / f_a * sp.diag(
        *axial_coefficients
    )
    theta_derivative = 2 * c_GA * derivative_a / f_a
    terms = {
        "V4": -V_4 * theta_derivative * sp.trace(vielbein),
        "O46": g_46 * theta_derivative * sp.trace(vielbein * 2 * mass_spurion),
        "O47": g_47 * theta_derivative * sp.trace(vielbein)
        * sp.trace(2 * mass_spurion),
    }
    contributions = {
        label: sp.simplify(sp.Matrix([
            sp.diff(term, derivative_a, derivative) for derivative in derivative_fields
        ]))
        for label, term in terms.items()
    }
    return theta_derivative, vielbein, contributions


def derive():
    theta_derivative, vielbein, contributions = source_kinetic_terms()
    show("Pure ALP singlet source derivative", theta_derivative)
    show("Linear CCWZ vielbein", vielbein)
    for label, contribution in contributions.items():
        show(f"{label} extra mixed kinetic column Delta in (3,8,0)", contribution)

    c_G, trace_axial, trace_kappa = sp.symbols("c_G tr_c_A tr_kappa")
    invariant_source = c_G + trace_axial
    transformed_invariant_source = invariant_source.subs(
        {c_G: c_G - trace_kappa, trace_axial: trace_axial + trace_kappa},
        simultaneous=True,
    )
    show("c_GA after the ALP-dependent axial transformation",
         sp.expand(transformed_invariant_source))

    dimension = sp.Symbol("n", integer=True, positive=True)
    kinetic = sp.MatrixSymbol("K", dimension, dimension)
    mass = sp.MatrixSymbol("M^2", dimension, dimension)
    axial_overlap = sp.MatrixSymbol("v", dimension, 1)
    extra_mixing = sp.MatrixSymbol("Delta", dimension, 1)
    field_shift = sp.MatrixSymbol("rho", dimension, 1)
    mixed_mass = sp.MatrixSymbol("M^2_mixed", dimension, 1)
    momentum_squared = sp.Symbol("s")
    mixed_kinetic = -kinetic * axial_overlap + extra_mixing
    show("Mixed kinetic column from the two CCWZ derivative structures", mixed_kinetic)
    show("Kinetic contribution to w", (-kinetic.I * mixed_kinetic).expand().doit())
    invariant_pole_source = (mass * kinetic.I * mixed_kinetic - mixed_mass).expand().doit()
    show("Invariant pole source", invariant_pole_source)
    overlap = (mass - momentum_squared * kinetic).I * (
        momentum_squared * mixed_kinetic - mixed_mass
    )
    show("Overlap from the LSZ equation", overlap)
    transformed_kinetic = mixed_kinetic - kinetic * field_shift
    transformed_mass = mixed_mass - mass * field_shift
    transformed_extra = (
        transformed_kinetic + kinetic * (axial_overlap + field_shift)
    ).expand().doit()
    show("Extra mixed kinetic column after the axial transformation", transformed_extra)
    overlap_change = ((mass - momentum_squared * kinetic).I * (
        momentum_squared * (transformed_kinetic - mixed_kinetic)
        - (transformed_mass - mixed_mass)
    )).expand().doit()
    show("Change of the overlap from the transformed mixed blocks", overlap_change)
    show("Factorized overlap change", (mass - momentum_squared * kinetic).I * (
        mass - momentum_squared * kinetic) * field_shift)


class KineticCovarianceTests(unittest.TestCase):
    def test_non_diagonal_meson_blocks(self):
        """Solve the LSZ equations before and after a field shift."""
        kinetic = sp.Matrix([[2, 1], [1, 3]])
        mass = sp.Matrix([[7, 2], [2, 11]])
        axial_overlap = sp.Matrix([sp.Rational(1, 7), sp.Rational(-2, 9)])
        extra_mixing = sp.Matrix([0, sp.Rational(1, 11)])
        mixed_kinetic = -kinetic * axial_overlap + extra_mixing
        mixed_mass = sp.Matrix([sp.Rational(3, 13), sp.Rational(-5, 17)])
        field_shift = sp.Matrix([sp.Rational(2, 19), sp.Rational(-1, 23)])
        momentum_squared = sp.Rational(2, 5)
        pole_operator = mass - momentum_squared * kinetic
        overlap = pole_operator.inv() * (momentum_squared * mixed_kinetic - mixed_mass)
        transformed_kinetic = mixed_kinetic - kinetic * field_shift
        transformed_mass = mixed_mass - mass * field_shift
        transformed_overlap = pole_operator.inv() * (
            momentum_squared * transformed_kinetic - transformed_mass
        )
        self.assertEqual(transformed_overlap - overlap, field_shift)
        self.assertEqual(transformed_kinetic + kinetic * (
            axial_overlap + field_shift), extra_mixing)


if __name__ == "__main__":
    if "--test" in sys.argv:
        unittest.main(argv=[sys.argv[0]])
    else:
        derive()
