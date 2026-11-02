#!/usr/bin/env python3
"""Derive the kinetic-block relation and its assumptions with SymPy.

Run: .venv/bin/python scripts/derive_alp_kinetic_alignment.py
Inputs: memory/nlo/quadratic-lagrangian.md and lsz-invariant-split.md.
Fields have mass dimension one; K and ALP mixing coefficients are dimensionless.
Only quadratic, local, two-derivative terms and first order in 1/f_a are used.
"""

import sympy as sp


def show(label, expression):
    print(f"{label}:\n{sp.latex(expression)}\n")


def derive():
    dimension = sp.Symbol("n", integer=True, positive=True)
    kinetic = sp.MatrixSymbol("K", dimension, dimension)
    meson_derivative = sp.MatrixSymbol("partial_P", dimension, 1)
    axial_projection = sp.MatrixSymbol("v", dimension, 1)
    alp_derivative = sp.Symbol("partial_a", real=True)
    derivative_combination = meson_derivative - axial_projection * alp_derivative
    quadratic = sp.Trace(derivative_combination.T * kinetic * derivative_combination) / 2
    show("Quadratic kinetic form; K is real, symmetric and invertible", quadratic)
    gradient = sp.diff(quadratic, meson_derivative)
    mixed = sp.diff(gradient, alp_derivative).subs(kinetic.T, kinetic).doit()
    show("Mixed Hessian K_meson,alp", mixed)
    show("Minus inverse meson kinetic matrix times mixed block", (-kinetic.I * mixed).doit())

    f_pi, f_a = sp.symbols("f_pi f_a", nonzero=True, real=True)
    L_4, L_5, L_18, Lambda_1 = sp.symbols("L_4 L_5 L_18 Lambda_1", real=True)
    derivative_field = sp.MatrixSymbol("partial_Phi", 3, 3)
    axial = sp.MatrixSymbol("c_q^A", 3, 3)
    mass = sp.MatrixSymbol("chi_0", 3, 3)
    meson_marker, alp_marker = sp.symbols("t b", real=True)
    derivative_operators = {
        "LO": sp.Trace(derivative_field * derivative_field) / 4,
        "L5": 2 * L_5 / f_pi**2 * sp.Trace(mass * derivative_field * derivative_field),
        "Lambda1": Lambda_1 / 12 * sp.Trace(derivative_field)**2,
        "L4": 2 * L_4 / f_pi**2 * sp.Trace(derivative_field * derivative_field) * sp.Trace(mass),
        "L18": 2 * L_18 / f_pi**2 * sp.Trace(mass * derivative_field) * sp.Trace(derivative_field),
    }
    shift = meson_marker * derivative_field - 2 * f_pi / f_a * axial * alp_marker
    show("Common derivative replacement from linear u_mu", shift)
    mixed_terms = {}
    for name, operator in derivative_operators.items():
        show(f"{name} meson quadratic form", operator)
        kinetic_operator = (2 * sp.diff(operator, derivative_field).T).doit()
        show(f"{name} kinetic operator applied to partial_Phi", kinetic_operator)
        # For L=tr(Y K[Y])/4, polarization gives the mixed term
        # -f_pi/f_a tr(Y K[c_q^A]); the first derivation above supplies
        # the same polarization for a general symmetric quadratic form.
        mixed_term = (-f_pi / f_a * sp.Trace(derivative_field * kinetic_operator.subs(derivative_field, axial))).expand().doit()
        mixed_terms[name] = mixed_term
        show(f"{name} mixed derivative coefficient", mixed_term)

    independent_mixed = sp.MatrixSymbol("Delta_k", dimension, 1)
    show("Additional independent mixed kinetic coefficient", independent_mixed)
    show("Result in its presence", (-kinetic.I * (mixed + independent_mixed)).expand().doit())
    scalar_kinetic, scalar_axial, independent = sp.symbols("K v Delta_k", nonzero=True)
    scalar_lagrangian = scalar_kinetic * (meson_marker - scalar_axial * alp_marker)**2 / 2 + independent * meson_marker * alp_marker
    scalar_mixed = sp.diff(scalar_lagrangian, meson_marker, alp_marker)
    show("One-meson counterexample mixed Hessian", scalar_mixed)
    show("One-meson counterexample kinetic overlap", sp.simplify(-scalar_mixed / sp.diff(scalar_lagrangian, meson_marker, 2)))
    return mixed_terms


if __name__ == "__main__":
    derive()
