#!/usr/bin/env python3
"""Derive quadratic ALP terms of O46, O47, O50 and O52 with SymPy.

Run: .venv/bin/python scripts/derive_theta_operator_quadratic_terms.py

Inputs: memory/nlo/quadratic-lagrangian.md and
memory/nlo/singlet-source-derivative-basis.md.
Fields and f_pi, f_a have mass dimension one; chi_0 has dimension two.
Lorentz indices are suppressed, with each derivative product contracted.
The matrix symbols are arbitrary noncommuting 3x3 matrices.
The expansion is about the real-mass vacuum and keeps quadratic field terms.
This script does not establish operator-level redundancy.
"""

import sympy as sp


def show(label, expression):
    print(f"{label}:\n{sp.latex(expression)}\n")


def derive():
    meson_marker, alp_marker = sp.symbols("t b", real=True)
    f_pi, f_a = sp.symbols("f_pi f_a", real=True, nonzero=True)
    c_G, c_GA = sp.symbols("c_G c_GA", real=True)
    alp_derivative, alp_second_derivative, alp = sp.symbols(
        "partial_a Box_a a", real=True
    )
    meson = sp.MatrixSymbol("Phi", 3, 3)
    meson_derivative = sp.MatrixSymbol("partial_Phi", 3, 3)
    axial = sp.MatrixSymbol("c_q^A", 3, 3)
    mass = sp.MatrixSymbol("chi_0", 3, 3)
    u_mu = (
        -meson_marker * meson_derivative / f_pi
        + 2 * alp_marker * alp_derivative * axial / f_a
    )
    chi_plus_vacuum = 2 * mass
    theta_derivative = 2 * c_GA * alp_marker * alp_derivative / f_a
    theta_second_derivative = (
        2 * c_GA * alp_marker * alp_second_derivative / f_a
    )
    show("Linear vielbein", u_mu)
    show("Vacuum chi_plus", chi_plus_vacuum)
    show("Covariant theta derivative in Alda notation", theta_derivative)
    operators = {
        46: theta_derivative * sp.Trace(u_mu * chi_plus_vacuum),
        47: theta_derivative * sp.Trace(u_mu) * sp.Trace(chi_plus_vacuum),
        50: -theta_derivative**2 * sp.Trace(chi_plus_vacuum),
        52: -theta_second_derivative * sp.Trace(chi_plus_vacuum),
    }
    results = {}
    for number, operator in operators.items():
        expanded = operator.expand().doit().expand()
        mixed = expanded.coeff(meson_marker, 1).coeff(alp_marker, 1)
        alp_quadratic = expanded.coeff(alp_marker, 2).subs(meson_marker, 0)
        results[number] = mixed
        show(f"O_{number} expanded about the vacuum", expanded)
        show(f"O_{number} mixed a-Phi coefficient", mixed)
        show(f"O_{number} quadratic ALP coefficient", alp_quadratic)

    # The operator O52 is parity odd before multiplication by its coefficient
    # function. Its first parity-even contribution may have a coefficient
    # linear in X_TC, rather than a nonzero constant at X_TC=0.
    coefficient_derivative = sp.Symbol("L_52_prime_0", real=True)
    topological_source = (
        -meson_marker * sp.Trace(meson) / f_pi
        - 2 * c_G * alp_marker * alp / f_a
    )
    coefficient_dependent_operator = (
        coefficient_derivative * topological_source * operators[52]
    ).expand().doit()
    mixed_52 = coefficient_dependent_operator.coeff(
        meson_marker, 1
    ).coeff(alp_marker, 1)
    show("Linear topological source X_TC", topological_source)
    show("O52 with coefficient L52'(0) X_TC, mixed term", mixed_52)

    coordinate = sp.Symbol("x", real=True)
    theta = sp.Function("theta")(coordinate)
    mass_trace = sp.Function("mass_trace")(coordinate)
    total_derivative = sp.diff(sp.diff(theta, coordinate) * mass_trace, coordinate)
    # Product-rule expansion makes the boundary term explicit.
    product_rule_remainder = sp.expand(
        total_derivative - sp.diff(theta, coordinate, 2) * mass_trace
    )
    show("Integration-by-parts boundary derivative", total_derivative)
    show("O52 remainder after the second derivative is removed", product_rule_remainder)
    mixed_52_after_parts = mixed_52.subs(
        alp_second_derivative * sp.Trace(meson),
        -alp_derivative * sp.Trace(meson_derivative),
    )
    show("O52 mixed term after integration by parts", mixed_52_after_parts)
    # At the operator level, integrate a coefficient linear in X_TC by parts.
    # Suppress Lorentz indices and retain the derivative of the coefficient:
    # -L52' X div(Dtheta) tr(chi_+) =
    # L52' Dtheta dX tr(chi_+) + L52' X Dtheta d tr(chi_+).
    source_derivative = sp.Symbol("D_theta", real=True)
    vielbein_trace = sp.Symbol("Q_mu", real=True)
    topological_derivative = sp.Symbol("partial_X_TC", real=True)
    mass_trace_symbol = sp.Symbol("trace_chi_plus", real=True)
    derivative_mass_trace = sp.Symbol("partial_trace_chi_plus", real=True)
    topological_symbol = sp.Symbol("X_TC", real=True)
    after_parts = coefficient_derivative * source_derivative * (
        topological_derivative * mass_trace_symbol
        + topological_symbol * derivative_mass_trace
    )
    show("Coefficient-dependent O52 after integration by parts", after_parts)
    source_identity_substitution = after_parts.subs(
        topological_derivative, vielbein_trace - source_derivative
    ).expand()
    show("O52 reshuffled using dX_TC = Q_mu - D_theta", source_identity_substitution)
    return results


if __name__ == "__main__":
    derive()
