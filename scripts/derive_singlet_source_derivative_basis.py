#!/usr/bin/env python3
"""Convert Kaiser:1998ds singlet derivative terms to Alda:2026nrd notation.

Run: .venv/bin/python scripts/derive_singlet_source_derivative_basis.py

Inputs: memory/extracts/Kaiser:1998ds/singlet-source-derivative-terms.md
and memory/nlo/lagrangian.md. Lorentz contractions are suppressed: each
polynomial product denotes the corresponding bilinear contracted with g_mu_nu.
The derivatives have mass dimension one, and V_3, V_4, V_5 have dimension two.
The algebra is exact at quadratic derivative order. The choice V_4=0 is the
field convention established in Kaiser:1998ds, rather than a new derivation
of the full nonlinear field redefinition.
"""

import sympy as sp


def show(label, expression):
    names = {
        sp.Symbol("vielbein_trace", real=True): r"\langle u_\mu\rangle",
        sp.Symbol("topological_derivative", real=True): r"\left(\partial_\mu X_{TC}\right)",
        sp.Symbol("psi_derivative", real=True): r"\left(D_\mu\psi\right)",
        sp.Symbol("theta_derivative", real=True): r"\left(D_\mu\theta\right)",
        sp.Symbol("alp_derivative", real=True): r"\left(\partial_\mu a\right)",
    }
    print(f"{label}:\n{sp.latex(expression, symbol_names=names)}\n")


def derive():
    vielbein_trace = sp.Symbol("vielbein_trace", real=True)
    topological_derivative = sp.Symbol("topological_derivative", real=True)
    psi_derivative = sp.Symbol("psi_derivative", real=True)
    theta_derivative = sp.Symbol("theta_derivative", real=True)
    V_3, V_4, V_5 = sp.symbols("V_3 V_4 V_5", real=True)
    replacements = {
        psi_derivative: -vielbein_trace,
        theta_derivative: vielbein_trace - topological_derivative,
    }
    for derivative, replacement in replacements.items():
        show("Source convention", sp.Eq(derivative, replacement))

    lagrangian = (
        V_3 * psi_derivative**2
        + V_4 * psi_derivative * theta_derivative
        + V_5 * theta_derivative**2
    )
    show("Kaiser singlet derivative terms", lagrangian)
    for coefficient in (V_3, V_4, V_5):
        converted_operator = (
            lagrangian.diff(coefficient).subs(replacements).expand()
        )
        show(f"Operator multiplying {coefficient}", converted_operator)
    converted = lagrangian.subs(replacements).expand()
    show(
        "Full expression in the Alda operator basis",
        sp.collect(converted, [vielbein_trace**2, vielbein_trace * topological_derivative]),
    )
    reduced = converted.subs(V_4, 0)
    show("Expression after choosing the V_4=0 field convention", reduced)
    show(
        "Remaining V_5 operator as an external-source square",
        sp.factor(reduced.diff(V_5)),
    )

    c_GA, f_a = sp.symbols("c_GA f_a", real=True, nonzero=True)
    alp_derivative = sp.Symbol("alp_derivative", real=True)
    source_difference = 2 * c_GA * alp_derivative / f_a
    show(
        "Alda Eq. (III.16) source identity",
        sp.Eq(vielbein_trace - topological_derivative, source_difference),
    )
    show(
        "V_5 term for the dynamical ALP source",
        (V_5 * theta_derivative**2).subs(theta_derivative, source_difference),
    )
    print(
        "The V_4 operator can be removed by the cited meson-field convention. "
        "The V_5 source square remains; its quadratic ALP term is of order 1/f_a^2."
    )


if __name__ == "__main__":
    derive()
