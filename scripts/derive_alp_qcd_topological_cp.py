"""Derive ALP-order anomaly shifts and CP parity in the QCD gluon sector.

Run: .venv/bin/python scripts/derive_alp_qcd_topological_cp.py --order 5

Natural units; x=a/f_a is dimensionless and the ALP is CP odd. The input
anomaly convention is Theta'(x)=Theta(x)-tr A(x), calibrated by the
Alda:2026nrd relation c_G'=c_G-tr kappa for A(x)=x kappa. This script
performs Taylor and CP-sign bookkeeping; it does not recalculate the QCD
anomaly normalization. A_n are abstract flavor matrices. Theta_n and
tr(A_n) are dimensionless. The selected finite series is local at x=0.
"""

import argparse

import sympy as sp


def show(label, expression):
    """Print a deterministic expression using the standard LaTeX printer."""
    print(f"{label}:\n{sp.latex(expression)}\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--order", type=int, default=5)
    arguments = parser.parse_args()
    if arguments.order < 1:
        parser.error("--order must be positive")

    alp_coordinate = sp.Symbol("x", real=True)
    flavor_dimension = sp.Symbol("N_f", integer=True, positive=True)
    topological_coefficients = sp.symbols(f"Theta_0:{arguments.order + 1}", real=True)
    axial_coefficients = tuple(
        sp.MatrixSymbol(f"A_{degree}", flavor_dimension, flavor_dimension)
        for degree in range(1, arguments.order + 1))
    topological_source = sp.Add(*(
        coefficient * alp_coordinate**degree / sp.factorial(degree)
        for degree, coefficient in enumerate(topological_coefficients)))
    axial_trace = sp.Add(*(
        sp.Trace(coefficient) * alp_coordinate**degree / sp.factorial(degree)
        for degree, coefficient in enumerate(axial_coefficients, start=1)))
    transformed_source = sp.expand(topological_source - axial_trace)

    print("Dimensionless x=a/f_a; A(0)=0; no flavor commutativity assumption.")
    print("Input anomaly convention: Theta'=Theta-tr A.")
    show("Input Theta(x)", topological_source)
    show("Input tr A(x)", axial_trace)
    show("Transformed Theta(x)", transformed_source)
    for degree in range(arguments.order + 1):
        transformed_derivative = sp.diff(transformed_source, alp_coordinate, degree).subs(
            alp_coordinate, 0)
        show(f"Transformed Taylor derivative {degree} at zero", transformed_derivative)

    print("CP inputs: x -> -x, G Gtilde -> -G Gtilde, G^2 -> +G^2.")
    for degree in range(arguments.order + 1):
        alp_monomial = alp_coordinate**degree
        transformed_monomial = alp_monomial.subs(alp_coordinate, -alp_coordinate)
        topological_sign = sp.simplify(-transformed_monomial / alp_monomial)
        scalar_sign = sp.simplify(transformed_monomial / alp_monomial)
        show(f"CP sign of x^{degree} G Gtilde", topological_sign)
        show(f"CP sign of x^{degree} G^2", scalar_sign)
    derivative_number = sp.Symbol("n", integer=True, nonnegative=True)
    show("General CP sign of x^n G Gtilde", -(-1)**derivative_number)
    show("General CP sign of x^n G^2", (-1)**derivative_number)
    print("CP-even topological monomials have odd n; CP-even scalar-gluon monomials have even n.")


if __name__ == "__main__":
    main()
