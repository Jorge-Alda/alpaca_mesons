"""Derive the ALP insertion recursion for partonic bilinear generators.

Run: .venv/bin/python scripts/derive_alp_partonic_operator_induction.py

Natural units. z records ALP insertion number and t records the mass
dimension of the normalized operator, including its explicit factors of
1/f_a. J_k denotes the open Lorentz tensor (partial^k a)/f_a, with grade
(1,k). S_r denotes the full breaking block (a/f_a)^r S_r, with grade
(r,source_dimension). These are formal insertion labels, not operator
counts. The ALP-free seeds M_{0,d} retain open Lorentz and gauge indices.

The script derives grade/dimension recursions by polynomial multiplication
and coefficient extraction. It does not contract Lorentz/gauge indices or
reduce integration-by-parts, equations-of-motion, or algebraic relations.
Every insertion position and compatible gauge representation must be
included when interpreting the formal tensor products as operators.
"""

import argparse

import sympy as sp


def show(label, expression):
    print(f"{label}:\n{sp.latex(expression)}\n")


def derive_recursion(maximum_order, maximum_dimension, source_dimension):
    """Extract a finite part of M=M_0+(J+S)M in two gradings.

    SymPy treats ALP derivative tensors as commuting bosonic factors and
    spurions/seeds as noncommuting matrices. Numerical multiplicities are
    repeated generation paths in a linear span, not independent operators.
    """
    insertion_grade, operator_dimension = sp.symbols("z t", commutative=True)
    derivative_labels = {
        degree: sp.Symbol(f"J_{{{degree}}}", commutative=True)
        for degree in range(1, maximum_dimension - 2)
    }
    breaking_labels = {
        degree: sp.Symbol(f"S_{{{degree}}}", commutative=False)
        for degree in range(1, maximum_order + 1)
    }
    seed_labels = {
        dimension: sp.Symbol(f"M_{{0,{dimension}}}", commutative=False)
        for dimension in range(3, maximum_dimension + 1)
    }
    derivative_insertions = sp.Add(*(
        insertion_grade * operator_dimension**degree * label
        for degree, label in derivative_labels.items()))
    breaking_insertions = sp.Add(*(
        insertion_grade**degree * operator_dimension**source_dimension * label
        for degree, label in breaking_labels.items()))
    seeds = sp.Add(*(
        operator_dimension**dimension * label
        for dimension, label in seed_labels.items()))
    all_labels = {
        (order, dimension): sp.Symbol(f"M_{{{order},{dimension}}}",
                                      commutative=False)
        for order in range(maximum_order + 1)
        for dimension in range(3, maximum_dimension + 1)
    }
    formal_generators = sp.Add(*(
        insertion_grade**order * operator_dimension**dimension * label
        for (order, dimension), label in all_labels.items()))
    expanded_right_side = sp.expand(
        seeds + (derivative_insertions + breaking_insertions) * formal_generators)
    recursions = {
        (order, dimension): sp.expand(expanded_right_side.coeff(insertion_grade, order)
                                      .coeff(operator_dimension, dimension))
        for order in range(maximum_order + 1)
        for dimension in range(3, maximum_dimension + 1)
    }
    # Each insertion has strictly positive ALP grade, so the following
    # substitution is an induction in order without an implicit equation.
    generated = {}
    substitutions = {}
    for order in range(maximum_order + 1):
        for dimension in range(3, maximum_dimension + 1):
            result = sp.expand(recursions[order, dimension].xreplace(substitutions))
            generated[order, dimension] = result
            substitutions[all_labels[order, dimension]] = result
    return derivative_insertions, breaking_insertions, seeds, recursions, generated


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--order", type=int, default=2)
    parser.add_argument("--dimension", type=int, default=6)
    parser.add_argument("--source-dimension", type=int, default=1)
    arguments = parser.parse_args()
    if arguments.order < 0 or arguments.dimension < 3 or arguments.source_dimension < 0:
        parser.error("require order >= 0, dimension >= 3, source-dimension >= 0")
    derivative_insertions, breaking_insertions, seeds, recursions, generated = derive_recursion(
        arguments.order, arguments.dimension, arguments.source_dimension)
    print("Partonic quark-bilinear tensor generators; natural units.")
    print("z grades ALP insertions; t grades normalized operator mass dimension.")
    print("J_k=(partial^k a)/f_a has grade (1,k); its unnormalized numerator")
    print("has canonical dimension k+1. S_r includes (a/f_a)^r and its spurion.")
    print("Open Lorentz and gauge indices, insertion positions, and gauge")
    print("representations remain symbolic. Coefficients above one count")
    print("generation paths in a span; they are not operator multiplicities.")
    show("Derivative-ALP insertion polynomial", derivative_insertions)
    show("Explicit-breaking insertion polynomial", breaking_insertions)
    show("ALP-free tensor seeds", seeds)
    print("Coefficient extraction from M=M_0+(J+S)M gives:")
    for (order, dimension), expression in recursions.items():
        if order and expression != 0:
            show(f"Recursion for M_{{{order},{dimension}}}", expression)
    print("Inductively substituted generators:")
    for (order, dimension), expression in generated.items():
        if expression != 0:
            show(f"Grade/dimension ({order},{dimension})", expression)


if __name__ == "__main__":
    main()
