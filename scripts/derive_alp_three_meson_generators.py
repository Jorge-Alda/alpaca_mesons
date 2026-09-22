"""Expand the one-ALP, three-meson Lagrangian in generic generators.

Invocation:
    .venv/bin/python scripts/derive_alp_three_meson_generators.py

The scalar meson fields commute, while their abstract flavour generators do
not. A symbol dphi_i denotes partial^mu phi_i; the common Lorentz contraction
with partial_mu a is suppressed in the symbolic flavour calculation.
"""

from itertools import product

import sympy as sp


def show(label: str, expression: sp.Expr) -> None:
    print(f"{label}:\n{sp.latex(expression)}\n")


def commutator(left: sp.MatrixExpr, right: sp.MatrixExpr) -> sp.MatrixExpr:
    return left * right - right * left


n = sp.Symbol("N_f", integer=True, positive=True)
f_pi, f_a, B_0 = sp.symbols("f_pi f_a B_0", positive=True)
a, da = sp.symbols("a partial_mu_a")
phi_i, phi_j, phi_k = sp.symbols("phi_i phi_j phi_k")
dphi_i, dphi_j, dphi_k = sp.symbols(
    "partial^mu_phi_i partial^mu_phi_j partial^mu_phi_k"
)

lambda_i = sp.MatrixSymbol("lambda_i", n, n)
lambda_j = sp.MatrixSymbol("lambda_j", n, n)
lambda_k = sp.MatrixSymbol("lambda_k", n, n)
cA = sp.MatrixSymbol("c_q^A", n, n)
P = sp.MatrixSymbol("P", n, n)

Phi = phi_i * lambda_i + phi_j * lambda_j + phi_k * lambda_k
dPhi = dphi_i * lambda_i + dphi_j * lambda_j + dphi_k * lambda_k

field_terms = (
    (phi_i, lambda_i),
    (phi_j, lambda_j),
    (phi_k, lambda_k),
)
derivative_terms = (
    (dphi_i, lambda_i),
    (dphi_j, lambda_j),
    (dphi_k, lambda_k),
)
scalar_variables = (phi_i, phi_j, phi_k, dphi_i, dphi_j, dphi_k)

# SymPy does not distribute powers of abstract MatrixExpr sums. Retain the
# finite distributive expansion as (commuting scalar, MatrixExpr) term pairs.
kinetic_terms = [
    (
        scalar_1 * scalar_2 * dscalar,
        cA * commutator(generator_1, commutator(generator_2, dgenerator)),
    )
    for (scalar_1, generator_1), (scalar_2, generator_2), (dscalar, dgenerator)
    in product(field_terms, field_terms, derivative_terms)
]
breaking_terms = [
    (
        scalar_1 * scalar_2 * scalar_3,
        generator_1 * generator_2 * generator_3 * P,
    )
    for (scalar_1, generator_1), (scalar_2, generator_2), (scalar_3, generator_3)
    in product(field_terms, repeat=3)
]


def coefficient(
    terms: list[tuple[sp.Expr, sp.MatrixExpr]], powers: dict[sp.Symbol, int]
) -> sp.MatrixExpr:
    """Extract a scalar monomial coefficient while preserving matrix order."""
    monomial = sp.prod(symbol**order for symbol, order in powers.items())
    result = sp.ZeroMatrix(n, n)
    for scalar, matrix in terms:
        scalar_coefficient = sp.Poly(scalar, *scalar_variables).coeff_monomial(monomial)
        result += scalar_coefficient * matrix
    return result

# Three identical mesons i,i,i.
kin_iii = sp.Trace(coefficient(kinetic_terms, {phi_i: 2, dphi_i: 1}))
break_iii = sp.Trace(coefficient(breaking_terms, {phi_i: 3}))

# Two identical mesons i,i and one distinct meson j. The two derivative
# placements are independent monomials in the Lagrangian.
kin_iij_dj = sp.Trace(coefficient(kinetic_terms, {phi_i: 2, dphi_j: 1}))
kin_iij_di = sp.Trace(
    coefficient(kinetic_terms, {phi_i: 1, phi_j: 1, dphi_i: 1})
)
break_iij = sp.Trace(coefficient(breaking_terms, {phi_i: 2, phi_j: 1}))

# Three distinct mesons i,j,k. Extract one term for each differentiated field.
kin_ijk_di = sp.Trace(
    coefficient(kinetic_terms, {phi_j: 1, phi_k: 1, dphi_i: 1})
)
kin_ijk_dj = sp.Trace(
    coefficient(kinetic_terms, {phi_i: 1, phi_k: 1, dphi_j: 1})
)
kin_ijk_dk = sp.Trace(
    coefficient(kinetic_terms, {phi_i: 1, phi_j: 1, dphi_k: 1})
)
break_ijk = sp.Trace(
    coefficient(breaking_terms, {phi_i: 1, phi_j: 1, phi_k: 1})
)

show("Input Phi", Phi)
show("Input partial^mu Phi", dPhi)
show("Three-identical kinetic generator coefficient", kin_iii)
show("Three-identical breaking generator coefficient", break_iii)
show("Two-identical coefficient multiplying phi_i^2 partial^mu phi_j", kin_iij_dj)
show(
    "Two-identical coefficient multiplying phi_i phi_j partial^mu phi_i",
    kin_iij_di,
)
show("Two-identical breaking generator coefficient", break_iij)
show("Three-distinct coefficient multiplying phi_j phi_k partial^mu phi_i", kin_ijk_di)
show("Three-distinct coefficient multiplying phi_i phi_k partial^mu phi_j", kin_ijk_dj)
show("Three-distinct coefficient multiplying phi_i phi_j partial^mu phi_k", kin_ijk_dk)
show("Three-distinct breaking generator coefficient", break_ijk)

# The extracted expressions must equal the compact ordered-generator forms.
K_iij = sp.Trace(cA * commutator(lambda_i, commutator(lambda_i, lambda_j)))
S_iij = sp.Trace(
    (lambda_i**2 * lambda_j + lambda_i * lambda_j * lambda_i + lambda_j * lambda_i**2)
    * P
)

def K_distinct(a: sp.MatrixExpr, b: sp.MatrixExpr, c: sp.MatrixExpr) -> sp.Expr:
    """Coefficient for phi_a phi_b partial(phi_c)."""
    return sp.Trace(
        cA * (commutator(a, commutator(b, c)) + commutator(b, commutator(a, c)))
    )


S_ijk = sp.Trace(
    (
        lambda_i * lambda_j * lambda_k
        + lambda_i * lambda_k * lambda_j
        + lambda_j * lambda_i * lambda_k
        + lambda_j * lambda_k * lambda_i
        + lambda_k * lambda_i * lambda_j
        + lambda_k * lambda_j * lambda_i
    )
    * P
)

K = 1 / (6 * f_a * f_pi)
M = -B_0 / (6 * f_a * f_pi)
L_iii = M * a * phi_i**3 * break_iii
L_iij = (
    K
    * da
    * K_iij
    * (phi_i**2 * dphi_j - phi_i * phi_j * dphi_i)
    + M * a * phi_i**2 * phi_j * S_iij
)
L_ijk = K * da * (
    phi_i * phi_j * dphi_k * K_distinct(lambda_i, lambda_j, lambda_k)
    + phi_i * phi_k * dphi_j * K_distinct(lambda_i, lambda_k, lambda_j)
    + phi_j * phi_k * dphi_i * K_distinct(lambda_j, lambda_k, lambda_i)
) + M * a * phi_i * phi_j * phi_k * S_ijk

show("Final Lagrangian for three identical mesons", L_iii)
show("Final Lagrangian for two identical and one different meson", L_iij)
show("Final Lagrangian for three different mesons", L_ijk)

# Deterministic matrices verify the compact forms without imposing relations
# on the abstract generators in the derivation.
matrix_values = {
    lambda_i: sp.Matrix([[1, 2], [3, 5]]),
    lambda_j: sp.Matrix([[7, 11], [13, 17]]),
    lambda_k: sp.Matrix([[19, 23], [29, 31]]),
    cA: sp.Matrix([[37, 41], [43, 47]]),
    P: sp.Matrix([[53, 59], [61, 67]]),
}


def numerical_value(expression: sp.Expr) -> sp.Expr:
    return expression.doit().subs(matrix_values).doit()


assert numerical_value(kin_iii) == 0
assert sp.simplify(numerical_value(kin_iij_dj - K_iij)) == 0
assert sp.simplify(numerical_value(kin_iij_di + K_iij)) == 0
assert sp.simplify(numerical_value(break_iij - S_iij)) == 0
assert sp.simplify(
    numerical_value(kin_ijk_di - K_distinct(lambda_j, lambda_k, lambda_i))
) == 0
assert sp.simplify(
    numerical_value(kin_ijk_dj - K_distinct(lambda_i, lambda_k, lambda_j))
) == 0
assert sp.simplify(
    numerical_value(kin_ijk_dk - K_distinct(lambda_i, lambda_j, lambda_k))
) == 0
assert sp.simplify(numerical_value(break_ijk - S_ijk)) == 0
