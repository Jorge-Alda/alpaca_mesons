"""Small shared helpers for the executable SymPy derivations.

Keep physics-specific conventions and expansions in their derivation scripts;
this module contains only algebraic operations that are identical across them.
"""

from __future__ import annotations

from itertools import permutations
from typing import Any

import sympy as sp


def show(label: str, expression: Any) -> None:
    """Print a labelled expression in deterministic LaTeX form."""
    print(f"{label}:\n{sp.latex(expression)}\n")


def commutator(left: Any, right: Any) -> Any:
    """Return the associative-algebra commutator ``[left, right]``."""
    return left * right - right * left


def anticommutator(left: Any, right: Any) -> Any:
    """Return the associative-algebra anticommutator ``{left, right}``."""
    return left * right + right * left


def map_power(linear_map, argument: Any, field: Any, order: int) -> Any:
    """Apply a linear map repeatedly, with order zero returning argument."""
    if order < 0:
        raise ValueError("order must be nonnegative")
    result = argument
    for _ in range(order):
        result = linear_map(result, field)
    return result


def scalar_coefficient(expression: Any, symbol: sp.Symbol, order: int) -> Any:
    """Extract the coefficient of ``symbol**order`` by differentiation."""
    if order < 0:
        raise ValueError("order must be nonnegative")
    derivative = expression
    for _ in range(order):
        derivative = sp.diff(derivative, symbol)
    return derivative.subs(symbol, 0) / sp.factorial(order)


def matrix_unit(row: int, column: int, dimension: int = 3) -> sp.ImmutableMatrix:
    """Return the flavor matrix unit with one at ``(row, column)``."""
    if not (0 <= row < dimension and 0 <= column < dimension):
        raise IndexError("matrix-unit indices must lie within the matrix dimension")
    result = sp.zeros(dimension)
    result[row, column] = 1
    return sp.ImmutableMatrix(result)


def u3_generators() -> tuple[sp.ImmutableMatrix, ...]:
    """Return Hermitian U(3) generators normalized by tr(lambda_i lambda_j)=2 delta_ij.

    The ordering is (lambda_0, lambda_1, ..., lambda_8): singlet, then the
    conventional Gell-Mann ordering. This matches the meson-field convention
    used in the Alda derivations.
    """
    i = sp.I
    sqrt = sp.sqrt
    matrices = (
        sqrt(sp.Rational(2, 3)) * sp.eye(3),
        sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]]),
        sp.Matrix([[0, -i, 0], [i, 0, 0], [0, 0, 0]]),
        sp.diag(1, -1, 0),
        sp.Matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]]),
        sp.Matrix([[0, 0, -i], [0, 0, 0], [i, 0, 0]]),
        sp.Matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]]),
        sp.Matrix([[0, 0, 0], [0, 0, -i], [0, i, 0]]),
        sp.diag(1, 1, -2) / sqrt(3),
    )
    return tuple(sp.ImmutableMatrix(matrix) for matrix in matrices)


def u3_completeness_projection(matrix: sp.MatrixBase) -> sp.Matrix:
    """Project a 3x3 matrix with the normalized U(3) completeness relation."""
    if matrix.shape != (3, 3):
        raise ValueError("U(3) completeness projection requires a 3x3 matrix")
    basis = u3_generators()
    return sp.simplify(
        sum((sp.trace(matrix * generator) * generator for generator in basis), sp.zeros(3))
        / 2
    )


MINKOWSKI_PLUS_MINUS = (1, -1, -1, -1)


def minkowski_dot(left: Any, right: Any) -> Any:
    """Contract contravariant four-vectors with metric diag(+1,-1,-1,-1)."""
    if len(left) != 4 or len(right) != 4:
        raise ValueError("Minkowski dot product requires two four-vectors")
    return sp.expand(sum(sign * left[index] * right[index]
                         for index, sign in enumerate(MINKOWSKI_PLUS_MINUS)))


def lower_minkowski(vector: Any) -> tuple[Any, ...]:
    """Lower a contravariant four-vector index with diag(+1,-1,-1,-1)."""
    if len(vector) != 4:
        raise ValueError("index lowering requires a four-vector")
    return tuple(sign * component for sign, component in zip(MINKOWSKI_PLUS_MINUS, vector))


def levi_civita_contract(a: Any, b: Any, c: Any, d: Any) -> Any:
    """Return epsilon^{mu nu rho sigma} a_mu b_nu c_rho d_sigma, epsilon^0123=+1.

    Entries are contravariant vector components; all four indices are lowered
    before contraction. Thus epsilon_0123=-1 for metric (+---).
    """
    return legacy_epsilon_determinant(*(lower_minkowski(v) for v in (a, b, c, d)))


def legacy_epsilon_determinant(a: Any, b: Any, c: Any, d: Any) -> Any:
    """Legacy lower0123=+1 contraction of contravariant components.

    Kept only for reproducing historical source calculations. New physical
    contractions must use levi_civita_contract (upper0123=+1).
    """
    if any(len(vector) != 4 for vector in (a, b, c, d)):
        raise ValueError("Levi-Civita contraction requires four-vectors")
    vectors = (a, b, c, d)
    total = 0
    for p in permutations(range(4)):
        inversions = sum(
            p[left] > p[right]
            for left in range(4)
            for right in range(left + 1, 4)
        )
        total += (-1) ** inversions * vectors[0][p[0]] * vectors[1][p[1]] * vectors[2][p[2]] * vectors[3][p[3]]
    return sp.expand(total)
