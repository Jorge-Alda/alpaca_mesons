#!/usr/bin/env python3
"""Derive the one-ALP photon coefficient in Alda conventions.

Run from the repository root with:
    .venv/bin/python scripts/derive_alp_photon_invariant.py

All Wilson coefficients are dimensionless.  The physical photon operator is
alpha_em/(4*pi) * C_agammagamma * (a/f_a) * F Ftilde.
"""

import sympy as sp


def show(label: str, expression: sp.Expr | sp.MatrixBase) -> None:
    """Print a deterministic labelled LaTeX expression."""
    print(f"{label}:\n{sp.latex(expression)}\n")


def matrix_symbols(prefix: str) -> sp.Matrix:
    """Return a fully general 3-by-3 flavour matrix."""
    entries = sp.symbols(
        " ".join(f"{prefix}_{row}{column}" for row in "uds" for column in "uds")
    )
    return sp.Matrix(3, 3, entries)


# Conventions and assumptions: N_f=3 U(3), tr(lambda_i lambda_j)=2 delta_ij;
# Q is diagonal, but c_q^A, M_a^{-1}(C), and kappa remain fully general.
I = sp.I
sqrt = sp.sqrt
zero = sp.Integer(0)
one = sp.Integer(1)

generators = (
    sqrt(sp.Rational(2, 3)) * sp.eye(3),
    sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]]),
    sp.Matrix([[0, -I, 0], [I, 0, 0], [0, 0, 0]]),
    sp.diag(1, -1, 0),
    sp.Matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]]),
    sp.Matrix([[0, 0, -I], [0, 0, 0], [I, 0, 0]]),
    sp.Matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]]),
    sp.Matrix([[0, 0, 0], [0, 0, -I], [0, I, 0]]),
    sp.diag(1, 1, -2) / sqrt(3),
)

c_q_A = matrix_symbols("cA")
m_a_inverse_C = matrix_symbols(r"\mathbb{M}_a^{-1}(C)")
kappa = matrix_symbols("kappa")
q_u, q_d, q_s = sp.symbols("q_u q_d q_s", real=True)
Q = sp.diag(q_u, q_d, q_s)
Q2 = Q**2
N_c, f_pi, f_a, c_gamma = sp.symbols(
    "N_c f_pi f_a c_gamma", nonzero=True
)


def projection(matrix: sp.Matrix) -> sp.Matrix:
    """Apply U(3) generator completeness to a flavour matrix."""
    return sp.simplify(
        sum(
            (sp.trace(matrix * generator) * generator for generator in generators),
            sp.zeros(3),
        )
        / 2
    )


# Alda Eq. (III.28): the spectral coefficients of M_a^{-1}(C) are
# tr[C lambda_i]/(m_i^2-m_a^2).  Thus Eq. (III.25) can be written as below.
completeness_c_q_A = projection(c_q_A)
completeness_m_a_inverse_C = projection(m_a_inverse_C)
w = tuple(
    sp.simplify(
        f_pi
        / f_a
        * sp.trace((c_q_A + m_a_inverse_C) * generator)
    )
    for generator in generators
)

# No-vector WZW meson-photon vertex:
# L_WZW = alpha_em/(4*pi) * (N_c/f_pi) tr(Q^2 Phi) F Ftilde.
# Replacing Phi_i by its physical-ALP component w_i a gives this dimensionless
# contribution to the coefficient multiplying a/f_a.
delta_wzw_mixing = sp.simplify(
    N_c
    * f_a
    / f_pi
    * sum(
        w_i * sp.trace(Q2 * generator)
        for w_i, generator in zip(w, generators)
    )
)
delta_expected = sp.simplify(
    2 * N_c * sp.trace(Q2 * (c_q_A + m_a_inverse_C))
)

# Derive the invariant under c_q^A -> c_q^A+kappa and
# c_gamma -> c_gamma-2*N_c*tr(Q^2*kappa).
c_gamma_inv = sp.simplify(c_gamma + 2 * N_c * sp.trace(Q2 * c_q_A))
c_gamma_inv_transformed = sp.simplify(
    c_gamma
    - 2 * N_c * sp.trace(Q2 * kappa)
    + 2 * N_c * sp.trace(Q2 * (c_q_A + kappa))
)
full_coefficient = sp.simplify(c_gamma + delta_wzw_mixing)
full_from_invariant = sp.simplify(
    c_gamma_inv + 2 * N_c * sp.trace(Q2 * m_a_inverse_C)
)

# Alda no-mixing choice kappa=-c_q^A-M_a^{-1}(C).
kappa_no_mixing = -(c_q_A + m_a_inverse_C)
hat_c_q_A = sp.simplify(c_q_A + kappa_no_mixing)
hat_c_gamma = sp.simplify(
    c_gamma - 2 * N_c * sp.trace(Q2 * kappa_no_mixing)
)

# N_c=3 and Q=diag(2/3,-1/3,-1/3).  General off-diagonal matrix entries
# remain present in the input and are eliminated only by the diagonal Q^2 trace.
three_flavour_substitution = {
    N_c: 3,
    q_u: sp.Rational(2, 3),
    q_d: -sp.Rational(1, 3),
    q_s: -sp.Rational(1, 3),
}
specialized_invariant = sp.simplify(c_gamma_inv.subs(three_flavour_substitution))
specialized_full = sp.simplify(full_coefficient.subs(three_flavour_substitution))

print("Conventions: q -> exp(+i a kappa gamma5/f_a) q for delta=0.")
print("Approximation: one ALP at O(1/f_a), leading no-vector WZW sector.")
print("Domain: M_a is invertible (away from ALP-meson degeneracies).\n")
show("U(3) completeness applied to c_q^A", completeness_c_q_A)
show("U(3) completeness applied to M_a^{-1}(C)", completeness_m_a_inverse_C)
show("Alda physical-ALP meson components w_i", sp.Matrix(w))
show("WZW/mixing photon contribution", delta_wzw_mixing)
show("Invariant explicit photon combination c_gamma^inv", c_gamma_inv)
show("Full one-ALP coefficient C_agammagamma(m_a)", full_coefficient)
show("No-mixing kappa", kappa_no_mixing)
show("No-mixing transformed hat c_q^A", hat_c_q_A)
show("No-mixing transformed hat c_gamma", hat_c_gamma)
show("N_c=3 diagonal-Q invariant", specialized_invariant)
show("N_c=3 diagonal-Q full coefficient", specialized_full)

# Diagnostic assertions are deliberately separate from the derivation output.
assert all(
    sp.simplify((completeness_c_q_A - c_q_A)[i, j]) == zero
    for i in range(3)
    for j in range(3)
)
assert all(
    sp.simplify((completeness_m_a_inverse_C - m_a_inverse_C)[i, j]) == zero
    for i in range(3)
    for j in range(3)
)
assert sp.simplify(delta_wzw_mixing - delta_expected) == zero
assert sp.simplify(c_gamma_inv_transformed - c_gamma_inv) == zero
assert sp.simplify(full_coefficient - full_from_invariant) == zero
assert all(
    sp.simplify((hat_c_q_A + m_a_inverse_C)[i, j]) == zero
    for i in range(3)
    for j in range(3)
)
assert sp.simplify(hat_c_gamma - full_coefficient) == zero

off_diagonal_symbols = [
    c_q_A[i, j]
    for i in range(3)
    for j in range(3)
    if i != j
] + [
    m_a_inverse_C[i, j]
    for i in range(3)
    for j in range(3)
    if i != j
]
assert not (specialized_full.free_symbols & set(off_diagonal_symbols))
print("All symbolic assertions passed.")
