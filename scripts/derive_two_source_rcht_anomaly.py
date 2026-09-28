#!/usr/bin/env python3
"""Derive two-source ALP anomaly and odd-RChT insertions.

Invocation from the repository root:
    .venv/bin/python scripts/derive_two_source_rcht_anomaly.py

Natural units are used.  Flavor coefficients and vector-source couplings are
dimensionless; B_0, m_q, derivatives, and resonance masses carry mass dimension
one.  The calculation is exact at linear order in a/f_a and at tree level.
"""

import sympy as sp


def show(label: str, expression: sp.Expr | sp.MatrixBase) -> None:
    """Print a deterministic labelled LaTeX expression."""
    print(f"{label}:\n{sp.latex(expression)}\n")


def matrix_symbols(prefix: str) -> sp.Matrix:
    """Return a fully general three-flavor matrix."""
    entries = sp.symbols(
        " ".join(f"{prefix}_{row}{column}" for row in "uds" for column in "uds")
    )
    return sp.Matrix(3, 3, entries)


I = sp.I
zero = sp.Integer(0)
generators = (
    sp.sqrt(sp.Rational(2, 3)) * sp.eye(3),
    sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]]),
    sp.Matrix([[0, -I, 0], [I, 0, 0], [0, 0, 0]]),
    sp.diag(1, -1, 0),
    sp.Matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]]),
    sp.Matrix([[0, 0, -I], [0, 0, 0], [I, 0, 0]]),
    sp.Matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]]),
    sp.Matrix([[0, 0, 0], [0, 0, -I], [0, I, 0]]),
    sp.diag(1, 1, -2) / sp.sqrt(3),
)


def projection(matrix: sp.Matrix) -> sp.Matrix:
    """Apply tr(lambda_i lambda_j)=2 delta_ij completeness."""
    return sp.simplify(
        sum(
            (sp.trace(matrix * generator) * generator for generator in generators),
            sp.zeros(3),
        )
        / 2
    )


# ---------------------------------------------------------------------------
# 1. Symmetric two-source anomaly seed and rotation invariance
# ---------------------------------------------------------------------------
T_A = matrix_symbols("TA")
T_B = matrix_symbols("TB")
c_q_A = matrix_symbols("cA")
m_a_inverse_C = matrix_symbols("m_a_inverse_C")
kappa = matrix_symbols("kappa")
N_c, c_AB, c_gamma = sp.symbols("N_c c_AB c_gamma")
delta_AB = sp.Symbol("delta_AB", integer=True, nonnegative=True)
s_AB = sp.Rational(2, 1) / (1 + delta_AB)

anticommutator_AB = T_A * T_B + T_B * T_A
seed_AB = sp.expand(
    c_AB
    + s_AB * N_c * sp.trace(anticommutator_AB * (c_q_A + m_a_inverse_C))
)
seed_AB_rotated = sp.expand(
    c_AB
    - s_AB * N_c * sp.trace(anticommutator_AB * kappa)
    + s_AB
    * N_c
    * sp.trace(anticommutator_AB * (c_q_A + kappa + m_a_inverse_C))
)
seed_invariance_remainder = sp.simplify(seed_AB_rotated - seed_AB)

q_u, q_d, q_s = sp.symbols("q_u q_d q_s", real=True)
Q = sp.diag(q_u, q_d, q_s)
photon_seed = sp.simplify(
    seed_AB.subs({c_AB: c_gamma, delta_AB: 1})
    .subs(dict(zip(T_A, Q)))
    .subs(dict(zip(T_B, Q)))
)
photon_expected = sp.simplify(
    c_gamma + 2 * N_c * sp.trace(Q**2 * (c_q_A + m_a_inverse_C))
)
distinct_seed = sp.simplify(seed_AB.subs(delta_AB, 0))
distinct_expected = sp.simplify(
    c_AB
    + 2 * N_c * sp.trace(anticommutator_AB * (c_q_A + m_a_inverse_C))
)

# Compact abstract expressions are printed; the component matrices above make
# SymPy perform and verify the identities without assuming flavor commutativity.
TA_abs = sp.MatrixSymbol("T_A", 3, 3)
TB_abs = sp.MatrixSymbol("T_B", 3, 3)
cA_abs = sp.MatrixSymbol("c_q^A", 3, 3)
m_a_inverse_C_abs = sp.MatrixSymbol(r"\mathbb M_a^{-1}(C)", 3, 3)
kappa_abs = sp.MatrixSymbol("kappa", 3, 3)
anti_abs = TA_abs * TB_abs + TB_abs * TA_abs
seed_abs = c_AB + s_AB * N_c * sp.Trace(
    anti_abs * (cA_abs + m_a_inverse_C_abs)
)
seed_rotated_abs = (
    c_AB
    - s_AB * N_c * sp.Trace(anti_abs * kappa_abs)
    + s_AB
    * N_c
    * sp.Trace(anti_abs * (cA_abs + kappa_abs + m_a_inverse_C_abs))
)
show("Unordered-pair multiplicity s_AB", s_AB)
show("Two-source symmetric anomaly seed C_seed_AB", seed_abs)
show("Seed after c_q^A and c_AB rotation", seed_rotated_abs)
show("Rotation-invariance remainder", seed_invariance_remainder)
Q_abs = sp.MatrixSymbol("Q", 3, 3)
show(
    "Equal-source photon reduction",
    c_gamma + 2 * N_c * sp.Trace(Q_abs**2 * (cA_abs + m_a_inverse_C_abs)),
)
show(
    "Distinct-source seed combining both ordered cross terms",
    c_AB + 2 * N_c * sp.Trace(anti_abs * (cA_abs + m_a_inverse_C_abs)),
)


# ---------------------------------------------------------------------------
# 2. Completeness derivation of physical-ALP CCWZ insertions
# ---------------------------------------------------------------------------
f_pi, f_a, B_0 = sp.symbols("f_pi f_a B_0", nonzero=True)
c_q_P = matrix_symbols("cP")
c_q_AP = c_q_A - c_q_P  # Alda convention c_q^{AP}=c_q^A-c_q^P.
m_q = sp.diag(*sp.symbols("m_u m_d m_s", real=True))
c_q_A_plus_m_a_inverse_C = c_q_A + m_a_inverse_C

w = tuple(
    sp.simplify(
        f_pi / f_a * sp.trace(c_q_A_plus_m_a_inverse_C * generator)
    )
    for generator in generators
)
Phi_physical_over_a = sp.simplify(
    sum((w_i * generator for w_i, generator in zip(w, generators)), sp.zeros(3))
)

# Coefficients below multiply partial_mu a/f_a or a/f_a, respectively.
u_direct = 2 * c_q_A
u_meson = sp.simplify(-f_a * Phi_physical_over_a / f_pi)
u_physical = sp.simplify(u_direct + u_meson)

# At u=1, chi_-=2 i chi_P.  Expanding u=exp(i Phi/(2 f_pi))
# also gives chi_-|meson=-2 i B_0 {Phi,m_q}/f_pi.
chi_direct = 4 * I * B_0 * (c_q_P * m_q + m_q * c_q_P)
chi_meson = sp.simplify(
    -2 * I * B_0 * f_a / f_pi
    * (Phi_physical_over_a * m_q + m_q * Phi_physical_over_a)
)
chi_physical = sp.simplify(chi_direct + chi_meson)
chi_expected = sp.simplify(
    -4 * I * B_0
    * (
        (c_q_AP + m_a_inverse_C) * m_q
        + m_q * (c_q_AP + m_a_inverse_C)
    )
)

cP_abs = sp.MatrixSymbol("c_q^P", 3, 3)
cAP_abs = cA_abs - cP_abs
mq_abs = sp.MatrixSymbol("m_q", 3, 3)
show(
    "Physical meson matrix Phi/a from LSZ completeness",
    2 * f_pi / f_a * (cA_abs + m_a_inverse_C_abs),
)
show("Direct u_mu coefficient multiplying partial_mu a/f_a", 2 * cA_abs)
show(
    "Meson u_mu coefficient multiplying partial_mu a/f_a",
    -2 * (cA_abs + m_a_inverse_C_abs),
)
show(
    "Physical u_mu coefficient multiplying partial_mu a/f_a",
    -2 * m_a_inverse_C_abs,
)
show(
    "Direct chi_- coefficient multiplying a/f_a",
    4 * I * B_0 * (cP_abs * mq_abs + mq_abs * cP_abs),
)
show(
    "Meson chi_- coefficient multiplying a/f_a",
    -4 * I * B_0
    * (
        (cA_abs + m_a_inverse_C_abs) * mq_abs
        + mq_abs * (cA_abs + m_a_inverse_C_abs)
    ),
)
show(
    "Physical chi_- coefficient multiplying a/f_a",
    -4 * I * B_0
    * (
        (cAP_abs + m_a_inverse_C_abs) * mq_abs
        + mq_abs * (cAP_abs + m_a_inverse_C_abs)
    ),
)


# The flavor factors multiplying the distinct d_i Lorentz structures.
F_m_a_inverse_C = sp.expand(
    sp.trace(anticommutator_AB * m_a_inverse_C)
)
F_chi = sp.expand(
    sp.trace(
        anticommutator_AB
        * (
            (c_q_AP + m_a_inverse_C) * m_q
            + m_q * (c_q_AP + m_a_inverse_C)
        )
    )
)
show(
    "Flavor factor for d_1, d_3, and d_4",
    sp.Trace(anti_abs * m_a_inverse_C_abs),
)
show(
    "Flavor factor for d_2",
    sp.Trace(
        anti_abs
        * (
            (cAP_abs + m_a_inverse_C_abs) * mq_abs
            + mq_abs * (cAP_abs + m_a_inverse_C_abs)
        )
    ),
)


# ---------------------------------------------------------------------------
# 3. Exact anomalous-bilinear congruence under photon-vector diagonalization
# ---------------------------------------------------------------------------
n = 3
Z = sp.Symbol("Z", positive=True)
lambda_gamma = sp.Symbol("lambda_gamma_gamma")
ell = sp.MatrixSymbol("ell", n, 1)
Lambda_V = sp.MatrixSymbol("Lambda_V", n, n)
epsilon = sp.MatrixSymbol("epsilon", n, 1)
rho = sp.MatrixSymbol("rho", n, 1)
H = sp.MatrixSymbol("H", n, n)

Lambda_old = sp.BlockMatrix(
    [[sp.Matrix([[lambda_gamma]]), ell.T], [ell, Lambda_V]]
)
S = sp.BlockMatrix(
    [
        [sp.Matrix([[1 / sp.sqrt(Z)]]), -rho.T],
        [epsilon / sp.sqrt(Z), H],
    ]
)
Lambda_phys = sp.block_collapse(S.T * Lambda_old * S)

s_gamma = sp.BlockMatrix([[sp.Matrix([[1 / sp.sqrt(Z)]])], [epsilon / sp.sqrt(Z)]])
s_W = sp.BlockMatrix([[-rho.T], [H]])
block_gg = sp.block_collapse(s_gamma.T * Lambda_old * s_gamma)
block_gW = sp.block_collapse(s_gamma.T * Lambda_old * s_W)
block_WW = sp.block_collapse(s_W.T * Lambda_old * s_W)

expected_gg_scalar = (
    lambda_gamma + 2 * (ell.T * epsilon)[0, 0]
    + (epsilon.T * Lambda_V * epsilon)[0, 0]
) / Z
expected_gW = (
    (ell.T + epsilon.T * Lambda_V) * H
    - (lambda_gamma + (epsilon.T * ell)[0, 0]) * rho.T
) / sp.sqrt(Z)
expected_WW = (
    H.T * Lambda_V * H
    - rho * ell.T * H
    - H.T * ell * rho.T
    + lambda_gamma * rho * rho.T
)
Lambda_phys_compact = sp.BlockMatrix(
    [[sp.Matrix([[expected_gg_scalar]]), expected_gW], [expected_gW.T, expected_WW]]
)

print(
    "Lambda convention: the anomalous bilinear is "
    "(1/2) X_old^T Lambda_old Xtilde_old; Lambda is independent of the "
    "earlier unordered-pair coefficient convention.\n"
)
show("Old-to-physical field-strength transformation S", S)
show("General symmetric old anomalous matrix Lambda_old", Lambda_old)
show("Derived Lambda_phys=S^T Lambda_old S", Lambda_phys_compact)
show("Physical gamma-gamma block", expected_gg_scalar)
show("Physical gamma-W block", expected_gW)
show("Physical W-W block", expected_WW)


# Diagnostic assertions are separate from the derivation output.
assert seed_invariance_remainder == zero
assert sp.simplify(photon_seed - photon_expected) == zero
assert sp.simplify(distinct_seed - distinct_expected) == zero
assert sp.simplify(s_AB.subs(delta_AB, 1) - 1) == zero
assert sp.simplify(s_AB.subs(delta_AB, 0) - 2) == zero
assert all(
    sp.simplify(
        (
            projection(c_q_A_plus_m_a_inverse_C)
            - c_q_A_plus_m_a_inverse_C
        )[i, j]
    )
    == zero
    for i in range(3)
    for j in range(3)
)
assert all(
    sp.simplify((u_physical + 2 * m_a_inverse_C)[i, j]) == zero
    for i in range(3)
    for j in range(3)
)
assert all(
    sp.simplify((chi_physical - chi_expected)[i, j]) == zero
    for i in range(3)
    for j in range(3)
)
# Independent component-level check of the block formulas for two vector modes.
l1, l2, e1, e2, r1, r2 = sp.symbols("l1 l2 e1 e2 r1 r2")
h11, h12, h21, h22 = sp.symbols("h11 h12 h21 h22")
v11, v12, v22 = sp.symbols("v11 v12 v22")
ell_test = sp.Matrix([l1, l2])
eps_test = sp.Matrix([e1, e2])
rho_test = sp.Matrix([r1, r2])
H_test = sp.Matrix([[h11, h12], [h21, h22]])
LV_test = sp.Matrix([[v11, v12], [v12, v22]])
S_test = sp.BlockMatrix(
    [
        [sp.Matrix([[1 / sp.sqrt(Z)]]), -rho_test.T],
        [eps_test / sp.sqrt(Z), H_test],
    ]
).as_explicit()
Lambda_test = sp.BlockMatrix(
    [[sp.Matrix([[lambda_gamma]]), ell_test.T], [ell_test, LV_test]]
).as_explicit()
expected_test = sp.BlockMatrix(
    [
        [
            sp.Matrix(
                [[
                    (
                        lambda_gamma
                        + 2 * (ell_test.T * eps_test)[0]
                        + (eps_test.T * LV_test * eps_test)[0]
                    )
                    / Z
                ]]
            ),
            (
                (ell_test.T + eps_test.T * LV_test) * H_test
                - (lambda_gamma + (eps_test.T * ell_test)[0]) * rho_test.T
            )
            / sp.sqrt(Z),
        ],
        [
            (
                (ell_test.T + eps_test.T * LV_test) * H_test
                - (lambda_gamma + (eps_test.T * ell_test)[0]) * rho_test.T
            ).T
            / sp.sqrt(Z),
            H_test.T * LV_test * H_test
            - rho_test * ell_test.T * H_test
            - H_test.T * ell_test * rho_test.T
            + lambda_gamma * rho_test * rho_test.T,
        ],
    ]
).as_explicit()
assert sp.simplify(S_test.T * Lambda_test * S_test - expected_test) == sp.zeros(3)
print("All symbolic assertions passed.")
