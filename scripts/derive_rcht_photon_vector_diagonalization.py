#!/usr/bin/env python3
"""Derive photon--neutral-vector diagonalization in antisymmetric-tensor RChT.

Invocation:
    .venv/bin/python scripts/derive_rcht_photon_vector_diagonalization.py

Natural units are used.  The three neutral flavor generators are arbitrary
orthonormal diagonal U(3) generators; no ideal-mixing assignment is made.
"""

import sympy as sp
try:
    from scripts.symbolic_utils import show
except ModuleNotFoundError:  # direct ``python scripts/...py`` invocation
    from symbolic_utils import show



# ---------------------------------------------------------------------------
# 1. Generic diagonal U(3) generators and electromagnetic overlaps
# ---------------------------------------------------------------------------
alpha, beta, gamma = sp.symbols("alpha beta gamma", real=True)

R12 = sp.Matrix(
    [
        [sp.cos(alpha), -sp.sin(alpha), 0],
        [sp.sin(alpha), sp.cos(alpha), 0],
        [0, 0, 1],
    ]
)
R13 = sp.Matrix(
    [
        [sp.cos(beta), 0, -sp.sin(beta)],
        [0, 1, 0],
        [sp.sin(beta), 0, sp.cos(beta)],
    ]
)
R23 = sp.Matrix(
    [
        [1, 0, 0],
        [0, sp.cos(gamma), -sp.sin(gamma)],
        [0, sp.sin(gamma), sp.cos(gamma)],
    ]
)

# The columns are the diagonal entries of T_i=diag(t_ui,t_di,t_si).
T = sp.simplify(R12 * R13 * R23)
q = sp.Matrix([sp.Rational(2, 3), -sp.Rational(1, 3), -sp.Rational(1, 3)])
c = sp.simplify(T.T * q)
c_norm = sp.trigsimp((c.T * c)[0])
c_singlet = sp.simplify((q.T * sp.ones(3, 1) / sp.sqrt(3))[0])

show("Generic orthogonal diagonal-generator matrix T", T)
show("Electromagnetic overlaps c_i=(T^T q)_i", c)
show("Basis invariant sum_i c_i^2", c_norm)
show("Electromagnetic overlap of the normalized U(3) singlet", c_singlet)

assert sp.simplify(T.T * T - sp.eye(3)) == sp.zeros(3)
assert c_norm == sp.Rational(2, 3)
assert c_singlet == 0


# ---------------------------------------------------------------------------
# 2. First-order elimination of one antisymmetric tensor
# ---------------------------------------------------------------------------
# Lower-case symbols stand for Lorentz-contracted commuting coefficients in
# the quadratic form: v=V_{mu nu}, h=P_{mu nu}, f=F_{mu nu}.  SymPy performs
# the Gaussian elimination; the Lorentz structures are unchanged.
M, g = sp.symbols("M g", positive=True, real=True)
v, h, f, p2 = sp.symbols("v h f p2", real=True)
L_first = sp.Rational(1, 2) * p2 - sp.Rational(1, 2) * h * v + M**2 * v**2 / 4 - g * v * f / 2
v_solution = sp.solve(sp.diff(L_first, v), v)[0]
L_eliminated = sp.factor(L_first.subs(v, v_solution))

show("First-order quadratic Lagrangian for one tensor mode", L_first)
show("Antisymmetric-tensor equation of motion", v_solution)
show("Lagrangian after eliminating the tensor", L_eliminated)

R2, Rf = sp.symbols("R2 Rf", real=True)
epsilon = sp.symbols("epsilon", real=True)
L_proca = sp.expand(
    L_eliminated.subs({p2: M**2 * R2, h: M * Rf, g: epsilon * M})
)
show("Exact Proca form after P_mu=M R_mu and epsilon=g/M", L_proca)

expected_proca = M**2 * R2 / 2 - Rf**2 / 4 - epsilon * Rf * f / 2 - epsilon**2 * f**2 / 4
assert sp.simplify(L_proca - expected_proca) == 0


# ---------------------------------------------------------------------------
# 3. Exact photon decoupling for three neutral vectors
# ---------------------------------------------------------------------------
eps1, eps2, eps3 = sp.symbols("epsilon_1 epsilon_2 epsilon_3", real=True)
eps = sp.Matrix([eps1, eps2, eps3])
eps_sq = (eps.T * eps)[0]
Z = 1 + eps_sq

# K is the gauge-field kinetic matrix in X=(A,R_1,R_2,R_3).
K = sp.eye(4)
K[0, 0] = Z
for i in range(3):
    K[0, i + 1] = eps[i]
    K[i + 1, 0] = eps[i]

# X=S Y, Y=(A_gamma,R_1,R_2,R_3), implements
# A=A_gamma/sqrt(Z)-epsilon^T R/Z.
S = sp.eye(4)
S[0, 0] = 1 / sp.sqrt(Z)
for i in range(3):
    S[0, i + 1] = -eps[i] / Z

K_transformed = sp.simplify(S.T * K * S)
K_V = sp.simplify(sp.eye(3) - eps * eps.T / Z)
expected_K_transformed = sp.diag(1, 1, 1, 1)
expected_K_transformed[1:4, 1:4] = K_V

show("Original photon-vector kinetic matrix K", K)
show("Field transformation X=S Y", S)
show("Decoupled kinetic matrix S^T K S", K_transformed)
show("Residual vector kinetic matrix K_V", K_V)

assert sp.simplify(K_transformed - expected_K_transformed) == sp.zeros(4)

# K_V^{-1/2} is explicit because K_V is a rank-one deformation of I.
K_inv_half = sp.simplify(sp.eye(3) + eps * eps.T / (sp.sqrt(Z) + 1))
canonical_check = sp.simplify(K_inv_half * K_V * K_inv_half)
show("Positive inverse square root K_V^{-1/2}", K_inv_half)
show("Canonical-normalization check K_V^{-1/2} K_V K_V^{-1/2}", canonical_check)
assert canonical_check == sp.eye(3)


# ---------------------------------------------------------------------------
# 4. Massive poles and the exact massless photon
# ---------------------------------------------------------------------------
s = sp.symbols("s", real=True)
M1, M2, M3 = sp.symbols("M_1 M_2 M_3", positive=True, real=True)
D = sp.diag(M1**2, M2**2, M3**2)
massive_determinant = sp.factor((s * K_V - D).det())
unmixed_product = (s - M1**2) * (s - M2**2) * (s - M3**2)
pole_bracket = sp.factor(massive_determinant / unmixed_product)

pole_equation_form_1 = sp.simplify(
    Z
    - s
    * (
        eps1**2 / (s - M1**2)
        + eps2**2 / (s - M2**2)
        + eps3**2 / (s - M3**2)
    )
)
pole_equation_form_2 = sp.simplify(
    1
    - eps1**2 * M1**2 / (s - M1**2)
    - eps2**2 * M2**2 / (s - M2**2)
    - eps3**2 * M3**2 / (s - M3**2)
)

show("Massive generalized characteristic determinant det(s K_V-D)", massive_determinant)
show("Matrix-determinant-lemma bracket", pole_bracket)
show("Generic rational massive-pole equation", pole_equation_form_2)
show("Transverse characteristic polynomial including the photon factor", sp.factor(s * massive_determinant))

decoupled_check = sp.factor(massive_determinant.subs(eps3, 0))
show("Decoupled-direction check at epsilon_3=0", decoupled_check)

assert sp.simplify(pole_equation_form_1 - pole_equation_form_2) == 0
assert sp.simplify(pole_bracket - pole_equation_form_1 / Z) == 0
assert sp.simplify(decoupled_check.subs(s, M3**2)) == 0


# ---------------------------------------------------------------------------
# 5. Electromagnetic normalization and leading-e expansion
# ---------------------------------------------------------------------------
e, F_V = sp.symbols("e F_V", real=True)
c1, c2, c3 = sp.symbols("c_1 c_2 c_3", real=True)
g_vector = sp.sqrt(2) * e * F_V * sp.Matrix([c1, c2, c3])
epsilon_vector = sp.diag(1 / M1, 1 / M2, 1 / M3) * g_vector
show("Tensor mixing coefficients g_i=sqrt(2) e F_V c_i", g_vector)
show("Dimensionless kinetic mixings epsilon_i=g_i/M_i", epsilon_vector)
show("Tree-level canonical photon charge e_gamma", e / sp.sqrt(Z))

eta = sp.symbols("eta", real=True)
leading_A_gamma = sp.Matrix(
    [
        sp.series(sp.sqrt(1 + eta**2 * eps_sq), eta, 0, 3).removeO(),
        eta * eps1,
        eta * eps2,
        eta * eps3,
    ]
)
show(
    "Small-mixing content of A_gamma (normalization, R_1, R_2, R_3)",
    leading_A_gamma,
)

print("Diagnostic assertions: all passed.")
print("Conventions: natural units; real orthonormal diagonal U(3) generators; bare e.")
