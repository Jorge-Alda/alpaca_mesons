#!/usr/bin/env python3
"""Derive diagonalized LO RChT and Goldstone-current flavor structures.

Invocation:
    .venv/bin/python scripts/derive_rcht_gv_flavor_structures.py

The calculation uses natural units, U=exp(i*sqrt(2)*Phi/F), generic
noncommuting flavor generators, and terms linear in the ALP.
"""

import sympy as sp


def show(label: str, expression: sp.Basic) -> None:
    """Print one deterministic, LaTeX-formatted derivation step."""
    print(f"{label}:\n{sp.latex(expression)}\n")


def bookkeeping_coefficient(expression: sp.Expr, *parameters: sp.Symbol) -> sp.Expr:
    """Extract a multilinear bookkeeping coefficient with SymPy."""
    result = expression
    for parameter in parameters:
        result = sp.diff(result, parameter)
    return sp.expand(result.subs({parameter: 0 for parameter in parameters}))


def oriented_antisymmetric_coefficient(
    expression: sp.Expr,
    first_mu: sp.Symbol,
    second_nu: sp.Symbol,
    first_nu: sp.Symbol,
    second_mu: sp.Symbol,
) -> sp.Expr:
    """Reduce X_mn V^mn to the displayed (mu,nu) derivative ordering."""
    direct = sp.diff(expression, first_mu, second_nu)
    exchanged = sp.diff(expression, first_nu, second_mu)
    zeroes = {
        first_mu: 0,
        second_nu: 0,
        first_nu: 0,
        second_mu: 0,
    }
    # The minus sign follows from V_{nu mu}=-V_{mu nu} after relabeling.
    return sp.expand((direct - exchanged).subs(zeroes))


# ---------------------------------------------------------------------------
# 1. Inputs and the required CCWZ expansion
# ---------------------------------------------------------------------------
F, f_a, G_V, F_V, e = sp.symbols("F f_a G_V F_V e", nonzero=True, real=True)
I = sp.I

# Noncommutative Symbols make SymPy expand and extract coefficients without
# assigning components or a matrix dimension.  Tr is applied only after the
# noncommutative flavor word has been derived.
T_A = sp.Symbol("T_A", commutative=False)
lambda_i = sp.Symbol("lambda_i", commutative=False)
lambda_j = sp.Symbol("lambda_j", commutative=False)
cA = sp.Symbol("c_q^A", commutative=False)
cV = sp.Symbol("c_q^V", commutative=False)
Tr = sp.Function("Tr")
Phi = sp.Symbol("Phi", commutative=False)
V_munu = sp.Symbol("V_{mu nu}", commutative=False)
u_mu_input = sp.Symbol("u^mu", commutative=False)
u_nu_input = sp.Symbol("u^nu", commutative=False)

phi_i, t, x_i, x_j, a_order = sp.symbols("phi_i t x_i x_j a_order", real=True)
dphi_i_mu, dphi_i_nu = sp.symbols("partial_mu_phi_i partial_nu_phi_i", real=True)
dphi_j_mu, dphi_j_nu = sp.symbols("partial_mu_phi_j partial_nu_phi_j", real=True)
da_mu, da_nu = sp.symbols("partial_mu_a partial_nu_a", real=True)

# U=u^2=exp(i*sqrt(2)*Phi/F) implies u=exp(i*varphi),
# with varphi=Phi/(sqrt(2)*F).  These are the Taylor coefficients of the
# established all-orders ALP vielbein through the orders needed below.
varphi_i = phi_i * lambda_i / (sp.sqrt(2) * F)
u_phi_mu = -sp.sqrt(2) * dphi_i_mu * lambda_i / F
u_phi_nu = -sp.sqrt(2) * dphi_i_nu * lambda_i / F
u_aA_mu = 2 * da_mu * cA / f_a
u_aA_nu = 2 * da_nu * cA / f_a
u_aV_mu = -2 * I * da_mu * (cV * varphi_i - varphi_i * cV) / f_a
u_aV_nu = -2 * I * da_nu * (cV * varphi_i - varphi_i * cV) / f_a

show("Input U convention", sp.exp(I * sp.sqrt(2) * Phi / F))
show("Input LO G_V operator", I * G_V * Tr(V_munu * u_mu_input * u_nu_input) / sp.sqrt(2))
show("Input dimensionless field varphi=Phi/(sqrt(2) F)", varphi_i)
show("Mesonic vielbein at O(phi_i)", u_phi_mu)
show("Axial-ALP vielbein at O(a phi^0)", u_aA_mu)
show("Vector-ALP vielbein at O(a phi_i)", u_aV_mu)

# Antisymmetry of V_{mu nu} projects u^mu u^nu onto [u^mu,u^nu]/2.
u_mu_one_phi = t * u_phi_mu + a_order * u_aA_mu + t * a_order * u_aV_mu
u_nu_one_phi = t * u_phi_nu + a_order * u_aA_nu + t * a_order * u_aV_nu
L_projected_one_phi = I * G_V / (2 * sp.sqrt(2)) * T_A * (
    u_mu_one_phi * u_nu_one_phi - u_nu_one_phi * u_mu_one_phi
)
L_order_phi_a = bookkeeping_coefficient(L_projected_one_phi, t, a_order)
show("Antisymmetric flavor word inside the trace at O(phi_i a)", L_order_phi_a)


# ---------------------------------------------------------------------------
# 2. V phi_i phi_j flavor structure
# ---------------------------------------------------------------------------
u_phi_pair_mu = -sp.sqrt(2) / F * (
    x_i * dphi_i_mu * lambda_i + x_j * dphi_j_mu * lambda_j
)
u_phi_pair_nu = -sp.sqrt(2) / F * (
    x_i * dphi_i_nu * lambda_i + x_j * dphi_j_nu * lambda_j
)
L_phi_pair = I * G_V / (2 * sp.sqrt(2)) * T_A * (
    u_phi_pair_mu * u_phi_pair_nu - u_phi_pair_nu * u_phi_pair_mu
)
L_phi_i_phi_j_raw = bookkeeping_coefficient(L_phi_pair, x_i, x_j)
L_V_phi_i_phi_j = oriented_antisymmetric_coefficient(
    L_phi_i_phi_j_raw,
    dphi_i_mu,
    dphi_j_nu,
    dphi_i_nu,
    dphi_j_mu,
) * dphi_i_mu * dphi_j_nu
expected_V_phi_i_phi_j = (
    sp.sqrt(2)
    * I
    * G_V
    / F**2
    * T_A
    * (lambda_i * lambda_j - lambda_j * lambda_i)
    * dphi_i_mu
    * dphi_j_nu
)

show("Flavor word extracted at O(phi_i phi_j) before Lorentz ordering", L_phi_i_phi_j_raw)
show("Original-tensor V_A phi_i phi_j coefficient", L_V_phi_i_phi_j)


# ---------------------------------------------------------------------------
# 3. V phi_i a flavor structure
# ---------------------------------------------------------------------------
L_V_phi_i_a = oriented_antisymmetric_coefficient(
    L_order_phi_a,
    dphi_i_mu,
    da_nu,
    dphi_i_nu,
    da_mu,
) * dphi_i_mu * da_nu
expected_V_phi_i_a = (
    -2
    * I
    * G_V
    / (F * f_a)
    * T_A
    * (lambda_i * cA - cA * lambda_i)
    * dphi_i_mu
    * da_nu
)
show("Original-tensor V_A phi_i a coefficient", L_V_phi_i_a)


# ---------------------------------------------------------------------------
# 4. Full photon-vector diagonalization of a linear tensor current
# ---------------------------------------------------------------------------
# First derive the tensor-to-Proca current map before canonical normalization.
M, g = sp.symbols("M g", positive=True, real=True)
v, h, F_gamma_munu, J = sp.symbols("v h F_gamma_munu J", real=True)
L_first_order_with_current = M**2 * v**2 / 4 - (h + g * F_gamma_munu) * v / 2 + J * v
v_solution = sp.solve(sp.diff(L_first_order_with_current, v), v)[0]
L_after_v = sp.expand(L_first_order_with_current.subs(v, v_solution))
linear_current = sp.diff(L_after_v, J).subs(J, 0)
photon_component = sp.diff(linear_current, F_gamma_munu)

show("First-order tensor Lagrangian including a generic interaction current J", L_first_order_with_current)
show("Tensor equation-of-motion solution", v_solution)
show("Term linear in the interaction current after tensor elimination", linear_current * J)
show("Photon component g/M^2 of the original tensor", photon_component)

M1, M2, M3 = sp.symbols("M_1 M_2 M_3", positive=True, real=True)
c1, c2, c3 = sp.symbols("c_1 c_2 c_3", real=True)
T1 = sp.Symbol("T_1", commutative=False)
T2 = sp.Symbol("T_2", commutative=False)
T3 = sp.Symbol("T_3", commutative=False)
g1, g2, g3 = [sp.sqrt(2) * e * F_V * c for c in (c1, c2, c3)]
eps1, eps2, eps3 = sp.symbols("epsilon_1 epsilon_2 epsilon_3", real=True)
eps = sp.Matrix([eps1, eps2, eps3])
eps_physical = sp.Matrix([g1 / M1, g2 / M2, g3 / M3])
Z = 1 + (eps.T * eps)[0]
Z_physical = sp.simplify(Z.subs(dict(zip(eps, eps_physical))))
K_V = sp.eye(3) - eps * eps.T / Z
K_half = sp.eye(3) - eps * eps.T / (sp.sqrt(Z) * (sp.sqrt(Z) + 1))
K_inv_half = sp.eye(3) + eps * eps.T / (sp.sqrt(Z) + 1)

# O is the real orthogonal matrix that diagonalizes
# K_V^{-1/2} diag(M_i^2) K_V^{-1/2}.  Its entries remain generic.
O = sp.MatrixSymbol("O", 3, 3)
W = sp.MatrixSymbol("W_{mu nu}", 3, 1)
mass_inverse = sp.diag(1 / M1, 1 / M2, 1 / M3)

# V=(R+epsilon F)/M, F=F_gamma/sqrt(Z)-epsilon^T R/Z, and
# R=K_V^{-1/2} O W imply the map requested in the review.
proca_numerator_map = K_V * K_inv_half * O * W + eps * F_gamma_munu / sp.sqrt(Z)
tensor_map = mass_inverse * proca_numerator_map
expected_tensor_map = mass_inverse * (K_half * O * W + eps * F_gamma_munu / sp.sqrt(Z))

show("Abstract dimensionless mixing vector epsilon", eps)
show("Physical relation epsilon_i=g_i/M_i", eps_physical)
show("Exact photon normalization Z", Z)
show("Residual vector kinetic matrix K_V", K_V)
show("Positive square root K_V^(1/2)", K_half)
show("Positive inverse square root K_V^(-1/2)", K_inv_half)
show("Derived fully canonical tensor map in matrix form", expected_tensor_map)

a = sp.symbols("a", integer=True)
H = sp.IndexedBase("H")
H_matrix = sp.MatrixSymbol("H", 3, 3)
show("Definition H=K_V^(1/2) O", sp.Eq(H_matrix, K_half * O, evaluate=False))
Tcal_a = sp.Add(
    T1 * H[0, a] / M1,
    T2 * H[1, a] / M2,
    T3 * H[2, a] / M3,
)
Tcal_gamma = sp.Add(
    T1 * eps[0] / (M1 * sp.sqrt(Z)),
    T2 * eps[1] / (M2 * sp.sqrt(Z)),
    T3 * eps[2] / (M3 * sp.sqrt(Z)),
)
Tcal_gamma_physical = Tcal_gamma.subs(dict(zip(eps, eps_physical)))
rho_row = eps.T * K_inv_half * O / Z
rho_row_expected = eps.T * O / sp.sqrt(Z)
rho = sp.IndexedBase("rho")
rho_a = rho[a]
show("Effective canonical-vector flavor matrix mathcal(T)_a", Tcal_a)
show("Effective mixing-induced photon flavor matrix mathcal(T)_gamma", Tcal_gamma)
show("Physical mathcal(T)_gamma displaying F_V", Tcal_gamma_physical)
show("Derived rho row=(epsilon^T K_V^(-1/2) O)/Z", rho_row)
show("Simplified rho_A definition", sp.Eq(rho_a, rho_row_expected[0, a], evaluate=False))

W_a_munu = sp.symbols("W_a_munu", real=True)
final_W_phi_i_phi_j = (
    sp.sqrt(2)
    * I
    * G_V
    / F**2
    * W_a_munu
    * dphi_i_mu
    * dphi_j_nu
    * Tr(Tcal_a * (lambda_i * lambda_j - lambda_j * lambda_i))
)
final_W_phi_i_a = (
    -2
    * I
    * G_V
    / (F * f_a)
    * W_a_munu
    * dphi_i_mu
    * da_nu
    * Tr(Tcal_a * (lambda_i * cA - cA * lambda_i))
)
mixing_gamma_phi_i_a = (
    -2
    * I
    * G_V
    / (F * f_a)
    * F_gamma_munu
    * dphi_i_mu
    * da_nu
    * Tr(Tcal_gamma_physical * (lambda_i * cA - cA * lambda_i))
)
mixing_gamma_phi_i_phi_j = (
    sp.sqrt(2)
    * I
    * G_V
    / F**2
    * F_gamma_munu
    * dphi_i_mu
    * dphi_j_nu
    * Tr(Tcal_gamma_physical * (lambda_i * lambda_j - lambda_j * lambda_i))
)

show("(1) Tensor-origin canonical W^a phi_i phi_j term (pair counted once)", final_W_phi_i_phi_j)
show("(1) Tensor-origin canonical W^a phi_i a term", final_W_phi_i_a)
show("(3) Mixing-induced canonical-photon phi_i phi_j term", mixing_gamma_phi_i_phi_j)
show("(3) Mixing-induced canonical-photon phi_i a term", mixing_gamma_phi_i_a)

# At one electromagnetic insertion, 1/sqrt(Z)=1+O(e^2).
Tcal_gamma_leading_e = sp.diff(Tcal_gamma_physical, e).subs(e, 0) * e
mixing_gamma_leading_e = (
    -2
    * I
    * G_V
    / (F * f_a)
    * F_gamma_munu
    * dphi_i_mu
    * da_nu
    * Tr(Tcal_gamma_leading_e * (lambda_i * cA - cA * lambda_i))
)
show("Mixing-induced photon phi_i a term at leading electromagnetic order", mixing_gamma_leading_e)


# ---------------------------------------------------------------------------
# 5. Goldstone-current terms under the exact field transformation
# ---------------------------------------------------------------------------
Q = sp.Symbol("Q", commutative=False)
A_old_mu, A_gamma_mu, W_a_mu = sp.symbols("A_old_mu A_gamma_mu W_a_mu", real=True)

# At O(Phi), gauging l_mu=r_mu=-e Q A_mu gives
# u_mu=-(sqrt(2)/F)(partial_mu Phi+i e A_mu[Q,Phi])+...
Phi_generic = sp.Symbol("Phi", commutative=False)
dPhi_generic = sp.Symbol("partial^mu_Phi", commutative=False)
u_phi_generic = -sp.sqrt(2) * dPhi_generic / F
u_em_phi_generic = (
    -sp.sqrt(2)
    * I
    * e
    * A_old_mu
    * (Q * Phi_generic - Phi_generic * Q)
    / F
)
goldstone_A_PP_word = sp.expand(
    F**2 / 4 * (u_phi_generic * u_em_phi_generic + u_em_phi_generic * u_phi_generic)
)
goldstone_A_PP = I * e * A_old_mu * Tr(
    Q * (Phi_generic * dPhi_generic - dPhi_generic * Phi_generic)
)

phi_j = sp.Symbol("phi_j", real=True)
pair_current = (phi_i * dphi_j_mu - phi_j * dphi_i_mu) * Tr(
    Q * (lambda_i * lambda_j - lambda_j * lambda_i)
)
goldstone_A_phi_i_phi_j = I * e * A_old_mu * pair_current

# Exact A_old=A_gamma/sqrt(Z)-rho_a W^a gives both direct photon and massive
# vector couplings from the same Goldstone current.
direct_gamma_phi_i_phi_j = goldstone_A_phi_i_phi_j.subs(A_old_mu, A_gamma_mu / sp.sqrt(Z))
goldstone_W_phi_i_phi_j = goldstone_A_phi_i_phi_j.subs(A_old_mu, -rho_a * W_a_mu)

# The direct A-phi-a contact is derived from the same gauged vielbein.
u_em_phi_mu = (
    -sp.sqrt(2)
    * I
    * e
    * A_old_mu
    * phi_i
    * (Q * lambda_i - lambda_i * Q)
    / F
)
u_axial_a_mu = 2 * da_mu * cA / f_a
goldstone_cross_word = sp.expand(
    F**2 / 4 * (u_em_phi_mu * u_axial_a_mu + u_axial_a_mu * u_em_phi_mu)
)
goldstone_A_phi_i_a = (
    sp.sqrt(2)
    * I
    * e
    * F
    / f_a
    * A_old_mu
    * da_mu
    * phi_i
    * Tr((Q * cA - cA * Q) * lambda_i)
)
direct_gamma_phi_i_a = goldstone_A_phi_i_a.subs(A_old_mu, A_gamma_mu / sp.sqrt(Z))
goldstone_W_phi_i_a = goldstone_A_phi_i_a.subs(A_old_mu, -rho_a * W_a_mu)

show("Generic gauged-Goldstone A_old Phi Phi flavor word before trace cyclicity", goldstone_A_PP_word)
show("Generic one-photon/two-pseudoscalar Goldstone term", goldstone_A_PP)
show("Generic component-pair A_old phi_i phi_j term (pair counted once)", goldstone_A_phi_i_phi_j)
show("Direct A_old phi_i a flavor word before trace cyclicity", goldstone_cross_word)
show("Direct A_old phi_i a contact", goldstone_A_phi_i_a)
show("(2) Goldstone-current induced W^a phi_i phi_j term", goldstone_W_phi_i_phi_j)
show("(2) Goldstone-current induced W^a phi_i a term", goldstone_W_phi_i_a)
show("(3) Direct canonical-photon phi_i phi_j Goldstone term", direct_gamma_phi_i_phi_j)
show("(3) Direct canonical-photon phi_i a Goldstone term", direct_gamma_phi_i_a)
show("Direct Pa terms under the EM-invariant spurion assumption [Q,c_q^A]=0", sp.Integer(0))


# ---------------------------------------------------------------------------
# 6. Independent diagnostics and generic component checks
# ---------------------------------------------------------------------------
assert sp.expand(L_V_phi_i_phi_j - expected_V_phi_i_phi_j) == 0
assert sp.expand(L_V_phi_i_a - expected_V_phi_i_a) == 0
assert not L_order_phi_a.has(cV), "c_q^V cannot contribute to V phi a at this order"
assert sp.simplify(photon_component - g / M**2) == 0
assert sp.simplify(K_half * K_half - K_V) == sp.zeros(3)
assert sp.simplify(K_half * K_inv_half - sp.eye(3)) == sp.zeros(3)
assert sp.simplify(K_V * K_inv_half - K_half) == sp.zeros(3)
assert sp.simplify(mass_inverse * K_V * K_inv_half - mass_inverse * K_half) == sp.zeros(3)
assert sp.simplify(eps.T * K_inv_half / Z - eps.T / sp.sqrt(Z)) == sp.zeros(1, 3)
for epsilon_component in (eps1, eps2, eps3):
    assert sp.simplify(K_half.diff(epsilon_component).subs({eps1: 0, eps2: 0, eps3: 0})) == sp.zeros(3)
assert sp.simplify(Z_physical.subs(e, 0) - 1) == 0
assert sp.simplify(sp.diff(Z_physical, e).subs(e, 0)) == 0
assert sp.simplify(Tcal_gamma_physical.subs(e, 0)) == 0
assert sp.expand(
    goldstone_W_phi_i_phi_j
    + I * e * rho_a * W_a_mu * pair_current
) == 0

# Verify the trace-cyclicity sign in the direct term with generic component
# matrices.  This supplements, rather than replaces, the noncommutative result.
component_symbols = sp.symbols("q11 q12 q21 q22 c11 c12 c21 c22 l11 l12 l21 l22")
q11, q12, q21, q22, ca11, ca12, ca21, ca22, l11, l12, l21, l22 = component_symbols
Q_component = sp.Matrix([[q11, q12], [q21, q22]])
cA_component = sp.Matrix([[ca11, ca12], [ca21, ca22]])
lambda_component = sp.Matrix([[l11, l12], [l21, l22]])
cyclic_left = sp.trace(cA_component * (Q_component * lambda_component - lambda_component * Q_component))
cyclic_right = -sp.trace((Q_component * cA_component - cA_component * Q_component) * lambda_component)
assert sp.expand(cyclic_left - cyclic_right) == 0

# Verify the sign and normalization of the generic one-photon/two-Goldstone
# term directly with generic symbolic matrix components.
p11, p12, p21, p22, dp11, dp12, dp21, dp22 = sp.symbols(
    "p11 p12 p21 p22 dp11 dp12 dp21 dp22"
)
Phi_component = sp.Matrix([[p11, p12], [p21, p22]])
dPhi_component = sp.Matrix([[dp11, dp12], [dp21, dp22]])
u_phi_component = -sp.sqrt(2) * dPhi_component / F
u_em_component = -sp.sqrt(2) * I * e * A_old_mu * (
    Q_component * Phi_component - Phi_component * Q_component
) / F
goldstone_APP_component = F**2 / 4 * sp.trace(
    u_phi_component * u_em_component + u_em_component * u_phi_component
)
expected_APP_component = I * e * A_old_mu * sp.trace(
    Q_component * (Phi_component * dPhi_component - dPhi_component * Phi_component)
)
assert sp.expand(goldstone_APP_component - expected_APP_component) == 0

# Generic component limits of the current map.
assert sp.simplify(Tcal_gamma_physical - sum(
    T_k * g_k / (M_k**2 * sp.sqrt(Z_physical))
    for T_k, g_k, M_k in ((T1, g1, M1), (T2, g2, M2), (T3, g3, M3))
)) == 0

# [Q,c_q^A]=0 does not force a differently weighted diagonal flavor matrix
# to commute through the mixing-induced trace.  A d-s component demonstrates
# the nonzero result whenever the effective d and s weights differ.
w_u, w_d, w_s = sp.symbols("w_u w_d w_s", real=True)
Q_em_component = sp.diag(sp.Rational(2, 3), -sp.Rational(1, 3), -sp.Rational(1, 3))
cA_ds_component = sp.zeros(3)
cA_ds_component[1, 2] = -I
cA_ds_component[2, 1] = I
lambda_ds_component = sp.zeros(3)
lambda_ds_component[1, 2] = 1
lambda_ds_component[2, 1] = 1
weighted_diagonal = sp.diag(w_u, w_d, w_s)
mixing_survival_check = sp.trace(
    weighted_diagonal * (lambda_ds_component * cA_ds_component - cA_ds_component * lambda_ds_component)
)
show("EM-invariant d-s spurion mixing trace for unequal effective weights", mixing_survival_check)
assert Q_em_component * cA_ds_component - cA_ds_component * Q_em_component == sp.zeros(3)
assert cA_ds_component.H == cA_ds_component
assert lambda_ds_component.H == lambda_ds_component
assert sp.simplify(mixing_survival_check - 2 * I * (w_d - w_s)) == 0

print("Assumptions and truncations:")
print("- Natural units; flavor brackets are traces; antisymmetric tensors obey V_{mu nu}=-V_{nu mu}.")
print("- Phi=sum_i phi_i lambda_i with generic noncommuting flavor generators.")
print("- Terms are LO in the RChT G_V operator and linear in the ALP.")
print("- O is real orthogonal and diagonalizes K_V^(-1/2) D K_V^(-1/2); no explicit mixing angles are chosen.")
print("- For nondegenerate masses and a continuously chosen O->I basis, H=K_V^(1/2)O=delta+O(e^2).")
print("- No ideal mixing, common vector mass, flavor commutativity, or short-distance relation is imposed.")
print("- The exact 1/sqrt(Z) photon normalization belongs to the truncated tree-level quadratic model.")
print("- The displayed leading-e result retains one electromagnetic mixing insertion.")
print("- [Q,c_q^A]=0 removes direct photon Pa and Goldstone-current WPa terms, but not generically mixing-induced Pa terms.")
print("Units: [F]=[f_a]=[G_V]=[F_V]=[M_k]=mass; e, c_k, and flavor generators are dimensionless.")
print("Diagnostic assertions: all passed.")
