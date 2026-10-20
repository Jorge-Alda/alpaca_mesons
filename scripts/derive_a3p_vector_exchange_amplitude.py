#!/usr/bin/env python3
"""Derive the canonical-vector pole contribution to a -> P_i P_j P_k.

Invocation:
    .venv/bin/python scripts/derive_a3p_vector_exchange_amplitude.py
"""

import sympy as sp
try:
    from scripts.symbolic_utils import show
except ModuleNotFoundError:  # direct ``python scripts/...py`` invocation
    from symbolic_utils import show



I = sp.I
sqrt2 = sp.sqrt(2)
G_V, f_pi, f_a, M_V, Gamma_V = sp.symbols(
    "G_V f_pi f_a M_V Gamma_V", positive=True, real=True
)
X, Y = sp.symbols("X Y")

print("Assumptions:")
print("- Natural units and metric dot products; all external pseudoscalars are on shell.")
print("- Alda:2026nrd convention: u=exp(i Phi_A/(2 f_pi)), hence U=u^2=exp(i Phi_A/f_pi).")
print("- RChT convention: U=exp(i sqrt(2) Phi_R/f_pi), with Phi_A=sqrt(2) Phi_R.")
print("- One canonical Proca pole with constant width; the induced representation-matching local term is derived separately below.")
print("- Pair labels are fixed as i<j and all momenta in the displayed decay kinematics are outgoing.\n")

# ---------------------------------------------------------------------------
# 1. Generator-normalization conversion
# ---------------------------------------------------------------------------
Lambda_i, Lambda_j, Lambda_k, T_V, c_q_A = sp.symbols(
    "Lambda_i Lambda_j Lambda_k T_V c_q^A", commutative=False
)
Tr = sp.Function("Tr")
lambda_i = Lambda_i / sqrt2
lambda_j = Lambda_j / sqrt2
lambda_k = Lambda_k / sqrt2

rcht_wpp_coefficient = sqrt2 * I * G_V / (f_pi**2 * M_V)
rcht_wpa_coefficient = -2 * I * G_V / (f_pi * f_a * M_V)
alda_2026nrd_wpp_coefficient = sp.simplify(rcht_wpp_coefficient / 2)
alda_2026nrd_wpa_coefficient = sp.simplify(rcht_wpa_coefficient / sqrt2)
normalization_product = sp.simplify(alda_2026nrd_wpp_coefficient * alda_2026nrd_wpa_coefficient)

show("Alda:2026nrd/RChT field relation Phi_A=sqrt(2) Phi_R", sp.Eq(sp.Symbol("Phi_A"), sqrt2 * sp.Symbol("Phi_R")))
show("RChT pseudoscalar generator lambda_i=Lambda_i/sqrt(2)", sp.Eq(sp.Symbol("lambda_i"), lambda_i))
show("Converted L_WPP coefficient", alda_2026nrd_wpp_coefficient)
show("Converted L_WPa coefficient", alda_2026nrd_wpa_coefficient)
show("Product of the two Lagrangian coefficients", normalization_product)
show("Positive normalization after i*(-i)=+1", G_V**2 / (f_pi**3 * f_a * M_V**2))

expected_wpp = I * G_V / (sqrt2 * f_pi**2 * M_V)
expected_wpa = -sqrt2 * I * G_V / (f_pi * f_a * M_V)
assert sp.simplify(alda_2026nrd_wpp_coefficient - expected_wpp) == 0
assert sp.simplify(alda_2026nrd_wpa_coefficient - expected_wpa) == 0
assert sp.simplify(normalization_product - G_V**2 / (f_pi**3 * f_a * M_V**2)) == 0

# Scalar i-factor audit.  With exp(-ip.x), each three-derivative vertex has
# (-i)^3=i.  Including the Feynman-rule i gives -i*a at WPP and +i*b at WPa.
a, b = sp.symbols("a b", positive=True, real=True)
derivative_factor = (-I) ** 3
vertex_wpp = sp.simplify(I * (I * a) * derivative_factor)
vertex_wpa = sp.simplify(I * (-I * b) * derivative_factor)
vertex_product = sp.simplify(vertex_wpp * vertex_wpa)
show("WPP scalar vertex i-factor", vertex_wpp)
show("WPa scalar vertex i-factor", vertex_wpa)
show("Product of scalar vertex i-factors", vertex_product)
assert vertex_wpp == -I * a
assert vertex_wpa == I * b
assert vertex_product == a * b

# ---------------------------------------------------------------------------
# 2. Derive currents and their Mandelstam contraction with a Gram matrix
# ---------------------------------------------------------------------------
m_a2, m_i2, m_j2, m_k2 = sp.symbols("m_a^2 m_i^2 m_j^2 m_k^2", real=True)
s, s_ik, s_jk = sp.symbols("s s_ik s_jk", real=True)

# Basis is (p_i,p_j,p_k).  The matrix entries are generated from on-shell
# conditions and s_ab=(p_a+p_b)^2.
gram = sp.Matrix(
    [
        [m_i2, (s - m_i2 - m_j2) / 2, (s_ik - m_i2 - m_k2) / 2],
        [(s - m_i2 - m_j2) / 2, m_j2, (s_jk - m_j2 - m_k2) / 2],
        [(s_ik - m_i2 - m_k2) / 2, (s_jk - m_j2 - m_k2) / 2, m_k2],
    ]
)

def dot(left: sp.Matrix, right: sp.Matrix) -> sp.Expr:
    """Contract coefficient vectors with the symbolic momentum Gram matrix."""
    return sp.expand((left.T * gram * right)[0])


p_i = sp.Matrix([1, 0, 0])
p_j = sp.Matrix([0, 1, 0])
p_k = sp.Matrix([0, 0, 1])
q = p_i + p_j
p_a = q + p_k
J_a = sp.simplify(dot(q, p_k)) * p_a - sp.simplify(dot(q, p_a)) * p_k
J_PP = sp.simplify(dot(q, p_i)) * p_j - sp.simplify(dot(q, p_j)) * p_i

q_dot_Ja = sp.factor(dot(q, J_a))
q_dot_JPP = sp.factor(dot(q, J_PP))
numerator_raw = sp.factor(dot(J_a, J_PP))
numerator_expected = sp.factor(
    s
    / 4
    * (
        (s + m_j2 - m_i2) * (s_ik - m_i2 - m_k2)
        - (s + m_i2 - m_j2) * (s_jk - m_j2 - m_k2)
    )
)

show("Momentum Gram matrix in basis (p_i,p_j,p_k)", gram)
show("q coefficient vector", q)
show("J_a coefficient vector", J_a)
show("J_PP coefficient vector", J_PP)
show("Transversality q.J_a", q_dot_Ja)
show("Transversality q.J_PP", q_dot_JPP)
show("Derived J_a.J_PP", numerator_raw)
show("Mandelstam numerator in implementation form", numerator_expected)

assert q_dot_Ja == 0
assert q_dot_JPP == 0
assert sp.simplify(numerator_raw - numerator_expected) == 0

# ---------------------------------------------------------------------------
# 3. Pole amplitude, pair orientation, and scaling checks
# ---------------------------------------------------------------------------
denominator = M_V**2 - s - I * M_V * Gamma_V
amplitude = sp.factor(
    G_V**2 / (f_pi**3 * f_a * M_V**2) * X * Y * numerator_expected / denominator
)
show("Final labeled canonical-vector pole amplitude", amplitude)

pair_reversed_numerator = numerator_expected.xreplace(
    {m_i2: m_j2, m_j2: m_i2, s_ik: s_jk, s_jk: s_ik}
)
show("Numerator after reversing pair i<->j", sp.factor(pair_reversed_numerator))
assert sp.simplify(pair_reversed_numerator + numerator_expected) == 0
# X=Tr(T_Vbar[Lambda_i,Lambda_j]) also changes sign; Y does not.
assert sp.simplify((-X) * pair_reversed_numerator - X * numerator_expected) == 0

z = sp.symbols("z", positive=True, real=True)
dimensionful_scaling = {
    G_V: z * G_V,
    f_pi: z * f_pi,
    f_a: z * f_a,
    M_V: z * M_V,
    Gamma_V: z * Gamma_V,
    s: z**2 * s,
    s_ik: z**2 * s_ik,
    s_jk: z**2 * s_jk,
    m_i2: z**2 * m_i2,
    m_j2: z**2 * m_j2,
    m_k2: z**2 * m_k2,
}
scaled_amplitude = sp.simplify(amplitude.xreplace(dimensionful_scaling))
show("Amplitude under a common mass-unit rescaling", scaled_amplitude)
assert sp.simplify(scaled_amplitude - amplitude) == 0
assert sp.simplify(sp.diff(amplitude, X) * X - amplitude) == 0
assert sp.simplify(sp.diff(amplitude, f_a) * f_a + amplitude) == 0

print("Sign audit: the executable scalar i-factor bookkeeping gives the requested positive pole contribution.")
print("Units: G_V, f_pi, f_a, M_V, Gamma_V have mass dimension 1; Mandelstam variables have dimension 2.")
print("Diagnostic assertions: normalization, sign bookkeeping, transversality, pair reversal, and scaling passed.")


# ---------------------------------------------------------------------------
# 4. Identical-particle Mandelstam dependence
# ---------------------------------------------------------------------------
Lambda_P = sp.MatrixSymbol("Lambda_P", 3, 3)
Lambda_Q = sp.MatrixSymbol("Lambda_Q", 3, 3)
T_Vbar = sp.MatrixSymbol("T_Vbar", 3, 3)

# If all final particles are P, every possible VPP pair is [Lambda_P,Lambda_P].
identical_vpp_commutator = sp.simplify(Lambda_P * Lambda_P - Lambda_P * Lambda_P)
identical_vpp_trace = sp.Trace(T_Vbar * identical_vpp_commutator).doit()
show("Three-identical VPP commutator [Lambda_P,Lambda_P]", identical_vpp_commutator)
show("Three-identical VPP trace for every spectator assignment", identical_vpp_trace)
assert identical_vpp_commutator == sp.ZeroMatrix(3, 3)
assert identical_vpp_trace == 0

# For P(p1) P(p2) Q(p3), use the implementation's fixed increasing pair
# order.  The Q-spectator graph has the identical pair (p1,p2), while the two
# P-spectator graphs have pair invariants s23 and s13 respectively.
m_P2, m_Q2 = sp.symbols("m_P^2 m_Q^2", nonnegative=True, real=True)
s12, s13, s23 = sp.symbols("s_12 s_13 s_23", real=True)
Sigma = sp.factor(m_a2 + 2 * m_P2 + m_Q2)

q_spectator_vpp = sp.simplify(Lambda_P * Lambda_P - Lambda_P * Lambda_P)
N_23 = sp.factor(
    numerator_expected.subs(
        {
            s: s23,
            m_i2: m_P2,
            m_j2: m_Q2,
            m_k2: m_P2,
            s_ik: s12,
            s_jk: s13,
        },
        simultaneous=True,
    )
)
N_13 = sp.factor(
    numerator_expected.subs(
        {
            s: s13,
            m_i2: m_P2,
            m_j2: m_Q2,
            m_k2: m_P2,
            s_ik: s12,
            s_jk: s23,
        },
        simultaneous=True,
    )
)

show("P(p1) P(p2) Q(p3): Q-spectator VPP commutator", q_spectator_vpp)
show("P(p1)-spectator pole numerator N_23", N_23)
show("P(p2)-spectator pole numerator N_13", N_13)
assert q_spectator_vpp == sp.ZeroMatrix(3, 3)
assert sp.simplify(N_23.xreplace({s13: s23, s23: s13}) - N_13) == 0

# Express the Bose-symmetric numerator sum in the elementary symmetric
# variables S=s13+s23 and P=s13*s23.  SymPy's standard symmetrize routine does
# this directly and also proves that no nonsymmetric remainder is present.
S_pair, P_pair = sp.symbols("S_pair P_pair", real=True)
numerator_symmetric, numerator_remainder, elementary = sp.symmetrize(
    sp.expand(N_13 + N_23), [s13, s23], formal=True
)
elementary_substitution = {elementary[0][0]: S_pair, elementary[1][0]: P_pair}
numerator_SP = sp.factor(numerator_symmetric.subs(elementary_substitution))
numerator_on_shell_SP = sp.factor(numerator_SP.subs(S_pair, Sigma - s12))
numerator_product_coefficient = sp.factor(sp.diff(numerator_on_shell_SP, P_pair))

show("Symmetric numerator remainder", numerator_remainder)
show("N_13+N_23 in S_pair=s13+s23 and P_pair=s13*s23", numerator_SP)
show("On-shell numerator sum after S_pair=Sigma-s12", numerator_on_shell_SP)
show("Coefficient of s13*s23 in the on-shell numerator sum", numerator_product_coefficient)
assert numerator_remainder == 0
assert sp.simplify(numerator_product_coefficient + (Sigma + s12) / 4) == 0

# Add the two finite pole denominators.  At fixed s12 their common denominator
# is itself affine in P_pair, so finite pole propagation adds further product
# dependence rather than restoring a one-variable contact form.
Delta_V = sp.Symbol("Delta_V")  # Delta_V=M_V^2-i M_V Gamma_V
pole_sum = sp.factor(N_13 / (Delta_V - s13) + N_23 / (Delta_V - s23))
pole_numerator, pole_denominator = sp.fraction(sp.together(pole_sum))
pole_num_symmetric, pole_num_remainder, pole_num_elementary = sp.symmetrize(
    sp.expand(pole_numerator), [s13, s23], formal=True
)
pole_den_symmetric, pole_den_remainder, pole_den_elementary = sp.symmetrize(
    sp.expand(pole_denominator), [s13, s23], formal=True
)
pole_num_SP = sp.factor(
    pole_num_symmetric.subs(
        {
            pole_num_elementary[0][0]: S_pair,
            pole_num_elementary[1][0]: P_pair,
        }
    )
)
pole_den_SP = sp.factor(
    pole_den_symmetric.subs(
        {
            pole_den_elementary[0][0]: S_pair,
            pole_den_elementary[1][0]: P_pair,
        }
    )
)
pole_on_shell_SP = sp.factor(
    (pole_num_SP / pole_den_SP).subs(S_pair, Sigma - s12)
)
pole_product_derivative = sp.factor(sp.diff(pole_on_shell_SP, P_pair))
pole_derivative_numerator = sp.factor(
    sp.fraction(sp.together(pole_product_derivative))[0]
)
pole_derivative_quadratic_coefficient = sp.factor(
    sp.Poly(sp.expand(pole_derivative_numerator), P_pair).coeff_monomial(P_pair**2)
)

show("Bose-symmetric finite-pole sum", pole_sum)
show("Finite-pole common numerator in (S_pair,P_pair)", pole_num_SP)
show("Finite-pole common denominator in (S_pair,P_pair)", pole_den_SP)
show("On-shell derivative of pole sum with respect to P_pair=s13*s23", pole_product_derivative)
show("Quadratic P_pair coefficient in the derivative numerator", pole_derivative_quadratic_coefficient)
assert pole_num_remainder == 0
assert pole_den_remainder == 0
assert pole_product_derivative != 0
assert pole_derivative_quadratic_coefficient != 0

# Exceptional reductions are exact flavor/coupling zeros, restriction to a
# one-dimensional kinematic locus, or the algebraic but nonphysical condition
# Sigma+s12=0 for the numerator sum.  Equal masses, zero width, or replacing
# both pole denominators by one constant do not remove the P_pair coefficient.
constant_denominator = sp.Symbol("D_0", nonzero=True)
constant_denominator_sum = sp.factor(numerator_on_shell_SP / constant_denominator)
show("Sum after replacing both pole denominators by one constant", constant_denominator_sum)
assert sp.simplify(
    sp.diff(constant_denominator_sum, P_pair)
    + (Sigma + s12) / (4 * constant_denominator)
) == 0
assert sp.factor(numerator_product_coefficient.subs(m_Q2, m_P2)) != 0

print("Identical-particle conclusion:")
print("- Three identical pseudoscalars give zero canonical vector-pole amplitude because every VPP commutator vanishes.")
print("- For P(p1)P(p2)Q(p3), the Q-spectator graph vanishes, but the P-spectator graphs probe s23 and s13.")
print("- Their Bose-symmetric sum remains a two-variable Dalitz amplitude; even its numerator sum depends on s13*s23.")
print("- Exact flavor/coupling cancellations can remove the amplitude; restricting to s13=s23 gives only a one-dimensional locus.")
print("- The nonphysical Sigma+s12=0 condition removes the numerator-product coefficient, but not finite-pole dependence.")
print("Identical-particle assertions: all passed.")


# ---------------------------------------------------------------------------
# 5. Representation-matching local current-current term
# ---------------------------------------------------------------------------
# First recover the original antisymmetric-tensor coefficients from the
# already-derived canonical Proca coefficients.  V_munu=W_munu/M_V+... at
# first order in the current, so multiplying the Proca coefficients by M_V
# must remove exactly one inverse mass.
tensor_wpp_coefficient = sp.simplify(M_V * alda_2026nrd_wpp_coefficient)
tensor_wpa_coefficient = sp.simplify(M_V * alda_2026nrd_wpa_coefficient)
expected_tensor_wpp = I * G_V / (sqrt2 * f_pi**2)
expected_tensor_wpa = -sqrt2 * I * G_V / (f_pi * f_a)
show("Recovered Alda:2026nrd-normalized original-tensor VPP coefficient", tensor_wpp_coefficient)
show("Recovered Alda:2026nrd-normalized original-tensor VPa coefficient", tensor_wpa_coefficient)
assert sp.simplify(tensor_wpp_coefficient - expected_tensor_wpp) == 0
assert sp.simplify(tensor_wpa_coefficient - expected_tensor_wpa) == 0

# Scalar placeholders encode the coefficients of the fully contracted tensor
# expression.  K_PP and K_Pa therefore mean the projected currents
# K_[mu nu]=(K_mu nu-K_nu mu)/2, not unprojected outer products.
V_component, W_component, K_PP, K_Pa = sp.symbols(
    "V_component W_component K_PP K_Pa"
)
K_total = K_PP + K_Pa
first_order_lagrangian = (
    M_V**2 * V_component**2 / 4
    - M_V * W_component * V_component / 2
    + V_component * K_total
)
V_solution = sp.solve(sp.diff(first_order_lagrangian, V_component), V_component)[0]
eliminated_lagrangian = sp.factor(first_order_lagrangian.subs(V_component, V_solution))
expected_eliminated = sp.factor(
    -W_component**2 / 4
    + W_component * K_total / M_V
    - K_total**2 / M_V**2
)
local_cross_term = sp.factor(
    sp.diff(eliminated_lagrangian, K_PP, K_Pa)
    .subs({K_PP: 0, K_Pa: 0})
    * K_PP
    * K_Pa
)

show("First-order tensor Lagrangian for K=K_PP+K_Pa", first_order_lagrangian)
show("Antisymmetric-tensor equation-of-motion solution", V_solution)
show("Lagrangian after eliminating V_munu", eliminated_lagrangian)
show("Extracted local K_PP-K_Pa cross term", local_cross_term)
assert sp.simplify(eliminated_lagrangian - expected_eliminated) == 0
assert sp.simplify(local_cross_term + 2 * K_PP * K_Pa / M_V**2) == 0

# Derive the factor from antisymmetrizing the raw derivative outer products.
# The coefficient-vector basis is (p_i,p_j,p_k,p_a), with an otherwise generic
# symmetric Gram matrix.  SymPy performs the full double Lorentz contraction.
d_ii, d_jj, d_kk, d_aa = sp.symbols("d_ii d_jj d_kk d_aa")
d_ij, d_ik, d_ia, d_jk, d_ja, d_ka = sp.symbols(
    "d_ij d_ik d_ia d_jk d_ja d_ka"
)
gram4 = sp.Matrix(
    [
        [d_ii, d_ij, d_ik, d_ia],
        [d_ij, d_jj, d_jk, d_ja],
        [d_ik, d_jk, d_kk, d_ka],
        [d_ia, d_ja, d_ka, d_aa],
    ]
)
e_i = sp.Matrix([1, 0, 0, 0])
e_j = sp.Matrix([0, 1, 0, 0])
e_k = sp.Matrix([0, 0, 1, 0])
e_a = sp.Matrix([0, 0, 0, 1])
K_PP_raw = e_i * e_j.T
K_Pa_raw = e_k * e_a.T
K_PP_antisym = (K_PP_raw - K_PP_raw.T) / 2
K_Pa_antisym = (K_Pa_raw - K_Pa_raw.T) / 2
antisymmetric_current_contraction = sp.factor(
    sum(
        K_PP_antisym[mu, nu]
        * K_Pa_antisym[rho, sigma]
        * gram4[mu, rho]
        * gram4[nu, sigma]
        for mu in range(4)
        for nu in range(4)
        for rho in range(4)
        for sigma in range(4)
    )
)
L_dot = sp.factor(d_ik * d_ja - d_jk * d_ia)
show("Contraction K_PP[mu nu] K_Pa^[mu nu] after explicit antisymmetrization", antisymmetric_current_contraction)
show("Definition L_(k;ij) in scalar products", L_dot)
assert sp.simplify(antisymmetric_current_contraction - L_dot / 2) == 0

# The two original-tensor coefficients multiply to a positive real coefficient.
# Combining -2 K_PP.K_Pa/M_V^2 with K_PP.K_Pa=L/2 gives a coordinate-space
# local coefficient -C*L.  For one incoming ALP and three outgoing mesons,
# the four derivative Fourier factors multiply to -1, hence i*M_match=+i*C*L.
tensor_coefficient_product = sp.simplify(
    tensor_wpp_coefficient * tensor_wpa_coefficient
)
matching_prefactor = G_V**2 * X * Y / (f_pi**3 * f_a * M_V**2)
local_lagrangian_coefficient = sp.factor(
    -tensor_coefficient_product * X * Y / M_V**2
)
four_derivative_fourier_factor = sp.simplify((I**3) * (-I))
contact_vertex = sp.simplify(
    I * local_lagrangian_coefficient * four_derivative_fourier_factor * L_dot
)
matching_amplitude = sp.simplify(contact_vertex / I)
show("Product of original-tensor VPP and VPa coefficients", tensor_coefficient_product)
show("Local matching Lagrangian coefficient multiplying L_(k;ij)", local_lagrangian_coefficient)
show("Fourier factor for three outgoing mesons and one incoming ALP", four_derivative_fourier_factor)
show("Derived S=1+iT matching amplitude", matching_amplitude)
assert sp.simplify(tensor_coefficient_product - G_V**2 / (f_pi**3 * f_a)) == 0
assert four_derivative_fourier_factor == -1
assert sp.simplify(matching_amplitude - matching_prefactor * L_dot) == 0

# This same exact elimination is the propagator decomposition.  Separate the
# source-dependent effective action from the source-free Proca term, then set
# W=0 to isolate the momentum-independent metric part of the full
# antisymmetric-tensor propagator.  Its PP-Pa cross coefficient must reproduce
# the matching amplitude; the remaining W.K/M terms generate the Proca pole.
source_dependent_action = sp.expand(
    eliminated_lagrangian
    - eliminated_lagrangian.subs({K_PP: 0, K_Pa: 0})
)
nonpole_source_action = sp.factor(source_dependent_action.subs(W_component, 0))
nonpole_cross_coefficient = sp.factor(
    sp.diff(nonpole_source_action, K_PP, K_Pa)
    .subs({K_PP: 0, K_Pa: 0})
)
nonpole_from_tensor_propagator = sp.simplify(
    nonpole_cross_coefficient
    * tensor_coefficient_product
    * X
    * Y
    * antisymmetric_current_contraction
    * four_derivative_fourier_factor
)
show("Source-dependent action after exact tensor elimination", source_dependent_action)
show("Momentum-independent non-pole source action", nonpole_source_action)
show("PP-Pa coefficient in the non-pole propagator piece", nonpole_cross_coefficient)
show("Non-pole part of the full antisymmetric-tensor exchange", nonpole_from_tensor_propagator)
assert sp.simplify(nonpole_cross_coefficient + 2 / M_V**2) == 0
assert sp.simplify(nonpole_from_tensor_propagator - matching_amplitude) == 0
assert not matching_amplitude.has(Gamma_V)
matching_dynamic_denominator = sp.factor(sp.denom(sp.together(matching_amplitude)))
show("Matching-amplitude denominator (no kinematic pole)", matching_dynamic_denominator)
assert not any(
    matching_dynamic_denominator.has(variable)
    for variable in (s, s_ik, s_jk, m_a2, m_i2, m_j2, m_k2, Gamma_V)
)

contact_unit_scaling = {
    G_V: z * G_V,
    f_pi: z * f_pi,
    f_a: z * f_a,
    M_V: z * M_V,
    d_ik: z**2 * d_ik,
    d_ja: z**2 * d_ja,
    d_jk: z**2 * d_jk,
    d_ia: z**2 * d_ia,
}
scaled_matching_amplitude = sp.simplify(
    matching_amplitude.xreplace(contact_unit_scaling)
)
show("Matching amplitude under a common mass-unit rescaling", scaled_matching_amplitude)
assert sp.simplify(scaled_matching_amplitude - matching_amplitude) == 0

# On-shell Mandelstam representation of L_(k;ij).
L_mandelstam = sp.factor(
    (
        (s_ik - m_i2 - m_k2) * (m_a2 + m_j2 - s_ik)
        - (s_jk - m_j2 - m_k2) * (m_a2 + m_i2 - s_jk)
    )
    / 4
)
show("L_(k;ij) in Mandelstam variables", L_mandelstam)
assert sp.simplify(
    L_dot.subs(
        {
            d_ik: (s_ik - m_i2 - m_k2) / 2,
            d_ja: (m_a2 + m_j2 - s_ik) / 2,
            d_jk: (s_jk - m_j2 - m_k2) / 2,
            d_ia: (m_a2 + m_i2 - s_jk) / 2,
        }
    )
    - L_mandelstam
) == 0
contact_pair_reversed = L_mandelstam.xreplace(
    {m_i2: m_j2, m_j2: m_i2, s_ik: s_jk, s_jk: s_ik}
)
assert sp.simplify(contact_pair_reversed + L_mandelstam) == 0

# Identical fields for the matching contact.  The Q-spectator flavor factor is
# zero.  The two P-spectator diagrams have the same X=Tr(Tbar[Lambda_P,Lambda_Q])
# and Y=Tr(T[Lambda_P,c_q^A]), so their kinematic factors add with a plus sign.
T_V_matrix = sp.MatrixSymbol("T_V", 3, 3)
c_q_A_matrix = sp.MatrixSymbol("c_q^A", 3, 3)
X_PQ = sp.Trace(
    T_Vbar * (Lambda_P * Lambda_Q - Lambda_Q * Lambda_P)
)
Y_P = sp.Trace(
    T_V_matrix * (Lambda_P * c_q_A_matrix - c_q_A_matrix * Lambda_P)
)
X_PP = sp.Trace(
    T_Vbar * (Lambda_P * Lambda_P - Lambda_P * Lambda_P)
).doit()
flavor_p1_spectator = X_PQ * Y_P
flavor_p2_spectator = X_PQ * Y_P
show("PPQ Q-spectator VPP flavor trace", X_PP)
show("PPQ p1-spectator flavor factor X_PQ Y_P", flavor_p1_spectator)
show("PPQ p2-spectator flavor factor X_PQ Y_P", flavor_p2_spectator)
assert X_PP == 0
assert sp.simplify(flavor_p1_spectator - flavor_p2_spectator) == 0

L_23 = sp.factor(
    L_mandelstam.subs(
        {
            m_i2: m_P2,
            m_j2: m_Q2,
            m_k2: m_P2,
            s_ik: s12,
            s_jk: s13,
        },
        simultaneous=True,
    )
)
L_13 = sp.factor(
    L_mandelstam.subs(
        {
            m_i2: m_P2,
            m_j2: m_Q2,
            m_k2: m_P2,
            s_ik: s12,
            s_jk: s23,
        },
        simultaneous=True,
    )
)
contact_sum = sp.factor(L_13 + L_23)
contact_symmetric, contact_remainder, contact_elementary = sp.symmetrize(
    sp.expand(contact_sum), [s13, s23], formal=True
)
contact_SP = sp.factor(
    contact_symmetric.subs(
        {
            contact_elementary[0][0]: S_pair,
            contact_elementary[1][0]: P_pair,
        }
    )
)
contact_on_shell_SP = sp.factor(contact_SP.subs(S_pair, Sigma - s12))
contact_product_coefficient = sp.factor(sp.diff(contact_on_shell_SP, P_pair))

show("Matching-contact factor for the p1 spectator (pair p2,Q)", L_23)
show("Matching-contact factor for the p2 spectator (pair p1,Q)", L_13)
show("Bose-symmetric matching-contact sum in (S_pair,P_pair)", contact_SP)
show("On-shell matching-contact sum", contact_on_shell_SP)
show("Coefficient of s13*s23 in the matching-contact sum", contact_product_coefficient)
assert contact_remainder == 0
assert sp.simplify(L_23.xreplace({s13: s23, s23: s13}) - L_13) == 0
assert contact_product_coefficient == -sp.Rational(1, 2)

print("Representation-matching contact conclusion:")
print("- K denotes the explicitly antisymmetrized current K_[mu nu]; tensor elimination gives W.K/M-K_[mu nu]K^[mu nu]/M^2.")
print("- The PP-Pa cross term yields +G_V^2 X Y L_(k;ij)/(f_pi^3 f_a M_V^2) in the S=1+iT amplitude.")
print("- L_(k;ij) has mass dimension four, is O(p^4), and contains no pole denominator or width.")
print("- Three identical fields vanish; for PPQ the equal flavor factors multiply L_13+L_23, which still depends on s13*s23.")
print("Representation-matching contact assertions: all passed.")
