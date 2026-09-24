#!/usr/bin/env python3
"""Derive the canonical-vector pole contribution to a -> P_i P_j P_k.

Invocation:
    .venv/bin/python scripts/derive_a3p_vector_exchange_amplitude.py
"""

import sympy as sp


def show(label: str, expression: sp.Basic) -> None:
    """Print one deterministic derivation step in LaTeX."""
    print(f"{label}:\n{sp.latex(expression)}\n")


I = sp.I
sqrt2 = sp.sqrt(2)
G_V, f_pi, f_a, M_V, Gamma_V = sp.symbols(
    "G_V f_pi f_a M_V Gamma_V", positive=True, real=True
)
X, Y = sp.symbols("X Y")

print("Assumptions:")
print("- Natural units and metric dot products; all external pseudoscalars are on shell.")
print("- Alda convention: u=exp(i Phi_A/(2 f_pi)), hence U=u^2=exp(i Phi_A/f_pi).")
print("- RChT convention: U=exp(i sqrt(2) Phi_R/f_pi), with Phi_A=sqrt(2) Phi_R.")
print("- One canonical Proca pole with constant width; no representation-matching local current-current term.")
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
alda_wpp_coefficient = sp.simplify(rcht_wpp_coefficient / 2)
alda_wpa_coefficient = sp.simplify(rcht_wpa_coefficient / sqrt2)
normalization_product = sp.simplify(alda_wpp_coefficient * alda_wpa_coefficient)

show("Alda/RChT field relation Phi_A=sqrt(2) Phi_R", sp.Eq(sp.Symbol("Phi_A"), sqrt2 * sp.Symbol("Phi_R")))
show("RChT pseudoscalar generator lambda_i=Lambda_i/sqrt(2)", sp.Eq(sp.Symbol("lambda_i"), lambda_i))
show("Converted L_WPP coefficient", alda_wpp_coefficient)
show("Converted L_WPa coefficient", alda_wpa_coefficient)
show("Product of the two Lagrangian coefficients", normalization_product)
show("Positive normalization after i*(-i)=+1", G_V**2 / (f_pi**3 * f_a * M_V**2))

expected_wpp = I * G_V / (sqrt2 * f_pi**2 * M_V)
expected_wpa = -sqrt2 * I * G_V / (f_pi * f_a * M_V)
assert sp.simplify(alda_wpp_coefficient - expected_wpp) == 0
assert sp.simplify(alda_wpa_coefficient - expected_wpa) == 0
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
