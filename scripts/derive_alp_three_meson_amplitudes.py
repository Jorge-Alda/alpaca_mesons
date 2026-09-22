"""Derive a -> three-meson amplitudes and Dalitz-plot decay widths.

Invocation:
    .venv/bin/python scripts/derive_alp_three_meson_amplitudes.py

Natural units and the mostly-minus metric are used. Momenta p, p_1, p_2,
p_3 describe a(p) -> phi_1(p_1) phi_2(p_2) phi_3(p_3), and
s_rs = (p_r + p_s)^2. The S-matrix convention is S = 1 + i T.
"""

import sympy as sp


def show(label: str, expression: sp.Expr) -> None:
    print(f"{label}:\n{sp.latex(expression)}\n")


def multilinear_coefficient(expression: sp.Expr, legs: tuple[sp.Symbol, ...]) -> sp.Expr:
    """Functionally differentiate once with respect to every external leg."""
    result = expression
    for leg in legs:
        result = sp.diff(result, leg)
    return sp.simplify(result.subs({leg: 0 for leg in legs}))


K, M = sp.symbols("K M", real=True)
m_a, m_i, m_j, m_k = sp.symbols("m_a m_i m_j m_k", positive=True)
s12, s13, s23 = sp.symbols("s_12 s_13 s_23", real=True)
d1, d2, d3 = sp.symbols("d_1 d_2 d_3", real=True)
x1, x2, x3 = sp.symbols("x_1 x_2 x_3")

T_iii_P, T_iij_A, T_iij_P, T_ijk_P = sp.symbols(
    "T_iii^P T_iij^A T_iij^P T_ijk^P"
)
T_ij_k_A, T_ik_j_A, T_jk_i_A = sp.symbols(
    "T_(ij;k)^A T_(ik;j)^A T_(jk;i)^A"
)

# Replace each field by the sum of its labeled external legs. Contracting a
# differentiated outgoing field with partial_mu a gives d_r = p . p_r.
phi_i_identical = x1 + x2 + x3
vertex_iii = M * T_iii_P * phi_i_identical**3
amplitude_iii = multilinear_coefficient(vertex_iii, (x1, x2, x3))

phi_i = x1 + x2
phi_j = x3
dphi_i = d1 * x1 + d2 * x2
dphi_j = d3 * x3
vertex_iij = (
    K * T_iij_A * (phi_i**2 * dphi_j - phi_i * phi_j * dphi_i)
    + M * T_iij_P * phi_i**2 * phi_j
)
amplitude_iij_dot = multilinear_coefficient(vertex_iij, (x1, x2, x3))

vertex_ijk = (
    K
    * (
        T_ij_k_A * x1 * x2 * d3 * x3
        + T_ik_j_A * x1 * x3 * d2 * x2
        + T_jk_i_A * x2 * x3 * d1 * x1
    )
    + M * T_ijk_P * x1 * x2 * x3
)
amplitude_ijk_dot = multilinear_coefficient(vertex_ijk, (x1, x2, x3))

# On-shell momentum conservation gives
# p.p_1=(m_a^2+m_1^2-s_23)/2 and cyclic permutations.
dot_iij = {
    d1: (m_a**2 + m_i**2 - s23) / 2,
    d2: (m_a**2 + m_i**2 - s13) / 2,
    d3: (m_a**2 + m_j**2 - s12) / 2,
}
dalitz_iij = {s13: m_a**2 + 2 * m_i**2 + m_j**2 - s12 - s23}
amplitude_iij = sp.factor(amplitude_iij_dot.subs(dot_iij).subs(dalitz_iij))

dot_ijk = {
    d1: (m_a**2 + m_i**2 - s23) / 2,
    d2: (m_a**2 + m_j**2 - s13) / 2,
    d3: (m_a**2 + m_k**2 - s12) / 2,
}
amplitude_ijk = sp.collect(sp.expand(amplitude_ijk_dot.subs(dot_ijk)), K)
amplitude_ijk_2d = sp.simplify(
    amplitude_ijk.subs(
        s13, m_a**2 + m_i**2 + m_j**2 + m_k**2 - s12 - s23
    )
)

show("Three-identical vertex polynomial", vertex_iii)
show("Three-identical amplitude", amplitude_iii)
show("Two-identical amplitude in momentum products", amplitude_iij_dot)
show("Two-identical amplitude in s_12", amplitude_iij)
show("Three-different amplitude in momentum products", amplitude_ijk_dot)
show("Three-different amplitude in Mandelstam variables", amplitude_ijk)
show("Three-different amplitude in independent variables s_12 and s_23", amplitude_ijk_2d)

# Derive the Jacobi-identity consequence before taking the flavour trace. It
# removes any common shift from all three p.p_r terms in the distinct amplitude.
lambda_i, lambda_j, lambda_k = sp.symbols(
    "lambda_i lambda_j lambda_k", commutative=False
)


def commutator(left: sp.Expr, right: sp.Expr) -> sp.Expr:
    return left * right - right * left


def kinetic_tensor(first: sp.Expr, second: sp.Expr, differentiated: sp.Expr) -> sp.Expr:
    return commutator(first, commutator(second, differentiated)) + commutator(
        second, commutator(first, differentiated)
    )


jacobi_matrix = sp.expand(
    kinetic_tensor(lambda_i, lambda_j, lambda_k)
    + kinetic_tensor(lambda_i, lambda_k, lambda_j)
    + kinetic_tensor(lambda_j, lambda_k, lambda_i)
)
show("Jacobi sum before the flavour trace", jacobi_matrix)


def kallen(x: sp.Expr, y: sp.Expr, z: sp.Expr) -> sp.Expr:
    return x**2 + y**2 + z**2 - 2 * x * y - 2 * x * z - 2 * y * z


m1, m2, m3 = sp.symbols("m_1 m_2 m_3", positive=True)
sigma = sp.Symbol("S", integer=True, positive=True)
amplitude = sp.Symbol(r"\mathcal{M}", complex=True)

s12_min = (m1 + m2) ** 2
s12_max = (m_a - m3) ** 2
s23_center = m2**2 + m3**2 + (
    (s12 + m2**2 - m1**2) * (m_a**2 - s12 - m3**2)
) / (2 * s12)
s23_half_width = sp.sqrt(
    kallen(s12, m1**2, m2**2) * kallen(m_a**2, s12, m3**2)
) / (2 * s12)
s23_minus = s23_center - s23_half_width
s23_plus = s23_center + s23_half_width
width_prefactor = 1 / (256 * sp.pi**3 * m_a**3 * sigma)

show("Lower s_12 boundary", s12_min)
show("Upper s_12 boundary", s12_max)
show("Lower s_23 boundary", s23_minus)
show("Upper s_23 boundary", s23_plus)
show("Dalitz differential width", width_prefactor * sp.Abs(amplitude) ** 2)

width_iii = width_prefactor.subs(sigma, sp.factorial(3)) * sp.Integral(
    sp.Abs(amplitude_iii) ** 2 * (s23_plus - s23_minus),
    (s12, s12_min, s12_max),
).subs({m1: m_i, m2: m_i, m3: m_i})
width_iij = width_prefactor.subs(sigma, sp.factorial(2)) * sp.Integral(
    sp.Abs(amplitude_iij) ** 2 * (s23_plus - s23_minus),
    (s12, s12_min, s12_max),
).subs({m1: m_i, m2: m_i, m3: m_j})
width_ijk = width_prefactor.subs(sigma, 1) * sp.Integral(
    sp.Abs(amplitude_ijk_2d) ** 2,
    (s23, s23_minus, s23_plus),
    (s12, s12_min, s12_max),
).subs({m1: m_i, m2: m_j, m3: m_k})

show("Width into three identical mesons", width_iii)
show("Width into two identical and one different meson", width_iij)
show("Width into three different mesons", width_ijk)

# Diagnostics for combinatorics, Mandelstam reduction, and phase-space factors.
assert sp.simplify(amplitude_iii - 6 * M * T_iii_P) == 0
assert sp.simplify(
    amplitude_iij_dot - (K * T_iij_A * (2 * d3 - d1 - d2) + 2 * M * T_iij_P)
) == 0
assert sp.simplify(
    amplitude_iij
    - (K * T_iij_A * (m_a**2 + 3 * m_j**2 - 3 * s12) / 2 + 2 * M * T_iij_P)
) == 0
assert sp.simplify(
    s12 + s13 + s23 - (m_a**2 + 2 * m_i**2 + m_j**2)
).subs(dalitz_iij) == 0
assert sp.simplify(s23_plus - s23_minus - 2 * s23_half_width) == 0
assert not amplitude_ijk_2d.has(s13)
assert jacobi_matrix == 0
