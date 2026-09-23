"""Derive charge-basis ingredients for isospin-broken K pi rescattering.

Invocation:
    .venv/bin/python scripts/derive_kpi_isospin_breaking.py

The Condon-Shortley convention is |pi+> = -|1,+1>. The charge basis is
(K+ pi0, K0 pi+) and the isospin basis is (I=1/2, I=3/2), both with I3=+1/2.
"""

import sympy as sp
from sympy.physics.wigner import clebsch_gordan


def show(label: str, expression: sp.Expr | sp.MatrixBase) -> None:
    print(f"{label}:\n{sp.latex(expression)}\n")


half = sp.Rational(1, 2)
sqrt2 = sp.sqrt(2)
sqrt3 = sp.sqrt(3)

# SymPy returns coefficients for abstract |K,m_K>|pi,m_pi>. The physical
# |pi+> state has an additional minus sign in the adopted convention.
cg_half_kp_pi0 = clebsch_gordan(half, 1, half, half, 0, half)
cg_half_k0_pip = -clebsch_gordan(half, 1, half, -half, 1, half)
cg_three_kp_pi0 = clebsch_gordan(sp.Rational(1, 2), 1, sp.Rational(3, 2), half, 0, half)
cg_three_k0_pip = -clebsch_gordan(
    sp.Rational(1, 2), 1, sp.Rational(3, 2), -half, 1, half
)

# Rows are charge states and columns are I=1/2,3/2 amplitudes.
C = sp.Matrix(
    [
        [cg_half_kp_pi0, cg_three_kp_pi0],
        [cg_half_k0_pip, cg_three_k0_pip],
    ]
)

T_half, T_three = sp.symbols("T^(1/2) T^(3/2)", complex=True)
T_isospin = sp.diag(T_half, T_three)
T_charge = sp.simplify(C * T_isospin * C.T)

rho_charged_neutral, rho_neutral_charged = sp.symbols(
    "rho_(K+pi0) rho_(K0pi+)", real=True, nonnegative=True
)
rho_charge = sp.diag(rho_charged_neutral, rho_neutral_charged)
rho_isospin = sp.simplify(C.T * rho_charge * C)

P_half, P_three = sp.symbols("P_(1/2) P_(3/2)", complex=True)
P_isospin = sp.Matrix([P_half, P_three])
P_charge = sp.simplify(C * P_isospin)

c_u, c_d = sp.symbols("c_u c_d", real=True)
c_zero = sp.simplify((c_u + c_d) / 2)
c_three = sp.simplify((c_u - c_d) / 2)

m_K_plus = sp.Rational(493677, 10**6)
m_K_zero = sp.Rational(497611, 10**6)
m_pi_plus = sp.Rational(13957039, 10**8)
m_pi_zero = sp.Rational(1349768, 10**7)
threshold_1 = m_K_plus + m_pi_zero
threshold_2 = m_K_zero + m_pi_plus

show("Physical Clebsch-Gordan transformation C", C)
show("C^T C", sp.simplify(C.T * C))
show("Charge-basis scattering matrix", T_charge)
show("Physical phase space rotated into the isospin basis", rho_isospin)
show("Charge-basis production vector", P_charge)
show("Isoscalar source coefficient c_0", c_zero)
show("Isovector source coefficient c_3", c_three)
show("K+ pi0 threshold in GeV", threshold_1)
show("K0 pi+ threshold in GeV", threshold_2)
show("Threshold splitting in GeV", threshold_2 - threshold_1)

expected_C = sp.Matrix(
    [[1 / sqrt3, sp.sqrt(sp.Rational(2, 3))], [sp.sqrt(sp.Rational(2, 3)), -1 / sqrt3]]
)
expected_T_charge = sp.Matrix(
    [
        [T_half / 3 + 2 * T_three / 3, sqrt2 * (T_half - T_three) / 3],
        [sqrt2 * (T_half - T_three) / 3, 2 * T_half / 3 + T_three / 3],
    ]
)
expected_rho_off_diagonal = sqrt2 * (
    rho_charged_neutral - rho_neutral_charged
) / 3

assert sp.simplify(C - expected_C) == sp.zeros(2)
assert sp.simplify(C.T * C - sp.eye(2)) == sp.zeros(2)
assert sp.simplify(T_charge - expected_T_charge) == sp.zeros(2)
assert sp.simplify(rho_isospin[0, 1] - expected_rho_off_diagonal) == 0
assert sp.simplify(
    rho_isospin.subs(rho_neutral_charged, rho_charged_neutral)
    - rho_charged_neutral * sp.eye(2)
) == sp.zeros(2)
