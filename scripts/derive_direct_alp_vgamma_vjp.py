#!/usr/bin/env python3
"""Derive direct physical-ALP VJP source insertions.

Invocation from the repository root:
    .venv/bin/python scripts/derive_direct_alp_vgamma_vjp.py

Natural units are used.  The calculation is exact at linear order in a/f_a,
at the standard electromagnetic-compatible Goldstone vacuum, and before any
integration by parts or on-shell reduction.  The sole on-shell check is the
external antisymmetric-tensor-vector contraction in O_7.
"""

import sympy as sp
try:
    from scripts.symbolic_utils import show
except ModuleNotFoundError:  # direct ``python scripts/...py`` invocation
    from symbolic_utils import show



def component_matrix(prefix: str) -> sp.Matrix:
    """Return a fully general 3x3 component matrix."""
    entries = sp.symbols(
        " ".join(f"{prefix}_{row}{column}" for row in "uds" for column in "uds")
    )
    return sp.Matrix(3, 3, entries)


def component_zero(expression: sp.Expr) -> bool:
    """Test a component-expanded scalar identity."""
    return sp.expand(expression) == 0


I = sp.I
e, f_a, B_0 = sp.symbols("e f_a B_0", nonzero=True)
M_V = sp.Symbol("M_V", positive=True)

# Fully general component matrices: no diagonality, Hermiticity, or mutual
# commutativity of the matrix products is assumed.  Their scalar entries commute.
T_V = component_matrix("TV")
Q = component_matrix("Q")
m_a_inverse_C = component_matrix("m_a_inverse_C")
c_q_AP = component_matrix("cAP")
m_q = component_matrix("mq")

f_plus = -2 * e * Q
u_physical = -2 * m_a_inverse_C / f_a  # coefficient multiplying partial_mu a
chi_physical = -4 * I * B_0 * (
    (c_q_AP + m_a_inverse_C) * m_q + m_q * (c_q_AP + m_a_inverse_C)
) / f_a  # multiplying a

# SymPy derives the flavor factors directly from the source substitutions.
derivative_raw = sp.trace((T_V * f_plus + f_plus * T_V) * u_physical)
derivative_expected = (
    4 * e / f_a * sp.trace((T_V * Q + Q * T_V) * m_a_inverse_C)
)
mass_raw = I * sp.trace((T_V * f_plus + f_plus * T_V) * chi_physical)
mass_expected = (
    -8
    * e
    * B_0
    / f_a
    * sp.trace(
        (T_V * Q + Q * T_V)
        * ((c_q_AP + m_a_inverse_C) * m_q + m_q * (c_q_AP + m_a_inverse_C))
    )
)

# Compact abstract forms for deterministic readable output.  Component matrices
# above independently establish each identity without assuming matrix commutation.
TVa = sp.MatrixSymbol("T_V", 3, 3)
Qa = sp.MatrixSymbol("Q", 3, 3)
m_a_inverse_C_abs = sp.MatrixSymbol(r"\mathbb M_a^{-1}(C)", 3, 3)
cAPa = sp.MatrixSymbol("c_q^{AP}", 3, 3)
mqa = sp.MatrixSymbol("m_q", 3, 3)
anti_TQ = TVa * Qa + Qa * TVa
anti_cAP_m_a_inverse_C_mq = (
    (cAPa + m_a_inverse_C_abs) * mqa
    + mqa * (cAPa + m_a_inverse_C_abs)
)

show(
    "Alda inverse mixing operator applied to C",
    m_a_inverse_C_abs,
)
show("Physical u_mu coefficient multiplying partial_mu a", -2 * m_a_inverse_C_abs / f_a)
show(
    "Physical chi_- coefficient multiplying a",
    -4 * I * B_0 * anti_cAP_m_a_inverse_C_mq / f_a,
)
show("Electromagnetic f_+ coefficient multiplying F^{mu nu}", -2 * e * Qa)
show(
    "Derived common O_1,O_2,O_5,O_6,O_7 flavor coefficient",
    4 * e * sp.Trace(anti_TQ * m_a_inverse_C_abs) / f_a,
)
show("O_3 scalar product i(-2e)(-4iB_0/f_a)", sp.simplify(I * (-2 * e) * (-4 * I * B_0 / f_a)))
show(
    "Derived O_3 mass flavor coefficient",
    -8 * e * B_0 * sp.Trace(anti_TQ * anti_cAP_m_a_inverse_C_mq) / f_a,
)

# Keep all five derivative Lorentz contractions distinct.  Indexed objects are
# formal component tensors; no equations of motion or integration by parts enter.
mu, nu, rho, sigma, alpha = sp.symbols("mu nu rho sigma alpha", integer=True)
eps4 = sp.IndexedBase("epsilon")
V = sp.IndexedBase("V")
F = sp.IndexedBase("F")
dda = sp.IndexedBase("partial partial a")
da = sp.IndexedBase("partial a")
dV = sp.IndexedBase("partial V")
a = sp.Symbol("a")

lorentz = {
    1: eps4[mu, nu, rho, sigma] * V[mu, nu] * F[rho, alpha] * dda[alpha, sigma],
    2: eps4[mu, nu, rho, sigma] * V[mu, alpha] * F[rho, sigma] * dda[alpha, nu],
    5: eps4[mu, nu, rho, sigma] * dV[alpha, mu, nu] * F[rho, alpha] * da[sigma],
    6: eps4[mu, nu, rho, sigma] * dV[alpha, mu, alpha] * F[rho, sigma] * da[nu],
    7: eps4[mu, nu, rho, sigma] * dV[sigma, mu, nu] * F[rho, alpha] * da[alpha],
}
common = 4 * e * sp.Trace(anti_TQ * m_a_inverse_C_abs) / f_a
for operator_index in (1, 2, 5, 6, 7):
    show(f"Inserted O_{operator_index} (distinct Lorentz contraction)", common * lorentz[operator_index])
show(
    "Inserted O_3 mass structure",
    -8
    * e
    * B_0
    * sp.Trace(anti_TQ * anti_cAP_m_a_inverse_C_mq)
    * eps4[mu, nu, rho, sigma]
    * V[mu, nu]
    * F[rho, sigma]
    * a
    / f_a,
)

# At the assumed electromagnetic-compatible vacuum f_-=0, hence O_4=0.
f_minus = sp.zeros(3)
chi_plus = component_matrix("chip")
O4_flavor = I * sp.trace(T_V * (f_minus * chi_plus - chi_plus * f_minus))
show("O_4 flavor factor at f_-=0", O4_flavor)

# Sole on-shell reduction: epsilon_{mu nu rho sigma} p^sigma V^{mu nu}=0
# for V^{mu nu}=(p^mu eps^nu-p^nu eps^mu)/M_V.  Compute every free-rho component.
p = sp.symbols("p0:4")
polarization = sp.symbols("varepsilon0:4")
c7_contraction = sp.Matrix(
    [
        sp.simplify(
            sum(
                sp.LeviCivita(i, j, r, k)
                * p[k]
                * (p[i] * polarization[j] - p[j] * polarization[i])
                / M_V
                for i in range(4)
                for j in range(4)
                for k in range(4)
            )
        )
        for r in range(4)
    ]
)
show("On-shell O_7 contraction for free rho=0,1,2,3", c7_contraction)

# Diagnostic assertions are separate from the derivation output.
assert component_zero(derivative_raw - derivative_expected)
assert component_zero(mass_raw - mass_expected)
assert O4_flavor == 0
assert c7_contraction == sp.zeros(4, 1)
print("All symbolic assertions passed.")
