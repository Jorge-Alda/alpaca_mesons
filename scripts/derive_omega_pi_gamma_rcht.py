#!/usr/bin/env python3
"""Reproduce the published on-shell omega -> pi0 gamma algebra.

Invocation from the repository root:
    .venv/bin/python scripts/derive_omega_pi_gamma_rcht.py

This is an exact SymPy reduction of the literature formulas in
Ruiz-Femenia:2003jdx, Eqs. (27)--(29), not a fresh operator derivation.
Natural units are used: masses and decay constants are in one common mass
unit, the scalar amplitude coefficient has inverse-mass units, and the width
has mass units.  Ideal mixing, a common operator-normalization mass M_V,
tree-level/NWA kinematics, and canonical external residue Z_omega=1 are the
published assumptions used here.
"""

import sympy as sp
try:
    from scripts.symbolic_utils import show
except ModuleNotFoundError:  # direct ``python scripts/...py`` invocation
    from symbolic_utils import show



e, F, F_V, M_V, M_omega, m_pi = sp.symbols(
    "e F F_V M_V M_omega m_pi", positive=True
)
N_C = sp.Symbol("N_C", positive=True, integer=True)
alpha = sp.Symbol("alpha", positive=True)
c1, c2, c3, c4, c5, c6, c7 = sp.symbols("c_1:8", real=True)
d1, d2, d3, d4 = sp.symbols("d_1:5", real=True)

# Published Eq. (27): coefficients of
# i epsilon_{alpha beta rho sigma} epsilon_omega^alpha
# epsilon_gamma^beta q^rho k^sigma.  The minus sign in rho exchange is part
# of the published scalar coefficient.
direct_mass_combination = (
    (c2 - c1 + c5 - 2 * c6) * M_omega**2
    + (c1 + c2 + 8 * c3 - c5) * m_pi**2
)
rho_mass_combination = d3 * M_omega**2 + (d1 + 8 * d2 - d3) * m_pi**2
G_direct = 2 * sp.sqrt(2) * e * direct_mass_combination / (M_omega * M_V * F)
G_rho = -4 * e * F_V * rho_mass_combination / (M_V**2 * M_omega * F)
G_uneliminated = G_direct + G_rho

show("Eq. 27 direct mass combination", direct_mass_combination)
show("Eq. 27 direct scalar coefficient", G_direct)
show("Eq. 27 rho-exchange mass combination", rho_mass_combination)
show("Eq. 27 rho-exchange scalar coefficient (including its sign)", G_rho)
show("Eq. 27 total uneliminated scalar coefficient", G_uneliminated)

# The five published odd-parity OPE equations, Eq. (19).  SymPy solves and
# applies all five; c4, c7, and d4 remain absent from this channel.
ope_equations = (
    sp.Eq(4 * c3 + c1, 0),
    sp.Eq(c1 - c2 + c5, 0),
    sp.Eq(c5 - c6, N_C * M_V / (64 * sp.pi**2 * sp.sqrt(2) * F_V)),
    sp.Eq(d1 + 8 * d2, -N_C * M_V**2 / (64 * sp.pi**2 * F_V**2) + F**2 / (4 * F_V**2)),
    sp.Eq(d3, -N_C * M_V**2 / (64 * sp.pi**2 * F_V**2) + F**2 / (8 * F_V**2)),
)
for index, equation in enumerate(ope_equations, start=1):
    show(f"OPE relation {index}", equation)

ope_solution = sp.solve(ope_equations, (c3, c2, c6, d1, d3), dict=True)[0]
for symbol in (c3, c2, c6, d1, d3):
    show(f"Solved OPE substitution for {sp.sstr(symbol)}", ope_solution[symbol])

direct_combination_ope = sp.simplify(direct_mass_combination.subs(ope_solution))
rho_combination_ope = sp.simplify(rho_mass_combination.subs(ope_solution))
G_direct_ope = sp.factor(G_direct.subs(ope_solution))
G_rho_ope = sp.factor(G_rho.subs(ope_solution))
G_constrained = sp.factor(G_uneliminated.subs(ope_solution))

show("Direct mass combination after all OPE substitutions", direct_combination_ope)
show("Rho mass combination after all OPE substitutions", rho_combination_ope)
show("Constrained direct contribution", G_direct_ope)
show("Constrained rho-exchange contribution", G_rho_ope)
show("Sum after OPE substitutions", G_constrained)

G_eq28 = e / F_V * (
    N_C * M_omega / (8 * sp.pi**2 * F)
    - F * M_omega / (2 * M_V**2) * (1 + m_pi**2 / M_omega**2)
)
show("Published constrained scalar coefficient, Eq. 28", G_eq28)

# Generic two-body formula for this scalar-amplitude normalization.  The even-
# parity short-distance relation F_V=sqrt(2)F and alpha=e^2/(4*pi) are needed
# only to put the result in the exact printed Eq. (29) form.
Gamma_generic = sp.Symbol("abs(G)**2", positive=True) * (M_omega**2 - m_pi**2) ** 3 / (
    96 * sp.pi * M_omega**3
)
show("Generic V -> P gamma width for scalar coefficient G", Gamma_generic)

Gamma_from_eq28 = sp.factor(
    G_eq28**2 * (M_omega**2 - m_pi**2) ** 3 / (96 * sp.pi * M_omega**3)
)
Gamma_fv_relation = sp.factor(Gamma_from_eq28.subs(F_V, sp.sqrt(2) * F))
Gamma_alpha = sp.factor(Gamma_fv_relation.subs(e**2, 4 * sp.pi * alpha))
show("Width after inserting Eq. 28", Gamma_from_eq28)
show("Width after F_V = sqrt(2) F", Gamma_fv_relation)
show("Width after alpha = e^2/(4 pi)", Gamma_alpha)

Gamma_eq29 = alpha * M_omega / 192 * (1 - m_pi**2 / M_omega**2) ** 3 * (
    N_C * M_omega**2 / (4 * sp.pi**2 * F**2)
    - M_omega**2 / M_V**2 * (1 + m_pi**2 / M_omega**2)
) ** 2
show("Published width, Eq. 29", Gamma_eq29)

# Assertions diagnose exact equality but do not replace the displayed steps.
assert sp.simplify(G_constrained - G_eq28) == 0
assert sp.simplify(Gamma_alpha - Gamma_eq29) == 0
assert not any(symbol in G_uneliminated.free_symbols for symbol in (c4, c7, d4))
print("All symbolic assertions passed.")
