---
title: Isoscalar VPPP projection onto the charged three-pion P wave
status: derived
created: 2026-10-27
updated: 2026-10-27
bibtex_keys: [GomezDumm:2010ra, Alda:2026nrd, Hoferichter:2025lcz]
provenance: analytical-derivation
---

# Assumptions and inputs

This derivation projects the five local $VPPP$ operators of
[@GomezDumm:2010ra] onto $V_i\to\pi^+\pi^-\pi^0$ at cubic pion order. It
uses exact isospin, equal on-shell pion masses, tree level, and the common
couplings $g_i$ of the published $SU(3)$ basis. The vector flavor generator is
$T_i=\operatorname{diag}(t_n,t_n,t_s)$, with the tensor leg projected as
$$
V_i^{\mu\nu}=\frac{q^\mu\varepsilon_i^\nu-q^\nu\varepsilon_i^\mu}{M_i},
\qquad q\cdot\varepsilon_i=0.
$$
Pion fields and normalization follow [@Alda:2026nrd]; the source operators
and literal coupling normalization follow the
[five-operator extract](../../extracts/GomezDumm:2010ra/vppp-basis-and-omega-three-pion.md).
The overall vector LSZ/source phase and photon--vector mixing factor are not
included; they multiply this local projection commonly and do not change its
coefficient rank.

# SymPy derivation

The executable
[projection script](../../../scripts/derive_vppp_isoscalar_three_pion.py)
expands all six pion assignments in the five operator words, contracts the
explicit $u,d$ flavor matrices, uses $p_i^2=m_\pi^2$ and momentum
conservation, and extracts the coefficients of $s=q^2$ and $m_\pi^2$.
The scalar divides out the common tensor
$i\epsilon(p_+,p_-,p_0,\varepsilon_i)$ and includes the source factor
$g_i/M_V$. It gives
$$
\mathcal F_{VPPP}^{(3\pi)}(s)=
\frac{16t_n}{M_i f_\pi^3M_V}
\left[(g_1-g_2-g_3)(s-3m_\pi^2)
+3m_\pi^2(2g_4+g_5)\right].
$$
The $s$ and $m_\pi^2$ coefficient rows, after removing the common
normalization, are
$$
\begin{pmatrix}
1&-1&-1&0&0\\
-3&3&3&6&3
\end{pmatrix}.
$$
This matrix has rank two. Thus the projected local amplitude uses exactly
the combinations $g_1-g_2-g_3$ and $2g_4+g_5$. Its three null directions
can be chosen as $(1,1,0,0,0)$, $(1,0,1,0,0)$, and $(0,0,0,1,-2)$.

At one fixed $s$, the coefficient row has rank one, so only the single
mass-weighted combination in the displayed amplitude is measured. At two
distinct energies and nonzero $m_\pi$, the local polynomial has rank two in
the five $g_i$. In the chiral limit the mass-suppressed operators disappear
from this projection and the rank returns to one. The contact term has no
Dalitz-angle dependence, so its $P$-wave projection equals the same scalar.

For ideal $\omega$ mixing, $t_n=1/\sqrt2$ and $M_i=M_\omega$. At
$s=M_\omega^2$ the result reproduces the local term in Appendix B.1 of
[@GomezDumm:2010ra], including its $8\sqrt2/(M_\omega f_\pi^3M_V)$
normalization. Ideal $\phi$ and exact-isospin $\rho^0$ flavor projections
vanish in this charged three-pion channel.

# Implication for the isoscalar TFF

The local $VPPP$ contribution to the isoscalar $\gamma^*\to3\pi$ P wave
therefore depends on only those two coefficient combinations. The on-shell
$\omega$ branching fraction used by [@GomezDumm:2010ra] supplies one
fixed-energy combination; its reported $2g_4+g_5$ value additionally uses
their short-distance relations and $G_V=F_V/F^2$ assumption. It is not a
second independent measurement from the same rate.

An off-shell energy-dependent amplitude could in principle separate the two
combinations through the distinct $s$ and mass terms. This is only the rank
of the isolated local polynomial. It does not establish that the fitted
dispersive normalization $a(q^2)$ separates them after the WZW subtraction,
vector-resonance pieces, and inelastic conformal terms are included. The
updated global analysis in [@Hoferichter:2025lcz] supplies a fitted
normalization and model variants but not a complex $f_1(s,q^2)$ grid or full
correlated parameter replicas; these remain prerequisites for a numerical
refit and a TFF spectral update, as detailed in the
[source audit](pion-tff-spectral-input-source-audit.md).

# Reproducibility

Run:

    .venv/bin/python scripts/derive_vppp_isoscalar_three_pion.py

The script prints the intermediate projections, the reduced amplitude,
coefficient matrix, rank and nullspace using SymPy's LaTeX printer. Its
focused checks reproduce the ideal-$\omega$ source expression, verify the
strange and isovector zeros, and establish the fixed-energy, two-energy and
chiral-limit ranks. This is a tree-level local-vertex result, not a fit to
the dispersive $e^+e^-\to3\pi$ data or a photon--vector source calculation.
