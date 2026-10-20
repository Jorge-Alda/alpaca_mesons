---
title: No-mixing ALP amplitudes to two neutral vector mesons
status: derived
created: 2026-10-20
updated: 2026-10-20
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx, Ecker:1988yg]
provenance: analytical-derivation
---

# Inputs and scope

Use the original seven VJP and four VVP operators and the canonical tensor
convention in [the preliminary fit](odd-preliminary-fit.md). The ALP is
specified directly in Alda:2026nrd's no-mixing basis by diagonal dimensionless
couplings $c_q^A$ and $c_q^P$. No field rotation is performed here. Let
$\chi=2B_0m_q$ be the diagonal quark-mass matrix, and let $M_V=0.770$ GeV
retain its role as the odd-operator normalization scale. The even-vector
source includes the fitted one-trace coefficient $z_V$; its signed residues
are the fitted $R_i$. The vector masses and mixing angle are the fitted
near-ideal values.

# ALP insertions in this basis

Define the ALP derivative-source matrix and ALP mass-source matrix by
$$
\mathcal U_a=\frac{2c_q^A}{f_a},\qquad
\mathcal Y_a=-\frac{2}{f_a}\{c_q^P,\chi\}.
$$
These follow directly from $c_q^A=-\mathbb M_a^{-1}(C)$ and
$c_q^P=-c_q^{AP}-\mathbb M_a^{-1}(C)$; they are not computed from UV
couplings. The no-mixing prescription already includes the meson pole
contribution. Do not add a second LSZ meson-pole term.

# Analytic amplitude

For a neutral vector pair define
$$
D_{ij}=\operatorname{tr}(T_iT_j\mathcal U_a),\qquad
P_{ij}=\operatorname{tr}(T_iT_j\mathcal Y_a),\qquad
q_{ij}=q_i\cdot q_j=\frac{m_a^2-s-t}{2},
$$
where $s=q_i^2$ and $t=q_j^2$. The original VJP operators give the ordered
scalar coefficient
$$
\mathcal V_{ij}(s,t)=\frac{4}{M_V}\left\{
D_{ij}\left[c_1(t+q_{ij})+c_2(s+q_{ij})-c_5q_{ij}-c_6s\right]
-2c_3P_{ij}\right\}.
$$
The original VVP operators, including both vector assignments, give
$$
\mathcal P_{ij}(s,t)=4\left[d_1m_a^2-d_3(m_a^2-s-t)\right]D_{ij}
-16d_2P_{ij}.
$$
For neutral diagonal generators, the $c_4,c_7,d_4$ contractions vanish.
The eleven coefficients are the canonical original fit coefficients; no
coefficient-coordinate conversion is needed.

The canonical amplitude is
$$
\mathcal M_{ij}=G_{ij}\,\epsilon_{\mu\nu\rho\sigma}
\varepsilon_i^\mu\varepsilon_j^\nu q_i^\rho q_j^\sigma,
\qquad
G_{ij}=G^{\rm VVP}_{ij}+G^{\rm conv}_{ij},
$$
where, with $s=M_i^2$ and $t=M_j^2$,
$$
G_{ij}=\frac{\mathcal P_{ij}(s,t)-2\sqrt2e^2\left[
R_j\mathcal V_{ij}(s,t)+R_i\mathcal V_{ji}(t,s)\right]}{M_iM_j}.
$$
Here $R_i$ is the signed photon-vector residue from the fitted even sector,
including its one-trace $z_V$ contribution. The conversion factor follows
from the tensor mixing vertex $\sqrt2eR_i$ and the external tensor
normalization $1/M_i$. It is a pole-residue insertion, not the unamputated
response $W(s)$ evaluated at its pole. Both assignments are required for
$i\ne j$; the expression is Bose symmetric. For identical vectors the two
terms coincide and are both retained.

Applying this expression to $(i,j)=(\omega,\omega),(\omega,\phi),
(\phi,\phi)$ gives the three channel coefficients. For identical vectors,
the two VJP conversion terms coincide and are both retained. For
$\omega\phi$, they use the two ordered generator/momentum assignments.
These expressions are analytic in the no-mixing ALP couplings, the eleven
fitted odd coefficients, and the fitted even-sector residues and masses. No
individual coefficient error is inferred from its minimum-norm fit
representative; the original-model null audit establishes fixed-input
blindness of the neutral $aVV$ family.

# Conventions and limitations

The result is tree level, uses narrow on-shell daughter vectors, the
near-ideal fitted mass branch, and the fitted single-trace even source.
Electromagnetic conversion is retained through two photon-vector mixing
insertions on the VJP amplitude. Loops, widths in internal propagators,
kinetic breaking and additional local odd operators are omitted. The
executable [`alp_vv_width.py`](../../../scripts/implementation/alp_vv_width.py)
reads the selected preliminary-fit JSON and evaluates the displayed trace
expressions directly.
