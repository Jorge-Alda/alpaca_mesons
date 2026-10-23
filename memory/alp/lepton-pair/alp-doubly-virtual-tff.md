---
title: Initial tree-level doubly virtual ALP transition form factor
status: derived
created: 2026-10-23
updated: 2026-10-23
bibtex_keys: [Alda:2026nrd, Bauer:2020jbp, Ecker:1988te, Ecker:1989yg, Messerli:2025rnv, Ruiz-Femenia:2003jdx]
provenance: analytical-derivation
---

# Inputs and scope

This starts the construction of the scalar TFF needed in the two-photon
contribution to $a\to\ell^+\ell^-$. It is a tree-level neutral-vector
resonance expression with the original eleven odd-parity couplings, not a
completed hadronic TFF. Use Messerli:2025rnv's unnormalized TFF convention
[@Messerli:2025rnv, Eq. (9)]:
$$
i\mathcal M=i e^2 R\,F_{a\gamma^*\gamma^*}(s_1,s_2),\qquad
R=\epsilon_{\mu\nu\rho\sigma}\varepsilon_1^\mu\varepsilon_2^\nu
q_1^\rho q_2^\sigma,
$$
with $s_i=q_i^2$, $p_a+q_1+q_2=0$, metric $(+---)$, and
$\epsilon^{0123}=+1$. The scalar $F$ has dimension $\mathrm{GeV}^{-1}$.
Messerli:2025rnv's normalized form factor is a separate quantity,
$\widetilde F_{a\gamma^*\gamma^*}=F_{a\gamma^*\gamma^*}/F_{a\gamma\gamma}$,
with $F_{a\gamma\gamma}=F_{a\gamma^*\gamma^*}(0,0)$ [@Messerli:2025rnv,
Eq. (7)].

The intermediate source-contraction calculation in the existing RChT scripts
uses $E=\epsilon(q_1,\varepsilon_1,q_2,\varepsilon_2)=-R$. Therefore its
scalar obeys $F=-G$; the final form below is converted to the $R$ tensor
convention before use.

The source matrices are already in Alda:2026nrd's no-mixing basis. Write
$S=\mathbb M_a^{-1}(\mathcal C)$, so that
$$
U_a=-\frac{2S}{f_a},\qquad
Y_a=\frac{2}{f_a}\{c_q^{AP}+S,\chi\},\qquad
\widehat c_\gamma=C_{a\gamma\gamma}(m_a).
$$
The source-EOM-compatible LO profile obeys
$Y_a+2m_a^2U_a-2Z_a\mathbf1=0$. These source and contact inputs follow the
[LO source/LSZ assembly](../../rcht/odd-parity/odd-alp-source-lsz-assembly.md)
and [photon invariant](../photon/alp-photon-invariant-from-mixing.md).
In the no-mixing basis the meson-pole WZW contribution is already included in
$\widehat c_\gamma$; do not add it again. The positive photon coefficient and
canonical epsilon convention are those in the
[photon/action phase record](../photon/alp-photon-odd-action-phase-matching.md)
[@Bauer:2020jbp].

The odd vertices use the canonical $c_i,d_i$ convention of
[the neutral SM TFF derivation](../../rcht/odd-parity/radiative-tff-amplitude.md),
with the ALP source pair $U_a,Y_a$ replacing the meson sources. Their
operator basis and tensor normalizations are from
[@Ruiz-Femenia:2003jdx; @Ecker:1988te; @Ecker:1989yg]. The vector response is
diagonal in the physical neutral-vector basis and is written
$$
W_i(s)=\frac{R_i}{M_i^2-s},
$$
where each signed residue $R_i$ is an input from the selected even-parity
photon source, as in the [vector photon-response record](../../rcht/odd-parity/vector-stage4-amplitudes.md).
Stable poles, tree level, one ALP, diagonal neutral sources,
and leading electromagnetic order are assumed. No vector widths or additional
local odd-parity terms are included.

# Ordered numerators and assembly

The executable
[`derive_alp_doubly_virtual_tff.py`](../../../scripts/derive_alp_doubly_virtual_tff.py)
uses the existing explicit Lorentz contractions and derives the source
numerators before inserting photon responses. For a vector on photon 1,
$$
\begin{aligned}
V_{i\gamma}(s_1,s_2)=\frac{2}{M_V}\Big[&
\operatorname{tr}(T_i Q U_a)\big((c_1+c_2-c_5)m_a^2
 +(-c_1+c_2+c_5-2c_6)s_1+(c_1-c_2+c_5)s_2\big)\\
&-4c_3\operatorname{tr}(T_i QY_a)\Big].
\end{aligned}
$$
For the two vector legs,
$$
P_{ij}(s_1,s_2)=4\Big[\operatorname{tr}(T_iT_jU_a)
\big((d_1-d_3)m_a^2+d_3(s_1+s_2)\big)
-4d_2\operatorname{tr}(T_iT_jY_a)\Big].
$$
The full retained scalar, in Messerli:2025rnv's $R$ tensor convention, is
$$
\boxed{\begin{aligned}
F_{a\gamma^*\gamma^*}(s_1,s_2)=&+\frac{\widehat c_\gamma}{4\pi^2f_a}\\
&+2\sqrt2\sum_i\left[W_i(s_1)V_{i\gamma}(s_1,s_2)
                  +W_i(s_2)V_{i\gamma}(s_2,s_1)\right]\\
&-2\sum_{i,j}W_i(s_1)W_j(s_2)P_{ij}(s_1,s_2).
\end{aligned}}
$$
The first two terms in brackets are the two VJP photon placements. The last
line includes both VVP vector assignments through the ordered sum. Keep the
local polynomial's photon labels attached until after selecting the matching
$W_i(s)$; swapping an already inserted pole would assign it to the wrong
photon. For diagonal neutral sources $P_{ij}=P_{ji}$ after simultaneous
exchange of $s_1,s_2$ and the vector labels, so the full expression is Bose
symmetric.

## Derivation and checks

Run
```text
.venv/bin/python scripts/derive_alp_doubly_virtual_tff.py
```
The script solves the tensor-field quadratic equation to obtain $W_i(s)$,
contracts the VJP and VVP structures in SymPy, strips their common tensor
$E$, and then performs the ordered response sums. The ordered real-photon
limits are
$$
V_{i\gamma}(0,0)=\frac{2}{M_V}\left[(c_1+c_2-c_5)m_a^2
\operatorname{tr}(T_iQU_a)-4c_3\operatorname{tr}(T_iQY_a)\right],
$$
$$
P_{ij}(0,0)=4\left[(d_1-d_3)m_a^2\operatorname{tr}(T_iT_jU_a)
-4d_2\operatorname{tr}(T_iT_jY_a)\right].
$$
The $Y_a$ term is left un-factored by $m_a^2$, so the expression remains
regular at $m_a=0$. The derivation keeps $U_a$ and $Y_a$ independent; no
equal-spurion identity is assumed.

# Limits and next inputs

This expression is a source-resolved resonance starting point, not the full
doubly virtual TFF. It omits hadronic cuts and widths, any enlarged odd-parity
operator families, and the source-specific short-distance coefficient and
matching prescription. In particular, the large-virtuality behavior of the
single-vector response must not be treated as the QCD asymptotic limit. The
[lepton-pair plan](alp-lepton-pair-calculation-plan.md) specifies the later
double-spectral/dispersive completion, matching, lepton-loop insertion, and
local ALP-lepton contribution. Numerical use first requires a chosen ALP
source profile, vector residues and pole prescription, and a declared
truncation uncertainty.
