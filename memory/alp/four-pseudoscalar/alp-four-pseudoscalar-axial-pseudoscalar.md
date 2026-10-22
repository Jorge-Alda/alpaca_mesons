---
title: Axial and pseudoscalar ALP contributions to four-pseudoscalar vector exchange
status: derived
created: 2026-10-23
updated: 2026-10-23
bibtex_keys: [Alda:2026nrd, Ecker:1988te, Ecker:1989yg, Ruiz-Femenia:2003jdx]
provenance: analytical-derivation
---

# Scope and inputs

This record derives the axial and pseudoscalar ALP contributions to
$a(p_a)\to P_1(p_1)P_2(p_2)P_3(p_3)P_4(p_4)$ in the vector-exchange
families catalogued in the [process overview](alp-four-pseudoscalar-vector-exchange.md).
The matrices $c_q^A,c_q^P$ are already in the no-ALP--meson-mixing basis and
are not rotated again [@Alda:2026nrd]. The retained interactions are the
original seven $VJP$ and four $VVP$ odd operators [@Ruiz-Femenia:2003jdx] and
the leading even-parity $G_V$ interaction [@Ecker:1988te], using the
preliminary-plus-$G_V$ fit point and its caveats in the [fit Lagrangian
record](../../rcht/odd-parity/preliminary-fit-lagrangian-and-null-directions.md).
We work at first order in the ALP sources, zero electromagnetic coupling,
tree level, natural units, metric $(+---)$, and $\epsilon^{0123}=+1$.
Pseudoscalar generators obey $\operatorname{tr}(\lambda_i\lambda_j)=2\delta_{ij}$;
vector generators obey $\operatorname{tr}(T_RT_S)=\delta_{RS}$.

The executable derivation is
[`scripts/derive_alp_four_pseudoscalar.py`](../../../scripts/derive_alp_four_pseudoscalar.py):

```text
.venv/bin/python scripts/derive_alp_four_pseudoscalar.py
.venv/bin/python scripts/derive_alp_four_pseudoscalar.py --test
```

# Absence of an $aPP$ vertex

Neither $c_q^A$ nor $c_q^P$ produces an $a\Phi^2$ vertex from the leading
chiral Lagrangian. The axial part of $u_\mu$ is even in $\Phi$, while the
source-free part is odd, so their cross term in $\langle u_\mu u^\mu\rangle$
has odd meson degree. The pseudoscalar source enters $\chi_+$ through the
sine of the anticommutator map and likewise has odd meson degree. The possible
$aP$ bilinear is removed by the prescribed no-mixing basis; that cancellation
is an input of the basis prescription, not an identity for arbitrary
independent source matrices. Therefore the $aPP$--pseudoscalar--$VPP$--$VPP$
family in the overview is a $c_q^{VS}$ contribution and is absent here.

# Odd $aVV$ vertex

For vector generators $T_R,T_S$ and virtualities $s_R=q_R^2$ and
$s_S=q_S^2$, set
$$
D^a_{RS}=\frac{1}{f_a}\operatorname{tr}\!\left(\{T_R,T_S\}c_q^A\right),
\qquad
P^a_{RS}=-\frac{1}{f_a}\operatorname{tr}\!\left(\{T_R,T_S\}\{c_q^P,\chi\}\right),
\qquad \chi=2B_0m_q.
$$
Including both vector assignments, the original fitted odd operators give
$$
\mathcal N^a_{RS}(s_R,s_S)
=4D^a_{RS}\left[d_1m_a^2-d_3(m_a^2-s_R-s_S)\right]-16d_2P^a_{RS}.
$$
The routed $d_4$ contribution vanishes. The numerator annihilates the four
original-fit coefficient-null directions; this is specific to this vertex
and does not determine a null-sensitive channel amplitude or its uncertainty.

# Axial $aVP$ vertex

The nonzero axial even-parity source with one vector and one pseudoscalar
comes from
$$
\mathcal L_{G_V}=\frac{iG_V}{\sqrt2}
\langle V_{\mu\nu}u^\mu u^\nu\rangle.
$$
At the required orders,
$u_\mu^{(a,A)}=2(\partial_\mu a)c_q^A/f_a$ and
$u_\mu^{(P)}=-(\partial_\mu\Phi)/f_\pi+O(\Phi^3)$. Antisymmetry of
$V_{\mu\nu}$ combines the two cross terms. For
$V_{\mu\nu}=T_RV_{R,\mu\nu}$ and $\Phi=\lambda_i^\dagger P_i$,
$$
\left.\mathcal L_{G_V}\right|_{aVP}=
-\frac{2iG_V}{\sqrt2 f_af_\pi}\,
V_R^{\mu\nu}\partial_\mu a\,\partial_\nu P_i\,
\operatorname{tr}\!\left(T_R[c_q^A,\lambda_i^\dagger]\right).
$$
This is a Lagrangian coefficient before Fourier and vertex factors.

# Sewn scalar amplitudes

The full tensor responses are sewn to the canonical odd tensor amplitude
$\mathcal N_{RS}\epsilon(\varepsilon_R,\varepsilon_S,q_R,q_S)/(M_RM_S)$.
Fourier signs use incoming $a$ waves $e^{-ip_ax}$ and outgoing meson waves
$e^{+ip_ix}$. Define $s_{ij}=(p_i+p_j)^2$ and
$$
\epsilon(p_1,p_2,p_3,p_4)=\epsilon_{\mu\nu\rho\sigma}
p_1^\mu p_2^\nu p_3^\rho p_4^\sigma.
$$
Each vector index is summed over the real orthonormal vector nonet. A
denominator means $q^2-M_R^2+i0$; no widths are inserted.

For the odd $aVV$ source followed by two $VPP$ vertices,
$$
\begin{aligned}
\mathcal M_{aVV}={}&\epsilon(p_1,p_2,p_3,p_4)
\sum_{(ij|kl)=(12|34),(13|24),(14|23)}\sigma_{ij|kl}\sum_{R,S}
\\[-2pt]&\times
\frac{-G_V^2\mathcal N^a_{RS}(s_{ij},s_{kl})
\operatorname{tr}(T_R[\lambda_i^\dagger,\lambda_j^\dagger])
\operatorname{tr}(T_S[\lambda_k^\dagger,\lambda_l^\dagger])}
{2f_\pi^4(s_{ij}-M_R^2+i0)(s_{kl}-M_S^2+i0)},
\\
&(\sigma_{12|34},\sigma_{13|24},\sigma_{14|23})=(+1,-1,+1).
\end{aligned}
$$
The odd numerator already includes both vector assignments at its vertex.

For the axial $aVP$ source followed by an odd $VVP$ vertex and a $VPP$
vertex, let $i,j$ be ordered distinct labels and $k<l$ the remaining two.
Then
$$
\mathcal M_{aVP}=\epsilon(p_1,p_2,p_3,p_4)
\sum_{\substack{i\ne j\\\{k,l\}=\{1,2,3,4\}\setminus\{i,j\}\\k<l}}
\operatorname{sgn}(i,j,k,l)\sum_{R,S}
\frac{G_V^2\mathcal N^j_{RS}((p_a-p_i)^2,s_{kl})
\operatorname{tr}(T_R[c_q^A,\lambda_i^\dagger])
\operatorname{tr}(T_S[\lambda_k^\dagger,\lambda_l^\dagger])}
{f_\pi^3f_a((p_a-p_i)^2-M_R^2+i0)(s_{kl}-M_S^2+i0)}.
$$
The sign is the parity of $(i,j,k,l)$ relative to $(1,2,3,4)$; the vector
momentum at the odd vertex is crossed, $q_R=-(p_a-p_i)$. The strong odd
numerator $\mathcal N^j_{RS}$ is the same VVP polynomial with source mass
squared $m_j^2$ and
$$
D^j_{RS}=-\frac{1}{2f_\pi}
\operatorname{tr}(\{T_R,T_S\}\lambda_j^\dagger),\qquad
P^j_{RS}=\frac{1}{2f_\pi}
\operatorname{tr}(\{T_R,T_S\}\{\lambda_j^\dagger,\chi\}).
$$
Thus the retained axial/pseudoscalar-source amplitude is
$\mathcal M_{aVV}+\mathcal M_{aVP}$; the first term contains both $c_q^A$
and $c_q^P$, the second only $c_q^A$. Flavor restrictions and admitted
final states are recorded separately in the [selection-rule note](alp-four-pseudoscalar-selection-rules.md).

These are generic labeled-channel amplitudes. Numerical predictions still
require a specified final channel, evaluation of the vector-flavor sums, and
four-body phase-space integration with the appropriate identical-particle
factor.
