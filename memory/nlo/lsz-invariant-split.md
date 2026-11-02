---
title: Invariant split of the ALP–meson LSZ overlap
status: derived
created: 2026-11-02
updated: 2026-11-02
bibtex_keys: [Alda:2026nrd, Kaiser:1998ds]
provenance: analytical-derivation
---

# Invariant split of the ALP–meson overlap

## Assumptions and inputs

Use the real neutral meson vector $P=(\pi_3,\eta_8,\eta_0)^T$ and the local quadratic convention

$$
\mathcal L_{\rm quad}=
\frac12\partial_\mu P^T K\partial^\mu P
-\frac12P^T M^2P
+\frac12(\partial a)^2-\frac12m_a^2a^2
+\partial_\mu a\,k^T\partial^\mu P-a\,\mu^TP.
$$

$K$ and $M^2$ are symmetric, $K$ and $M^2-sK$ are invertible at the pole under consideration, and only first order in the ALP–meson mixing is retained. The local NLO contributions to the meson matrices are those in [the neutral meson matrix record](neutral-meson-matrices.md); the local ALP–meson terms are to be projected from [the quadratic NLO Lagrangian](quadratic-lagrangian.md). The operator basis and ALP source conventions are those of [@Kaiser:1998ds; @Alda:2026nrd]. Loops and the displayed higher-order remainder in the NLO operator basis are omitted.

## Derivation

Let $s=p^2$ and define the meson and mixed inverse-propagator blocks

$$
A(s)=sK-M^2,\qquad B(s)=sk-\mu.
$$

To first order in the off-diagonal block, the ALP-pole meson overlap is

$$
x(s)=-A(s)^{-1}B(s)=-(sK-M^2)^{-1}(sk-\mu).
$$

The combination

$$
\mathcal I=M^2K^{-1}k-\mu
$$

does not depend on $s$. The identity

$$
sk-\mu=\mathcal I-(sK-M^2)K^{-1}k
$$

then gives the split

$$
\boxed{x(s)=-(sK-M^2)^{-1}\mathcal I-K^{-1}k.}
$$

Now make a general ALP-dependent meson field shift $P=P'-r a$. Direct substitution into the quadratic form gives

$$
k'=k-Kr,\qquad \mu'=\mu-M^2r.
$$

Therefore

$$
M^2K^{-1}k'-\mu'
=M^2K^{-1}(k-Kr)-(\mu-M^2r)
=M^2K^{-1}k-\mu=\mathcal I.
$$

Thus $\mathcal I$ is invariant under this field shift. The other term in the split transforms as

$$
-(K^{-1}k')=-K^{-1}k+r,
$$

so $x'=x+r$, as required by $P=P'-ra$. The resolvent term is unchanged. Under a regular meson coordinate change $P=R Q$, $\mathcal I$ transforms as a source covector, while $x$ and $K^{-1}k$ transform as field vectors; the expression is consequently basis covariant.

For an arbitrary source linearized as $\mathcal J=V_a a+V^TP$, the same field shift gives $V_a'=V_a-V^Tr$. The direct-plus-pole coefficient is

$$
V_a+V^Tx(s),
$$

and is unchanged because $V_a'+V^Tx'=V_a+V^Tx$. In the split form,

$$
V_a+V^Tx(s)
=\underbrace{V_a-V^TK^{-1}k}_{\text{shift invariant}}
-V^T(sK-M^2)^{-1}\mathcal I.
$$

The first bracket is invariant since both $V_a$ and $K^{-1}k$ shift by the compensating amount. Equivalently, if the source is expressed in the shifted no-mixing fields, its direct coefficient is $V_a'=V_a-V^Tr$; it is the bare direct coefficient by itself that is basis dependent. The invariant source $\mathcal I$ and the pole resolvent give the compact NLO generalization of the LO LSZ organization. Its explicit U(3) flavor form follows from the general-matrix projection of the quadratic ALP–meson Lagrangian below.

## U(3) source through the retained local NLO terms

Use a complete normalized $U(3)$ generator basis and arbitrary Hermitian flavor matrices $c_q^A$ and $c_q^P$. No commutativity with the diagonal quark-mass spurion is assumed. Define

$$
(r_A)_i=\frac{f_\pi}{f_a}\operatorname{tr}(\lambda_i c_q^A).
$$

Projection of every derivative bilinear gives $k=-Kr_A$, including the $L_4,L_5,L_{18},\Lambda_1$ terms. Hence $-K^{-1}k=r_A$. With the established source combinations $c_q^{AP}=c_q^A-c_q^P$ and $c_{GA}=c_G+\operatorname{tr}c_q^A$, the invariant source is

$$
\mathcal I_i=\frac{f_\pi}{f_a}\operatorname{tr}\!\left[\lambda_i\left(C+\Delta C_{\rm NLO}\right)\right],
\qquad C=-\frac12\{\chi_0,c_q^{AP}\}-\frac{m_0^2}{3}c_{GA}\mathbf 1,
$$

and, to first order in each retained local LEC,

$$
\begin{aligned}
\Delta C_{\rm NLO}={}&-
\frac{4L_8}{f_\pi^2}\{\chi_0,\{\chi_0,c_q^{AP}\}\}
-\frac{8L_6}{f_\pi^2}\operatorname{tr}(\chi_0)\{\chi_0,c_q^{AP}\}
-\frac{16L_7}{f_\pi^2}\operatorname{tr}(\chi_0c_q^{AP})\chi_0
\\
&-\frac{\Lambda_2}{3}\left[
c_{GA}\chi_0+\operatorname{tr}(\chi_0c_q^{AP})\mathbf 1
\right]
+\frac{8L_{25}}{f_\pi^2}\left[
c_{GA}\chi_0^2+\operatorname{tr}(\chi_0^2c_q^{AP})\mathbf 1
\right].
\end{aligned}
$$

For diagonal spurions, the anticommutators reduce to the corresponding products. For general $c_q^{AP}$, the LO term is not $-\chi_0c_q^{AP}$; the $L_8$ term is the nested anticommutator shown; and the $L_6$ term also requires an anticommutator. The $L_7$, $\Lambda_2$, and $L_{25}$ trace structures retain the displayed orderings. The axial matrix cancels from $\mathcal I_i$ after combining $M^2K^{-1}k$ with $\mu$. The derivative LECs $L_4,L_5,L_{18},\Lambda_1$ affect $K$ and $k$ together and do not enter this invariant source.

The ALP-pole overlap is therefore

$$
\langle0|P|a\rangle
=r_A-(m_a^2K-M^2)^{-1}\mathcal I
=r_A+(M^2-m_a^2K)^{-1}\mathcal I.
$$

This matrix result applies to the retained local operators without a diagonal-coupling assumption. It does not include loops or omitted higher-order operators.

## Reproducibility

The general field-shift derivation is [`scripts/derive_nlo_lsz_invariant_split.py`](../../scripts/derive_nlo_lsz_invariant_split.py); the general noncommuting flavor-matrix derivation is [`scripts/derive_nlo_lsz_source_projection.py`](../../scripts/derive_nlo_lsz_source_projection.py). Run

```text
.venv/bin/python scripts/derive_nlo_lsz_invariant_split.py
.venv/bin/python scripts/derive_nlo_lsz_source_projection.py
```

The first script constructs the inverse-propagator block matrix and derives the invariant split and direct-plus-pole cancellation. The second differentiates the LO and NLO trace operators using noncommuting $3\times3$ matrix symbols and tests the ordered formulas on a non-diagonal Hermitian example. Its tests establish $k=-Kr_A$, cancellation of $c_q^A$, recovery of the general LO $C$, and the NLO spurion structures.
