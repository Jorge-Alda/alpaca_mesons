---
title: Consolidated ALP charged-pion production amplitude for the isovector TFF
status: mixed
created: 2026-10-28
updated: 2026-10-28
bibtex_keys: [Alda:2026nrd, Bai:2025fvl, Wess:1971yu, Ecker:1988te, Ecker:1989yg, Ruiz-Femenia:2003jdx, GomezDumm:2010ra, Messerli:2025rnv]
provenance: analytical-derivation
---

# Scope and conventions

This record collects the completed leading WZW, VJP/VVP, and local VPPP
pieces for the charged-pion production wave entering the isovector ALP TFF.
It organizes results derived in the linked records; it does not replace their
operator-by-operator derivations. Use natural units, metric $(+---)$,
$\epsilon^{0123}=+1$, and Alda:2026nrd's no-mixing ALP sources
[@Alda:2026nrd]. Set
$$
s=(p_++p_-)^2,\qquad t=(k+p_-)^2,\qquad u=(k+p_+)^2,
\qquad p_a=k+p_++p_-,\qquad k^2=y,
$$
where $k$ is the outgoing photon momentum (or the crossed virtual-current
momentum). Define the scalar by
$$
\mathcal M(a\to\gamma^*(k)\pi^+\pi^-)
=\epsilon_{\mu\nu\rho\sigma}\varepsilon^{*\mu}
p_+^\nu p_-^\rho k^\sigma\,
\mathcal A_{a\pi\pi\gamma^*}(s,t,u;y).
$$
The physical charged-pion mass is used in phase space. The leading odd
resonance derivations use strong isospin and diagonal ALP sources; their
stated source and vector-sector assumptions remain in force below.

# Coherent production amplitude

At the retained order, organize the scalar as
$$
\mathcal A_{a\pi\pi\gamma^*}
=\mathcal A_{\rm WZW}
+\mathcal A_{\rm VJP/VVP}
+e\,\mathcal F_{\rm VPPP}.
$$
These are distinct contributions in the declared operator basis and must be
added coherently. The first is fixed by the gauged anomaly and the source/LSZ
convention. The latter two contain independent homogeneous resonance
couplings and are not fixed by WZW matching [@Ruiz-Femenia:2003jdx;
@GomezDumm:2010ra].

## WZW contact: matrix derivation

Use Alda:2026nrd's no-mixing axial coupling matrix $\widehat c_q^A$ and the
full mass-dependent $\chi$ from the linked
[Alda:2026nrd source record](../../../extracts/Alda:2026nrd/appendix-a-expansion.md).
Keep its light-flavor entries independent; there is no isospin approximation
in the WZW trace below.

1. **Start from the full gauged functional.** Use Eq. (A11) of
   [the full-functional extract](../../extracts/Bai:2025fvl/background-fields-wzw.md)
   [@Bai:2025fvl]. In the no-mixing basis, the photon and ALP occur only in
   the external left/right currents $A_L,A_R$; the pions occur only through
   the $U$-dependent terms. The counterterm $\Gamma_c$ is source-only and
   has no pion fields. Retain the terms of Eq. (A11) quadratic in external
   currents and quadratic in pion fields; the source-linear term that gives
   $\gamma3P$ alone is not sufficient for this vertex.

2. **Insert the photon and ALP through the external currents.** In the
   no-mixing basis the physical ALP has no pseudoscalar-meson overlap.
   The currents in Eq. (A11) are
   $$
   A_L=-eQA+\widehat c_q^A\frac{da}{f_a},\qquad
   A_R=-eQA-\widehat c_q^A\frac{da}{f_a}.
   $$
   This source assignment follows Alda:2026nrd's
   $u_\mu|_{U=1}=r_\mu-l_\mu=2\widehat c_q^A\partial_\mu a/f_a$.
   Keep the pion dependence only in the $U$ factors and expand to second
   order in the charged-pion fields.

3. **Expand the four two-current groups.** From Eq. (A11), the terms
   $-\frac{i}{2}[(A_L\alpha)^2-(A_R\beta)^2]$ and
   $i(A_LUA_RU^\dagger\alpha^2-A_RU^\dagger A_LU\beta^2)$ give the first
   two groups below. The two terms containing $dA_{L/R}$ and one
   Maurer--Cartan form give the other groups. Fourier transformation is
   included. With the common factor $N_Ce/(48\pi^2 f_af_\pi^2)$ removed,
   the four contributions are
   $$
   \begin{aligned}
   \text{first group: }&-2\operatorname{tr}(Q\lambda_+^\dagger\widehat c_q^A
   \lambda_-^\dagger)
   +2\operatorname{tr}(Q\lambda_-^\dagger\widehat c_q^A\lambda_+^\dagger),\\
   \text{second group: }&+2\operatorname{tr}(Q\widehat c_q^A\lambda_+^\dagger
   \lambda_-^\dagger)
   -2\operatorname{tr}(Q\widehat c_q^A\lambda_-^\dagger\lambda_+^\dagger)\\
   &+2\operatorname{tr}(Q\lambda_+^\dagger\lambda_-^\dagger
   \widehat c_q^A)
   -2\operatorname{tr}(Q\lambda_-^\dagger\lambda_+^\dagger
   \widehat c_q^A),\\
   \text{third group: }&-2\operatorname{tr}(Q\lambda_+^\dagger\widehat c_q^A
   \lambda_-^\dagger)
   +2\operatorname{tr}(Q\lambda_-^\dagger\widehat c_q^A\lambda_+^\dagger),\\
   \text{fourth group: }&+2\operatorname{tr}(Q\widehat c_q^A\lambda_+^\dagger
   \lambda_-^\dagger)
   -2\operatorname{tr}(Q\widehat c_q^A\lambda_-^\dagger\lambda_+^\dagger)\\
   &+2\operatorname{tr}(Q\lambda_+^\dagger\lambda_-^\dagger
   \widehat c_q^A)
   -2\operatorname{tr}(Q\lambda_-^\dagger\lambda_+^\dagger
   \widehat c_q^A).
   \end{aligned}
   $$
   The first and third traces are equal, and the second and fourth are
   equal. Their sum is four times the six-term trace.

4. **Restore the common factor.** The $N_C/(48\pi^2)$ coefficient in
   Eq. (A11) combines with the external-current and quadratic-pion factors
   $e/(f_af_\pi^2)$. The result is
   $$
   \mathcal A_{\rm WZW}(s,t,u;y)
   =-\frac{N_C e}{12\pi^2 f_a f_\pi^2}
   \operatorname{tr}\!\left[Q\left(
   \widehat c_q^A\lambda_+^\dagger\lambda_-^\dagger
   -\widehat c_q^A\lambda_-^\dagger\lambda_+^\dagger-\lambda_+^\dagger\widehat c_q^A\lambda_-^\dagger
   +\lambda_+^\dagger\lambda_-^\dagger\widehat c_q^A
   +\lambda_-^\dagger\widehat c_q^A\lambda_+^\dagger
   -\lambda_-^\dagger\lambda_+^\dagger\widehat c_q^A
   \right)\right].
   $$
   The common sign is fixed by matching the same Eq. (A11) to the signed
   $\pi^0\to\gamma\gamma$ WZW scalar in the [canonical TFF record](../../rcht/odd-parity/radiative-tff-amplitude.md).
   In the script's exterior basis, the raw source scalar is positive in
   $E=\epsilon(q_1,\varepsilon_1,q_2,\varepsilon_2)$; the tensor dictionary
   gives $E=-R$, and the indexed external-pion phase is applied before
   comparison with the canonical $R$ scalar. This comparison fixes the common
   A11 multiplier to $-1$. It is not a free orientation in the resulting
   project amplitude. This sign is also required for the WZW/mixing term to
   combine with the photon Jacobian into the invariant $a\to\gamma\gamma$
   coefficient in the [photon record](../../photon/alp-photon-invariant-from-mixing.md).
   Bai:2025fvl does not state a numerical $\epsilon^{0123}$ convention and
   notes a sign difference among cited WZW conventions; that is why the
   signed $\pi^0$ match is part of the convention conversion
   [@Bai:2025fvl; @Alda:2026nrd].

   The decay-constant conversion is $f_\pi^{\rm Bai}=\sqrt2 f_\pi$ in the
   project convention, obtained by matching the kinetic prefactors in
   Eqs. (A10) and Alda:2026nrd. The charged-pion exponent coefficients then
   agree, so there is no extra $\sqrt2$ in the boxed amplitude. The tensor
   obeys the photon Ward identity also at $k^2=y\ne0$. This local coefficient
   is independent of $s,t,u,y$, and its charge-stripped production wave is
   $\mathcal F_{\rm WZW}=\mathcal A_{\rm WZW}/e$.

## VJP/VVP resonance contribution

The retained seven VJP and four VVP operators, with the leading $G_V$ VPP
vertex, give the neutral-rho exchange production scalar. In the notation of
[the source-resolved partial-wave derivation](alp-pion-vector-and-production-form-factors.md),
$$
\mathcal A_{\rm VJP/VVP}(s,t,u;y)
=\frac{N_s(s,y)}{s-M_\rho^2+i0}
+\frac{N_t(t,y)}{t-M_\rho^2+i0}
+\frac{N_u(u,y)}{u-M_\rho^2+i0}.
$$
The numerators $N_j$ are the displayed contractions of the VJP coefficients
$c_i$, VVP coefficients $d_i$, ALP source matrices, and the full neutral
photon response $W_\gamma(y)$ in that record. Retain both crossed pion
assignments. The $c_4$ trace vanishes for the charged channel in the exact
strong-isospin limit. The real-photon expression in
[alp-photon-pp-amplitude.md](../photon-two-pseudoscalar/alp-photon-pp-amplitude.md)
is a separate $y=0$ limit; the formula above is the required virtual-current
amplitude.

## Local VPPP resonance contribution

The five local $VPPP$ operators carry independent coefficients $g_i$ and are
converted to the photon current with
$$
W_\gamma(y)=\sum_i\frac{R_iT_i}{M_i^2-y-i0}.
$$
Their charged-channel contribution to the charge-stripped production wave is
$$
\mathcal F_{\rm VPPP}(s,y)
=\mathcal F^{VPPP}_{a\pi\pi\gamma^*}(s,y),
$$
where the full trace formula and its angular projection are given in
[alp-vppp-isovector-tff.md](alp-vppp-isovector-tff.md)
[@GomezDumm:2010ra; @Ecker:1988te; @Ecker:1989yg]. This piece is local at
the odd $VPPP$ vertex but has the photon-vector response dependence on $y$.
It is not the WZW contact and must not be used as a substitute for it.

# Production P wave and isovector cut

For each contribution, use the same reduced $P$-wave projector,
$$
f_1(s,y)=\frac34\int_{-1}^{1}dz\,(1-z^2)
\mathcal A_{a\pi\pi\gamma^*}(s,t(z),u(z);y),
\qquad
\mathcal F_{a\pi\pi\gamma^*}(s,y)=\frac{f_1(s,y)}{e}.
$$
The WZW scalar is angle independent, so its projection is unchanged:
$$
f_1^{\rm WZW}(s,y)=\mathcal A_{\rm WZW},\qquad
\mathcal F^{\rm WZW}_{a\pi\pi\gamma^*}(s,y)
=-\frac{N_C}{12\pi^2 f_a f_\pi^2}
\operatorname{tr}\!\left[Q\left(\widehat c_q^A\lambda_+^\dagger
\lambda_-^\dagger-\widehat c_q^A\lambda_-^\dagger\lambda_+^\dagger
-\lambda_+^\dagger\widehat c_q^A\lambda_-^\dagger
+\lambda_+^\dagger\lambda_-^\dagger\widehat c_q^A
+\lambda_-^\dagger\widehat c_q^A\lambda_+^\dagger
-\lambda_-^\dagger\lambda_+^\dagger\widehat c_q^A\right)\right].
$$
In the raw decay scalar convention, the VPPP term is therefore
$e\mathcal F^{VPPP}$. The complete charge-stripped projected wave is the sum
of WZW, VJP/VVP, and VPPP projected waves. Use the same phase and
normalization for the pion vector form factor
$F_\pi^V$, with $F_\pi^V(0)=1$ [@Messerli:2025rnv].

The single-cut relation is
$$
\operatorname{Im}_s F_{a\gamma^*\gamma^*}(s,y)
=\frac{s\sigma_\pi^3(s)}{96\pi}
[F_\pi^V(s)]^*\mathcal F_{a\pi\pi\gamma^*}(s,y),
\qquad
\sigma_\pi(s)=\sqrt{1-\frac{4m_{\pi^\pm}^2}{s}}.
$$
For the per-ordering double-spectral convention, the density is
$$
\rho_a^{I=1}(s,y)=\frac{s\sigma_\pi^3(s)}{192\pi}
\operatorname{Im}\!\left\{[F_\pi^V(s)]^*
\mathcal F_{a\pi\pi\gamma^*}(s,y)\right\}.
$$
Insert this density in the cutoff double dispersion relation in
[alp-isovector-tff-spectral-construction.md](alp-isovector-tff-spectral-construction.md),
then add the exchanged photon ordering once. The factor of two between the
single-cut and per-ordering density conventions is already accounted for
there; do not add a second pion or photon symmetry factor.

The WZW coefficient above is the leading local production input. A physical
spectral density requires phase-compatible pion rescattering and the chosen
spectral/subtraction prescription. The contact term alone is not a complete
doubly-virtual TFF, and the stable rho-pole approximation in
[alp-isovector-doubly-virtual-tff.md](alp-isovector-doubly-virtual-tff.md)
must not be added to the same rho strength without a matching prescription.
No numerical TFF is fixed here because the phase-compatible amplitudes,
vector responses over the selected range, and cutoff/subtraction inputs have
not been selected together.

# Calculation and input records

- [Full gauged WZW functional](../../extracts/Bai:2025fvl/background-fields-wzw.md) — Eq. (A11), external-current prescription, and source-only counterterm.
- [Executable WZW expansion](../../../scripts/derive_wzw_external_currents.py) — SymPy expansion of the four two-current groups in Eq. (A11). Run `.venv/bin/python scripts/derive_wzw_external_currents.py`.
- [VJP/VVP production wave](alp-pion-vector-and-production-form-factors.md) — symbolic source contractions, neutral-rho exchange, and $P$-wave projection.
- [VPPP contribution](alp-vppp-isovector-tff.md) — five local operators, photon response, and its pion-cut insertion.
- [Pion vector form factor and cut normalization](alp-pion-vector-and-production-form-factors.md) — $F_\pi^V$, angular weight, and single-cut normalization.
- [Double-spectral construction](alp-isovector-tff-spectral-construction.md) — per-ordering density, dispersion relation, and numerical-input requirements.
