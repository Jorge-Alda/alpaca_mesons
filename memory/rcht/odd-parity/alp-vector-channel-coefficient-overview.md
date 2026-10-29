---
title: Even- and odd-sector coefficient directions for ALP vector channels
status: derived
created: 2026-10-28
updated: 2026-10-28
bibtex_keys: [Alda:2026nrd, Bauer:2020jbp, Ecker:1988te, Ecker:1989yg, GomezDumm:2010ra, Messerli:2025rnv, Ruiz-Femenia:2003jdx]
provenance: analytical-derivation
---

# Scope and conventions

This overview maps the coefficient directions entering the retained vector-meson amplitudes for $a\to3P$, $a\to\gamma PP$, $a\to VV$, $a\to VPP$, $a\to\gamma\gamma$, and the doubly virtual $a\gamma^*\gamma^*$ transition form factor (TFF). “Enters” means that the direction occurs in the displayed amplitude at the stated truncation. It does not mean that the direction is uniquely determined by a fit or that the channel prediction is complete.

ALP axial and pseudoscalar sources remain independent no-mixing inputs in Alda:2026nrd's convention [@Alda:2026nrd]. The odd-parity vector operators include $VJP$, $VVP$, and $VPPP$, with coefficients $c_1,\ldots,c_7$, $d_1,\ldots,d_4$, and $g_1,\ldots,g_5$, respectively [@Ruiz-Femenia:2003jdx; @GomezDumm:2010ra]. No numerical prescription for the $g_i$ is selected here. The even-sector coupling $G_V$ controls the leading strong $VPP$ and $VPa$ vertices; the photon-vector response is specified by residues $R_i$, which combine the fitted source coefficients $F_V,z_V,y_V$. These even inputs are independent of the odd-parity coefficients.

# Channel map

The cells show coefficient combinations only; common source traces, pole responses, and flavor projections are suppressed. For $a\to VV$, the combinations are read from the on-shell width amplitude with $s=M_i^2$, $t=M_j^2$, and $q_{ij}=(m_a^2-M_i^2-M_j^2)/2$. Since $m_a$ is free, its kinematic weight does not identify a fixed coefficient combination. For the two-pion production terms, $j=s,t,u$ labels the ALP, crossed-pion, and crossed-pion exchange, respectively. In the even column, $R_i=\operatorname{tr}(T_iS_\gamma)$, with $S_\gamma=F_VQ+z_V\{\chi_+,Q\}+y_V\operatorname{tr}(\chi_+Q)\mathbf1$.

| Channel or contribution | Even operators | VVP | VJP | VPPP |
|---|---|---|---|---|
| $a\to3P$ | $G_V^2$ | $0$ | $0$ | $0$ |
| $a\to\gamma PP$: $s,t,u$ channels | $G_V,\ R_i$ | $d_1-d_3,\ d_3,\ d_2$ | $c_1+c_2-c_5,\ -c_1+c_2+c_5-2c_6,\ c_3,\ c_4$ | $0$ |
| $a\to\gamma PP$: VPPP | $R_i$ | $0$ | $0$ | $g_1-g_3,\ g_2,\ g_4,\ g_5$ |
| $a\to VV$ | $R_i$ | $d_1-d_3,\ d_3,\ d_2$ | $c_1+c_2-c_5,\ -c_1+c_2+c_5-2c_6,\ c_1-c_2+c_5,\ c_3$ | $0$ |
| $a\to VPP$: $s,t,u$ channels | $G_V$ | $d_1-d_3,\ d_3,\ d_2$ | $0$ | $0$ |
| $a\to VPP$: VPPP | $0$ | $0$ | $0$ | $g_1-g_3,\ g_2,\ g_4,\ g_5$ |
| $a\to\gamma\gamma$ (pole/contact) | $R_iR_j$ | $d_1-d_3,\ d_2$ | $c_1+c_2-c_5,\ c_3$ | $0$ |
| TFF: $\omega,\phi$ VJP| $R_i$ | $0$ | $c_1+c_2-c_5,\ -c_1+c_2+c_5-2c_6,\ c_1-c_2+c_5,\ c_3$ | $0$ |
| TFF: $\omega,\phi$ VVP poles | $R_i R_j$ | $d_1-d_3,\ d_3,\ d_2$ | $0$ | $0$ |
| TFF: mixed $\rho^0$--$\omega/\phi$ VVP | $R_\rho R_i$ | $d_1-d_3,\ d_3,\ d_2$ | $0$ | $0$ |
| TFF: isovector $2\pi$ cut | $G_V,\ R_i$ | $d_1-d_3,\ d_3,\ d_2$ | $c_1+c_2-c_5,\ -c_1+c_2+c_5-2c_6,\ c_1-c_2+c_5,\ c_3$ | $g_1-g_3,\ g_2,\ g_4,\ g_5$ |
| TFF: isoscalar $3\pi$ cut (* outside provisional scope) | $G_V,\ R_i$ | $d_1-d_3,\ d_3,\ d_2$ | $c_1+c_2-c_5,\ -c_1+c_2+c_5-2c_6,\ c_1-c_2+c_5,\ c_3$ | $g_1-g_2-g_3,\ 2g_4+g_5$ |

The $a\to\gamma PP$ $c_4$ term is in crossed pseudoscalar-source vertices and vanishes for the charged-pion projection in exact strong isospin. The $a\to VV$ width uses both VJP assignments; $c_4,c_7,d_4$ vanish in the neutral on-shell projection. The VPPP entries for $a\to VPP$ and $a\to\gamma PP$ follow the source-dependent $aVPP$ vertex, with the latter followed by photon-vector conversion. The provisional TFF scope uses the isovector $2\pi$ component, narrow $\omega,\phi$ poles, and mixed $\rho^0$--$\omega/\phi$ VVP terms with the $\rho$ response dressed; it does not use a separate stable-$\rho$ pole. The isoscalar $3\pi$ row is retained for completeness but remains outside that scope. Its VPPP combinations are the isoscalar vector-current projection. Pole and cut contributions require spectral matching before combination [@Messerli:2025rnv]. See the [channel derivations](../../alp/three-pseudoscalar/a3p-vector-exchange-amplitude.md), [photon-PP amplitude](../../alp/photon-two-pseudoscalar/alp-photon-pp-amplitude.md), [vector-PP amplitude](../../alp/vector-two-pseudoscalar/alp-vpp-amplitude.md), [two-vector amplitude](alp-vv-amplitude.md), [two-vector width](alp-vv-width.md), [isovector TFF construction](../../alp/lepton-pair/alp-isovector-tff-spectral-construction.md), [isoscalar pole decomposition](../../alp/lepton-pair/alp-isoscalar-tff-source-decomposition.md), [mixed $\rho I$ term](../../alp/lepton-pair/alp-tff-mixed-rho-isoscalar-vvp.md), and [isoscalar VPPP projection](isoscalar-vppp-three-pion-projection.md).

# Ranks and simpler spanning bases

Taking the union of the combinations displayed in the channel map, the coefficient-linear-form ranks are $5$ for VJP, $3$ for VVP, and $4$ for VPPP. Equivalent bases with simpler combinations are
$$
\begin{aligned}
\text{VJP: }&\{c_1,\ c_2-c_5,\ c_2-c_6,\ c_3,\ c_4\},\\
\text{VVP: }&\{d_1,\ d_2,\ d_3\},\\
\text{VPPP: }&\{g_1-g_3,\ g_2,\ g_4,\ g_5\}.
\end{aligned}
$$
For VJP, the first three elements span the three independent combinations in $c_1,c_2,c_5,c_6$: the two table combinations $c_1+c_2-c_5$ and $-c_1+c_2+c_5-2c_6$ are recovered by adding/subtracting the displayed basis elements, and $c_1-c_2+c_5$ is also in their span. For VVP, replacing $d_1-d_3,d_3$ by $d_1,d_3$ is an invertible basis change. The isoscalar VPPP combinations $g_1-g_2-g_3$ and $2g_4+g_5$ are already in the listed VPPP span, so they do not raise its rank. These are ranks of the coefficient combinations in the overview, not channel-specific measurement ranks.

# Coverage by the current vector fits

The preliminary low-energy fit covers only the original VJP and VVP coefficients $c_i,d_i$; at its fixed even/source inputs its measured amplitude image has rank seven in the eleven-dimensional coefficient space. Resolving that rank against the overview basis, the fit covers four VJP directions and all three listed VVP directions. The missing eighth direction in the preliminary fit's combined VJP/VVP overview span is
$$
c_1-c_2+c_5.
$$
The measured VJP combinations include $c_1+c_2-c_5$, $-c_1+c_2+c_5-2c_6$, $c_3$, and $c_4$; the fit null vectors in $(c_1,c_2,c_5,c_6)$ annihilate the first two but not the displayed missing combination. Later finite-virtuality $\omega\to\pi^0\ell^+\ell^-$ and pion-TFF fits raise the original odd-coefficient rank to eight by probing $c_1-c_2+c_5$. Their retained tree-level description is poor, so this is an algebraic sensitivity result rather than an accepted coefficient calibration. The other null directions include $c_7$ and $d_4$, outside the overview's rank-eight span. Thus the exported central $c_i,d_i$ values remain a representative, not a determination of the null coordinates. The fit has no $g_i$ or VPPP observable. The separate four-row $G_V$ fit constrains only the leading even-parity $VPP$ coupling, conditionally on the preliminary fit inputs; it supplies no constraint on $c_i,d_i$ or $g_i$.

The VPPP fit is still at the preparation stage, with no numerical result. Its proposed isoscalar $3\pi$ analysis targets only
$$
A_{3\pi}=g_1-g_2-g_3,\qquad B_{3\pi}=2g_4+g_5,
$$
and hence has rank two within the four-dimensional overview VPPP span. The two missing independent VPPP directions can be chosen as $g_2$ and $g_4$ alongside the proposed pair: from $g_1-g_2-g_3$ and $g_2$ one recovers $g_1-g_3$, and from $2g_4+g_5$ and $g_4$ one recovers $g_5$. The cited on-shell $\omega$ rate gives one fixed-energy combination, while the quoted conditional $B_{3\pi}$ value in [@GomezDumm:2010ra] depends on imposed short-distance relations and $G_V=F_V/F^2$; it is not an independent fit constraint. The validated dispersive response and full likelihood are still missing, so even the two intended directions are not yet numerically fitted. In particular, none of the four VPPP combinations needed for the isovector $2\pi$ cut and local $a\to VPP$ and $a\to\gamma PP$ terms has a repository fit result.

# Mesonic observables for the uncovered combinations

The finite-virtuality $\omega\to\pi^0\ell^+\ell^-$ spectrum and pion transition form factor already provide shape sensitivity to $c_1-c_2+c_5$, the direction absent from the preliminary real-photon fit. A renewed determination needs a finite-width, dispersive description that gives an acceptable account of both data sets before its coefficient profile can be treated as a calibration. The $KK\pi$ vector-current form factors also contain this VJP combination multiplying the total invariant mass squared, so their energy dependence offers a separate cross-check [@GomezDumm:2010ra].

The most direct mesonic candidates for the two VPPP directions not targeted by the proposed isoscalar fit are the vector-current Dalitz distributions in
$$
\tau^-\to K^-K^0\pi^0\nu_\tau,\qquad
\tau^-\to K^+K^-\pi^-\nu_\tau.
$$
The first channel contains the local combination $(g_1+2g_2-g_3)(x-y)$, which is $(A_{3\pi}+3g_2)(x-y)$; together with the isoscalar $A_{3\pi}$ projection it can resolve $g_2$. The second channel contains a different mass polynomial in $g_2,g_4,g_5$. Once $A_{3\pi}$, $g_2$, and $B_{3\pi}=2g_4+g_5$ are constrained, its remaining $g_4$ dependence is proportional to $m_\pi^2-m_K^2$, allowing $g_4$ in principle to be separated. The source form factors give this operator-level dependence, but a new fit must recompute the full coefficient Jacobian in the conventions used here and include resonance, rescattering, and vector/axial-current uncertainties. Isovector $e^+e^-\to K\bar K\pi$ distributions related by CVC can add statistics and a cross-check when the isospin component is separated [@GomezDumm:2010ra].

Use energy spectra and Dalitz distributions jointly, retaining the vector-current structure functions rather than only integrated branching fractions. An amplitude term odd under exchange of two Dalitz variables does not by itself imply an odd event-rate asymmetry. These channels are proposed measurements and projections, not fit results; short-distance relations used in the cited analysis must not be imposed as if they were independent experimental constraints.

# Doubly-virtual TFF sensitivity to the common shift

Consider the coefficient displacement
$$
\delta c_2=\delta c_5=\delta c_6=\lambda,
\qquad \delta c_1=\delta c_3=\delta c_4=\delta c_7=0,
\qquad \delta d_i=0.
$$
For the VJP source polynomial attached to a vector on photon 1, define
$$
L(s_1,s_2)=(c_1+c_2-c_5)m_a^2
+(-c_1+c_2+c_5-2c_6)s_1
+(c_1-c_2+c_5)s_2.
$$
The three coefficient variations are respectively
$$
\delta(c_1+c_2-c_5)=0,\qquad
\delta(-c_1+c_2+c_5-2c_6)=0,\qquad
\delta(c_1-c_2+c_5)=0,
$$
so $\delta L(s_1,s_2)=0$ for arbitrary virtualities. The independent mass-source term proportional to $c_3\operatorname{tr}(T_iQY_a)$ does not vary. Consequently each VJP numerator is unchanged, including after exchanging the two photon placements. The VVP numerator contains only $d_1,d_2,d_3$, so it is unchanged as well. The local photon coefficient and the even vector responses are independent of this odd-coefficient displacement. Therefore
$$
\boxed{\delta F_{a\gamma^*\gamma^*}(s_1,s_2)=0}
$$
for arbitrary $s_1,s_2$, independent $U_a,Y_a$, and arbitrary retained vector residues and pole masses. This is exact within the displayed tree-level TFF, not merely a zero at the real-photon point or an equal-spurion result.

This zero is specific to the common $c_2,c_5,c_6$ displacement. The $g_i$ combinations in the VPPP column can change the TFF's hadronic spectral input even though they are absent from the neutral-pole formula.

The result is a coefficient-direction sensitivity check, not a statement that every combination involving $c_2,c_5,c_6$ is invisible: other displacements, including correlated changes involving $c_1$, need not cancel. It also does not promote the pole TFF into a complete hadronic form factor.

The executable check is [the symbolic audit](../../../scripts/audit_alp_tff_c256_direction.py); run it with:
```text
.venv/bin/python scripts/audit_alp_tff_c256_direction.py
```
It varies the displayed VJP and VVP expressions symbolically and asserts zero without assigning numerical sources, residues, or pole parameters. The input TFF expression and its conventions are documented in [the TFF derivation](../../alp/lepton-pair/alp-doubly-virtual-tff.md). The channel-specific tree-zero input is [the three-pseudoscalar audit](../../alp/three-pseudoscalar/odd-parity-three-pseudoscalar-audit.md).
