---
title: One-loop corrections to mixed two-point functions and ALP-meson LSZ mixing
status: derived
created: 2026-11-05
updated: 2026-11-05
bibtex_keys: [Alda:2026nrd, Kaiser:1998ds, Herrera-Siklody:1996tqr]
provenance: analytical-derivation
---

# One-loop corrections to mixed two-point functions and ALP-meson LSZ mixing

## Assumptions and inputs

Use the fields and block notation of the [LSZ invariant-split record](lsz-invariant-split.md): $\Phi_I=(P_i,a)$, with meson blocks $K_{ij},M^2_{ij}$ and mixed blocks $K_{i\hat a},M^2_{i\hat a}$. The local tree-level quadratic Lagrangian and its first-order mixing solution are derived there. The local NLO operators are expanded in the [quadratic NLO record](quadratic-lagrangian.md), while [the singlet-source derivative basis](singlet-source-derivative-basis.md) records the additional $O_{47}$ source structure. The leading ALP-dependent vertices follow the field-redefinition-invariant setup of [@Alda:2026nrd]. The one-loop meson effective action and $U(3)$ nonet loop content follow [@Kaiser:1998ds; @Herrera-Siklody:1996tqr].

Assume perturbative one-loop corrections and first order in ALP-meson mixing. No specific ALP-source matrices or meson masses are specialized. Thus this note derives the general block structure, but not the source-dependent loop integrals. Loop effects on the mixed two-point function are $O(1/f_a)$ when one ALP-dependent vertex is inserted.

## Mixed loop self-energy

Let $H[\bar P,\bar a]$ be the Euclidean quadratic fluctuation operator for meson fields about background fields $\bar P_i,\bar a$. The one-loop meson effective action is

$$
\Gamma_{\rm loop}[\bar P,\bar a]=\frac12\operatorname{Tr}\log H[\bar P,\bar a].
$$

The one-loop mixed two-point function is the second variation

$$
\Sigma_{i\hat a}^{\rm loop}(p)=
\left.\frac{\delta^2\Gamma_{\rm loop}}
{\delta \bar P_i(-p)\,\delta \bar a(p)}\right|_{\bar P=\bar a=0}.
$$

Writing $G=H^{-1}$ and denoting background derivatives of $H$ by subscripts, the Hessian identity gives

$$
\Sigma_{i\hat a}^{\rm loop}
=\frac12\operatorname{Tr}\!\left(GH_{i\hat a}-GH_{\hat a}GH_i\right),
$$

with momentum routing and internal indices understood. The $GH_{i\hat a}$ term is the tadpole from a quartic vertex with external $P_i,a$ and two internal meson legs. The $GH_{\hat a}GH_i$ term is a possible bubble formed from separate cubic vertices, one attached to each external leg. Its presence depends on the vertices. In the ordinary parity-even real-mass vacuum, the leading purely mesonic chiral Lagrangian has no three-pseudoscalar vertex, so this bubble is absent for the standard leading meson-loop calculation unless additional sources or interactions provide the required cubic vertex. The mixed tadpole remains. It gives mass-dependent chiral logarithms and polynomial dependence on $s$; it does not by itself generate a two-meson threshold cut. The analogous second variation with respect to two meson backgrounds gives $\Sigma_{ij}^{\rm loop}(p)$.

The loop functions contain local ultraviolet poles and finite mass- and momentum-dependent terms. The poles are canceled by local counterterms. Chiral logarithms remain after subtraction. Threshold dependence and absorptive parts occur when the relevant bubble topologies are allowed; they should not be presumed for the leading mixed two-point function in the parity-even vacuum. The propagating $U(3)$ singlet contributes when retained. The one-loop effective-action structure follows [@Kaiser:1998ds; @Herrera-Siklody:1996tqr]; evaluating $\Sigma_{i\hat a}^{\rm loop}$ requires specifying the ALP source vertices and checking which cubic vertices survive.

## Full inverse propagator and LSZ eigenvector

Define the sign of the self-energy blocks by writing the momentum-space inverse propagator as

$$
\Gamma(s)=
\begin{pmatrix}
A(s)&b(s)\\
b^T(s)&d(s)
\end{pmatrix},\qquad s=p^2,
$$

with

$$
\begin{aligned}
A_{ij}(s)&=sK_{ij}-M^2_{ij}+\Sigma_{ij}^{\rm loop}(s),\\
b_i(s)&=sK_{i\hat a}-M^2_{i\hat a}+\Sigma_{i\hat a}^{\rm loop}(s).
\end{aligned}
$$

Here $K,M^2$ include the local terms to the chosen order, including the renormalized NLO LECs. If the ALP-like pole eigenvector is $u^{(a)}=(w_i^{(a)},1)^T$ at first order in mixing, the meson row of $\Gamma(s_a)u^{(a)}=0$ gives

$$
w^{(a)}(s_a)=-A^{-1}(s_a)b(s_a).
$$

When the loop terms are removed, this reduces to the LSZ equation and solution in [the LSZ invariant-split record](lsz-invariant-split.md). With loops, $A(s)$ and $b(s)$ depend non-polynomially on $s$ in general, so a constant shift of $K$ and $M^2$ cannot represent the full correction.

For a strict NLO expansion, write $A=A_0+\delta A$ and $b=b_0+\delta b$, where the subscript zero denotes the LO blocks and $\delta A,\delta b$ include local NLO and one-loop pieces. At a fixed pole value $s$, the LO and NLO admixtures are

$$
w_0^{(a)}(s)=-A_0^{-1}(s)b_0(s),\qquad
\delta w^{(a)}(s)=-A_0^{-1}(s)\left[\delta b(s)+\delta A(s)w_0^{(a)}(s)\right].
$$

In the project notation, the block corrections are

$$
\begin{aligned}
\delta A_{ij}(s)&=s\,\delta K_{ij}^{\rm local}-\delta M^2_{ij,\rm local}+\Sigma_{ij}^{\rm loop}(s),\\
\delta b_i(s)&=s\,\delta K_{i\hat a}^{\rm local}-\delta M^2_{i\hat a,\rm local}+\Sigma_{i\hat a}^{\rm loop}(s).
\end{aligned}
$$

Thus $\Sigma_{i\hat a}^{\rm loop}$ directly corrects the mixed source column, and $\Sigma_{ij}^{\rm loop}w_{0j}^{(a)}$ dresses the existing LO mixing through the meson propagator. If the pole position itself shifts by $\delta s_a$, the evaluation adds

$$
\delta s_a\,\partial_s w_0^{(a)}(s_a^{(0)})
=-\delta s_a\,[A_0^{-1}(s_a^{(0)})]_{ij}
\left(K_{j k}^{(0)}w_{0k}^{(a)}(s_a^{(0)})+K_{j\hat a}^{(0)}\right)
$$

at this order. Whether that term is retained depends on the ALP power counting and whether the pole mass is an input or a prediction.

The [local invariant source](lsz-invariant-split.md) $\mathcal C_i^{\rm pole}=M^2_{ij}(K_{\rm meson}^{-1})_{jk}K_{k\hat a}-M^2_{i\hat a}$ organizes the local, momentum-independent quadratic problem. The one-loop extension is most directly and safely expressed by the full block equation above: a momentum-dependent $\Sigma^{\rm loop}(s)$ cannot generally be folded into that same constant source. Under a regular change of meson coordinates the full inverse propagator transforms as a quadratic form, while its pole and physical LSZ amplitudes remain invariant when the direct vertices are transformed with it.

## Pole residue and amplitudes

For a symmetric inverse propagator and a pole eigenvector $u^{(a)}$, its pole normalization is governed by

$$
Z_a^{-1}=u^{(a)T}\left.\frac{\partial\Gamma(s)}{\partial s}\right|_{s=s_a}u^{(a)},
$$

under the chosen eigenvector normalization. Therefore the loop effects enter an amplitude both through the changed meson components $w_i^{(a)}$ and through the pole residue $Z_a^{1/2}$. The meson-pole contribution combines with direct ALP vertices; neither term separately is invariant under a change of meson field coordinates [@Alda:2026nrd].

## Reproducible block expansion

The script [derive_loop_corrected_alp_lsz.py](../../scripts/derive_loop_corrected_alp_lsz.py) uses noncommuting SymPy matrix symbols to expand $-(A_0+\delta A)^{-1}(b_0+\delta b)$ through first order in the NLO corrections, and prints the intermediate and final expressions in LaTeX. Run it with `.venv/bin/python scripts/derive_loop_corrected_alp_lsz.py`. The derivation assumes $A_0$ is invertible at the ALP pole. It does not evaluate the tadpole or bubble integrals; those require the chosen ALP sources and the explicit leading vertices.
