---
title: Basics of Resonance Chiral Theory with the light vector nonet
status: mixed
created: 2026-09-23
updated: 2026-09-23
bibtex_keys: [Ecker:1988te, Ecker:1989yg]
provenance: literature
---

# Basics of Resonance Chiral Theory

## Scope and terminology

Resonance Chiral Theory (usually abbreviated **RChT**, although “RChPT” is also used informally) extends chiral perturbation theory by retaining explicit meson resonances while preserving the nonlinear realization of $SU(3)_L\times SU(3)_R$. It is intended for momenta near resonance poles, where Goldstone-only ChPT is insufficient. The standard working approximations are leading order in $1/N_c$ (tree-level, narrow resonances) and truncation to the lightest resonance multiplet in each channel. Neither approximation follows from chiral symmetry. [@Ecker:1988te]

The pseudo-Goldstone field is

$$
U=u^2=\exp\!\left(\frac{i\sqrt2\,\Phi}{F}\right),
$$

and transforms as $U\to g_RUg_L^\dagger$. Its square root transforms nonlinearly,

$$
u\longrightarrow g_Ruh^\dagger=hug_L^\dagger,
$$

where $h(g_L,g_R,U)\in SU(3)_V$ is the CCWZ compensator. Resonance multiplets transform homogeneously, $R\to hRh^\dagger$. [@Ecker:1988te, Sec. 2]

## Vector nonet

The light $J^{PC}=1^{--}$ states are collected into a $U(3)$ nonet represented here by an antisymmetric Lorentz tensor,

$$
V_{\mu\nu}=\frac{1}{\sqrt2}\sum_{a=1}^{8}\lambda_a V^a_{\mu\nu}
+\frac{1}{\sqrt3}\,\mathbf 1\,V^0_{\mu\nu},
\qquad
V_{\mu\nu}\longrightarrow hV_{\mu\nu}h^\dagger,
$$

with $\langle\lambda_a\lambda_b\rangle=2\delta_{ab}$. The singlet is included by extending the octet construction of the foundational $SU(3)$ theory to $U(3)$, as motivated at leading order in $1/N_c$. This nonet assignment and single-multiplet truncation are modeling assumptions. [@Ecker:1988te]

Before specifying neutral-vector mixing,

$$
V_{\mu\nu}=
\begin{pmatrix}
\dfrac{\rho^0}{\sqrt2}+\dfrac{\omega_8}{\sqrt6}+\dfrac{\omega_0}{\sqrt3}
&\rho^+&K^{*+}\\
\rho^-&-\dfrac{\rho^0}{\sqrt2}+\dfrac{\omega_8}{\sqrt6}+\dfrac{\omega_0}{\sqrt3}
&K^{*0}\\
K^{*-}&\bar K^{*0}&-\dfrac{2\omega_8}{\sqrt6}+\dfrac{\omega_0}{\sqrt3}
\end{pmatrix}_{\!\mu\nu}.
$$

In the ideal-mixing convention $\omega=(u\bar u+d\bar d)/\sqrt2$ and $\phi=s\bar s$, this becomes

$$
V_{\mu\nu}=
\begin{pmatrix}
\dfrac{\rho^0+\omega}{\sqrt2}&\rho^+&K^{*+}\\
\rho^-&\dfrac{-\rho^0+\omega}{\sqrt2}&K^{*0}\\
K^{*-}&\bar K^{*0}&\phi
\end{pmatrix}_{\!\mu\nu}.
$$

Ideal $\omega$--$\phi$ mixing is a conventionally useful approximation, not part of the definition of RChT; realistic applications may use a fitted nonideal mixing angle and $SU(3)$-breaking masses. Overall signs of $\phi$ and charged fields are phase conventions if all interactions are changed consistently. [@Ecker:1988te]

## Chiral building blocks

For external left- and right-handed sources $\ell_\mu,r_\mu$ and $\chi=2B_0(s+ip)$, define

$$
\begin{aligned}
u_\mu&=i\left[u^\dagger(\partial_\mu-ir_\mu)u-u(\partial_\mu-i\ell_\mu)u^\dagger\right],\\
\Gamma_\mu&=\frac12\left[u^\dagger(\partial_\mu-ir_\mu)u+u(\partial_\mu-i\ell_\mu)u^\dagger\right],\\
\nabla_\mu V_{\alpha\beta}&=\partial_\mu V_{\alpha\beta}+[\Gamma_\mu,V_{\alpha\beta}],\\
f_\pm^{\mu\nu}&=uF_L^{\mu\nu}u^\dagger\pm u^\dagger F_R^{\mu\nu}u,\\
\chi_\pm&=u^\dagger\chi u^\dagger\pm u\chi^\dagger u.
\end{aligned}
$$

Except for the connection $\Gamma_\mu$, these blocks transform homogeneously as $X\to hXh^\dagger$, so flavor traces of their products are chiral invariants. [@Ecker:1988te, Sec. 2]

## Lowest-order Lagrangian

With $\langle\cdots\rangle$ denoting the flavor trace, the lowest-order Goldstone Lagrangian is

$$
\mathcal L_2^\chi=\frac{F^2}{4}\left\langle u_\mu u^\mu+\chi_+\right\rangle.
$$

For one antisymmetric-tensor vector nonet with a common chiral-limit mass $M_V$, the free term is

$$
\mathcal L_{\mathrm{kin}}^V
=-\frac12\left\langle
\nabla^\lambda V_{\lambda\mu}\nabla_\nu V^{\nu\mu}
-\frac{M_V^2}{2}V_{\mu\nu}V^{\mu\nu}
\right\rangle.
$$

The leading even-intrinsic-parity interactions linear in $V_{\mu\nu}$ are

$$
\mathcal L_2^V
=\frac{F_V}{2\sqrt2}\left\langle V_{\mu\nu}f_+^{\mu\nu}\right\rangle
+\frac{iG_V}{\sqrt2}\left\langle V_{\mu\nu}u^\mu u^\nu\right\rangle.
$$

Thus the minimal leading action for the Goldstones and light vectors is

$$
\boxed{\mathcal L_{\mathrm{LO}}^{\mathrm{even}}
=\mathcal L_2^\chi+\mathcal L_{\mathrm{kin}}^V+\mathcal L_2^V.}
$$

$F_V$ couples vectors to external vector sources (and therefore produces neutral-vector--photon transitions when electromagnetism is switched on), while $G_V$ generates the leading $VPP$ vertex. The subscript “2” labels two-derivative/$O(p^2)$ chiral tensors; RChT does not possess the same complete homogeneous momentum counting as Goldstone-only ChPT because $M_V$ remains nonzero in the chiral limit. [@Ecker:1988te, Sec. 2]

## What is not contained in this minimal expression

The boxed Lagrangian is only the lowest-order **even-parity, one-vector** sector. A phenomenological theory generally also needs explicit $SU(3)$-breaking operators, nonideal neutral mixing, resonance widths or loops, operators with two or more resonance fields, and the odd-intrinsic-parity sector (Wess--Zumino--Witten plus independent resonance operators). Their coefficients are not fixed by chiral symmetry alone. [@Ecker:1988te]

Antisymmetric tensors are a convenient representation rather than a physical requirement. Proca and other spin-one formulations are equivalent only after all local operators induced by the field redefinition are retained; truncated Lagrangians can otherwise differ. Separately, QCD short-distance constraints are needed to obtain the representation-independent $O(p^4)$ resonance predictions analyzed in the foundational comparison. [@Ecker:1989yg]

## Evidence status

- **Literature established:** nonlinear transformation laws, covariant building blocks, antisymmetric-tensor kinetic term, and the two leading even-parity vector couplings. [@Ecker:1988te]
- **Explicit assumptions:** leading $1/N_c$, one light vector nonet, a common chiral-limit mass, universal octet/singlet couplings $F_V$ and $G_V$ with no independent double-trace terms, and (only in the simplified physical-state matrix) ideal $\omega$--$\phi$ mixing.
- **Convention dependent:** field phases, generator normalization, and neutral-vector mixing parametrization. Physical amplitudes are unchanged when conventions are used consistently.
