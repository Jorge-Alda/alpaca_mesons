---
title: Flavor algebra of the vector radiative channels
status: derived
created: 2026-10-04
updated: 2026-10-04
bibtex_keys: [Ruiz-Femenia:2003jdx, Aloni:2018vki]
provenance: analytical-derivation
---

# Flavor Algebra of the Vector Radiative Channels

## Scope

This record establishes which ALP spurion combinations the on-shell radiative
channels $V\to P\gamma$ can reach. It covers the flavor half only: the
Lorentz-level mapping from operator coefficients to the scalar $G$ of
$\mathcal M=i\epsilon_{\alpha\beta\rho\sigma}\epsilon_V^\alpha
\epsilon_\gamma^\beta q^\rho k^\sigma G$ is not derived here, so the coupling
combination entering each channel remains open and the published
$\omega\to\pi^0\gamma$ combination $c_2-c_1+c_5-2c_6$, $c_1+c_2+8c_3-c_5$ is not
reproduced.

Reproduce with
[`scripts/derive_radiative_flavor_algebra.py`](../../../scripts/derive_radiative_flavor_algebra.py),
which runs in SymPy and prints every table below.

## Method

With source matrices $T_A,T_B$, the kernel of
[the numerical-kernel record](odd-parity-rcht-numerical-kernel.md) gives

$$
D=\mathrm{tr}\bigl(\{T_A,T_B\}\,\mathbb M_a^{-1}(C)\bigr),\qquad
M=\mathrm{tr}\bigl(\{T_A,T_B\}\,\{c_q^{AP}+\mathbb M_a^{-1}(C),B_0m_q\}\bigr),
$$

$$
(c_1,\dots,c_7)\propto(4eD,\,4eD,\,-8eM,\,0,\,4eD,\,4eD,\,4eD),\qquad
(d_1,\dots,d_4)\propto(-2D,\,4M,\,-2D,\,-2D).
$$

Thus $c_1,c_2,c_5,c_6,c_7$ share $D$ with $d_1,d_3,d_4$; $c_3$ and $d_2$ share
$M$; $c_4$ vanishes identically. No channel can reach more than **two**
independent spurion combinations, whatever its meson content.

For diagonal sources $\{T_V,T_P\}=\mathrm{diag}(s_i)$ with
$s_i=2(T_V)_{ii}(T_P)_{ii}$, so

$$
D=\sum_i s_iC_{ii},\qquad M=2\sum_i s_i\bigl(C_{ii}+q_{ii}\bigr)m_i .
$$

This reduction is stated for diagonal sources only. The script evaluates
$D$ and $M$ as $\mathrm{tr}(\{T_A,T_B\}\mathbb M_a^{-1}(C))$ and
$\mathrm{tr}(\{T_A,T_B\}\{c_q^{AP}+\mathbb M_a^{-1}(C),B_0m_q\})$ with the
**full** anticommutator, so off-diagonal entries of $\mathbb M_a^{-1}(C)$ would
be picked up whenever the anticommutator has off-diagonal entries. Projecting
the anticommutator onto its diagonal before tracing would discard them and is
not done.

Channels with diagonal sources therefore differ **only** through the weight
pattern $s$. Because the flavor algebra involves only the anticommutator, it is
blind to crossing: $X\to Y\gamma$ and $Y\to X\gamma$ give identical $D$ and $M$.

## Conventions

Meson sources follow the
[on-shell radiative extract](../../extracts/Ruiz-Femenia:2003jdx/on-shell-radiative-amplitudes.md):
octet fields enter as $\Phi=\lambda^a\phi^a/\sqrt2$, the singlet is added as
$\omega_1/\sqrt3$ with $\omega_1=\mathrm{diag}(1,1,1)$. Every state is normalized
to unit $\mathrm{Tr}(V^\dagger V)$, which is the nonet requirement; the script asserts
this for each source. In the $(u,d,s)$ basis,

$$
\pi^0=\frac{\mathrm{diag}(1,-1,0)}{\sqrt2},\quad
\eta=\frac{\mathrm{diag}(1,1,-1)}{\sqrt3},\quad
\rho=\frac{\mathrm{diag}(1,-1,0)}{\sqrt2},\quad
\omega=\frac{\mathrm{diag}(1,1,0)}{\sqrt2},\quad
\phi=\mathrm{diag}(0,0,-1).
$$

No separate factor is attached to $\phi$: the ideal-mixing definition already
fixes it. With $\omega_8=\mathrm{diag}(1,1,-2)/\sqrt6$ the two combinations

$$
\omega=\sqrt{\tfrac23}\,\omega_1+\sqrt{\tfrac13}\,\omega_8=\tfrac{\mathrm{diag}(1,1,0)}{\sqrt2},
\qquad
\phi=-\sqrt{\tfrac13}\,\omega_1+\sqrt{\tfrac23}\,\omega_8=\mathrm{diag}(0,0,-1)
$$

arrive at unit $\mathrm{Tr}(V^2)$ on their own, which is the sense in which the
normalization is a nonet requirement rather than a free choice.

The pseudoscalar textures are the $\eta$–$\eta'$ shapes
$\mathrm{diag}(1,1,-1)$ and $\mathrm{diag}(1,1,2)$, rescaled to unit trace:

$$
\eta=\frac{\mathrm{diag}(1,1,-1)}{\sqrt3},\qquad
\eta'=\frac{\mathrm{diag}(1,1,2)}{\sqrt6}.
$$

They are exactly orthogonal, as $\mathrm{tr}(\eta\eta')\propto 1+1-2=0$
confirms. The texture corresponds to $|\sin\theta|=1/3$,
$\cos\theta=2\sqrt2/3$ in the $\eta_0$–$\eta_8$ system with
$\hat\eta_8=\mathrm{diag}(1,1,-2)/\sqrt6$ and $\hat\eta_0=I/\sqrt3$; the unit
trace factors above, not Aloni's, set the overall scale.

Kaon sources are **single off-diagonal elements**, one distinct field each. Only
the two electric-charge conserving channels are used, $K^{*0}\to K^0\gamma$ and
$K^{*+}\to K^+\gamma$:

$$
K^0=\begin{pmatrix}0&0&0\\0&0&1\\0&0&0\end{pmatrix},
\qquad
K^+=\begin{pmatrix}0&0&1\\0&0&0\\0&0&0\end{pmatrix},
$$

with $K^*$ sharing the flavor shape of $K$.

Charged pions follow the same single-element texture, with negative-charge
states given by the transpose, $\pi^+\propto E_{12}$ and $\pi^-\propto E_{21}$.

**No isospin limit is imposed.** $B_0m_q=\mathrm{diag}(m_u,m_d,m_s)$ is fully
general, and both $\mathbb M_a^{-1}(C)$ and $c_q^{AP}$ are general Hermitian
matrices.

## Results

| channel | weights $s$ | $D$ | $M$ |
|---|---|---|---|
| $\omega\to\pi^0\gamma$ | $(1,-1,0)$ | $C_{00}-C_{11}$ | $2(C_{00}m_u-C_{11}m_d)+2(m_uq_{00}-m_dq_{11})$ |
| $\rho\to\eta\gamma$ | $\tfrac{\sqrt6}{3}(1,-1,0)$ | $\tfrac{\sqrt6}{3}(C_{00}-C_{11})$ | $\tfrac{\sqrt6}{3}\times$ the above |
| $\eta'\to\rho^0\gamma$ | $\tfrac1{\sqrt3}(1,-1,0)$ | $\tfrac1{\sqrt3}(C_{00}-C_{11})$ | $\tfrac1{\sqrt3}\times$ the above |
| $\rho\to\pi^0\gamma$ | $(1,1,0)$ | $C_{00}+C_{11}$ | $2(C_{00}m_u+C_{11}m_d)+2(m_uq_{00}+m_dq_{11})$ |
| $\omega\to\eta\gamma$ | $\tfrac{\sqrt6}{3}(1,1,0)$ | $\tfrac{\sqrt6}{3}(C_{00}+C_{11})$ | $\tfrac{\sqrt6}{3}\times$ the above |
| $\eta'\to\omega\gamma$ | $\tfrac1{\sqrt3}(1,1,0)$ | $\tfrac1{\sqrt3}(C_{00}+C_{11})$ | $\tfrac1{\sqrt3}\times$ the above |
| $\phi\to\eta\gamma$ | $(0,0,\tfrac{2\sqrt3}{3})$ | $\tfrac{2\sqrt3}{3}C_{22}$ | $\tfrac{4\sqrt3}{3}(C_{22}+q_{22})m_s$ |
| $K^{*0}\to K^0\gamma$ | $(0,1,1)$ | $C_{11}+C_{22}$ | $2(C_{11}m_d+C_{22}m_s)+2(m_dq_{11}+m_sq_{22})$ |
| $K^{*+}\to K^+\gamma$ | $(1,0,1)$ | $C_{00}+C_{22}$ | $2(C_{00}m_u+C_{22}m_s)+2(m_uq_{00}+m_sq_{22})$ |
| $\rho^+\to\pi^+\gamma$ | $(1,1,0)$ | $C_{00}+C_{11}$ | $2(C_{00}m_u+C_{11}m_d)+2(m_uq_{00}+m_dq_{11})$ |

**Ten channels share five weight patterns:**

- $(1,-1,0)$: $\omega\to\pi^0\gamma$, $\rho\to\eta\gamma$, $\eta'\to\rho^0\gamma$
- $(1,1,0)$: $\rho\to\pi^0\gamma$, $\rho^+\to\pi^+\gamma$, $\omega\to\eta\gamma$, $\eta'\to\omega\gamma$
- $(0,0,1)$: $\phi\to\eta\gamma$, alone in reaching the $s$ component
- $(0,1,1)$: $K^{*0}\to K^0\gamma$
- $(1,0,1)$: $K^{*+}\to K^+\gamma$

Dropping the isospin limit makes $D$ and $M$ independent: $D$ depends on
$\mathbb M_a^{-1}(C)$ alone, while $M$ also carries $c_q^{AP}$ weighted by the
individual quark masses. Under $m_u=m_d=m_s$ they become proportional and the
redundancy is total.

## Reach of the off-diagonal spurion

Every channel above reaches only the **diagonal** of $\mathbb M_a^{-1}(C)$. For
the channels with diagonal meson sources the anticommutator is itself diagonal,
so $\mathrm{tr}$ of it against a Hermitian $\mathbb M_a^{-1}(C)$ keeps only the
diagonal. The kernels contract the **full** anticommutator, not a
diagonal projection, so this is a consequence of the sources rather than an
assumption built into the evaluation.

The kaon channels are the one place where the anticommutator is not diagonal in
the naive sense. The quantity to evaluate is $\{T^*,T\}$, and $T^*$ is the
transpose of $T$ by crossing symmetry. Substituting the same matrix on both
sides, $\{T,T\}=2T^2$, is not an equivalent shortcut: for a single-element matrix
it vanishes identically, which would make the channel vanish while its crossed
partners do not, contradicting crossing symmetry. With the transpose,

$$
\{E_{13}^T,E_{13}\}=\mathrm{diag}(1,0,1),\qquad
\{E_{23}^T,E_{23}\}=\mathrm{diag}(0,1,1),
$$

so both channels are non-empty and probe $C_{00}+C_{22}$ and $C_{11}+C_{22}$
respectively.

Note that $K^{*+}\to K^+\gamma$ and $K^{*0}\to K^0\gamma$ are $\Delta S=0$.

## Consequences for the fit

- Flavor redundancy is not observable redundancy. The width carries
  $(M^{\rm ext\,2}-m_P^2)^3$, so channels sharing a weight pattern can weight the
  coupling combination very differently. The crossed $\eta'$ channels are the
  clearest case: they duplicate a pattern but carry a large kinematic lever.
- Within this class the reachable spurion content is $C_{00},C_{11}$ in two sign
  combinations plus $C_{22}$, with $c_q^{AP}$ entering only through $M$. The two
  kaon channels add the $C_{00}+C_{22}$ and $C_{11}+C_{22}$ combinations, which
  no other channel in this class reaches, so they are the informative additions
  here. Ten channels carry five patterns, and $\rho^+\to\pi^+\gamma$ is redundant
  with the $(1,1,0)$ group, adding nothing new: its pattern involves $u,d$ only.
- Together with the observed saturation of $\omega\to\pi^0\gamma$ (normalized
  residual $\approx0.045$), the discriminating power of the remaining radiative
  channels looks limited. The Jacobian measurement of
  [the saturation probe](radiative-channel-saturation-probe.md) is still the
  deciding test and has not been performed.

## Accuracy caveats

- The $\sin\theta=1/3$ pseudoscalar texture is approximate; the final fit should
  substitute accurate nonet matrices for both pseudoscalar and vector states from
  one common parametrization, which may split the degeneracies above.
- Absolute normalizations of $D$ and $M$ are convention dependent; only the
  weight patterns are physical here.
- Results hold for the flavor level only. Equal weights mean equal $D$ and $M$,
  not equal observables, since channels differ in external masses and in which
  operators contribute.
- The kaon and charged rho results use $T^\star$ from crossing symmetry. Evaluating $\{T,T\}$
  instead would report both channels as vanishing, which is wrong; the two
  differ whenever $T^2\neq (T^\star)T$.