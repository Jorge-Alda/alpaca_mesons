---
title: Up-only sources - leading large-epsilon size estimate of the one-loop sigma to pi+ pi- amplitude
status: derived
created: 2026-11-27
updated: 2026-11-27
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: analytical-derivation
---

# Up-only sources: $M_1^{\rm est}/M_0$ at large $\epsilon$

OpenSpec change `estimate-sigma-pipi-one-loop-murayama` (task 4.6). This record derives the leading large-$\epsilon$ behaviour of the size estimate of the one-loop $\sigma\to\pi^+\pi^-$ amplitude for the up-only sources:
$$
\frac{M_1^{\rm est}}{M_0}\to-\frac{21B_0m_u}{64\pi^2f_\pi^2}\,\epsilon^4,\qquad R=\frac{|M_1^{\rm est}|}{|M_0|}\to\frac{21B_0m_u}{64\pi^2f_\pi^2}\,\epsilon^4 .
$$
$R$ is a size ratio, not a computed correction. No loop integral, loop function or loop amplitude is evaluated. Units: GeV for $B_0$, $m_u$, $f_\pi$ and the source letters; $\epsilon=f_\pi/f_a$ is dimensionless.

## Inputs

- $M_1^{\rm est}$, $X_\Gamma^{\max}$ and $\delta Z^{\max}$ (size rule, per-term loop scales, termwise routing maximum): [estimate record](../one-loop-estimate/murayama-sigma-pipi-loop-estimate.md).
- Large-$\epsilon$ assumptions (fixed sources, fixed $\lambda>0$, ALP mass entry $\epsilon^2A$ with $A=-2B_0\langle\widehat S_{0,2}\rangle$, $m_\sigma^2=2\lambda f_\pi^2/\epsilon^2+\epsilon^2\kappa_\sigma$ with $\kappa_\sigma=B_0(\langle\widehat S_{1,0}\rangle-2\langle\widehat S_{2,0}\rangle)$, $m_P$ fixed, routing choice at large $\epsilon$ and its regions of $\kappa_\sigma/A$): [large-$\epsilon$ input record](murayama-sigma-pipi-loop-large-epsilon-inputs.md).
- Leading $\epsilon^4$ terms of $M_1^{\rm est}/M_0$ for generic fixed sources, graph by graph: [large-$\epsilon$ graph record](murayama-sigma-pipi-loop-large-epsilon-graphs.md).
- Exact tree $M_0=-\frac{B_0\epsilon}{2f_\pi}\langle\widehat S_{1,0}\{\lambda_{\pi^+},\lambda_{\pi^-}\}\rangle$ and the up-only relations of the pion matrices: [tree and legs record](../one-loop-estimate/murayama-sigma-pipi-loop-tree-and-legs.md).
- Up-only pole data ($A=B_0m_u$, $m_a^2=\epsilon^2B_0m_u+2B_0m_u+\dots$): [large-$\epsilon$ record](murayama-sigma-pipi-loop-large-epsilon-scaling.md), Step 4.

## Assumptions

1. Up-only sources (model assumption): $\widehat S_{r,s}+i\widehat P_{r,s}=m_u(-i)^sP_u/s!$ for $r\le1$, $(r,s)\neq(0,0)$, and zero for $r\ge2$. $P_u$ is a symbolic projector: $P_u^2=P_u$, $\widehat S_{0,0}P_u=P_u\widehat S_{0,0}=m_uP_u$, $\langle P_u\rangle=1$. No matrix components are used.
2. Pion matrices ($\pi^+=u\bar d$): $\lambda_{\pi^+}\lambda_{\pi^-}=2P_u$, $\lambda_{\pi^-}\lambda_{\pi^+}=2P_d$, $P_d^2=P_d$, $P_uP_d=P_dP_u=0$, $\langle P_d\rangle=1$, $\lambda_{\pi^+}=P_u\lambda_{\pi^+}P_d$, $\lambda_{\pi^-}=P_d\lambda_{\pi^-}P_u$, $\lambda_{\pi^\pm}^2=0$. The last three relations are needed only for the lower powers of graphs 3, 14, 15, 16, 22 and π1–π3; the $\epsilon^4$ terms do not use them.
3. All assumptions of the large-$\epsilon$ input record and the limits of the [limits record](../one-loop-estimate/murayama-sigma-pipi-loop-limits.md) apply. In particular, the on-shell cancellations between the kinetic $\sigma aa$ vertex and the $\widehat S_{1,2}$ source vertex, which the size rule cannot resolve, concern exactly graphs 2, 13 and σ4, which give the whole result below, and the up-only model has $m_\sigma\simeq m_a$.

## Derivation

**Step 1. Up-only letters.** Assumption 1 gives
$$
\widehat S_{1,0}=m_uP_u,\quad \widehat S_{0,2}=\widehat S_{1,2}=-\frac{m_u}{2}P_u,\quad \widehat P_{0,1}=\widehat P_{1,1}=-m_uP_u,
$$
$$
\widehat S_{0,1}=\widehat S_{1,1}=\widehat S_{2,0}=\widehat S_{2,1}=\widehat S_{3,0}=\widehat P_{2,0}=0,
$$
and $\widehat S_{0,0}$ stays a letter with $\widehat S_{0,0}P_u=P_u\widehat S_{0,0}=m_uP_u$. These are all the source letters of the kept graphs.

**Step 2. Derived inputs.** Inserting Step 1 and the relations of assumptions 1 and 2:
$$
A=-2B_0\langle\widehat S_{0,2}\rangle=B_0m_u,\qquad \kappa_\sigma=B_0\big(\langle\widehat S_{1,0}\rangle-2\langle\widehat S_{2,0}\rangle\big)=B_0m_u,\qquad M_0=-\frac{B_0m_u\epsilon}{f_\pi},
$$
$$
\langle\widehat S_{1,0}\{\lambda_{\pi^+},\lambda_{\pi^-}\}\rangle=2m_u,\qquad \frac{\langle\widehat S_{0,2}\{\lambda_{\pi^+},\lambda_{\pi^-}\}\rangle}{\langle\widehat S_{1,0}\{\lambda_{\pi^+},\lambda_{\pi^-}\}\rangle}=\frac{\langle\widehat S_{1,2}\{\lambda_{\pi^+},\lambda_{\pi^-}\}\rangle}{\langle\widehat S_{1,0}\{\lambda_{\pi^+},\lambda_{\pi^-}\}\rangle}=-\frac12,\qquad \langle\widehat S_{1,2}\rangle=-\frac{m_u}{2},
$$
and the ratios with $\widehat S_{1,1}$, $\widehat S_{2,0}$, $\widehat S_{3,0}$, as well as $\langle\widehat S_{2,1}\rangle$ and $\langle\widehat S_{3,0}\rangle$, vanish.

**Step 3. Region.** $A=\kappa_\sigma$, so the up-only model sits on the boundary $\kappa_\sigma=A$ of the regions of the routing choice. The $\epsilon^2$ terms of $m_a^2-m_\sigma^2$ cancel, and the comparison is decided at order $\epsilon^0$:
$$
m_a^2-m_\sigma^2=2B_0m_u+\Big(\frac{2m_1^2}{3}-2\lambda f_\pi^2\Big)\epsilon^{-2}+O(\epsilon^{-4})\to2B_0m_u>0,
$$
so $\max(m_\sigma,m_a)=m_a$. The script evaluates the up-only limit on both sides of the boundary, $\kappa_\sigma<A$ and $A<\kappa_\sigma<2A$, and finds the same result.

**Step 4. Generic $\epsilon^4$ terms and their up-only values.** For generic fixed sources the graphs that reach $\epsilon^4$ in $M_1^{\rm est}/M_0$ are 1, 2, 10, 11, 13, σ1, σ2 and σ4 ([large-$\epsilon$ graph record](murayama-sigma-pipi-loop-large-epsilon-graphs.md)). The contribution of a 1PI graph is $X_\Gamma^{\max}/M_0$ and that of a $\sigma$ self-energy is $\frac12\delta Z^{\max}$. With $\{\cdot\}=\{\lambda_{\pi^+},\lambda_{\pi^-}\}$, the coefficients of $\epsilon^4$ in region $\kappa_\sigma<A$ are:

| graph | generic coefficient of $\epsilon^4$ | up-only value |
|---|---|---|
| 1 | $\dfrac{3\kappa_\sigma\langle\widehat S_{3,0}\{\cdot\}\rangle}{16\pi^2f_\pi^2\langle\widehat S_{1,0}\{\cdot\}\rangle}$ | $0$ |
| 2 | $\dfrac{A\langle\widehat S_{1,2}\{\cdot\}\rangle}{16\pi^2f_\pi^2\langle\widehat S_{1,0}\{\cdot\}\rangle}$ | $-\dfrac{B_0m_u}{32\pi^2f_\pi^2}$ |
| 10 | $\dfrac{3B_0\langle\widehat S_{3,0}\rangle\langle\widehat S_{2,0}\{\cdot\}\rangle}{8\pi^2f_\pi^2\langle\widehat S_{1,0}\{\cdot\}\rangle}$ | $0$ |
| 11 | $\dfrac{B_0\langle\widehat S_{2,1}\rangle\langle\widehat S_{1,1}\{\cdot\}\rangle}{8\pi^2f_\pi^2\langle\widehat S_{1,0}\{\cdot\}\rangle}$ | $0$ |
| 13 | $\dfrac{\big(B_0\langle\widehat S_{1,2}\rangle+3A\big)\langle\widehat S_{0,2}\{\cdot\}\rangle}{8\pi^2f_\pi^2\langle\widehat S_{1,0}\{\cdot\}\rangle}$ | $-\dfrac{5B_0m_u}{32\pi^2f_\pi^2}$ |
| σ1 | $-\dfrac{9B_0^2\langle\widehat S_{3,0}\rangle^2}{16\pi^2f_\pi^2\kappa_\sigma}$ | $0$ |
| σ2 | $-\dfrac{B_0^2\langle\widehat S_{2,1}\rangle^2}{8\pi^2f_\pi^2A}$ ($A\to\kappa_\sigma$ for $A<\kappa_\sigma<2A$) | $0$ |
| σ4 | $-\dfrac{A}{4\pi^2f_\pi^2}-\dfrac{B_0\langle\widehat S_{1,2}\rangle}{4\pi^2f_\pi^2}-\dfrac{B_0^2\langle\widehat S_{1,2}\rangle^2}{16\pi^2f_\pi^2A}$ | $-\dfrac{9B_0m_u}{64\pi^2f_\pi^2}$ |

The up-only values are obtained by inserting the values of Step 2 into the generic coefficients. The coefficients of graphs 1, 10, 11, σ1 and σ2 vanish because each contains $\widehat S_{3,0}$, $\widehat S_{2,0}$, $\widehat S_{1,1}$ or $\widehat S_{2,1}$; graphs 1, 10, 11 and σ2 vanish exactly, and σ1 starts at $\epsilon^{-4}$. In region $A<\kappa_\sigma<2A$ only the generic coefficient of σ2 changes, and the up-only values are the same.

**Step 5. The other kept graphs.** For the up-only sources every other kept graph either vanishes exactly or stays below $\epsilon^4$:

- exactly zero (each product of traces contains a vanishing letter of Step 1): 4–9, 12, 17–21, 23–26, σ3, and π4, π5 on each charged leg;
- $\epsilon^2$: 3, 14, 15, σ5, σ6;
- $\epsilon^0$: 22, and π1, π2, π3 on each charged leg;
- $\epsilon^{-2}$: 16.

The leading powers are the same in both regions.

**Step 6. Sum.** The three nonzero coefficients of Step 4 ($-\frac{1}{32}$, $-\frac{5}{32}$ and $-\frac{9}{64}$, in units of $B_0m_u/(\pi^2f_\pi^2)$) add up to
$$
\frac{M_1^{\rm est}}{M_0}\to-\frac{21B_0m_u}{64\pi^2f_\pi^2}\,\epsilon^4,\qquad R\to\frac{21B_0m_u}{64\pi^2f_\pi^2}\,\epsilon^4\quad(B_0m_u>0),
$$
on both sides of $\kappa_\sigma=A$.

## Verification

The script prints every item of Steps 1–6 as exact sympy expressions; the values above are copied from its output. It checks that the sum of Step 6 equals the up-only limit of the full large-$\epsilon$ total of $M_1^{\rm est}/M_0$ (sources with $r\le1$) in both regions, and the tests check each coefficient of Steps 2–6 exactly. No numerical matrix is used.

## Reproduction

```sh
PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_loop_up_only_estimate.py
PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_loop_up_only_estimate
```

[Script](../../../../../../scripts/derive_sigma_pipi_loop_up_only_estimate.py), [tests](../../../../../../scripts/test_derive_sigma_pipi_loop_up_only_estimate.py), built on [`derive_sigma_pipi_loop_large_epsilon.py`](../../../../../../scripts/derive_sigma_pipi_loop_large_epsilon.py).
