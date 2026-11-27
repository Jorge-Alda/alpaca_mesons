---
title: Up-only parameters of the radial-mode loop-to-tree crossing
status: derived
created: 2026-11-23
updated: 2026-11-25
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: analytical-derivation
---

# Up-only parameters of the radial-mode loop-to-tree crossing

This record gives the up-only parameters of the crossing criterion of [the general record](../../radial-mode/finite-mixing/loop-estimates/power-counting/murayama-radial-crossing-implicit.md). The estimate and its numbers are in [the estimate record](../../radial-mode/couplings-and-decays/finite-mixing/murayama-radial-up-only-estimate-result.md).

## Parameters

With $f_a=\sqrt2\,m_u/\kappa_u$,
$$
\epsilon=\frac{f_\pi}{f_a}=\frac{f_\pi\,\kappa_u}{\sqrt2\,m_u}.
$$
The radial mass is
$$
m_\sigma^2=\frac{B_0 f_\pi^2\kappa_u^2}{2m_u}+\frac{4\lambda m_u^2}{\kappa_u^2},
$$
and the ratio of the $\lambda$ term to the first term is
$$
\frac{4\lambda m_u^2/\kappa_u^2}{B_0 f_\pi^2\kappa_u^2/(2m_u)}=\frac{8\lambda m_u^3}{B_0 f_\pi^2\kappa_u^4}.
$$
The value $\lambda=10^{-3}$ is the default of the numerical evaluation.

## Sources

- $\widehat S_{00}=\operatorname{diag}(m_u,m_d,m_s)$.
- $\widehat S_{10}=m_u E_{11}$ and $\widehat P_{01}=-m_uE_{11}$, where $E_{11}=\operatorname{diag}(1,0,0)$.
- $\widehat S_{02}=\widehat S_{12}=-\tfrac{m_u}{2}E_{11}$.
- $\widehat S_{20}=\widehat S_{21}=\widehat P_{11}=\widehat P_{21}=\widehat S_{22}=0$.
- $\mathcal C=B_0\widehat P_{01}=-B_0m_uE_{11}$.

## ALP mass

The ALP-like pole $m_a^2$ is the root of $F(s)=0$ with the largest ALP component, where
$$
F(s)=s+2B_0\epsilon^2\langle\widehat S_{02}\rangle+2\epsilon^2\big\langle\mathcal C\,\mathbb M^{-1}_{s}(\mathcal C)\big\rangle,
$$
$$
\mathbb M_s(X)=B_0\{\widehat S_{00},X\}+\tfrac{m_1^2}{3}\operatorname{tr}(X)\,\mathbf 1-sX,\qquad N(s)=1+2\epsilon^2\operatorname{tr}\big[(\mathbb M^{-1}_{s}(\mathcal C))^2\big]=F'(s).
$$
The root is found numerically from $F(s)=0$; no closed form is used.

## Estimate

The one-loop estimate in these parameters is the one defined in [the estimate record](../../radial-mode/couplings-and-decays/finite-mixing/murayama-radial-up-only-estimate-result.md). All dependence on $\kappa_u$ enters through $\epsilon$, $f_a$, $m_\sigma^2$, and the pole data $m_a^2$, $m_{\pi^0}^2$, $N$ and $\mathbb M^{-1}_{s}(\mathcal C)$.
