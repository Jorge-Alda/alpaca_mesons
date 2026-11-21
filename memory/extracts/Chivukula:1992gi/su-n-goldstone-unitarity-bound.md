---
title: SU(N) Goldstone scattering at leading order, largest S-wave channel and unitarity bound (Chivukula:1992gi)
status: literature
created: 2026-11-21
updated: 2026-11-21
bibtex_keys: [Chivukula:1992gi, Lee:1977eg, Manohar:1983md]
provenance: literature
---

# SU(N) Goldstone scattering: largest S-wave channel and unitarity bound

Source: Chivukula:1992gi, arXiv hep-ph/9206222 TeX source read directly. Equation numbers below follow the default section-wise numbering of the harvmac macros used by the source (Section 1: (1.1)-(1.3); Section 2: (2.1)-(2.5)); they are inferred from the TeX order, not read from a rendered PDF.

## Setup

Symmetry breaking $SU(N)_L\times SU(N)_R\to SU(N)_V$, $\Sigma=\exp(2i\pi^aT^a/f)$ (Eq. (1.1)), $\mathrm{tr}\,T^aT^b=\delta^{ab}/2$,
$$\mathcal{L}_2=\frac{f^2}{4}\mathrm{tr}\,\partial^\mu\Sigma^\dagger\partial_\mu\Sigma+\frac{m^2f^2}{4}\mathrm{tr}(\Sigma+\Sigma^\dagger)\quad\text{(Eqs. (1.2), (1.3))},$$
with a common mass $m$ for all Goldstone bosons.

## Amplitude and channels

General amplitude (Eq. (2.2)), with $d^{abc}=2\,\mathrm{tr}\{T^a,T^b\}T^c$:
$$a^{ab;cd}=\delta^{ab}\delta^{cd}A(s,t,u)+\delta^{ac}\delta^{bd}A(t,s,u)+\delta^{ad}\delta^{bc}A(u,t,s)+d^{abe}d^{cde}B(s,t,u)+d^{ace}d^{bde}B(t,s,u)+d^{ade}d^{bce}B(u,t,s).$$
Leading order (text after Eq. (2.3)): $A=\frac{2}{N}\frac{s-m^2}{f^2}$, $B=\frac{s-m^2}{f^2}$.

Channels: adjoint $\otimes$ adjoint $=\Delta\oplus F\oplus D\oplus Y\oplus X\oplus T\oplus\bar T$ ($\Delta$ singlet; for $N=2$ only $\Delta,F,X$, i.e. isospin $0,1,2$). Singlet channel (Eq. (2.3)):
$$a_\Delta=(N^2-1)A(s,t,u)+A(t,s,u)+A(u,t,s)+\frac{N^2-4}{N}\big[B(t,s,u)+B(u,t,s)\big].$$

## Partial-wave normalization (Eq. (2.4))

$$a_{I\ell}(s)=\frac{1}{64\pi}\int_{-1}^{1}a_I(s,\cos\theta)P_\ell(\cos\theta)\,d\cos\theta.$$

## Leading-order partial waves (Eq. (2.5))

$$\begin{aligned}
a_{\Delta0}&=\frac{Ns}{32\pi f^2}-\frac{m^2}{16\pi Nf^2},&
a_{F1}&=\frac{Ns}{192\pi f^2}-\frac{Nm^2}{48\pi f^2},\\
a_{D0}&=\frac{Ns}{64\pi f^2}-\frac{m^2}{8\pi N},&
a_{Y0}&=\frac{s}{32\pi f^2}-\frac{m^2}{16\pi f^2},\\
a_{X0}&=-\frac{s}{32\pi f^2}+\frac{m^2}{16\pi f^2}.
\end{aligned}$$
All other partial waves vanish at this order, including $a_{T1}$ and $a_{\bar T1}$. The source prints $a_{D0}$ with $-m^2/(8\pi N)$ and no $f^2$; this is as printed in the TeX (probable typo in the source, not corrected here). Only $a_{\Delta0}$, $a_{F1}$, $a_{D0}$ grow with $N$.

## Largest S-wave eigenvalue and unitarity bound

- The largest S-wave channel is the $SU(N)_V$ singlet, $a_{\Delta0}$ of Eq. (2.5). For small $m^2$: $a_{\Delta0}\simeq Ns/(32\pi f^2)$.
- Source statement (Section 2, after Eq. (2.5)): $a_{\Delta0}$ at tree level is real and, for small $m^2$, "would exceed 1 when $\sqrt{s}>4\pi f/\sqrt{N}$"; a physical amplitude must lie on or inside the Argand circle, so loop and higher-order terms must then be as large as the two-derivative term. Hence $\Lambda_\chi\lesssim 4\pi f/\sqrt{N}$, independent of any $1/N$ expansion.
- Normalization used: the criterion in the source is that $a_{I\ell}$ not exceed 1 ("unitarity limit", Argand circle); the source does not write $|\mathrm{Re}\,a_0|\le1/2$.
- Arithmetic note (derived by me, not in the source): setting $Ns/(32\pi f^2)=1$ gives $\sqrt{s}=\sqrt{32\pi}\,f/\sqrt{N}\approx 10.0\,f/\sqrt{N}$, slightly below $4\pi f/\sqrt{N}\approx12.6\,f/\sqrt{N}$, where $a_{\Delta0}=\pi/2$. So $4\pi f/\sqrt{N}$ is an order-of-magnitude scale, not the exact crossing of $|a|=1$. With a $|\mathrm{Re}\,a_0|\le1/2$ criterion the crossing would be at $\sqrt{s}=\sqrt{16\pi}f/\sqrt{N}\approx7.1\,f/\sqrt{N}$ (also derived).
- One-loop corrections are of relative size $Ns/(16\pi^2f^2)$ and the leading-$N$ $L$-loop corrections scale as $\big(\sqrt{N}p/(4\pi f)\big)^{2L}$ (Eq. (3.5)-type relation, source label `typical`), supporting the same scale. For $N=2$ this reproduces the naive-dimensional-analysis scale $4\pi f$ cited from Manohar:1983md; the argument is stated to generalize that of Lee:1977eg for longitudinal gauge bosons (the source cites the Phys. Rev. Lett. 38 (1977) 883 paper, not the Phys. Rev. D 16 paper, for that reference).
