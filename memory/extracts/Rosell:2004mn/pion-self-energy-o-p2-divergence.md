---
title: One-loop pion self-energy in RChT has an O(p^2) divergence from resonance loops
status: literature
created: 2026-11-05
updated: 2026-11-05
bibtex_keys: [Rosell:2004mn]
provenance: literature
---

# Pion self-energy: O(p^2) divergence from $V\pi$ and $S\pi$ loops

Source: `memory/sources/Rosell:2004mn/1LVFF.tex`, subsection "Pion self-energy" (lines 660-715), with setup at lines 300-306 and 650-656. [@Rosell:2004mn]

## Setup (assumptions of the source)

- Two flavours, $U(2)_L\otimes U(2)_R$, massless quarks (chiral limit), so Goldstone tadpoles vanish (lines 300-306, 654-656). The $U(1)_A$ effects are neglected.
- Single-resonance approximation, couplings linear in resonance fields: $G_V$ ($\rho\pi\pi$-type, $V\pi$ loop) and $c_d$ ($S\pi$ loop). The $F_V$ coupling does not enter the self-energy. The kinetic terms $\mathcal{L}_{2Z}(R)$ give resonance-tadpole topologies that vanish identically (lines 672-676). Couplings $\lambda^{RR}$ of two-resonance operators are not included.
- Scheme: $\overline{MS}-1$ with $\lambda_\infty$; on-shell condition for the pion self-energy (lines 640-653).

## Result

Divergences of $\mathcal{O}(p^2)$ are absorbed in the pion wave function, $\pi^{(B)}=(1+\delta Z_\pi)^{1/2}\pi^{(r)}$ (lines 677-680). In the on-shell scheme (lines 683-690):

$$
\delta Z_\pi = -\frac{2G_V^2}{F^2}\,\frac{3M_V^2}{16\pi^2F^2}\left\{\lambda_\infty+\ln\frac{M_V^2}{\mu^2}+\frac16\right\}
+\frac{4c_d^2}{F^2}\,\frac{M_S^2}{16\pi^2F^2}\left\{\lambda_\infty+\ln\frac{M_S^2}{\mu^2}-\frac12\right\}.
$$

Divergences of $\mathcal{O}(p^4)$ renormalize the coupling $\widetilde{\ell}_{12}$ of $\widetilde{\mathcal{L}}_{4\chi}$ (lines 692-700):
$\delta\widetilde\ell_{12}=-\frac{G_V^2+2c_d^2}{F^2}\frac{\lambda_\infty}{32\pi^2}$.

The renormalized self-energy (lines 702-712) is $-i\Sigma_\pi^r(p^2)=-i\frac{p^4}{16\pi^2F^2}\{64\pi^2\widetilde\ell_{12}^r(\mu)+\frac{2G_V^2}{F^2}[\ln\frac{M_V^2}{\mu^2}+\phi(p^2/M_V^2)]+\frac{4c_d^2}{F^2}[\ln\frac{M_S^2}{\mu^2}+\phi(p^2/M_S^2)]\}$: after on-shell subtraction of the $\mathcal{O}(p^2)$ piece the remainder starts at $\mathcal{O}(p^4)$.

## Related statements in the same paper

- $\rho$ self-energy: only an $\mathcal{O}(p^4)$ divergence, $\delta X_Z$; $\delta M_V^2=\delta Z_V=0$ (lines 735-756).
- $\langle v V\rangle$ vertex: only $\mathcal{O}(p^4)$ divergence; "no divergences of $\mathcal{O}(p^2)$, $\delta F_V=0$" (lines 800-808).
- $\langle V\pi\pi\rangle$ vertex: $\mathcal{O}(p^2)$ divergence renormalizing $G_V$ (lines 840-852): $\delta G_V=G_V[3M_V^2(2G_V^2/F^2-1/2)-4c_d^2M_S^2/F^2]\lambda_\infty/(16\pi^2F^2)$.

## Verdict

- O(p^2) affected: **yes**. The $V\pi$ and $S\pi$ loops give an $\mathcal{O}(p^2)$ ultraviolet divergence in the pion two-point function, of size $G_V^2M_V^2/(16\pi^2F^4)$ and $c_d^2M_S^2/(16\pi^2F^4)$ relative to the tree kinetic term, i.e. suppressed by $1/N_C$ since $G_V,c_d\sim\sqrt{N_C}$, $F\sim\sqrt{N_C}$. This confirms the estimate $\propto G_V^2M_V^2/(16\pi^2F^4)$.
- Here it is absorbed in $\delta Z_\pi$ (wave function). The source does not write it as a shift of $F$ or $B_0$; it works at $m_q=0$, so $B_0$ cannot be discussed. The equivalent chirally covariant statement (renormalization of the coefficient of $\langle u_\mu u^\mu\rangle$, i.e. $F^2$) is in the extract of Rosell:2005ai (see `../Rosell:2005ai/o-p2-counterterms-one-loop-generating-functional.md`). That equivalence is an interpretation, not a statement of this source.
- Condition: no short-distance constraint is needed for the divergence to appear; the finite parts and the $G_V$, $c_d$ values are where constraints would enter.
- The term $c_d^2M_S^2$ agrees with the $N F^{-2}c_d^2M_S^2$ term of $\gamma_1$ in Rosell:2005ai.
- Thresholds: this source has massless pions and only the $\pi\pi$ cut in the vector channel; it does not discuss pseudoscalar-channel thresholds of $R\pi$ loops. The pion self-energy cut $V\pi$ starts at $p^2=M_V^2$ in the chiral limit, by kinematics (not stated explicitly in the source).
