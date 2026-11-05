---
title: Rpi and RR' cuts and divergence structure of the one-loop RChT vector form factor
status: literature
created: 2026-11-05
updated: 2026-11-05
bibtex_keys: [Pich:2010sm]
provenance: literature
---

# Divergence structure and thresholds of resonance loops

Source: `memory/sources/Pich:2010sm/L9.tex`. [@Pich:2010sm]

## Setup (assumptions of the source)

- $U(3)$ multiplets $V,A,S,P$, chiral limit throughout, only operators with $\mathcal{O}(p^2)$ chiral tensors (lines 335-352). Matching $U(3)\to SU(3)$ gives identical $L_9$, $C_{88}-C_{90}$ at LO and NLO in $1/N_C$ (lines 943-946).
- Absorptive cuts kept: $\pi\pi$ and $R\pi$. Two-resonance channels $RR'$ have thresholds at $(M_R+M_{R'})^2\gtrsim2$ GeV$^2$ and "their impact on the LEC determination is expected to be negligible" (lines 655-658; line 407 states the same for higher thresholds). The condition for the simplified dispersive formula is a threshold above $M_V$ (lines 1225-1230).

## Divergence structure (Appendix "Description in terms of Feynman diagrams", lines 1627-1705)

One-loop part: $\sum\mathcal{F}(s)|_{m_1,m_2}=\sum\mathcal{F}^{1\ell}+\frac{2\hat\delta_2}{F^2}s+\hat\delta_0\frac{s}{M_V^2-s}+\hat\delta_{-2}\frac{s}{(M_V^2-s)^2}$ (Eq. at lines 1646-1656). The $\hat\delta_i$ carry the ultraviolet divergences; $\widetilde L_9$ ($\mathcal{O}(p^4)$) renormalizes $\hat\delta_2$, the product $F_VG_V$ ($\mathcal{O}(p^2)$ resonance couplings) cancels $\hat\delta_0$, and $\delta M_V^2$ cancels $\hat\delta_{-2}$ (on-shell scheme) (lines 1665-1673). Short-distance matching: $F_VG_V/F^2+\hat\delta_0=1+\delta_0$, $\widetilde L_9+\hat\delta_2=0$ (lines 1685-1688). The source gives these only for the form factor; it does not state the renormalization of $F$ or $B_0$, and does not treat the pion self-energy.

## Verdict

- O(p^2) affected: **yes, for the resonance couplings**: the $\pi\pi$ and $R\pi$ loops shift the product $F_VG_V$ at NLO, $F_VG_V/F^2=1+\delta_{\rm NLO}$ (line 618), an $\mathcal{O}(p^2)$-order coupling. Nothing about $F$ or $B_0$ is stated here (see Rosell:2004mn and Rosell:2005ai extracts for $F$).
- Conditions: this needs the short-distance matching (form factor vanishing at $s\to\infty$), and is for $V\pi$/$A\pi$ cuts in the vector channel; $\delta_{\rm NLO}$ is dominated by the $\pi\pi$ and $A\pi$ cuts (line 1066).
- Thresholds: $R\pi$ thresholds at $(M_R+m_\pi)^2\approx M_R^2$ (chiral limit); the $RR'$ ones are $\gtrsim2$ GeV$^2$. The source gives no pseudoscalar-channel thresholds below 1 GeV.
