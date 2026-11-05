---
title: Two-particle cuts in the NLO Pi_VV - Pi_AA correlator: included, neglected, and why
status: literature
created: 2026-11-05
updated: 2026-11-05
bibtex_keys: [Pich:2008jm]
provenance: literature
---

# Two-particle cuts of $\Pi(t)=\Pi_{VV}(t)-\Pi_{AA}(t)$ at NLO in $1/N_C$

Source: `memory/sources/Pich:2008jm/src/L10.tex`. [@Pich:2008jm] Chiral limit, single-resonance approximation (SRA) with $U(3)$ multiplets $V,A,S,P$, only $\mathcal{O}(p^2)$ chiral tensors (lines 204-208).

## Cuts kept

- Spectral functions, Eq. `eq:ImSS` (lines 488-515): $\frac1\pi{\rm Im}\Pi_V=2F_V^{r2}\delta(t-M_V^{r2})+\frac{n_f}{2}\frac{\rho_V}{24\pi^2}$, $\frac1\pi{\rm Im}\Pi_A=2F^2\delta(t)+2F_A^{r2}\delta(t-M_A^{r2})+\frac{n_f}{2}\frac{\rho_A}{24\pi^2}$.
- Vector channel: $\pi\pi$ (threshold $0$), $A\pi$ (threshold $M_A^2$, form factors $\mathcal{F}^v_{A\pi},\mathcal{G}^v_{A\pi}$), $P\pi$ (threshold $M_P^2$, $\mathcal{F}^v_{P\pi}$) (lines 501-508).
- Axial channel: $V\pi$ (threshold $M_V^2$, $\mathcal{F}^a_{V\pi},\mathcal{G}^a_{V\pi}$), $S\pi$ (threshold $M_S^2$, $\mathcal{F}^a_{S\pi}$) (lines 510-514).
- The $V\pi$ threshold ($M_V\simeq0.77$ GeV, line 1186) lies below 1 GeV; the $A\pi$ and $S\pi$ thresholds are at $M_A^r=1.0\pm0.05$ GeV and $M_S=1.09\pm0.11$ GeV (lines 1245); $P\pi$ at $M_P=1.3\pm0.1$ GeV (line 1246). In the chiral limit thresholds are $M_R^2$.
- The text lists the order by threshold as $V\pi$, $S\pi$, $A\pi$, $P\pi$ (line 681).

## Cuts neglected and the stated reasons

- $RR'$ (two resonances): "kinematically suppressed" in the region of interest (lines 521-523, 683-686, 888-889). Justification in Appendix C.2 (lines 2927-2936): for a cut with threshold $\Lambda_{th}^2>M_R^2$ and $q^2\ll\Lambda_{th}^2$, $\Pi(q^2)|_\rho=\sum_n (q^2/\Lambda_{th}^2)^n\int_1^\infty dx\,\frac1\pi{\rm Im}\Pi(x\Lambda_{th}^2)/x^{n+1}$; for $Q^2\gg\Lambda_{th}^2$, $|\Pi|_\rho|\le Q^{-2}\int_{\Lambda_{th}^2}^\infty dt\,\frac1\pi|{\rm Im}\Pi|$. Both shrink as $\Lambda_{th}^2$ grows. This needs a spectral function that vanishes at large $t$ (line 2937). No numerical size for $RR'$ is given.
- Cuts with the extra multiplets $V',A'$ ("heavier thresholds", lines 888-889) and with extra $S',P'$ (lines 965-968, "high threshold"): neglected.
- "Higher two-meson cuts" are the dots in Eq. `eq:ImSS` (line 517). The explicit $\pi\pi$-type channels with kaons/$\eta$ are encoded only by the flavour factor $n_f/2$ multiplying two-flavour results (lines 687-689, 2432-2436); this is a counting statement of the source, not a separate kaon calculation.

## Short-distance reason for the structure

- Each cut must have a spectral function vanishing at large $t$ (lines 372-390, 446-447), which fixes the form factors (Eqs. `eq:VFF`, `eq:AFF`, lines 604-625). The correlator then has no polynomial or $t^m\ln(-t)$ pieces (lines 539-541, 2858-2862).

## SRA result for the cuts

- $\mathcal{F}^v_{P\pi}=0$ and $\mathcal{G}^v_{A\pi}=0$ (lines 611-613, 643), so $\Pi|_{P\pi}=0$ (Eq. `totc3b`).
- $\Pi|_{\pi\pi}$ is Eq. `totc1` (lines 655-658); $\Pi|_{V\pi}$ and $\Pi|_{A\pi}$ are proportional to $(1-M_V^2/M_A^2)$ and vanish for the NLO mass condition $M_A=M_V$, $M_P=\sqrt2M_S$ (Eq. `masses`, lines 702-705, 727, 743); $\Pi|_{S\pi}$ remains and cancels the leading large-$t$ term of $\Pi|_{\pi\pi}$ (line 743).
- With $V',A'$, all $R\pi$ cuts are non-zero; $V\pi$ is unchanged and $S\pi$ has $M_A\to M_{A'}$ (lines 961-963, 2535-2538). The conditions fix $M_{V'}^2=M_A^2$ and $M_{A'}$ (Eqs. lines 975-991).
