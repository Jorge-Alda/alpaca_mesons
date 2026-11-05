---
title: Resonance-Goldstone cuts in the NLO S-P correlator and their thresholds
status: literature
created: 2026-11-05
updated: 2026-11-05
bibtex_keys: [Rosell:2006dt]
provenance: literature
---

# $R\pi$ loops in $\Pi_S-\Pi_P$ at NLO in $1/N_C$

Source: `memory/sources/Rosell:2006dt/L8.tex`. [@Rosell:2006dt]

## Setup (assumptions of the source)

- Chiral limit, single-resonance approximation, only operators with $\mathcal{O}(p^2)$ chiral tensors; $U(3)$ matched to $SU(3)$ (lines 350-356, 497-510).
- The loop part is obtained from unitarity and analyticity (dispersion relation) rather than from Feynman diagrams, "avoiding all technicalities associated with the renormalization procedure" (lines 179-182). The undetermined subtraction constants are absorbed in $c_m,d_m,M_S,M_P$ (lines 366-372).
- Only the lowest two-particle cuts are kept: $\pi\pi$ and $R\pi$ ($R=S,P,V,A$). "Exchanges of two heavy resonances or higher multiplicity states are kinematically suppressed" (lines 312-318); no numerical threshold is quoted for $RR'$ here.

## Spectral functions and thresholds (lines 295-315, 320-340)

$\rho_S(t)=\frac{\theta(t)}{2}|F_S^{\pi\pi}|^2+\theta(t-M_P^2)(1-\frac{M_P^2}{t})|F_S^{P\pi}|^2+\theta(t-M_A^2)\frac{t^2}{4M_A^4}(1-\frac{M_A^2}{t})^3|F_S^{A\pi}|^2$,
$\rho_P(t)=\theta(t-M_V^2)\frac{t^2}{4M_V^4}(1-\frac{M_V^2}{t})^3|F_P^{V\pi}|^2+\theta(t-M_S^2)(1-\frac{M_S^2}{t})|F_P^{S\pi}|^2$.

In the chiral limit the $R\pi$ thresholds sit at $t=M_R^2$: $V\pi$ and $S\pi$ in the pseudoscalar channel at $M_V^2$ and $M_S^2$; $P\pi$ and $A\pi$ in the scalar channel at $M_P^2$ and $M_A^2$. With $M_V\simeq0.77$ GeV the $V\pi$ threshold ($\sqrt t=M_V$) is below 1 GeV (numerical value of $M_V$ is my input, not given at these lines); the others are not stated numerically here. $F_S^{A\pi}=0$ within the single-resonance approximation (lines 351-356).

## Effect on the $\mathcal{O}(p^2)$ term

- Eq. `eq:Pi_SP` (line 373): $\Pi(t)=2B_0^2\{\frac{8c_m^{r2}}{M_S^{r2}-t}-\frac{8d_m^{r2}}{M_P^{r2}-t}+\frac{F^2}{t}+\Delta\Pi(t)|_\rho\}$. The Goldstone pole residue $F^2/t$ carries no explicit NLO correction. $\Delta\Pi|_\rho$ is a finite, once-subtracted dispersive function given in appendix A (lines 700-725: $\eta\pi$, $V\pi$, $A\pi$, $S\pi$, $P\pi$ with prefactor $\frac{n_f}{2}\frac{1}{16\pi^2}$). At large $t$, $\Delta\Pi|_\rho\to\frac{F^2}{t}\delta^{(1)}_{\rm NLO}+\frac{F^2M_S^2}{t^2}(\delta^{(2)}_{\rm NLO}+\tilde\delta^{(2)}_{\rm NLO}\ln\frac{-t}{M_S^2})$ (line 385), with $\delta^{(m)}_{\rm NLO}=\frac{3M_S^2}{32\pi^2F^2}\{\dots\}$ (line 410).
- The $\mathcal{O}(p^2)$-relevant constants are traded through the short-distance conditions Eq. `eq:NLO_rel` (lines 401-405): $F^2(1+\delta^{(1)}_{\rm NLO})-8c_m^{r2}+8d_m^{r2}=0$, hence $c_m^r,d_m^r$ (Eq. `dmr`, lines 445-446) absorb the $R\pi$ loops. This is an $R\pi$ loop effect of order $\frac{M_S^2}{32\pi^2F^2}$ on $c_m,d_m$ ($\mathcal{O}(p^2)$ resonance couplings), which are $1/N_C$ corrections. The $\mathcal{O}(p^4)$ LEC $L_8$ gets $\Gamma_8$-type logs (Eq. `eq:Pi_chpt`, lines 199-203).

## Verdict

- O(p^2) affected: **yes**, in the sense that short-distance matching shifts the leading-order resonance couplings $c_m,d_m$ at NLO (a renormalization of the $\mathcal{O}(p^2)$ resonance Lagrangian), by $R\pi$ loops with $R=S,P,V$ ($A$ vanishes). The source does not renormalize $F$ or $B_0$ explicitly (chiral limit, $F^2/t$ pole kept as is, $B_0^2$ prefactor).
- Conditions: requires the short-distance constraints that fix the form factors and Eq. `eq:NLO_rel`; the result relies on $\mathcal{O}(p^2)$-only operators, single resonance approximation, and $U(3)\to SU(3)$ matching.
- Caveat: dispersive derivation; divergences are not shown, so the UV structure at $\mathcal{O}(p^2)$ must be read from Rosell:2004mn and Rosell:2005ai.
