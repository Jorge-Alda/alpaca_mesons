---
title: NLO Weinberg sum rules and effect of loops on F, F_V, F_A, M_V, M_A
status: literature
created: 2026-11-05
updated: 2026-11-05
bibtex_keys: [Pich:2008jm]
provenance: literature
---

# Loop effects on $F$, $F_V$, $F_A$, masses and the NLO Weinberg sum rules

Source: `memory/sources/Pich:2008jm/src/L10.tex`. [@Pich:2008jm]

## LO relations (SRA)

- $\Pi=\frac{2F_V^2}{M_V^2-t}-\frac{2F_A^2}{M_A^2-t}+\frac{2F^2}{t}$ (Eq. `eq:LO`, lines 474-476). Requiring $\Pi\sim1/t^3$ gives the Weinberg sum rules $F_V^2-F_A^2=F^2$, $F_V^2M_V^2-F_A^2M_A^2=0$ (Eq. `eq:WSR-LO`, lines 561-566), hence $F_V^2=F^2\frac{M_A^2}{M_A^2-M_V^2}$, $F_A^2=F^2\frac{M_V^2}{M_A^2-M_V^2}$ (Eq. `eq:FV_FA`), $F_VG_V=F^2$ (Eq. `eq:FG=F2`, line 451), $G_V^2=F^2\frac{M_A^2-M_V^2}{M_A^2}$ (Eq. `eq:GV`).

## NLO relations

- Large-$t$ expansion of the loop part (lines 693-700): $\Pi|_\rho=\frac{2F^2}{t}(\delta^{(1)}_{\rm NLO}+\tilde\delta^{(1)}_{\rm NLO}\ln\frac{-t}{M_V^2})+\frac{2F^2M_V^2}{t^2}(\delta^{(2)}_{\rm NLO}+\tilde\delta^{(2)}_{\rm NLO}\ln\frac{-t}{M_V^2})+\mathcal{O}(t^{-3})$.
- Log terms must vanish: $\tilde\delta^{(1)}_{\rm NLO}=\tilde\delta^{(2)}_{\rm NLO}=0$ (Eq. `masses`, line 705).
- NLO Weinberg sum rules, Eq. `eq:NLO_rel` (lines 711-716): $F^2(1+\delta^{(1)}_{\rm NLO})-F_V^{r2}+F_A^{r2}=0$ and $F^2M_V^2\delta^{(2)}_{\rm NLO}-F_V^{r2}M_V^{r2}+F_A^{r2}M_A^{r2}=0$.
- Solution, Eq. `dmr` (lines 718-725): $F_V^{r2}=F^2\frac{M_A^{r2}}{M_A^{r2}-M_V^{r2}}(1+\delta^{(1)}_{\rm NLO}-\frac{M_V^2}{M_A^2}\delta^{(2)}_{\rm NLO})$, $F_A^{r2}=F^2\frac{M_V^{r2}}{M_A^{r2}-M_V^{r2}}(1+\delta^{(1)}_{\rm NLO}-\delta^{(2)}_{\rm NLO})$.
- With $V',A'$ (Eq. `dmrbis`, lines 1011-1019): $\delta^{(m)}_{\rm NLO}\to\epsilon_m+\delta^{(m)}_{\rm NLO}$ for $m=1,2$, where $\epsilon_1=\frac{F_{A'}^2-F_{V'}^2}{F^2}$, $\epsilon_2=\frac{F_{A'}^2M_{A'}^2-F_{V'}^2M_{V'}^2}{F^2M_V^2}$ (Eq. `FVprime`, lines 865-870).
- $\delta^{(m)}_{\rm NLO}$ in the SRA: Eqs. lines 744-768, prefactor $\frac{n_f}{2}\frac{M_V^2}{48\pi^2F^2}$ (depends on $M_V$, $M_S$ only after $M_A=M_V$, $M_P=\sqrt2M_S$). With $V',A'$: lines 995-1005 (sum of $\xi^{(m)}_{A\pi},\xi^{(m)}_{P\pi},\xi^{(m)}_{V\pi},\xi^{(m)}_{S\pi}$, Appendix B).

## What is and is not renormalized

- $F$ (chiral-limit pion decay constant) appears without correction in the Goldstone pole $2F^2/t$ of the correlator and in the $2F^2\delta(t)$ term of ${\rm Im}\Pi_A$ (Eqs. `eq:ImSS`, `eq:Pi_structure`, lines 493, 531). The loop part has no pole at $t=0$ beyond what is shown. The source states no pion wave-function renormalization and no explicit shift of $F$ from resonance loops; $F$ is an input, $F=(89\pm2)$ MeV (line 1187). It is not discussed whether $F$ in the pole equals the LO $\mathcal{O}(p^2)$ constant.
- The resonance parameters are renormalized: $F_V=F_V^r+\delta F_V$, $M_V^2=M_V^{r2}+\delta M_V^2$ (lines 2882-2890, 2907-2913). The $A\pi$/$V\pi$/$S\pi$ cuts therefore shift $F_A^r$, $M_A^r$ (lines 1204-1209) through the NLO sum rules above, not $F$.
- The axial-channel loop ($V\pi$, $S\pi$) enters $F_A^r$ only through $\delta^{(1,2)}_{\rm NLO}$ in Eq. `dmr`; the $\mathcal{O}(p^4)$ operators $\tilde L_i$ of the resonance Lagrangian are not renormalized (lines 2866-2872, citing the one-loop generating functional analysis `RPP:05`).
- In the SRA the NLO constraints force $M_A=M_V$ (Eq. `masses`, line 727), so $V\pi$ and $A\pi$ vanish; this is a statement of the single-resonance approximation, relieved by $V',A'$ (Section "Numerical impact of heavier vector and axial multiplets", lines 772-1019). The paper takes $M_A^r=1.0\pm0.05$ GeV as an input instead (lines 1204-1209).
