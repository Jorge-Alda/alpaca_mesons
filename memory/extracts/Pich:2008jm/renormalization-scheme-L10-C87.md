---
title: Renormalization scheme and mu dependence of L_10 and C_87
status: literature
created: 2026-11-05
updated: 2026-11-05
bibtex_keys: [Pich:2008jm]
provenance: literature
---

# Renormalization conventions for $L_{10}^r(\mu)$, $C_{87}^r(\mu)$

Source: `memory/sources/Pich:2008jm/src/L10.tex`. [@Pich:2008jm]

## Low-energy expansion (Eq. `eq:Pi_chpt`, lines 1070-1079)

$\Pi(t)=\frac{2F^2}{t}-8L_{10}^r(\mu)-\frac{\Gamma_{10}}{4\pi^2}(\frac53-\ln\frac{-t}{\mu^2})+\frac{t}{F^2}[16C_{87}^r(\mu)-\frac{\Gamma_{87}^{(L)}}{2\pi^2}(\frac53-\ln\frac{-t}{\mu^2})+\mathcal{O}(N_C^0)]+\mathcal{O}(t^2)$, with $\Gamma_{10}=-1/4$ and $\Gamma_{87}^{(L)}=-L_9/2$. $F^2$, $L_{10}$, $C_{87}/F^2$ are $\mathcal{O}(N_C)$; $\Gamma_{10}$, $\Gamma_{87}^{(L)}/F^2$ are $\mathcal{O}(N_C^0)$, i.e. NLO. The $\chi$PT scheme is the standard one, with $\frac53-\ln(-t/\mu^2)$ as the finite part (no extra constants). $\mu_0=770$ MeV (line 1218).

## How $\mu$ dependence arises and cancels

- No Feynman-diagram renormalization is done in the main text: the loop part is $\sum\Pi|_{m_1,m_2}$ from unsubtracted dispersion relations, Eq. `eq.master` (lines 2820-2824) and Eq. `DeltaPi` (lines 2826-2842). The only freedom is the real part $D(M_R^2)$ and $D'(M_R^2)$, i.e. a prescription for $F_V^r$ and $M_V^r$ (lines 2851-2854).
- The $\pi\pi$ cut reproduces the $\chi$PT $\ln(-t)$ exactly at low energies; the $R\pi$ cuts have thresholds and give analytic terms at low $t$ that go into the LECs (comment block lines 1158-1163, which is commented out in the source, so it is the authors' unpublished remark).
- The explicit $\ln(M_V^2/\mu^2)$ in Eqs. `eq:L10_SRA` (line 1176), `eq:C87_SRA` (line 1278), `eq:L10_2Res` (line 1236), `eq:C87_2Res` (line 1294) is the only $\mu$ dependence. The tree-level $F_V^r$, $F_A^r$, $M^r$ and $\delta^{(m)}_{\rm NLO}$ are $\mu$ independent. It reproduces $\partial L_{10}^r/\partial\ln\mu^2=-\Gamma_{10}/(32\pi^2)=1/(128\pi^2)\approx0.8\cdot10^{-3}$ (line 1125), since the SRA expression contains $-\frac{1}{128\pi^2}\ln\frac{M_V^2}{\mu^2}$ (my check).
- Diagrammatic version (lines 2876-2918): UV divergences $\gamma_i\lambda_\infty$ of $\Pi_{VV}$ cancel by $\delta F_V$, $\delta M_V^2$ (lines 2907-2913). Polynomial divergences $\gamma_{-m}t^m\lambda_\infty$ cancel identically because of short-distance behaviour (lines 2902-2905). Setting $c_1=c_2=0$ makes $F_V^r$, $M_V^r$ scale independent (lines 2918-2920). $\mathcal{O}(p^4)$ $\tilde L_i$ of the resonance Lagrangian do not run (lines 2866-2872).
- Running of $L_{10}^r$ and $C_{87}^r$ is thus generated only by the $\pi\pi$ loop with chiral logs ($\Gamma_{10},\Gamma_{87}^{(L)}$), while resonance loops give $\mu$-independent shifts.
- The LECs are independent of quark masses, so the chiral limit loses no information (lines 205-206).
