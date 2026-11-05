---
title: Numerical size of the NLO-in-1/N_C corrections to L_10 and C_87
status: literature
created: 2026-11-05
updated: 2026-11-05
bibtex_keys: [Pich:2008jm]
provenance: literature
---

# Size of the NLO corrections to $L_{10}^r$, $C_{87}^r$

Source: `memory/sources/Pich:2008jm/src/L10.tex`. [@Pich:2008jm] All numbers at $\mu_0=770$ MeV unless stated.

## Inputs (lines 1186-1187, 1245-1246)

$M_V=770\pm5$ MeV, $M_S=1090\pm110$ MeV, $F=89\pm2$ MeV, $M_P=1300\pm100$ MeV, $M_A^r=1000\pm50$ MeV. $M_A$ from other sources: $938\pm13$ MeV (WSR with $F_V$ from $\Gamma(\rho\to e^+e^-)$), $1160\pm40$ MeV (from $G_V$ with Eq. `eq:GV`), $998\pm49$ MeV (`jorge_vicent`, $\pi\to e\nu\gamma$) (lines 1199-1203). $F_V=155.7\pm1.5$, $G_V=66.7\pm0.9$ MeV (line 1197). In the 2-multiplet case $\epsilon_i=0.0\pm0.1$, $\epsilon_3\approx0.007$ (lines 873-878, 1254).

## Results

| Quantity | LO large-$N_C$ (Eq. `eq:L10-LO`, lines 1101-1103) | NLO, SRA | NLO, with $V',A'$ | Final |
|---|---|---|---|---|
| $L_{10}^r(\mu_0)$ | $-5.3\cdot10^{-3}$ | $(-5.2\pm0.4)\cdot10^{-3}$ (Eq. `eq:L10_SRA_num`, line 1216) | $(-3.6\pm0.9\pm0.3)\cdot10^{-3}$ (line 1258) | $(-4.4\pm0.9)\cdot10^{-3}$ (line 1386) |
| $C_{87}^r(\mu_0)$ | $4.3\cdot10^{-5}$ | $(3.9\pm0.6)\cdot10^{-5}$ (line 1282) | $(2.2\pm1.0\pm0.4)\cdot10^{-5}$ (line 1305) | $(3.1\pm1.1)\cdot10^{-5}$ (line 1388) |

The second error is from $\epsilon_{1,2}$. Final values are the midpoint of SRA and two-multiplet results with the larger errors (lines 1391-1392). Phenomenology: $L_{10}^r(\mu_0)=(-5.5\pm0.7)\cdot10^{-3}$, $(-5.13\pm0.19)\cdot10^{-3}$ (ALEPH fit); $C_{87}=(4.5\pm0.4)\cdot10^{-5}$ (Pad\'e, `C87-Peris`) (lines 1106-1112).

## Statements useful for the truncation estimate

- The SRA error "does not include the systematic uncertainties due to neglecting higher resonance effects and the inconsistencies between form-factor constraints" (lines 1219-1222).
- Adding $V',A'$ shifts $L_{10}^r$ by about $+1.6\cdot10^{-3}$ and $C_{87}^r$ by about $-1.7\cdot10^{-5}$ (differences of the central values above; my arithmetic). The source calls the $\epsilon_i$ corrections "tiny" but notes the $V\pi$/$A\pi$ contributions cancel against the tree-level $V,A$ exchange (lines 1262-1265).
- NLO size estimated from the $\mu$ dependence: $\partial L_{10}^r/\partial\ln\mu^2=-\Gamma_{10}/(32\pi^2)=0.8\cdot10^{-3}$, $\partial C_{87}^r/\partial\ln\mu^2=\Gamma_{87}^{(L)}/(32\pi^2)=-1.1\cdot10^{-5}$ (lines 1124-1129); LO and NLO agree for $\mu\sim0.5$ GeV (lines 1397-1400, Fig. `run-LECs`).
- The source does not give numerical values of $\delta^{(1,2)}_{\rm NLO}$, of the separate $V\pi$/$A\pi$/$S\pi$/$\pi\pi$ contributions, or of the shift of $M_A$, $F_A$ at NLO. Only formulas (lines 744-768, 995-1005, 1011-1019) and the quoted range $M_A^r=1000\pm50$ MeV are given. The "$\delta_{\rm NLO}$ dominated by..." statement is not present in this paper.
- Width: $\Gamma_{V\to\pi\pi}=\frac{n_f}{2}\frac{G_V^2M_V^3}{48\pi F^4}$ (line 676); the $1/N_C$ suppression is stated as the factor $\frac1\pi\frac{\Gamma}{M_V}$ in Eq. `eq.Gamma` (lines 668-673), i.e. $\Gamma_\rho/(\pi M_\rho)\approx0.06$ for $\Gamma=149.4$ MeV (my arithmetic).
