---
title: Validity range and truncation uncertainty of the U(3) NLO theory with explicit resonances
status: mixed
created: 2026-11-05
updated: 2026-11-05
bibtex_keys: [Kaiser:1998ds, Kaiser:2000gs, Pich:2008jm, Rosell:2004mn, Rosell:2005ai, ParticleDataGroup:2024cfk]
provenance: assumption
---

# Validity range and truncation uncertainty

## Assumptions and inputs

- The declared calculation is $\mathcal O(\delta^2)$ in the Goldstone/ALP sector ([power-counting.md](power-counting.md)). The omitted terms are $\mathcal O(\delta^3)$, suppressed by one more factor of $\delta$ relative to the $\delta^2$ terms.
- The excluded loop classes, and the conditions under which their exclusion is consistent, are in [loop-classes.md](loop-classes.md). The inputs to be fixed from data are listed there and in [counterterm-basis.md](counterterm-basis.md).
- Inputs of the estimate:
  - $N_c=3$;
  - $f_\pi=0.0921\,$GeV, the value of $F$ at this order in the project convention (PDG $130.2\,$MeV$/\sqrt2$ rounded);
  - $M_K=0.493677\,$GeV ($K^\pm$);
  - the $\rho\pi$ threshold $0.91483\,$GeV, from $M_\rho=0.77526\,$GeV and $m_\pi=0.13957\,$GeV [@ParticleDataGroup:2024cfk].
- **Assumption:** $\epsilon_R=0.30$ is the size of the excluded resonance-loop effect relative to the NLO correction. It is taken from the NLO-in-$1/N_c$ shifts of $L_{10}^r$ and $C_{87}^r$ of about 17% and 28% in Pich:2008jm, in which the second multiplets $V'$ and $A'$ account for about 30–40% of the correction ([extract](../../extracts/Pich:2008jm/nlo-size-L10-C87-numerics.md)) [@Pich:2008jm]. It is a stated working value, not a result.
- **Assumption:** the truncation uncertainty and the excluded-loop uncertainty are independent, so they are added in quadrature.

The executable estimate is [`scripts/derive_nlo_u3_resonance_counting.py`](../../../scripts/derive_nlo_u3_resonance_counting.py):

```bash
.venv/bin/python scripts/derive_nlo_u3_resonance_counting.py --truncation 0.5 1.0
```

## Validity range

The framework is declared valid for external momenta and ALP masses with $\sqrt s\lesssim1\,$GeV. The reasons are as follows.

- The chiral scale is $\Lambda_\chi=4\pi f_\pi=1.157\,$GeV. The expansion parameter $s/\Lambda_\chi^2$ reaches $0.75$ at $1\,$GeV and exceeds $1$ for $\sqrt s>\Lambda_\chi$. Convergence is poor in the upper part of the range, and the framework is not expected to hold above $\Lambda_\chi$.
- Resonance-resonance cuts lie at $(M_R+M_{R'})^2\gtrsim2\,$GeV$^2$, outside the range. The lowest resonance–Goldstone cut, $\rho\pi$, opens at $0.915\,$GeV, inside it ([loop-classes.md](loop-classes.md)).
- The $\rho$, $\omega$ and $K^*$ poles, and the $\phi$ pole at about $1.02\,$GeV, lie inside or at the edge of the range. See the next section.

## Truncation uncertainty

For a generic observable whose NLO correction (the relative shift from the $\delta$ and $\delta^2$ terms) is $\Delta$, supplied by the downstream calculation, the omitted $\mathcal O(\delta^3)$ terms are estimated as

$$
u_{\rm trunc}=|\Delta|\,\delta_{\rm est}(s),\qquad
\delta_{\rm est}(s)=\max\!\left(\frac1{N_c},\ \frac{s}{\Lambda_\chi^2},\ \frac{M_K^2}{\Lambda_\chi^2}\right).
$$

The three terms are the three expansion parameters of the counting ($1/N_c$, $p^2$ and $m_q$). Numerically:

- $M_K^2/\Lambda_\chi^2=0.182<1/N_c=0.333$, so the quark-mass term never dominates.
- $s/\Lambda_\chi^2$ overtakes $1/N_c$ at $\sqrt s=\Lambda_\chi/\sqrt3=0.668\,$GeV.

The excluded resonance loops (class R-int, which starts at $\delta$) are added in quadrature,

$$
u_{\rm tot}=|\Delta|\sqrt{\delta_{\rm est}^2(s)+\epsilon_R^2\,\theta},
$$

where $\theta$ is set by the inputs:

| Inputs | $\theta$ |
|---|---|
| fixed from data, with the short-distance constraints imposed (condition of [loop-classes.md](loop-classes.md)) | $1$ for $\sqrt s\geq0.91483\,$GeV, otherwise $0$ |
| from large $N_c$, resonance saturation or LO short-distance estimates | $1$ at all $\sqrt s$ |

In the first row, below the $\rho\pi$ threshold the effect of the resonance loops is analytic and is absorbed in the constants fixed from data. Above it the non-analytic part is open and is not included. In the second row, the analytic part is not absorbed either.

## Result (script output, per unit $|\Delta|$)

| $\sqrt s$ (GeV) | $\delta_{\rm est}$ | dominant term | $u_{\rm tot}/|\Delta|$, inputs from data | $u_{\rm tot}/|\Delta|$, inputs not from data |
|---|---|---|---|---|
| 0.3 | 0.333 | $1/N_c$ | 0.333 | 0.448 |
| 0.5 | 0.333 | $1/N_c$ | 0.333 | 0.448 |
| 0.7 | 0.366 | $s/\Lambda_\chi^2$ | 0.366 | 0.473 |
| 0.9 | 0.605 | $s/\Lambda_\chi^2$ | 0.605 | 0.675 |
| 1.0 | 0.747 | $s/\Lambda_\chi^2$ | 0.805 | 0.805 |
| 1.2 | 1.075 | $s/\Lambda_\chi^2$ | 1.116 | 1.116 |

Remarks:

- The estimate is a rule of thumb. It does not give a probability.
- At $1\,$GeV the uncertainty is about 80% of the NLO correction, and above $\Lambda_\chi$ (the $1.2\,$GeV row) it exceeds the correction. The NLO result is then not quantitatively reliable.
- The step at the $\rho\pi$ threshold, between $0.9$ and $1.0\,$GeV in the table, is a conservative simplification. The non-analytic part of the loops turns on smoothly above the threshold.
- The estimate does not replace the $\mu$ dependence: varying $\mu\in[0.5,1.0]\,$GeV changes $L_i^r(\mu)$ by the amounts of [subtraction-scheme-and-scales.md](subtraction-scheme-and-scales.md). That variation tests only the loop part of the calculation.

## Where the counting does not apply

- **Above the range.** For $\sqrt s>1\,$GeV the framework is outside its declared validity. For $\sqrt s\geq\Lambda_\chi$, $\delta_{\rm est}\geq1$ and the expansion parameter is not small.
- **Near a resonance pole.** The tree propagator $1/(M_R^2-s)$ is used without a width. When $|M_R^2-s|$ is of the order of $M_R\Gamma_R$ it fails, and the strict counting is not applicable. The line-shape prescription of [vector-width-lineshape-prescriptions.md](../../rcht/foundations/vector-width-lineshape-prescriptions.md) replaces it there, as a resummation outside the strict counting ([loop-classes.md](loop-classes.md)). The windows are set by the downstream calculation with the PDG widths. Results in these windows are reported together with the strict-order result.
- **Near thresholds.** At the $\rho\pi$ and $\omega\pi$ thresholds (0.915 and $0.922\,$GeV) the excluded non-analytic part of R-int becomes large. The estimate above covers it only as the constant $\epsilon_R$.

## Verification

The script asserts the following, and all runs exit with code 0 (`--truncation`, `--truncation 0.5 1.0`, `--orders`, `--loops`, `--scheme`, the default run and `--test`):

- $\delta_{\rm est}$ does not decrease with $s$;
- the maximum is $1/N_c$ below the crossover $0.668\,$GeV and $s/\Lambda_\chi^2$ above it, and the $M_K$ term never dominates;
- $\delta_{\rm est}(1.0\,{\rm GeV})$ lies between 0.7 and 0.8, and $\delta_{\rm est}<1$ exactly for $\sqrt s<\Lambda_\chi$;
- the excluded-loop term is zero below $0.91483\,$GeV and nonzero at and above it for inputs fixed from data;
- the quadrature total is at least as large as each term.

The numbers in the tables are copied from the script output.
