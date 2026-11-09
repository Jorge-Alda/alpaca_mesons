---
title: Adequacy verdict for the joint vector-sector fit
status: mixed
created: 2026-11-09
updated: 2026-11-09
bibtex_keys: [Alda:2026nrd, ParticleDataGroup:2024cfk, BaBar:2007ceh]
provenance: numerical-computation
---

# Adequacy verdict

## Criterion (assumption)

The criterion is fixed in `openspec/changes/assess-joint-vector-fit-adequacy/design.md`. The model is adequate when the nominal total $p>0.01$ and no block has nominal $p<10^{-3}$. The nominal $p$-values are not corrected for bounded or nonlinear parameters.

## Verdict (numerical)

**Inadequate.** This holds for the baseline model, in which the single-plus-double-trace source term $y\,\mathrm{tr}(\chi_+Q)\mathbb{1}$ is absent. That variation was deferred by user decision and was not tested.

- The total gives $\chi^2=225.73$ for 43 dof, $p=1.3\times10^{-26}$, at all four accepted minima ([goodness of fit](joint-vector-fit-goodness-of-fit.md)).
- Three blocks fail the block condition:
  - the shared total widths, $p=4\times10^{-20}$;
  - $\phi\to K\bar K$, $p=1.3\times10^{-4}$;
  - the odd block, $p=4.4\times10^{-4}$.
- None of the 47 variations that were run reaches $p>10^{-20}$, whether accepted or not ([variation results](joint-vector-fit-sensitivity-variation-results.md)).
- The BaBar $KK\pi$ block is well described, with $\chi^2=52.8$ for 45 bins and $p=0.20$. Its residuals show no coherent energy trend.

**Consequence:** the parameter estimates at the accepted minima are not interpretable. They are not passed to the identifiability study (alp-uzh.7).

## Responsible rows (diagnostic, from block exclusions)

The evidence comes from the [block-exclusion refits](joint-vector-fit-block-exclusion-results.md). Several of those refits are conditional or did not converge.

1. **$K^*$ total widths.** $\Gamma_{K^{*+}}$ (pull $-8.27$, $\chi^2=68.4$) and $\Gamma_{K^{*0}}$ (pull $-4.39$, $\chi^2=19.3$) together give 87.7 of the 225.7.
   - The $K^*$ widths cannot be fitted while the vector masses are constrained. Removing the mass block lowers the total by about 150.
   - The mechanism is not identified. The charged $K^*$ mass pull of $-3.20$ suggests the $K^*$ mass is involved, but the mass block also constrains $m_\rho$, $\Delta_V$ and the $\omega,\phi$ masses.
2. **$\phi\to K\bar K$ and $\Gamma_\phi$.** The whitened $\phi\to K^+K^-$ component has pull $+3.50$ and $\Gamma_\phi$ has $+2.83$. Removing the $\phi\to K\bar K$ block lowers the total by about 104.
3. **Odd block.** It is driven by $\rho^0\to\eta\gamma$ (pull $-3.27$), together with $K^{*0}\to K^0\gamma$ and $\eta'\to\gamma\gamma$. $\Gamma_{\eta'}$ has a pull of $-3.31$ in the width block. Removing the odd block does not relieve the total (166.3 for 31 dof, conditional refit).

Possible physical causes, and the extensions for the follow-up issue, are recorded separately as conjectures in [candidate causes](joint-vector-fit-inadequacy-candidate-causes.md).
