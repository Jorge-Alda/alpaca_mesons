---
title: Block-exclusion and sensitivity refit procedure for the joint vector fit
status: numerical
created: 2026-11-06
updated: 2026-11-06
bibtex_keys: [Alda:2026nrd, ParticleDataGroup:2024cfk, BaBar:2007ceh, Qin:2024ulb]
provenance: numerical-computation
---

# Block-exclusion and sensitivity refit procedure

## Inputs

- Best accepted endpoint `basin_probes/even_shared_widths`, the other three accepted endpoints and the fixed optimizer scales from [the endpoint comparison](joint-vector-fit-endpoint-comparison.md).
- Baseline local errors and block $\chi^2$ from [the adequacy assessment](joint-vector-fit-adequacy-assessment-method.md).
- Likelihood options from [the variation record](joint-vector-likelihood-variation-options.md).

## Method

1. Refit with scipy `least_squares` (trf, `max_nfev=5000`, `ftol=1e-8`, `xtol=1e-9`, `gtol=1e-7`, fixed optimizer scales, existing bounds) at quadrature $(2,16,16)$ with four-bin blocks, warm-started from the best accepted endpoint.
2. Re-evaluate the endpoint at $(4,128,128)$; report $\chi^2$, degrees of freedom (retained rows minus free coordinates), nominal $p$, per-block $\chi^2$ on retained rows, and the pulls of removed rows.
3. Block exclusion: one of the six blocks removed at a time. Near-singular directions of the fit-order Jacobian (times optimizer scales) with relative singular value below $10^{-6}$ are listed; such refits are labelled conditional. When a removed block is the only one sensitive to some coordinates (for example the BaBar block for the odd couplings entering only $e^+e^-\to KK\pi$ and the BaBar normalization and excited-pole residues), those Jacobian columns vanish and the default count can be zero or negative; the nominal $p$ is then not given. An alternative count, retained rows minus free coordinates with non-vanishing columns, is reported next to it (assumption: unconstrained coordinates do not absorb $\chi^2$).
4. Variations: starts are moved into the variation domain where needed (stage-3 $F_V,z_V$ for the source hypotheses; stage-3 second root and negative-source $F_V,z_V$ for the second branch; tensor products set to zero). Each endpoint receives the quadrature ladder and acceptance rule of the endpoint comparison. If the scaled distance to the best endpoint exceeds 1, the variation is refit from the other three accepted endpoints and the lowest accepted high-order $\chi^2$ is selected.
5. Parameter shifts are given in baseline local errors after converting even-sector coordinates to the baseline reference point.

## Output and verification

Executable: [`scripts/refit_joint_vector_sensitivity.py`](../../../scripts/refit_joint_vector_sensitivity.py). Ledgers: [`joint_vector_block_exclusion.json`](../../../data/vector_sector/joint_vector_block_exclusion.json) and [`joint_vector_sensitivity_variations.json`](../../../data/vector_sector/joint_vector_sensitivity_variations.json); unsuccessful refits and not-implemented variations are kept.

```text
.venv/bin/python scripts/refit_joint_vector_sensitivity.py --study exclusions --output data/vector_sector/joint_vector_block_exclusion.json
.venv/bin/python scripts/refit_joint_vector_sensitivity.py --study variations --output data/vector_sector/joint_vector_sensitivity_variations.json
```

[`scripts/test_refit_joint_vector_sensitivity.py`](../../../scripts/test_refit_joint_vector_sensitivity.py) runs a reduced exclusion ($\phi\to K\bar K$ block) and a tensor-zero variation at cheap quadrature, and checks the singular-direction detection and shift conversion.
