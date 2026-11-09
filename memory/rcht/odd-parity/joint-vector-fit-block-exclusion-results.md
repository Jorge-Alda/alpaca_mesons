---
title: Block-exclusion refits of the joint vector-sector fit
status: numerical
created: 2026-11-09
updated: 2026-11-09
bibtex_keys: [Alda:2026nrd, ParticleDataGroup:2024cfk, BaBar:2007ceh]
provenance: numerical-computation
---

# Block-exclusion refits

## Inputs and provenance

- Procedure, warm start, optimizer settings and gates: [refit procedure](joint-vector-block-exclusion-and-variation-refits.md).
- Ledger: `data/vector_sector/joint_vector_block_exclusion.json`, written by `scripts/refit_joint_vector_sensitivity.py --study exclusions`.
- Table printed by `.venv/bin/python scripts/summarize_joint_vector_fit_adequacy.py`.
- Each refit removes one block, starts from `basin_probes/even_shared_widths`, runs at $(2,16,16)$ with `max_nfev=5000`, and is re-evaluated at $(4,128,128)$.
- Provenance qualification: this ledger records earlier SHA-256 hashes of `joint_vector_likelihood_options.py` and `refit_joint_vector_sensitivity.py` than the current files. The pole options were added after these refits were run. The baseline residuals at default options are tested to be bitwise identical.

## Result

| excluded block | optimizer | remaining $\chi^2$ / dof / $p$ | excluded-block $\chi^2$ (rows) | conditional |
|---|---|---|---|---|
| odd | stopped at max_nfev | 166.3 / 31 / $1.9\times10^{-20}$ | 86952 (12) | yes |
| masses and leptonic widths | converged | 75.2 / 36 / $1.4\times10^{-4}$ | $1.5\times10^{6}$ (7) | no |
| shared total widths | stopped at max_nfev | 57.5 / 36 / $1.3\times10^{-2}$ | 15585 (7) | no |
| $\phi\to K\bar K$ | converged | 122.2 / 41 / $5.3\times10^{-10}$ | 776 (2) | no |
| BaBar $KK\pi$ | converged | 130.1 / $-2$ (9 without zero Jacobian columns) / undefined ($1.1\times10^{-23}$ with 9) | $1.5\times10^{15}$ (45) | yes |
| normalization constraints | stopped at max_nfev | 181.9 / 41 / $1.2\times10^{-19}$ | 598 (2) | yes |

"Conditional" means that near-singular Jacobian columns were found, so some coordinates are not constrained by the remaining blocks.

## Interpretation

- Removing the masses and leptonic widths lowers the total $\chi^2$ by about 150, although that block contributes only 12.9 in the baseline. The vector masses therefore prevent the other blocks, mainly the shared widths, from being fitted. The excluded mass rows then sit far from their data ($\chi^2\sim10^6$).
- Removing the shared total widths gives the lowest remaining $\chi^2$ (57.5 for 36 dof). This refit did not converge, so 57.5 is an upper bound on its minimum. Five width-pull coordinates lose their direct constraint, so the degrees-of-freedom count is not well defined. The result is qualitative: without the width rows the remaining data are not grossly inconsistent.
- Removing $\phi\to K\bar K$ lowers the $\chi^2$ by about 104, so these two rows are also in tension with the rest.
- Removing the odd block or the normalization constraints does not relieve the tension.
