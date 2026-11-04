# Proposal

## Why

The stored joint vector-sector fit has a large even-sector residual contribution, while the conditional four-row $G_V$ fit is explicitly correlated with the provisional fit inputs. The current joint-fit diagnostic already shows substantial charged/neutral $K^*$ tension and a smaller $phi\to K\bar K$ tension, but it does not provide a complete, reproducible audit of data provenance, mass and width prescriptions, duplicated measurements, or covariance assumptions needed to interpret those numbers.

## What Changes

- Add a reproducible blockwise audit of the shared even-sector and conditional $G_V$ constraints, including baseline block contributions, tensions that can be isolated analytically, and their effect on the joint likelihood.
- Trace each input row to its source and distinguish measured physical masses used in phase space from fitted tensor masses and averaged masses used in vertices or propagators.
- Audit shared total widths, the two charged/neutral $K^*$ partial-width summaries, $K^*\to K\pi$ saturation, and correlated $phi\to K\bar K$ branching fractions so overlapping information is not counted twice.
- Compare only alternatives justified by the source data or the declared model, and state the limits from unavailable cross-covariances, frozen-parameter uncertainty, and optimizer nonconvergence.
- Add an indexed research record and update the relevant memory indexes with the numerical diagnosis and retained assumptions.

## Capabilities

### New Capabilities

- `vector-sector-fit-audit`: Reproducibly diagnose likelihood blocks with shared even-sector and conditional decay constraints, including input provenance, overlap handling, covariance assumptions, and justified alternatives.

### Modified Capabilities

None.

## Impact

The work will use `scripts/fit_babar_kkpi_global_vector_sector.py`, `scripts/fit_vector_gv.py`, and the corresponding JSON/CSV inputs under `data/vector_sector/` and `data/odd_parity/`. It will add an executable audit calculation under `scripts/`, a numerical research record under `memory/rcht/`, and the required parent index links. The current global fit snapshot reports $\chi^2=307.985$ for 43 degrees of freedom and did not converge within 1500 function evaluations; any likelihood change must therefore be reported as a comparison against that stored snapshot, not as a converged refit.
