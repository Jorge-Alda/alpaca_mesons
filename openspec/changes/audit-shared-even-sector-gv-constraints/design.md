# Design

## Context

The current stored global fit is produced by `scripts/fit_babar_kkpi_global_vector_sector.py`. Its reported block contributions include $\chi^2=102.455$ from seven shared-width rows; the stored fit has $\chi^2=307.985$ for 43 degrees of freedom and terminated at its 1500-evaluation limit. A separately stored diagnostic extracts 12 even-sector and $G_V$ residuals, with $\chi^2=125.029$; its exact common-coupling subfits give minimum $\chi^2=22.211$ for the charged/neutral $K^*$ width pair and $4.924$ for the two $\phi\to K\bar K$ branching fractions. These are diagnostic subsets and must not be added to the joint total.

The conditional input record documents the two $K^*$ partial-width summaries, the correlated $\phi$ branching fractions, the source-level provenance, and missing cross-correlations. `memory/rcht/foundations/vector-gv-fit.md` also explains that the conditional four-row fit overlaps existing provisional constraints. The global fit additionally adopts $K^*\to K\pi$ saturation, shared on-shell widths, charge-mean $K^*$ propagator width, averaged vertex masses, and physical phase-space masses. The audit must read these conventions from current input artifacts instead of treating a prior conditional fit as independent data.

## Goals / Non-Goals

**Goals:** Produce a deterministic blockwise and row-level diagnosis, quantify supported comparisons against the stored baseline, and record why each retained data and model assumption is used.

**Non-Goals:** Claim a converged global minimum, combine overlapping conditional and joint-fit $\chi^2$ values, infer unavailable covariance entries, or introduce new fit parameters solely to reduce residuals.

## Decisions

- Keep the existing global fit output as the baseline snapshot. Read its block values, shared-even diagnostic, input checksums, optimizer status, and assumptions; label any comparison with that snapshot and do not recast it as a converged fit.
- Build a separate executable diagnostic, rather than changing the fit objective. This preserves the baseline and makes analytic pair minima and alternative prescriptions easy to inspect without silently changing the accepted fit inputs.
- Use `gv_measurements.json` and its source metadata as the authority for the four conditional rows and their covariance: independent quoted errors for the charged and neutral $K^*$ summaries, and the stated $-0.72$ correlation for the two $\phi$ branching fractions. Keep missing cross-covariances explicitly zero only as the adopted approximation. Do not add the conditional four-row $\chi^2$ to the joint likelihood.
- Separate mass roles in all predictions: PDG physical vector and daughter masses for phase space, fitted tensor masses for the on-shell vertex normalization, and the declared charge-mean mass/width inputs for the shared propagator. Compare charge-specific and isospin-averaged prescriptions only where the source data and implemented model define both consistently.
- Treat $K^*\to K\pi$ saturation as a model assumption, not a measured identity. Report the consequence of the charged and neutral width rows under that assumption; evaluate a non-saturated alternative only if sourced channel fractions and their covariance are available. Otherwise record it as unavailable instead of inventing a constraint.
- Recompute diagnostic residuals and subset minima with an explicitly stated covariance and fixed/free parameter set. Distinguish raw block contributions, projected pair minima, and changes from any alternative model prescription so quantities with different fit conditioning are not compared as if they were interchangeable.
- Store machine-readable results beside the existing vector-sector data, and pair them with an indexed memory record that names input paths, hashes, commands, assumptions, and limitations. Do not overwrite source inputs or the baseline fit JSON.

## Risks / Trade-offs

- [The baseline optimizer did not converge] → Label its values as the current stored snapshot and avoid claims that the local optimum is a global minimum.
- [The PDG summary covariance omits cross-block correlations] → Preserve the documented covariance, expose the zero-cross-covariance approximation, and avoid presenting nominal fit probabilities as fully calibrated.
- [The conditional $G_V$ fit shares information with the joint fit] → Keep its result as a standalone diagnostic and never add its $\chi^2$ to the joint total.
- [An alternative saturation or mass prescription may lack input support] → Evaluate only source-backed alternatives; list unsupported alternatives and the missing inputs needed to test them.

## Migration Plan

Add the diagnostic script and a new result artifact without changing existing fit outputs. Add the memory record and update its containing and parent indexes. Reproduce the result from the documented command and compare the recorded input hashes before interpreting numerical differences. Rollback consists of removing only the newly added diagnostic, result, memory record, and index links.
