# Design

## Context

See `proposal.md`. The adequacy criterion is fixed in `openspec/changes/assess-joint-vector-fit-adequacy/design.md`: nominal total $p>0.01$ and no block with nominal $p<10^{-3}$. The baseline assumes one shared $K^*$ width coordinate with $K\pi$ saturation, a negative single-trace source, no double-trace source, and the isospin-symmetric even sector (see `babar-kkpi-global-vector-sector-fit-baseline.md`). Evidence on responsible rows is in `joint-vector-fit-block-exclusion-results.md`; several of those refits were conditional or did not converge.

## Goals / Non-Goals

**Goals:**
- Separate the candidate causes with diagnostics that need no new operator.
- Add only those extensions that the diagnostics or the issue's priority list justify, each switchable and exactly reducing to the baseline.
- Reach the adequacy criterion, or document why it cannot be reached.

**Non-Goals:**
- Loop corrections, new data sets, or error inflation without justification.
- Profile intervals and identifiability (alp-uzh.7).
- Changing the adequacy criterion to obtain a pass.

## Decisions

- **Two stages, with a gate.** Stage A runs diagnostics without new operators. Stage B adds extensions in the issue's priority order. After each extension the criterion is re-evaluated; later extensions are run only if it is not met. Extensions 1–2 are tried regardless of the diagnostics' outcome, since the issue ranks them first and they are the cheapest.
- **Diagnostics are scripts.** Each diagnostic is a deterministic script from the `coder` agent, using the accepted endpoint `basin_probes/even_shared_widths` and quadrature $(4,128,128)$ for reporting, $(2,16,16)$ for refits, as in alp-uzh.6.
- **Extensions as likelihood options.** Each extension is an option with a default that reproduces baseline residuals bitwise; a regression test enforces this. Added coordinates get explicit bounds and starting values recorded before running.
- **Analytic input first.** For extensions that add operator contributions (double-trace source, SU(3)-breaking VPP, $\rho$–$\omega$ mixing), the amplitude is derived with a SymPy script and recorded in `memory/` before it enters the likelihood. Whether the double-trace source enters the vector mass spectrum is unconfirmed (candidate-causes item 3) and is settled first by that derivation; if it does not, the extension is recorded as not applicable.
- **Isospin-breaking inputs.** Charged and neutral $K^*$ masses and widths are taken from the PDG 2024 values already in the likelihood rows; no new data are needed. Whether the PDG neutral/charged width difference is itself significant is checked in Stage A.
- **Acceptance of refits.** Same gates as the baseline (optimizer status, certified quadrature, multi-start from the accepted endpoints, local-error extraction). Failed or non-qualified refits are kept in the ledger.
- **Parameter-count honesty.** Degrees of freedom are recounted for every extension; $p$-values are nominal, as in alp-uzh.6, and an extension that adds a coordinate is not preferred on $\chi^2$ alone but on the criterion and on $\Delta\chi^2$ per added coordinate.
- **Notation.** Any new coordinate (e.g. the $K^*K\pi/\rho\pi\pi$ ratio, $y$) needs a descriptive name that follows Alda:2026nrd; the user approves it before it enters the memory record.

## Risks / Trade-offs

- [The $K^*$ width tension is an experimental inconsistency, not a model failure] → The irreducible-$\chi^2$ diagnostic bounds what any model with a common width can achieve; if it already exceeds the block threshold, the result is documented as a data-level limit.
- [Several extensions together over-fit 43 dof] → Extensions are added one at a time, and the final record reports every combination tried.
- [Double-trace source requires a new derivation] → It is scheduled after the cheaper extensions and gated on the criterion not yet being met.
- [Mass-block removal refits were non-convergent before] → Diagnostics report them as conditional, as in alp-uzh.6.

## Migration Plan

Add options, drivers and ledgers; existing outputs are untouched. Rollback removes the new options and files.

## Open Questions

- Whether the user wants every extension run even after the criterion is met (default here: stop at the first adequate model).
