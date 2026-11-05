# Design

## Context

See `proposal.md` for motivation and `specs/joint-vector-fit-minima/spec.md` for the observable contract. The fit driver already uses scaled even-sector coordinates, analytic and channel-specific Jacobian paths, shared propagator evaluation, cached channel terms, and energy-bin blocks. The optimized end-to-end run reduced wall time by about a factor of 2.5 but reproduced the prior endpoint and still terminated unsuccessfully at `max_nfev=1500`. The four-start campaign also reached that cap for every start. Its lowest-objective valid endpoint is `stored_negative_c4`, with $\chi^2=226.194810$ and unsuccessful status. The high-order quadrature record certifies $(2,64,64)\to(4,128,128)$ for five historical candidate points, not for a converged minimum. The even-sector audit fixes the shared $G_V$ definitions and physical branch assumptions.

## Goals / Non-Goals

**Goals:**

- Find and retain candidate endpoints from deterministic, physically permitted initial points for the unchanged 32-coordinate, 75-residual fit.
- Make each attempt reproducible, including scaling, bounds, branch, quadrature, optimizer settings, residuals, input/source hashes, and termination diagnostics.
- Carry the existing runtime improvements through default and high-order calculations so multiple starts can be run without weakening numerical precision.
- Establish a specific acceptance check combining optimizer success, domain validity, residual finiteness, and high-order quadrature stability.
- Continue the best valid endpoint with a larger evaluation cap before adding new starts, and retain enough progress information to diagnose a second cap hit.
- If continuation converges, check stationarity under the certified quadrature and probe a small, deterministic set of additional basins.

**Non-Goals:**

- Alter the measured summaries, likelihood, covariance assumptions, physical branch, bounds, or parameter meanings to make convergence easier.
- Claim that multi-start sampling proves a global minimum or that local quadratic errors are calibrated confidence intervals.
- Accept a candidate merely because it improves the objective or reaches a small step while the optimizer reports failure.
- Infer that no minimum exists solely because the planned continuation remains unsuccessful.

## Decisions

- **Preserve the current likelihood and numerical path.** Reuse the existing propagator calculation reuse, cached channel evaluations, channel-specific Jacobian finite differences, and bounded energy-bin blocks. Do not introduce a second evaluator or lower quadrature orders as a runtime shortcut. Compare fixed-candidate residuals, predictions, and Jacobians against the existing route within the established tolerances.
- **Use named deterministic starts and retain the complete attempt ledger.** Include the existing stored candidate and starts that vary the allowed fitted coordinates while remaining inside the same physical branch and bounds. Persist start identifiers and coordinates, full final vectors, residual vectors or equivalent block diagnostics, objective and block contributions, optimizer status/message, `nfev`/`njev`, optimality/gradient, active bounds, and timing/software/input hashes. Preserve attempts that fail or converge to distinct endpoints. Deduplicate only for presentation using a declared coordinate/objective tolerance; never delete the underlying attempts.
- **Keep scaling explicit.** Retain the script's existing `EVEN_SCALES` for even-sector parameters, record the scale for every coordinate and the optimizer's `x_scale`/transformation convention, and set documented characteristic scales for other coordinates from their declared units/bounds or local Jacobian. Compare alternatives only through a recorded, deterministic setup. Do not silently change scale conventions between starts.
- **Separate optimizer termination from physical validity.** A candidate must satisfy the existing `unpack` branch/domain checks, bounds, and finite 75-component residual vector. Record both SciPy's success and termination message as returned; independently report any near-bound coordinates and physical checks. A low objective cannot override a failed status.
- **Continue the best endpoint before extending the start set.** Choose `stored_negative_c4`, the lowest-objective valid endpoint in the saved four-start ledger, and use its exact 32-coordinate vector as the initial point for a new optimizer call. Preserve the original ledger, use the same likelihood, bounds, fixed scale vector, and $(2,16,16)$ quadrature with four-bin blocks, and set `max_nfev=5000`. Record the source ledger and seed name. Treat the new call's evaluation count as per-run and report any cumulative count separately.
- **Capture continuation progress.** Write periodic optimizer callback checkpoints containing the current coordinates, objective, and evaluation count; include optimality when the optimizer exposes it. Keep the final complete endpoint diagnostics and termination fields in a separate output ledger. This supports diagnosis if the 5000-evaluation cap is reached again.
- **Check a successful continuation at higher quadrature.** Starting from a successful lower-order endpoint, refit sequentially at $(2,64,64)$ and $(4,128,128)$ with four-bin blocks and the same bounds and coordinate scales. Require successful termination and domain validity at each stage. At the final candidate's fixed coordinates, compare predictions, residuals, scaled-Jacobian columns, and $\chi^2$ across the certified orders; do not accept a candidate that fails any gate.
- **Probe other basins only after a stable continuation.** Define three deterministic, named, physically permitted starts by perturbing separate coordinate blocks of the converged continuation point using the recorded scale vector. Validate each vector before fitting. Run them at the established default quadrature, preserve every result, and refine any successful candidate that competes with or improves on the continuation result.
- **Report nonconvergence without overclaiming.** If the 5000-evaluation continuation or a required higher-order fit remains unsuccessful, use saved objective/optimality checkpoints, coordinate movement, active bounds, and final residual blocks to describe the failure. Record that no minimum was accepted and that the run did not establish whether a minimum exists. Changes to the likelihood or a different optimizer require a separate reviewed plan.
- **Refine only promising candidates with the certified ladder.** Re-evaluate at $(2,64,64)$ and $(4,128,128)$ using four-bin blocks. Compare predictions in experimental standard deviations, standardized residuals, scaled-Jacobian columns, and total $\chi^2$ against the existing tolerances (0.01, 0.01, $10^{-3}$ relative for non-negligible columns, 0.01). If any comparison fails, do not label the candidate accepted; refine the quadrature or report the failure. This extends the check to each candidate rather than presuming the five historical points cover new minima.
- **Store a separate result artifact.** Keep historical baseline and optimization outputs untouched. Add machine-readable start configuration, per-attempt diagnostics, and accepted-candidate summary under `data/vector_sector/`, with a focused memory record linked from the odd-parity index. Record the exact command and environment. This keeps the existing fit output schema and consumers stable.

## Risks / Trade-offs

- [The likelihood may have several basins or flat directions] → Retain distinct endpoints and state that the finite set of starts does not establish global optimality.
- [High-order Jacobians multiply runtime] → Reuse the validated bounded-block and channel-specific paths, run the high-order check only for successful promising candidates, and report timing by stage.
- [Tiny endpoint differences can be numerically amplified by ill-conditioned directions] → Record scaled coordinate differences, objective/residual differences, singular information, and explicit deduplication tolerances; retain raw attempts regardless.
- [A SciPy success flag may coexist with physical or numerical problems] → Apply independent domain and finite-residual checks and report all gates separately.
- [Existing quadrature certification covers only five old points] → Evaluate refinement at each candidate endpoint and preserve a failed result without promoting it.
- [The continuation may reach 5000 evaluations without success] → Retain progress checkpoints and report the observed termination and diagnostics as no accepted minimum; do not infer nonexistence.
- [Higher-order refitting may move to a different basin] → Keep the lower-order continuation and every higher-order endpoint as separate ledger entries; compare them with the declared scaled-coordinate and objective criteria.
- [Three additional starts cannot establish global optimality] → Use them only as a bounded basin check and state that the finite search does not prove a global minimum.

## Migration Plan

Add support for restarting from a saved endpoint and retaining progress checkpoints, then write continuation, higher-order, and additional-start results to new ledgers and an indexed memory record. Do not overwrite existing baseline, optimization, or four-start artifacts. Reproduce the documented commands and confirm input/source hashes before interpreting results. If continuation and higher-order fitting remain unsuccessful, preserve the failure diagnosis and leave the minimum unaccepted. Rollback removes only the newly added continuation workflow and its outputs; existing fit records remain available.
