# Design

## Context

See [proposal.md](./proposal.md) for motivation. The joint-fit implementation builds a bin-averaged Dalitz grid with separate energy, pair, and angle orders and pole-aware maps around the ϕ and K* structures. Its current optional audit varies pair and angle order together, uses one fit point, and omits residual and Jacobian comparisons. The joint residual Jacobian combines analytic odd-coefficient derivatives with finite differences for other parameters.

## Goals / Non-Goals

**Goals:**
- Make refinement independent in all three integration dimensions and compare complete likelihood outputs.
- Sample reproducible points that cover central, alternative, and less constrained parts of the fitted parameter region.
- Persist pass/fail evidence and enough input provenance to repeat it.

**Non-Goals:**
- Changing the physics model, experimental inputs, likelihood definition, or optimizer.
- Proving that an optimizer found a global minimum or that the model is adequate.
- Treating a quadrature pass as a fit-quality or parameter-uncertainty certification.

## Decisions

1. **Refine one integration dimension at a time, then verify the combined order.** Start from the fit's recorded orders. Increase energy, pair, and angle orders independently while holding the other two fixed, then compare the combined refined result. This exposes which integral controls a discrepancy; a single diagonal sequence can hide cancellation. Keep the physical Dalitz limits and existing ϕ/K* pole mappings unchanged.

2. **Use uncertainty-scaled output tolerances.** Require the largest per-bin prediction change divided by that bin's experimental uncertainty to be at most $0.01$, and the largest change in standardized residual to be at most $0.01$. Require each nonzero scaled-parameter Jacobian column to change by at most $10^{-3}$ in relative Euclidean norm (and report absolute changes for zero or numerically negligible columns). Require the total χ² change to be at most $0.01$. These thresholds are deliberately smaller than the experimental resolution and are recorded in every result so later reviews can revise them explicitly.

3. **Keep finite-difference settings fixed across refinements.** Evaluate the residual Jacobian with the same parameter scaling and finite-difference steps at every quadrature level; otherwise derivative changes would mix integration error with differentiation error. Record the derivative method and steps. Include the analytic odd-coefficient columns and the numerical remaining columns.

4. **Define the point set from saved fit information, without refitting.** Include the reference candidate, every available competing minimum, and deterministic valid displacements along representative well-constrained and weakly constrained directions inferred from the scaled fit Jacobian/covariance. Record each complete parameter vector and why it was selected. Do not call a point a minimum unless its source fit says so. If required covariance information is unavailable, use documented bounded parameter displacements and mark that basis in the output.

5. **Store a machine-readable certification beside the fit data.** Record software revision, hashes of all input data, point vectors and sources, orders at each comparison, pole-map description, tolerances, prediction/residual/Jacobian/χ² changes, per-point outcomes, and exact command. An invalid point or a maximum-order tolerance failure remains visible as unevaluated or uncertified; neither is silently dropped.

## Risks / Trade-offs

- [Risk] Full Jacobians at many refinement levels are computationally expensive → Use a small documented point set that includes both narrow-pole-sensitive and weakly constrained directions, and cache parameter-independent integration geometry where safe.
- [Risk] A finite point set cannot prove uniform convergence everywhere in the fitted region → Report the sampled region and point-selection rule; describe the result as certification at those points only.
- [Risk] Tight tolerances can fail at practical maximum orders → Preserve the failure data and allow higher orders to be requested; do not relax tolerances implicitly.
- [Risk] A local covariance can be unreliable near bounds or flat directions → Include competing minima and bounded displacements, and label the source of each point rather than treating covariance directions as exact confidence intervals.

## Migration Plan

No migration is required. Existing fit outputs remain readable. New certification output is additive and records its own settings; the ordinary likelihood and default fit behavior remain unchanged.
