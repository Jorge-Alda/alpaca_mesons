---
title: Guidelines for optimizing decay-width integrands
status: mixed
created: 2026-10-22
updated: 2026-10-22
bibtex_keys: []
provenance: numerical-computation
---

# Guidelines for optimizing decay-width integrands

## Scope

Use these steps when adding or improving numerical widths in `alpaca_chiral`. They summarize the validated scalar-algebra work for $a\to\gamma PP$, $a\to3P$, and $a\to VPP$. They are a workflow, not an assumption that another channel has the same amplitude identities or bottlenecks. Channel-specific transformations require their own derivation and numerical checks.

## 1. Freeze the reference calculation

Before editing, identify the amplitude convention, invariant definitions, masses, units, propagator widths, phase-space boundaries, symmetry factor, pole subdivisions, quadrature method, and tolerances. Record representative open, closed, threshold-adjacent, smooth, and pole-accessible cases. Preserve an executable reference calculation during development.

Do not change the amplitude model, integration domain, tolerances, or pole prescription in the same step as an algebraic optimization. This keeps discrepancies attributable to one change. Measure the reference full-width time and pointwise integrand time before optimizing.

## 2. Separate preparation from pointwise evaluation

Construct a prepared amplitude object or closure once per width call. Hoist channel-independent work out of the repeated Dalitz evaluation, including:

- flavour traces, source combinations, pole groupings, and cancellation checks;
- squared masses, total invariant sums, fixed numerator coefficients, and normalization factors;
- complex pole locations that combine mass and constant width;
- physical-limit formula coefficients and pole locations needed for integration splits.

At fixed outer invariant, compute outer-invariant contributions and kinematic combinations once, then reuse them across the inner integration. Keep the scalar evaluator small, avoid temporary arrays, dictionaries, matrices, generators, and repeated calls that rebuild closures where scalar arithmetic suffices. Retain coherent sums and interference terms before taking the absolute square.

## 3. Derive before simplifying

Write down the original numerator, denominators, polarization contractions, and on-shell invariant relations. Reduce them algebraically with a reproducible SymPy script when a nontrivial identity is used. Keep finite widths, crossed terms, and matching pieces in the derivation; do not infer a cancellation from a zero-width or narrow-width limit.

Useful validated examples are documented in [three-pseudoscalar-scalar-optimization.md](three-pseudoscalar-scalar-optimization.md) and [vpp-scalar-optimization.md](vpp-scalar-optimization.md). For polarization sums, a factorized Gram determinant must retain its complete dependence on every integration variable. In $a\to VPP$, the exact form is $s(t-t_-)(t_+-t)/4$; replacing it by a function of $s$ alone would remove the angular dependence and is incorrect.

## 4. Preserve phase space and error control

Use the same invariant definitions and physical boundaries as the reference. Check threshold limits, collapsed boundaries, and positivity at physical points. Retain identical-particle factors and channel-specific neutral-state projections. If the integrand has poles inside the physical region, keep the reference subdivisions or prove and test a replacement. Do not weaken tolerances to make a candidate look faster.

Relative accuracy alone is insufficient near destructive-interference zeros. Use a nonzero absolute tolerance in pointwise comparisons and report when a reference width is only a roundoff residual; do not call such a residual a physical nonzero prediction.

## 5. Compare pointwise and integrated results

Make a deterministic test/benchmark script that contains a slow reference evaluator or an otherwise independent check. Cover:

- interior points and points approaching each Dalitz boundary;
- points near every accessible pole band and crossed-pole subdivision;
- equal and unequal masses, repeated-particle symmetry factors, and coherent neutral states where implemented;
- thresholds, forbidden channels, exact zeros, and representative complex couplings;
- full widths across channels with smooth and resonance-enhanced regions.

Compare complex amplitudes before squared amplitudes, then compare polarization sums or $|\mathcal M|^2$, and finally compare integrated widths. Also verify any derived scalar identities symbolically. Keep the original integration tolerances, finite-width prescription, and symmetry factors the same in both width calculations.

## 6. Benchmark total cost and point cost separately

Report both repeated integrand-evaluation speed and complete-width speed. Complete-width timing should include amplitude preparation, because callers pay that cost. Use repeated runs and report the median, channel and mass point, coupling input, machine/runtime context when relevant, and numerical discrepancy alongside speedup. A fast integrand does not guarantee a fast width if adaptive subdivision or preparation dominates.

Do not generalize a measured speedup from one channel or environment. Profile the next channel independently. Larger changes such as fixed grids, compilation, caching, interpolation, or analytic integration need their own convergence and error study.

## 7. Promote only after validation

Initially expose the candidate through a separate optimized function so it can be compared directly with the current API. Once the algebra, pointwise behavior, widths, edge cases, and performance benefit are checked, promote it to the canonical function and remove the duplicate slow production evaluator if requested. Keep a slow, test-only reference in the benchmark script when it remains useful for regression checks; do not retain two production paths indefinitely without a documented need.

Record the derivation, exact inputs, invocation, benchmark cases, tolerances, measured differences, timings, and limitations in a focused `memory/` record, and update the containing `INDEX.md` in the same change.

## Provenance

These guidelines are distilled from the verified implementations and comparisons in [photon-pp-scalar-optimization.md](photon-pp-scalar-optimization.md), [three-pseudoscalar-scalar-optimization.md](three-pseudoscalar-scalar-optimization.md), and [vpp-scalar-optimization.md](vpp-scalar-optimization.md). The channel-independent workflow is recommended practice; each future algebraic reduction remains unverified until it receives its own derivation and tests.
