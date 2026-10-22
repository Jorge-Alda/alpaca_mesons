---
title: Strategies for faster phase-space integration
status: conjectured
created: 2026-10-22
updated: 2026-10-22
bibtex_keys: []
provenance: conjecture
---

# Strategies for Faster Phase-Space Integration

## Scope and current bottleneck

This note collects options for improving width-evaluation speed while preserving numerical accuracy. It is a discussion record, not a selected implementation plan. The current three-body resonance calculators use nested adaptive quadrature over two Dalitz variables. In practice, the cost can include both repeated adaptive-integrator overhead and work inside each integrand evaluation. No profiling or accuracy benchmark has yet established which cost dominates for each channel.

The leading-order chiral-only implementation has analytic reductions for repeated-particle channels because its amplitude depends on a single pair invariant. That simplification does not generally extend to amplitudes with vector exchange: Bose symmetry alone does not remove the second invariant. Any reduction must be justified for the full amplitude being integrated.

## Candidate improvements

### Fixed quadrature on mapped phase space

Map the outer invariant to a fixed interval and the inner invariant to $z\in[-1,1]$ using the Dalitz boundaries. Evaluate tensor-product Gauss--Legendre nodes in batches. For smooth integrands, this may remove substantial nested-adaptive overhead. Increase quadrature order until the integrated width stabilizes, and compare against independently converged adaptive results.

Resonances and threshold boundaries can spoil convergence of a global fixed grid. Partition the domain around accessible pole bands or choose coordinates that resolve the resonant structure; retain an adaptive fallback for difficult points. Existing pole subdivision in a given calculator should be accounted for before claiming this as a new benefit.

### Cheaper integrand evaluation

Prepare flavor traces, pole groupings, masses, mass-square combinations, and invariant-independent coefficients once per width call. Profile amplitude construction separately from quadrature. The inspected vector--pseudoscalar--pseudoscalar implementation performs small NumPy allocations and a determinant in its pointwise polarization contraction; the photon--pseudoscalar--pseudoscalar path recomputes quantities such as Dalitz limits inside the inner evaluation. These are candidates for algebraically equivalent scalar expressions, hoisted constants, batched evaluation, or compilation. Their benefit should be measured before expanding the implementation complexity.

### Analytic inner integration

For the present constant-width propagator prescriptions, investigate whether the complete coherent inner integral can be performed analytically at fixed outer invariant. The polarization contractions are polynomial in invariants, while amplitudes contain rational pole terms. A derivation must retain crossed exchanges and interference and handle coincident poles, narrow widths, collapsed intervals near threshold, and cancellation among primitive terms. Keep numerical integration as a fallback if stable formulas are unavailable. Running widths or added non-rational amplitude components could change this conclusion.

### Reuse across coupling scans

At fixed kinematics and fixed hadronic inputs, amplitudes linear in source couplings imply widths quadratic in their independent real components. One can integrate basis interference terms once and evaluate later coupling points through a quadratic form. This is attractive for scans over couplings at fixed $m_a$, but its gain is less clear when most scan points vary $m_a$.

For mass scans, possible compromises include a bounded cache of exact-mass results or integrated coefficient matrices, and interpolation of those matrices on an adaptively refined $m_a$ grid. Any cache needs a clear key over all kinematic and hadronic inputs, a memory limit, and a policy for nearby but unequal floating-point masses. Interpolation needs explicit threshold handling, additional resolution when resonance bands enter phase space, and validation at withheld masses. Its error adds to the quadrature error. Whether reuse pays off depends on scan density, channel count, memory budget, and the cost of evaluating each point; this remains open for specific scan workloads.

## Error control and benchmark design

The existing nested routines use relative tolerances, but returned inner and outer error estimates and evaluation counts should be recorded before changing tolerances. Near zeros caused by destructive interference, a relative error alone is not an adequate criterion; use an absolute accuracy floor as well. Tightened adaptive calculations and a second numerical method can serve as references for candidate methods.

A useful benchmark set includes the reported $a\to\omega\pi^+\pi^-$ and $a\to\gamma K_LK_S$ examples, smooth regions, masses near threshold, masses where a resonance band first becomes accessible, and coupling choices with destructive interference. Report wall time, amplitude-evaluation count, convergence failures, and width discrepancies. Do not infer numerical error from the size of omitted physical effects, and do not use a narrow-width approximation without assessing its threshold and interference errors.

## Four-body outlook

For a generic scalar $1\to4$ decay, removing overall orientation leaves five independent kinematic variables. A three-variable integral is possible only if the amplitude and phase-space structure permit two further integrations to be carried out or eliminated. Establish the variable count for the intended channel before choosing an integrator. If several dimensions remain, benchmark deterministic cubature against randomized quasi-Monte Carlo with resonance-aware phase-space mappings and independent repetitions; the preferred method is an empirical choice.

## Inputs and provenance

This discussion is based on the repository's numerical width implementations and the derivations linked in this directory's index. Proposed speedups, cache designs, and interpolation schemes are conjectures until profiled and compared with reference integrations. No code changes or numerical benchmarks are included in this record.
