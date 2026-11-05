# Design

## Context

See `proposal.md` for motivation. The indexed records select a source-specific two-pion dispersive contribution for the broad isovector $\rho$ and narrow $\omega/\phi$ poles for the pure isoscalar resonance terms. They also select a mixed $\rho I$ VVP term, but its overlap with the two-pion reconstruction is unresolved. The $a\gamma^*\to\pi\pi$ partial wave is an input to the two-pion unitarity cut; its retained tree expression has rho exchange in the dipion, crossed-channel, and photon-virtuality variables. The pion vector form factor carries its own rho dynamics. A generic finite-width self-energy derivation exists, but no rho-specific running-width or subtraction convention has been selected. Standalone three-pion and kaon cuts, their complex inputs, and a complete higher-energy continuation remain later extensions. The current vector-sector likelihood is nonconverged and is not a calibrated input.

## Goals / Non-Goals

**Goals:**

- Define the analytic component content and rho line-shape conventions before acquiring or combining numerical inputs.
- Decide whether the mixed $\rho I$ term is consistent in the first approximation, using an explicit overlap analysis with the two-pion discontinuity.
- Assemble selected channel responses with the independent no-mixing ALP sources while preserving source, current, phase, uncertainty, and completeness status.
- Add three-pion, kaon, and higher-energy components only after their required inputs and pole/continuum matching are defined.

**Non-Goals:**

- Claim a complete TFF at unbounded virtuality or a complete prediction for $a\to\ell^+\ell^-$.
- Infer source-resolved ALP responses from $η$ or $η'$ amplitudes, or use the nonconverged vector fit as validated calibration.
- Apply one energy-dependent width indiscriminately to rho propagators in different invariants or on top of rho dynamics already represented in $F_\pi^V$.
- Treat the production amplitude entering a unitarity cut as an independent TFF term.

## Decisions

1. **Use the two-pion spectrum for the broad isovector rho.** Do not also add a stable neutral-rho pole for the same source response. Keep rho exchange inside the $a\gamma^*\to\pi\pi$ production partial wave as part of the cut input; it is not an extra TFF contribution by itself.

2. **Map each rho occurrence before choosing its width.** Distinguish the dipion invariant in the right-hand cut, crossed $t/u$ poles in the production partial wave, and the photon virtuality in $W_\gamma(y)$. Identify the rho dynamics already present in $F_\pi^V$ and its phase. Do not replace all occurrences with one Breit-Wigner denominator or duplicate pion rescattering.

3. **Derive the dressed rho response with explicit subtractions.** For a rho propagator represented by a self-energy, derive the channel-dependent $\Gamma_\rho(s)$ and use $\operatorname{Im}\Delta_\rho^{\mathrm{self}}(s+i0)=+\sqrt{s}\,\Gamma_\rho(s)$; determine the real part through a subtracted dispersion relation with stated mass and residue conditions. The existing formal identity does not fix the numerical width or subtraction constants.

4. **Resolve the mixed-term overlap at the TFF level.** Compare the single and double discontinuities and subtraction terms generated from the two-pion production/current amplitudes with the mixed $F^{(\rho I)}$ term, including both photon orderings. Derive a replacement/subtraction prescription or establish that the selected terms are disjoint before numerical addition. Until this is done, the mixed term is conditional and excluded from a numerical first-pass sum.

5. **Make the first isoscalar approximation explicit.** Use the pure narrow $\omega/\phi$ poles as the baseline isoscalar resonance terms. The indexed memory currently retains the mixed $\rho I$ term; task 2.4 decides whether it can be included consistently or must be deferred. Standalone three-pion and kaon continua are later extensions and do not enter the baseline without compatible amplitudes and $\phi$-region matching.

6. **Require an input manifest before numerical assembly.** For every included channel record thresholds, kinematic variables, current and production conventions, phases, cut normalization, source basis, data or derivation provenance, and uncertainty representation. An unavailable item is a named blocker, not a default value.

7. **Separate uncertainty sources and gate fit-dependent numbers.** Propagate correlated input variations when supplied; label cutoff, subtraction, line-shape, pole/continuum, and omitted-channel changes as model variations. Keep resonance calibration and fit-driven uncertainty propagation conditional on a converged, identifiable `alp-uzh.8` result.

## Risks / Trade-offs

- [A width inserted in several rho-dependent terms can count the same pion rescattering twice] → Map each invariant and the rho content of $F_\pi^V$ before selecting the dressed responses.
- [The $F^{(\rho I)}$ overlap is not determined by the presence of a rho in the production wave alone] → Compare the TFF discontinuities and subtraction terms; keep the mixed term out of numerical sums until the audit is complete.
- [Missing complex amplitudes or relative phases prevent numerical interference] → Keep affected channel terms symbolic or mark them blocked; do not replace phases with magnitudes.
- [The $\phi$ pole overlaps its three-pion and kaon cuts] → Add those continua only with a common response or documented replacement that represents the shared strength once.
- [A low-energy spectrum may be mistaken for a full loop input] → Store and report the virtuality domain and missing higher-energy matching with each result.
- [The vector fit remains unresolved] → Keep calibration tasks conditional and keep the fit out of the accepted input manifest.

## Migration Plan

No existing numerical interface is replaced. Add the source-resolved assembly and its input manifest alongside existing derivation scripts, link new records from `memory/alp/lepton-pair/INDEX.md`, and preserve the current conditional and symbolic results. A rollback can remove the new assembly artifacts without changing existing derivations.
