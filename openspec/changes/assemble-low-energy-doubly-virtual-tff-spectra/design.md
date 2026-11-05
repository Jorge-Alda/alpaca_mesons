# Design

## Context

See `proposal.md` for motivation. The indexed records already provide formal source-resolved constructions for two-pion, three-pion, and kaon channels, a narrow-pole $ω/φ$ approximation, and a dressed mixed $ρ I$ VVP term. They do not yet provide all compatible complex production and current amplitudes, the relative three-pion phase, a numerical mixed-term/continuum match, or a complete higher-energy continuation. The current vector-sector likelihood is nonconverged and is not a calibrated input.

## Goals / Non-Goals

**Goals:**

- Define one bookkeeping layer that combines channel responses with the independent no-mixing ALP sources while preserving source, channel, current, phase, and uncertainty provenance.
- Resolve or explicitly gate each shared pole/continuum contribution, especially the dressed $ρ$ overlap and the common $φ$ strength in three-pion and kaon channels.
- Make partial coverage machine-readable so an incomplete low-energy contribution cannot be mistaken for the full TFF needed by the lepton loop.

**Non-Goals:**

- Claim a complete TFF at unbounded virtuality or a complete prediction for $a\to\ell^+\ell^-$.
- Infer source-resolved ALP responses from $η$ or $η'$ amplitudes, or use the nonconverged vector fit as validated calibration.
- Replace the selected two-pion treatment by a standalone stable-$ρ$ approximation.

## Decisions

1. **Keep channel terms explicit until overlap handling is applied.** Represent the low-energy contribution as separately identified source responses for the two-pion double-spectral term, three-pion and kaon cuts, narrow $ω/φ$ poles, and mixed dressed $ρ I$ terms. Combine only compatible responses; do not merge the poles and continua into an opaque sum before their shared strength is resolved. This follows the existing source decomposition and exposes overlap that a single effective pole fit would hide.

2. **Require an input manifest before numerical assembly.** For every channel record thresholds, kinematic variables, current and production conventions, complex phases, cut normalization, source basis, data or derivation provenance, and uncertainty representation. An unavailable item is a named blocker, not a default value. This is preferred to filling gaps with stable-pole residues or independent marginal errors that are not phase- or covariance-compatible.

3. **Treat pole/continuum matching as channel-level replacement or subtraction.** For the isovector sector retain the two-pion dispersive contribution and match the mixed dressed $ρ I$ term so the common $ρ$ strength appears once. For the $ω/φ$ terms, retain narrow poles only where no included continuum already represents that strength; in the $φ$ region coordinate three-pion and kaon cuts through a shared response or a documented replacement. Do not select a numerical matching function until compatible complex inputs and subtraction conventions are available.

4. **Separate uncertainty sources.** Propagate correlated experimental/input variations together when supplied; keep cutoff, subtraction, pole/continuum matching, and omitted-channel variations labeled as model uncertainty. If no defensible bound exists for higher isoscalar strength, report the omission and limit the stated domain instead of inventing an error estimate.

5. **Gate fit-driven numbers.** Structural derivations and input acquisition proceed independently. Resonance calibration and fit-driven error propagation wait for `alp-uzh.8` and its convergence/identifiability criteria; until then, outputs remain conditional and carry no fitted coefficient interval.

## Risks / Trade-offs

- [Missing complex amplitudes or relative phases prevent numerical interference] → Keep affected channel terms symbolic or mark them blocked; do not replace phases with magnitudes.
- [The $ρ$ or $φ$ overlap prescription changes the predicted spectrum] → Record the matching convention and its variation as a model uncertainty once the required amplitudes exist.
- [A low-energy spectrum may be mistaken for a full loop input] → Store and report the virtuality domain and missing higher-energy matching with each result.
- [The vector fit remains unresolved] → Keep calibration tasks conditional and keep the fit out of the accepted input manifest.

## Migration Plan

No existing numerical interface is replaced. Add the source-resolved assembly and its input manifest alongside existing derivation scripts, link new records from `memory/alp/lepton-pair/INDEX.md`, and preserve the current conditional and symbolic results. A rollback can remove the new assembly artifacts without changing existing derivations.

## Open Questions

No user decision is required before implementation. The input audit must determine which complex three-pion and charged/neutral kaon datasets, subtraction conditions, and higher-isoscalar bounds are supported. Anything that cannot be established from compatible inputs remains an explicit blocker or bounded omission, as required by the spec.
