---
title: Plan for a nonideal vector sector from an explicit action
status: mixed
created: 2026-10-13
updated: 2026-10-13
bibtex_keys: [Ecker:1988te, Ecker:1989yg, Moussallam:2021dpk, Kucukarslan:2006wk, Alda:2026nrd]
provenance: assumption
---

# Vector-sector execution plan

## Objective and source policy

Replace the ideal omega–phi assignment and independently inserted pole masses
with generators, masses and source residues derived from a declared quadratic
antisymmetric-tensor action. Preserve the repository's positive-strange phi
phase, generator normalization and electromagnetic-source convention.
Use the arXiv version identified by the Alda BibTeX entry for subsequent neutral
pseudoscalar inputs; the local manuscript is not an authoritative substitute.
This policy supersedes the local-manuscript input selection in the previous
[next-steps record](../odd-parity/next-steps.md).

## Branch preference

User instruction on 2026-10-13: retain both mass-compatible branches in the
notes and calculations. When all other adequacy, source-sign and truncation
considerations are equivalent, prefer the branch closer to ideal mixing.
This is a tie-breaker, not evidence that the near-ideal branch is measured,
and does not pool charged/neutral scenarios or select a residue-sign hypothesis.
The stage-4 implementation exposes both branches and uses the near-ideal one
as its explicit default selection.

## Indexed inputs

- [RChT baseline](rcht-basics.md) — tensor kinetic and leading source couplings.
- [Photon diagonalization](rcht-photon-vector-diagonalization.md) — generic
  neutral generators and representation-induced local terms.
- [Vector-mixing source](../../extracts/Moussallam:2021dpk/vector-mixing.md) —
  explicit flavor/nonet breaking example and limitations of mass-only calibration.
- [Radiative action](../odd-parity/radiative-sm-amplitude.md) and
  [TFF action](../odd-parity/radiative-tff-amplitude.md) — later integration targets.

## Stages and completion conditions

1. **Derive the quadratic tensor action.** Retain canonical kinetic terms,
   a leading quark-spurion mass operator and a singlet mass term as a declared
   minimal extension. Derive its mass matrix, physical generators, mixing,
   photon overlaps and matrix response in an executable SymPy script.
   Check the ideal/common-mass limits, orthogonality and basis invariance.
   Identify assumptions and any representation matching still absent.
2. **Audit the tensor literature and operator truncation.** Compare a direct
   tensor mixing treatment with the chosen extension. Separate mass operators
   from symmetry-breaking source/interaction operators; do not imply that
   quadratic matching establishes full odd-parity matching.
3. **Calibrate masses and electromagnetic residues.** Establish a fixed data
   snapshot and likelihood using masses and leptonic widths. Test whether the
   minimal action can describe the inputs. If not, quantify the residual and
   introduce only explicitly justified additional operators. Keep model
   uncertainty separate from input errors. No decay-fitted angle is an
   independent prior when the same decays enter the likelihood.
4. **Propagate into radiative and two-photon amplitudes.** Rotate every vector
   slot and source consistently; use the action's matrix propagator and
   canonical external-state normalization. Recover the existing baseline,
   recompute rank and state which symmetry-breaking vertices remain omitted.
   Physical eta inputs require the arXiv-verified Alda generators and spurion.
5. **Joint decay calibration, then assess dynamical mixing.** Include radiative
   channels and the separately planned dispersive three-pion treatment with
   documented correlations. Add momentum-dependent self-energies only after
   deciding the precision target and subtraction prescription.

## Execution log

- Plan saved before completing the first derivation.
- Stage 1: [tensor mass-action derivation](vector-nonideal-mass-action.md)
  completed and its executable script passed. It establishes the quadratic
  mass matrix, physical generators, photon overlaps and common-pole basis
  invariance within the explicitly assumed minimal action.
- Stage 2: [direct tensor source audit](../../extracts/Kucukarslan:2006wk/tensor-mixing.md)
  saved. Strong versus electromagnetic mixing and off-shell field definitions
  are separated; photon exchange must not be counted twice.
- Alda arXiv v1 source verified and archived: recorded central inputs agree,
  and local-manuscript provenance is superseded.
- Stage 3: [mass/source calibration](vector-stage3-calibration.md) completed
  using checksum-pinned PDG 2024 summaries. The minimal mass action and its
  universal source are inadequate at quoted errors. An explicitly derived
  additional mass invariant gives two branches; source/sign sensitivity and
  the loss of independent mixing information with three residue coefficients
  are recorded. Conditional input errors remain separate from model limits.
- Stage 4: [action amplitude integration](../odd-parity/vector-stage4-amplitudes.md)
  and [conditional rank diagnostics](../odd-parity/vector-stage4-rank.md)
  completed. Both mass branches and all source hypotheses remain selectable.
  ArXiv Alda generators are reconstructed, and the new kernels connect to the
  existing TFF/omega observations with their likelihood policy unchanged.
  Baseline, photon Ward, Bose, phase/basis, crossing and data-regression checks
  pass. All 64 retained diagnostic scenarios give rank 8 with inherited OPE
  rows; those rows have not been rematched for the enlarged action.
- No unique precision physical-angle determination or full odd-parity
  split-pole matching result is established.

## Next executable task

The user requested an ALP-target audit before advancing this plan.
[The audit and completion plans](../odd-parity/alp-target-amplitude-audit.md)
are now saved. The checked leading projections do not access the current
three blind directions, but normalized physical resonance-corrected ALP
predictors remain incomplete. Resolve the selected target-completion scope
before advancing the joint-fit stage below; this is a task-order change,
not an additional physical fit constraint.

The three blind odd-coefficient directions are explicitly **not a blocker**
for the desired tree-level ALP accuracy: the checked target contributions
do not depend on them. They could enter at higher loops or orders, so this
decision must be revisited if the truncation or target accuracy is extended.
Missing normalized predictors and the other matching/data requirements
remain separate from this blind-direction decision.

Stage 5 begins with the joint-fit specification: state the interaction-breaking
truncation, retained branch/sign hypotheses and source-parameter treatment;
declare the covariance/independence assumptions for the selected radiative,
TFF and leptonic data. The inherited OPE targets cannot become an enlarged-action
precision prior without renewed short-distance matching; the single-virtual
high-energy behavior remains a separate blocker. Three-pion inputs require the
already scoped dispersive treatment. Preserve both mass branches, preferring
the nearer ideal solution only when other adequacy/truncation criteria are equal.
No physical fit conclusion follows from the stage-4 rank alone.
