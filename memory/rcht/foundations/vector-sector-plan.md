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
- No unique precision physical-angle determination or full odd-parity
  split-pole matching result is established.

## Next executable task

Stage 4: integrate the declared mass and source operators into the complete
radiative and two-photon amplitudes. Retain the separately calibrated mass
scenarios, mass branches and relative source-sign hypotheses. Verify canonical
external normalization, every vector-slot rotation, Ward identities and the
baseline limit before repeating coupling-rank diagnostics. Additional amplitude
data and an explicit interaction-breaking truncation must select any physical
prescription; neither an exact residue interpolation nor a small quoted
conditional angle error is sufficient.
