# Proposal

## Why

The axial-vector sector of the ALP programme (epic alp-ct8) has no Lagrangian yet: `memory/plan.md` §4 is marked future work, and the existing RChT records cover only the vector nonet and the odd-parity operators. Every later axial task (alp-ct8.3 coefficient mapping, alp-ct8.4 decays with axial final states, alp-ct8.5 axial exchange, alp-ct8.6 code) needs one agreed field representation, flavor basis, set of even-parity operators and set of vertices. Issue alp-ct8.1 provides exactly that.

## What Changes

- Fix the axial-meson field representation (antisymmetric tensor $A_{\mu\nu}$, consistent with the $V_{\mu\nu}$ convention of Alda:2026nrd and Ecker:1988te), the nonet flavor basis, phase conventions and nonideal mixing (the $f_1(1285)$–$f_1(1420)$ pair) in memory records.
- Write the kinetic and mass terms, including flavor and singlet breaking, in the same form as the vector mass action.
- Construct the complete set of even-parity operators with one axial field and at most the chiral order needed for the $A P$, $A V$, $A P P$, $A V P$, $A a$ (ALP source) and electromagnetic-source vertices; separate coefficients fixed by symmetry from independent ones.
- Derive, with executable `sympy` scripts, the $A$–$P$, $A$–$V$, $A$–Goldstone and external-source vertices, including the field redefinitions and induced local terms needed so that the tensor and Proca representations agree.
- State the power counting and the width treatment (stable versus finite-width external states, line shape in exchange).
- Record every independent coefficient and every induced contact term in a table consumed by alp-ct8.3–8.6.
- Odd-parity axial operators and anomaly matching are **not** included (alp-ct8.2).

## Capabilities

### New Capabilities
- `axial-meson-rcht-lagrangian`: the even-parity axial nonet field content, conventions, quadratic terms, nonideal mixing, complete operator basis, independent-versus-fixed coefficients, power counting and width treatment.
- `axial-meson-rcht-vertices`: the executable derivation of the $A$–$P$, $A$–$V$, $A$–Goldstone and source vertices, representation-equivalence contacts, and the coefficient/contact table for the ALP processes.

### Modified Capabilities

None (`openspec/specs/` is empty).

## Impact

- New records under `memory/rcht/axial-even-parity/` (new topic folder) with `INDEX.md`, plus updates to `memory/rcht/INDEX.md`, `memory/plan.md` §4 and `bibliography.bib`.
- New derivation scripts under `scripts/` (written by the `coder` agent only).
- Literature extracts under `memory/extracts/<bibtex-key>/` for the axial-RChT references.
- Unblocks alp-ct8.3, 8.4, 8.5 and 8.6. No change to existing vector or odd-parity results.
