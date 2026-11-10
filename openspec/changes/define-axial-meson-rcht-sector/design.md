# Design

## Context

Existing RChT records use an antisymmetric-tensor vector nonet with $\langle T_iT_j\rangle=\delta_{ij}$, a minimal mass action with $\chi_+$ and singlet shift, and ALP sources in the no-ALP–meson-mixing basis of Alda:2026nrd (`memory/rcht/foundations/`). The $a\to4P$ records already use axial *ALP sources* $c_q^A$, but no axial *resonance*. `memory/plan.md` §4 is empty. Odd-parity axial operators are alp-ct8.2. Scripts may only be edited by the `coder` agent; memory records follow the AGENTS.md front-matter and index rules.

## Goals / Non-Goals

**Goals:**
- One convention-complete even-parity axial Lagrangian and verified vertices usable by alp-ct8.3–8.6.

**Non-Goals:**
- Odd-parity (Levi-Civita) operators and anomaly matching (alp-ct8.2).
- Numerical values of coefficients, fits, or widths (alp-ct8.3).
- Loops and NLO corrections.

## Decisions

1. **Antisymmetric tensor $A_{\mu\nu}$ as the primary representation.** It matches the vector nonet and Ecker:1988te, so the shared $\nabla_\mu$, $f_\pm^{\mu\nu}$ building blocks and mass-action technique carry over. Alternative: Proca $A_\mu$ (as in Jamin:2000wn-type and Roig:2013baa-type analyses). Kept as a cross-check: the equivalence map and induced contacts are derived, so results can be compared with literature written in either form.
2. **Nonet with nonideal mixing for the neutral isoscalars.** Same structure as the $\omega$–$\phi$ record: one single-trace $\chi_+$ operator plus a singlet shift, mixing angle derived symbolically. Other double-trace operators are listed as omitted assumptions.
3. **Operators classified by chiral order, using the standard building blocks, with redundancy removal.** Basis built by enumerating $\langle A\,X\rangle$ with $X$ from $u_\mu$, $f_\pm$, $\chi_\pm$ and derivatives, then reducing with Schouten, Bianchi and leading-order equations of motion. Counting is done by a script, mirroring `count_odd_leading_covariant.py`. Alternative (cite a literature basis only) rejected: acceptance requires completeness for our channels.
4. **Vertices from a single symbolic expansion.** One script expands the Lagrangian in $\Phi$, $V$, photon and ALP sources, extracts vertices, diagonalizes quadratic mixing, and reports induced contacts; checks (Ward identity, representation equivalence, ideal-limit) are separate test functions.
5. **Width treatment:** narrow-width stable state only where $\Gamma/M$ justifies it (decide per state in the record); otherwise treat $a_1\to\rho\pi\to3\pi$ as a decay chain. Exchange propagators use a constant width as a baseline, as in the vector records, consistent with `vector-width-lineshape-prescriptions.md`.
6. **Process inventory:** no explicit axial inventory exists yet. Task 1 builds it from the ALP processes already in `memory/alp/` plus new axial-state channels, and lists exclusions.
7. **Delegation:** the `bibliography-manager` obtains and digests the references; the `coder` writes scripts; memory records are written by the main agent/orchestrator.

## Risks / Trade-offs

- [Tensor-versus-Proca differences hide in contacts] → derive equivalence explicitly and test amplitudes.
- [Operator basis larger than the data constrain] → record all, mark which are redundant for ALP channels; reduction belongs to alp-ct8.3.
- [Broad $a_1$ is not an asymptotic state] → decay-chain default; narrow-width is opt-in with stated validity.
- [Reference conventions differ from Alda:2026nrd] → translate and record the dictionary; new notation requires user approval.
- [Strange-axial mixing $K_{1A}$–$K_{1B}$ (nonet of two $1^{++}$/$1^{+-}$ states)] → keep $K_1(1270)$/$K_1(1400)$ mixing out of the quadratic sector unless an operator requires it; note as an open assumption.

## Open Questions

- Whether the $1^{+-}$ nonet ($b_1$, $h_1$, $K_{1B}$) is part of the even-parity sector or deferred; default: only the $1^{++}$ nonet here.
