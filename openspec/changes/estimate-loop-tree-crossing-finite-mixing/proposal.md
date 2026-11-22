# Proposal

## Why

The finite-mixing treatment of the ALP and the pseudoscalar mesons (`lo-finite-mixing-spectrum`) is exact in $\epsilon=f_\pi/f_a$ at tree level. It is therefore not clear for which $f_\pi/f_a$ it stops being a good description once one-loop corrections are included. ALP lines carry $1/f_a$, so loops with ALP lines are not suppressed by $1/f_\pi$ alone, and for $f_\pi/f_a\gg1$ they may exceed tree-level terms. Earlier attempts to answer this gave no defensible result:

- The tadpole-to-tree ratio of an off-shell mass-matrix entry is not invariant under chiral redefinitions (it depends on $c_q^A$ and $c_q^P$ separately).
- The pole-mass shift compared with $m_n^2$ is not a reliability criterion, and compared with the tree-level mixing shift it is $\epsilon$-independent at leading order, so it cannot locate any value of $f_\pi/f_a$.
- Tree-only conditions (large mixing angle, partial-wave unitarity) do not compare loops with tree and are excluded here.

This change derives an analytical, frame-invariant estimate of the value of $f_\pi/f_a$ at which the one-loop correction to an ALP-sourced quantity becomes as large as its tree-level value, after choosing and justifying the observable.

## What Changes

Scope is the LO $U(3)$ Lagrangian of Alda:2026nrd in the dressed frame (the ALP enters only through $\widehat M(x)$ and $\widehat\Theta(x)$ with $x=a/f_a$), tree level against one loop, at finite mixing. The deliverable is closed formulas in the redefinition invariants $\widehat S_n$, $\widehat P_n$, $\widehat\Theta_n$, $\mathcal C=B_0\widehat P_1-\frac{m_0^2}{N_f}\widehat\Theta_1\mathbf 1$ and $\mathbb M_{m_a^2}^{-1}$, with symbolic proofs. Scripts only verify. No scan over parameters, no assumed values or representatives for the ALP coefficients, and no numerical inputs except as a clearly labelled check.

- **Choice of observable, with arguments.** A written comparison of candidate observables (meson pole mass, ALP-meson mixing entry, ALP pole residue, on-shell ALP-sourced amplitude) against explicit criteria: frame invariance, whether the ALP enters at tree level at first order in $1/f_a$ (so that loop and tree have different powers of $\epsilon$ and the ratio can depend on $\epsilon$), insensitivity to the subtraction scheme, and independence of any fitted number. The selected observable and the reasons are recorded.
- **Tree against loop, exact in $\epsilon$.** The tree quantity is the exact finite-mixing pole expression (not its $\epsilon$ expansion). The loop quantity is the one-loop correction split by loop function and by the number of ALP lines. The ratio is derived analytically, with the $\epsilon$ dependence of the exact pole masses and vectors kept.
- **Crossing estimate.** The closed formula for the $f_\pi/f_a$ at which the loop-to-tree ratio of the chosen observable reaches order one, written in the invariants, with the stated order-one criterion and its dependence on the subtraction scale made explicit.
- **Counterterm handling.** Which loop contributions are absorbed by the local NLO operators and which are not, so that the estimate uses only the non-absorbable part.
- **Validity statement**: the conditions (non-degenerate poles, leading order in the loop expansion, CP-conserving or not) under which the formula holds, and what replaces it near a degeneracy.
- Memory records (one concept each, front matter and indexes per AGENTS.md); scripts written by the `coder` agent; draft text only through the `drafter` agent if the user asks.
- **Out of scope**: tree-level-only bounds (mixing-angle and partial-wave unitarity conditions), numerical scans, fits, NLO low-energy constants and resonances, the ALP pole itself in the near-degenerate regime, and loops with photons or weak sources.

## Capabilities

### New Capabilities
- `finite-mixing-loop-reliability`: invariant analytical criterion for the value of $f_\pi/f_a$ at which one-loop corrections to an ALP-sourced observable become comparable to its exact tree-level value in the finite-mixing LO theory.

### Modified Capabilities
<!-- None. `lo-finite-mixing-spectrum` is consumed unchanged. -->

## Impact

- New records under `memory/alp/shared/` (observable choice, tree and loop ratio, crossing formula, validity), with `INDEX.md` there and in parents updated.
- New scripts under `scripts/` (coder agent), reusing `derive_alp_invariant_loop_crossing.py`, `derive_alp_loop_power_counting.py`, `derive_alp_dressed_frame_unitarity.py` and `solve_alp_meson_spectrum.py`.
- Supersedes the representative-based crossing numbers in `alp-loop-power-counting.md` and `alp-loop-crossing-physical-inputs.md`; those records get a pointer, not a deletion.
- Beads: link to the existing loop issues under alp-gii; no commit unless asked.
