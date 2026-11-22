# Proposal

## Why

The finite-mixing treatment of the ALP and the pseudoscalar mesons (`lo-finite-mixing-spectrum`) is exact in $\epsilon=f_\pi/f_a$ at tree level. It is therefore not clear for which $f_\pi/f_a$ it stops being a good description once one-loop corrections are included. ALP lines carry $1/f_a$, so loops with ALP lines are not suppressed by $1/f_\pi$ alone, and for $f_\pi/f_a\gg1$ they may exceed tree-level terms. Earlier attempts to answer this gave no defensible result:

- The tadpole-to-tree ratio of an off-shell mass-matrix entry is not invariant under chiral redefinitions (it depends on $c_q^A$ and $c_q^P$ separately).
- The pole-mass shift compared with $m_n^2$ is not a reliability criterion, and compared with the tree-level mixing shift it is $\epsilon$-independent at leading order, so it cannot locate any value of $f_\pi/f_a$.
- Tree-only conditions (large mixing angle, partial-wave unitarity) do not compare loops with tree and are excluded here.
- Computing a one-loop process (for example the bubbles of $A\to3\pi^0$) defeats the purpose: the goal is an analytical estimate of where the tree-level approximation fails, not a loop calculation.

This change derives an analytical, frame-invariant power-counting estimate of the value of $f_\pi/f_a$ at which the one-loop correction to an ALP-sourced amplitude becomes as large as its tree-level value.

## What Changes

Scope is the LO $U(3)$ Lagrangian of Alda:2026nrd in the dressed frame (the ALP enters only through $\widehat M(x)$ and $\widehat\Theta(x)$ with $x=a/f_a$), tree level against the size of the one-loop term, at finite mixing. The deliverable is closed formulas in the redefinition invariants $\widehat S_n$, $\widehat P_n$, $\widehat\Theta_n$, $\mathcal C=B_0\widehat P_1-\frac{m_0^2}{N_f}\widehat\Theta_1\mathbf 1$ and $\mathbb M_{m_a^2}^{-1}$, in matrix form. Scripts only verify the vertex scalings and identities used. No scan over parameters, no assumed values, representatives or relations among the ALP coefficients or the quark masses, and no numerical inputs except as a clearly labelled check.

- **Observable class, with arguments.** The comparison of candidates already delivered (task 1.1) stands. The class used for the counting is an on-shell amplitude in which the ALP enters at first order in $1/f_a$ at tree level, so that ALP-line loops carry extra powers of $\epsilon$ relative to the tree; the reasons are recorded.
- **Vertex power counting in the dressed frame.** For each loop class, the factor of $\epsilon$ per closed ALP line, the loop factor $m^2/(16\pi^2f_\pi^2)$, and the dimensionless matrix contraction of the invariants that multiplies them, written without flavour indices.
- **Tree quantity exact in $\epsilon$.** The tree amplitude uses the exact finite-mixing pole solution, not its expansion.
- **Crossing estimate.** The closed formula for the $f_\pi/f_a$ at which the ratio of the leading ALP-line class to the tree reaches an order-one criterion $R_{\max}$, in the invariants.
- **Naturalness assumption.** The neglected next-order term is taken to have the size of the loop, i.e. no accidental cancellation; the cancellations that this excludes are listed.
- **Validity statement**: non-degenerate poles, one-loop truncation, behaviour near an ALP-meson degeneracy, and the limit $f_a\to\infty$.
- Memory records (one concept each, front matter and indexes per AGENTS.md); scripts written by the `coder` agent; draft text only through the `drafter` agent if the user asks.
- **Out of scope**: tree-only conditions (mixing-angle and partial-wave unitarity), numerical scans, fits, the one-loop evaluation of any process or loop function, the analysis of what counterterms and $L_i$ absorb, assumptions such as $m_u=m_d$ or $\langle\lambda_3\widehat P_1\rangle=0$, NLO constants and resonances, the ALP pole in the near-degenerate regime, and loops with photons or weak sources.

## Capabilities

### New Capabilities
- `finite-mixing-loop-reliability`: invariant analytical power-counting estimate of the value of $f_\pi/f_a$ at which the one-loop size of an ALP-sourced amplitude becomes comparable to its exact tree-level value in the finite-mixing LO theory.

### Modified Capabilities
<!-- None. `lo-finite-mixing-spectrum` is consumed unchanged. -->

## Impact

- New records under `memory/alp/shared/` (observable class, power counting, crossing formula, validity), with `INDEX.md` there and in parents updated. The task 1.1 and 1.2 records already written stay; the Dalitz and absorption analysis of the 1.2 record is not used by the estimate and gets a pointer saying so.
- New scripts under `scripts/` (coder agent), reusing `derive_alp_invariant_loop_crossing.py`, `derive_alp_loop_power_counting.py`, `derive_alp_dressed_frame_unitarity.py` and `solve_alp_meson_spectrum.py`.
- Supersedes the representative-based crossing numbers in `alp-loop-power-counting.md` and `alp-loop-crossing-physical-inputs.md`; those records get a pointer, not a deletion.
- Beads: link to the existing loop issues under alp-gii; no commit unless asked.
