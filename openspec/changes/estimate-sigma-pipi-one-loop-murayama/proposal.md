# Proposal

## Why

The radial mode $\sigma$ of the Murayama model decays at tree level to two pseudoscalars through $\sigma P_iP_j$ couplings fixed by the source coefficients $\widehat S_{r,s}$ and $\widehat P_{r,s}$ (`murayama-sigma-two-pseudoscalar-decays.md`, `radial-decay-polar-cartesian-general-sources.md`). The same sources also generate loops in which the radial mode, the ALP-like state and the mesons run, and for the small decay constants $f_a$ of interest the couplings $\epsilon=f_\pi/f_a$ are not small. The change `phi-dependent-loop-tree-crossing-murayama` already sizes these loops for $a\to\pi^+\pi^-\pi^0$ in the up-only model. There is no size estimate for the one-loop correction to $\sigma\to\pi^+\pi^-$, and none for generic sources.

## What Changes

- Estimate the one-loop contribution to $\sigma\to\pi^+\pi^-$ in the generalized Murayama model, written for generic source coefficients $\widehat S_{r,s}$ and $\widehat P_{r,s}$ (no restriction to $r\le1$ or to the up-only sources).
- Keep exact everything that depends on the invariants: the flavour traces of $\widehat S_{r,s}$, $\widehat P_{r,s}$, the external pion matrices and the pole vectors, with their powers of $\epsilon$ and of the pole normalization $N$.
- Estimate everything else, namely the loop integrals, with the rule already used for $a\to\pi^+\pi^-\pi^0$: each scalar moment with $r$ denominators is replaced by $(m^2)^{2-r}/(16\pi^2)$, with $m$ the scale of the loop (`alp-loop-momentum-reduction.md`, `murayama-radial-up-only-estimate-result.md`).
- **No loop integral, scalar loop function, loop coefficient or loop amplitude is computed, in any form or frame.** Only the graphs, their vertices, symmetry factors, traces and the size rule are written down.
- Exclude only the meson-only chiral background: the term with the bare $U(3)$ generator sum on every internal line and a single $\sigma$ end, which contains no letter of the ALP sector and whose ratio to the tree depends on $\epsilon$ only through $p_\sigma^2=m_\sigma^2$ (D4). Its scale is stated as a limit.
- Reuse the vertex rules, the generator-sum treatment of internal meson lines and the external-leg residues of the $a\to\pi^+\pi^-\pi^0$ estimate. Carry out a check of the result under a constant axial redefinition of the sources.
- Write all flavour traces in compact matrix form (traces, commutators, anticommutators, sums over permutations or cyclic permutations), with no component notation or $E_{ij}$, and make every test an exact symbolic identity (no random or numerical matrices).
- Compare with the exact tree $\sigma\to\pi^+\pi^-$ coupling through the ratio $R=|M_1^{\rm est}|/|M_0|$, with $M_0$ the exact tree amplitude. No numerical scan is part of this change, and thresholds are stated as limits.

## Capabilities

### New Capabilities

- `radial-mode-pipi-loop-estimate`: analytical size estimate of the one-loop $\sigma\to\pi^+\pi^-$ amplitude in generic sources, with the loop integrals replaced by the scalar-moment size and the source traces kept exact.

### Modified Capabilities

None. The $a\to\pi^+\pi^-\pi^0$ requirements belong to `phi-dependent-loop-tree-crossing-murayama` and are reused, not modified.

## Impact

- New records under `memory/alp/shared/` (graph list, vertices, traces, estimate, limits) with indexes updated.
- New symbolic and bookkeeping scripts under `scripts/`, written by the coder agent only.
- No change to the draft in `latex/` (the drafter agent is not involved), to existing scripts, or to the main specs.
- No new bibliography entries expected; sources are Alda:2026nrd and DiLuzio:2026poh.
