# Proposal

## Why

Beads issue alp-sjp.1 (GitHub issue 22) asks for the all-orders ALP-dependent chiral redefinitions needed when $f_a$ is close to $f_\pi$ and the ALP cannot be treated as an unmixed state. The record [alp-general-chiral-redefinition-invariants.md](../../../memory/alp/shared/alp-general-chiral-redefinition-invariants.md) already gives the partonic transformation laws, the transport matrices $W_{R,L}$ and the dressed mass function $\widehat M(x)$. The record [finite-alp-meson-mixing.md](../../../memory/alp/shared/finite-alp-meson-mixing.md) aligns the fields with one pole, using a linear shift only. Two things are missing: (i) the hadronic transformation of the sources, anomaly and nonlinear vertices for general $g_{R,L}(x)$, and (ii) a statement of the combined kinetic-plus-mass problem as a coordinate change of one inverse propagator. This is a prerequisite for alp-sjp.2, .3, .5 and .6.

## What Changes

- Derive, from the restricted partonic Lagrangian, how $g_{R,L}(x)$ act on the full CCWZ Lagrangian: transformed $r_\mu,\ell_\mu$, $\chi$, the anomalous $\Theta(x)$ shift and the Wess-Zumino term, with the transport-matrix form $\widehat U=W_R^\dagger U W_L$.
- Derive the transformation of the full quadratic form $sK-M^2$ in the field space $(P_i,a)$ under a general invertible field redefinition, showing how the kinetic and mass mixing blocks reshuffle while poles, eigenvectors, residues and amplitudes are invariant.
- Derive the pole-aligned coordinates for one selected eigenstate, and the conditions under which one transformation aligns several poles.
- Derive the nonlinear terms of the coordinate map needed for multi-meson vertices, and the transformation of amputated vertices.
- Add executable SymPy verification and a memory record for each derivation, as required by the repository protocol.
- No numerical spectrum or decay width is produced here; that belongs to alp-sjp.2, .3 and .5.

## Capabilities

### New Capabilities
- `alp-chiral-redefinition-covariance`: all-orders transformation of the ALP-dependent sources, mass spurion, anomaly and nonlinear meson vertices under $g_{R,L}(x)$, within the restricted partonic Lagrangian.
- `finite-mixing-pole-alignment`: transformation of the combined inverse propagator under field redefinitions, invariance of the physical poles, residues and amplitudes, and pole-aligned coordinate choices.

### Modified Capabilities
<!-- None: openspec/specs/ is empty. -->

## Impact

- New memory records under `memory/alp/shared/` (with `INDEX.md` updated), written following the protocol in AGENTS.md.
- New scripts under `scripts/` (written by the `coder` agent), reusing `derive_alp_invariant_coefficient_functions.py`, `derive_alp_finite_mixing.py`, `derive_alp_lo_finite_mixing.py` and `derive_alp_ccwz_arbitrary_order.py`.
- The "Verification status" and pole-alignment sections of `finite-alp-meson-mixing.md` may need cross-references; their content is not altered.
- Unblocks alp-sjp.2, alp-sjp.3, alp-sjp.5 and alp-sjp.6.
