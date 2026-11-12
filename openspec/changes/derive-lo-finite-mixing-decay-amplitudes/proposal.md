# Proposal

## Why

Beads issue alp-sjp.3 (GitHub issue 20) asks for ALP decay and production amplitudes built on the physical mixed poles of the combined neutral meson-plus-ALP system, rather than on an unmixed external ALP. Its prerequisites are in place: alp-sjp.1 (chiral redefinitions, closed) and the LO part of alp-sjp.2 (archived change `solve-lo-finite-mixing-spectrum`, capability `lo-finite-mixing-spectrum`), which returns every physical pole $m_n^2$, its $K$-orthonormal vector $v^{(n)}$ and the residues of $D^{-1}(s)$ from the full $6\times6$ neutral system. What exists so far is channel-by-channel and incomplete: $a\to3\pi$ ([finite-mixing three-pion record](../../../memory/alp/three-pseudoscalar/alp-finite-mixing-three-pion-amplitudes.md)), one two-meson vertex ([charged $U(3)$ record](../../../memory/alp/shared/alp-charged-u3-vector-redefinition.md)) and the baryon one-ALP vertex ([baryon record](../../../memory/alp/baryon/finite-mixing-baryon-couplings.md)). There is no process inventory, no common amplitude construction that also projects the outgoing mesons onto their physical poles, and no statement of which channels open or close at finite mixing. alp-sjp.5 (implementation) depends on this. The work is purely analytical and does not depend on the numerical inputs of alp-sjp.4.

## What Changes

Scope is **leading order** only, as for the spectrum: the $O(p^2)$ even-parity chiral Lagrangian with the full ALP-dependent sources $M(x)$, $C_{R,L}(x)$, $\Theta(x)$ (no linear-in-$a/f_a$ assumption; the linear case of Alda:2026nrd is a special case), the leading anomalous (WZW) sector for odd-intrinsic-parity channels, and the LO baryon CCWZ Lagrangian. Tree level.

- **Process inventory**: list every channel with a physical-pole ALP parent at LO: $\gamma\gamma$; one, two and three pseudoscalars; $\gamma+$ pseudoscalars; baryon-antibaryon; and the channels that need resonances or loops (vector mesons, $\ell^+\ell^-$, resonance final states). For each, give the LO status (derived here, already derived elsewhere, or deferred with the reason), the quantum numbers and exact selection rules, the threshold condition in terms of physical masses $m_n$, and whether the channel exists at all in the unmixed limit.
- **Amplitude construction**: tree amplitude with one decaying physical pole and every outgoing meson leg also projected onto its own physical pole by $v^{(n)}$. The vertex tensors come from the joint expansion of the LO Lagrangian in $(P_i,a)$. Exchange diagrams use the exact $D^{-1}(s)$, not only the pole part. Residues give the LSZ factors.
- **Complete transformed vertices**: derive the same vertices in the pole-aligned coordinates, including every nonlinear contact term of the exact exponential chiral map, and show that the two routes agree.
- **Generic amplitudes**: pseudoscalar channels are derived for generic pseudoscalar indices, for example $a\to P_iP_jP_k$ and $a\to P_iP_j$ with $P_i$ the generators of the $U(3)$ nonet, not channel by channel ($a\to\eta\pi^+\pi^-$, $a\to\eta K^+K^-$, and so on). The amplitude is a function of the flavour indices, the pole vectors and the momenta. Specific channels are obtained afterwards by assigning the indices, and the existing specific records ($a\to3\pi$ and others) are reproduced only as checks.
- **Consistency checks** per channel: Ward identities (photon channels), crossing, Bose symmetry and symmetry factors, invariance under constant field redefinitions, weak-mixing limit against Alda:2026nrd, and agreement with the three existing finite-mixing records.
- **Channel-appearance rules**: exact statements of which decays are allowed or forbidden by mixing (for example, CP-odd versus CP-even content of the poles; $K^0$ pole with no ALP component), and which channels exist only at finite mixing.
- **Stable versus broad states**: LO treats all final mesons as stable poles; the inventory records which final states are broad resonances and states the rule for building rates on physical daughters, without implementing widths.
- Memory records (one concept each, front matter and indexes per AGENTS.md) and scripts written by the `coder` agent.
- **Out of scope** (to be tracked in a new follow-up bead): $a\to4P$ (excluded from this change, existing no-mixing records are unchanged); vector-meson exchange and $VPP$, $V\gamma$-type resonance channels; $a\to\ell^+\ell^-$ (loop-induced, needs the transition form factor); NLO local terms, loop self-energies and resonance mixing in the poles (alp-sjp.7); complex poles and widths; numerical evaluation; phase-space integration and performance (alp-sjp.5). alp-sjp.3 is therefore only partly resolved and stays open.

## Capabilities

### New Capabilities
- `lo-finite-mixing-process-inventory`: the LO process list with status, quantum numbers, exact selection rules, thresholds and channel-appearance rules.
- `lo-finite-mixing-decay-amplitudes`: the LO tree amplitudes with physical-pole projections on all legs, complete transformed vertices and the symmetry checks.

### Modified Capabilities
<!-- None. `lo-finite-mixing-spectrum` is consumed as a given and its requirements do not change. -->

## Impact

- New records under `memory/alp/shared/` (construction, inventory, selection rules) and in the existing channel folders (`two-pseudoscalar`, `three-pseudoscalar`, `photon`, `photon-two-pseudoscalar`, `baryon`); every affected `INDEX.md` updated.
- New scripts under `scripts/` (coder agent), reusing `solve_alp_meson_spectrum.py`, `derive_alp_ccwz_arbitrary_order.py`, `derive_alp_ccwz_redefinitions.py`, `derive_alp_three_pion_vector_redefinition.py`, `derive_alp_two_meson_vector_redefinition.py` and the photon and baryon derivation scripts.
- Beads: note on alp-sjp.3 (LO part done, remainder listed), a follow-up bead for the beyond-LO channels linked to alp-sjp.3. Do not close alp-sjp.3.
- Provides the amplitude definitions used by alp-sjp.5.
