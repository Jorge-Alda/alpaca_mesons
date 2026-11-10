# Proposal

## Why

Beads issue alp-gii.2 (GitHub #30) asks for the one-loop two-point functions of the $U(3)$ theory. alp-gii.1 fixed the calculated order, the included loop classes and the subtraction scheme, and it is closed. The existing record [loop-corrections-to-alp-meson-lsz-mixing.md](../../../memory/nlo/loop-corrections-to-alp-meson-lsz-mixing.md) gives only the block structure and the LSZ formula, with $\Sigma^{\rm loop}$ left symbolic. alp-gii.4 (poles and LSZ residues), alp-gii.5 (numerical inputs) and alp-gii.8 (decay vertices) are blocked until $\Sigma_{ij}(s)$, $\Sigma_{i\hat a}(s)$ and their derivatives exist as explicit, renormalized expressions.

## What Changes

- Derive the one-loop meson self-energy matrix $\Sigma_{ij}(s)$ and the mixed ALP–meson function $\Sigma_{i\hat a}(s)$ from the classes G-LO and G-ALP of alp-gii.1: Goldstone-only loops (including $\eta'$) with LO $U(3)$ vertices and one LO ALP-source vertex at first order in $1/f_a$. Resonance loops are not calculated (excluded by alp-gii.1); they enter only the truncation uncertainty.
- Enumerate which cubic and quartic vertices survive for the retained sources. Decide from this which topologies (tadpole, bubble) contribute, and hence which absorptive thresholds exist. The project has no CP-violating sources and no constant $\theta$ term, so there is no $\Phi^3$ vertex. Cubic $a\Phi\Phi$ vertices exist for flavour-violating vector and scalar ALP couplings, and their effect on each two-point function is stated explicitly.
- Give the UV pole, the finite local part and the non-analytic part of each function in the $\overline{\rm MS}-1$ scheme of alp-gii.1, with the flavour structure for the full $U(3)$ nonet in the LO mass eigenbasis (including $\pi^0$–$\eta$–$\eta'$ mixing).
- Show how the poles and $\mu_\chi$ dependence are removed using the local terms built from $L_i^r(\mu)$ and $B^r(\mu)$, together with the divergent meson field redefinition used for the meson self-energy. Apply that redefinition to the mixed block as well. Compare the combinations of $\Gamma_i$ constrained by the calculation with those of alp-gii.1.
- Give the wave-function derivatives $\partial_s\Sigma$ and the pole-residue contributions, and the perturbative NLO pole shift. State the prescription for energies near a threshold or a resonance pole, and label any resummation as outside the strict order.
- Add an executable script and step-by-step memory records, and update the memory indexes.

## Capabilities

### New Capabilities

- `nlo-u3-loop-self-energies`: the renormalized one-loop meson and mixed ALP–meson two-point functions of the $U(3)$ NLO theory, their divergence cancellation, threshold structure, derivatives, pole-residue terms and the NLO/resummation prescription.

### Modified Capabilities

None. `openspec/specs/` is empty; the capability of the unarchived change `fix-u3-nlo-resonance-counting` is only referenced, not modified.

## Impact

- New records under `memory/nlo/loops/` (new topic folder, with `INDEX.md`); `memory/nlo/INDEX.md` and `memory/README.md` updated. [loop-corrections-to-alp-meson-lsz-mixing.md](../../../memory/nlo/loop-corrections-to-alp-meson-lsz-mixing.md) gets a link to the explicit results; its content is not changed.
- New scripts under `scripts/`, written by the `coder` agent. The existing `derive_loop_corrected_alp_lsz.py` stays symbolic and is reused for the final pole-residue formula.
- `bibliography.bib`: possible additions through `bibliography-manager` if a further INSPIRE record is consulted (for example a $U(3)$ one-loop reference for $\eta'$ tadpoles).
- Unblocks alp-gii.4, alp-gii.5 and alp-gii.8. No changes to `latex/`.
