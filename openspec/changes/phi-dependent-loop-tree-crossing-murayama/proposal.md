# Proposal

## Why

The implicit crossing equation for the loop-to-tree ratio of an ALP-sourced three-meson amplitude is derived for the ALP alone (`estimate-loop-tree-crossing-finite-mixing`). Its exact tree denominator is valid for any $\epsilon=f_\pi/f_a$, but its loop numerator is truncated at the leading mixing-dependent order. In the Murayama model (and its DiLuzio generalization) the radial mode $\sigma$ couples to the pseudoscalars at the same order in $1/f_a$ as the ALP, so it enters the tree exchange and the one-loop classes. The crossing must therefore be restated with the radial mode included, with the loop numerator covering all classes and not truncated, in the regime of large $\kappa_u$ (small $f_a$, $\epsilon$ of order $1$ to $30$) that the model is intended to explore.

## What Changes

- Derive the general implicit crossing equation in $\varphi$-valued invariants: the dressed sources are functions of the complex radial-phase field $\varphi$, and the ALP and radial mode enter through its polar components. This is more general than the $n_q$ parametrization of DiLuzio:2026poh, which becomes a special case.
- Replace the truncated loop numerator by a defensible loop-size estimate covering all classes (meson-only loops and loops with the radial mode), with each class's $\epsilon$ power, reference scale and invariant contraction stated. **No loop contribution is computed**; only the estimate is proposed.
- Keep the tree denominator exact, with the $\sigma$ exchange included.
- State the naturalness assumption of the ALP-only change for the neglected next-order term, including its renormalization-scale-dependent and counterterm parts, without analysing what counterterms absorb.
- Keep $R_{\max}$ and the fixed independent Dalitz invariants of the existing derivation.
- Specialize the general result to the up-only model with $f_a=\sqrt2\,m_u/\kappa_u$, keeping $\kappa_u$ and $\lambda$ as analytical parameters.
- Reuse the ALP-only results (exact pole solver, dressed-frame invariants, loop power-counting and naturalness records, implicit crossing form) wherever their hypotheses still hold, and state explicitly where they do not.
- Caveats, stated but not computed: the near-degeneracy of the ALP-like state with $\eta'$, and the validity of leading-order chiral perturbation theory at $m_a\approx 1$–$2.3$ GeV.

## Capabilities

### New Capabilities

- `radial-mode-loop-tree-crossing`: implicit, analytical loop-to-tree crossing for the ALP-sourced three-meson amplitude in the Murayama model with the radial mode included, in $\varphi$-valued invariants, and its up-only specialization.

### Modified Capabilities

None. The ALP-only crossing requirements are in change `estimate-loop-tree-crossing-finite-mixing` (capability `finite-mixing-loop-reliability`, not yet in the main specs) and are reused by reference, not modified.

## Impact

- Records and derivations under `memory/alp/shared/` and `memory/alp/three-pseudoscalar/` (new radial-mode vertex, one-loop amplitude and crossing records).
- Symbolic scripts under `scripts/` (written by the coder agent only) for the $\varphi$-valued mass relation, the one-loop vertices and the amplitude with counterterms, and the up-only specialization.
- Bibliography entries for the scalar one-loop functions only (`bibliography.bib`).
- No change to existing code or to the main specs.
