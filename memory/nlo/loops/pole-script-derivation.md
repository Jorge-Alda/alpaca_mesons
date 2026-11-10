---
title: Derivation implemented for one-loop pole corrections
status: derived
created: 2026-11-10
updated: 2026-11-10
bibtex_keys: [Alda:2026nrd, Gasser:1983yg, Kaiser:1998ds]
provenance: analytical-derivation
---

# One-loop pole calculation in the script

## Inputs and assumptions

The input self-energies are the finite tadpoles of the [meson](meson-self-energy.md) and [mixed](mixed-alp-meson-self-energy.md) calculations, with the subtraction of [loop functions](loop-functions-and-subtraction.md) [@Gasser:1983yg; @Kaiser:1998ds]. The inverse propagator is $Q=sK-M^2+\Sigma$ and the LSZ block conventions and first-order expansion are those of [the LSZ loop record](../loop-corrections-to-alp-meson-lsz-mixing.md) [@Alda:2026nrd]. Units are GeV for masses and GeV$^2$ for $s$ and $\Sigma$.

Assume nondegenerate LO meson poles, an invertible LO meson block at the ALP pole, and first order in the loop terms and ALP mixing. The ALP self-energy and the mixing-induced ALP pole shift are of order $1/f_a^2$ and are excluded. An independently supplied ALP pole shift may still be included through the pole-shift term in the LSZ record. The [truncation estimate](../counting/validity-range-and-truncation.md) is a working estimate, not a probability.

## Derivation

1. The script differentiates the meson tadpole expression with respect to $s$. The mixed axial contribution is proportional to $s$, so its derivative is the axial trace coefficient printed for each nonet field. The pseudoscalar mixed contribution is independent of $s$; the vector, scalar and gluonic mixed tadpoles vanish under the stated assumptions. Both derivatives are independent of $s$.
2. For the ALP meson components at fixed $s$, the first-order expansion of the block equation yields

   $$
   w_0^{(a)}=-A_0^{-1}b_0,\qquad
   \delta w^{(a)}=-A_0^{-1}\left[\Sigma_{i\hat a}+\Sigma_{ij}w_{0j}^{(a)}\right],
   $$

   where the indices in the second expression denote the mixed column and the meson matrix. The terms from local NLO corrections and from an independent pole shift are given in the linked LSZ record and are printed by the script.
3. For a normalized LO meson eigenvector and nondegenerate LO pole, implicit differentiation of $\det Q=0$ gives the first-order shift as minus the expectation value of $\Sigma$ at that pole. The inverse residue is one plus the expectation value of $\partial_s\Sigma$ there; the LO kinetic matrix is the identity. The script verifies the shift by differentiating an explicit two-meson determinant.
4. Decision 8 compares every entry through $|\Sigma_{ij}(s)|\leq\delta_{\rm est}(s)|Q^{(0)}_{ij}(s)|$. Equality passes. A zero LO entry passes only if its loop entry also vanishes. The comparison depends on field basis and is a heuristic; at an LO pole it can be triggered simply by a vanishing diagonal entry. If it fails, the script's example solves $\det[Q^{(0)}+\Sigma]=0$ and uses the inverse of the same loop-corrected meson block. This resummation is outside strict order and does not include excluded resonance loops.

## Verification

Run `.venv/bin/python scripts/derive_nlo_loop_self_energies.py --poles`. The illustrative two-meson example has a lower LO pole at $0.25\,\mathrm{GeV}^2$ and a first-order shift of $-13/800\,\mathrm{GeV}^2$. Its lower first-order and resummed poles are $0.23375$ and $0.23413447\,\mathrm{GeV}^2$, respectively. The script derives the first-order overlap from the block equation and checks that its exact inverse agrees through first order and reduces to the tree-level overlap when loops vanish. The unit tests also check the determinant through first order, the quadratic scaling of the difference between the resummed and first-order poles, the criterion at equality and for zero entries, and the truncation estimate at a stated point. The example is illustrative and is not a fit.

```bash
PYTHONPATH=scripts .venv/bin/python -m unittest derive_nlo_loop_self_energies.LoopSelfEnergyTests.test_pole_expansion_and_tree_limit derive_nlo_loop_self_energies.LoopSelfEnergyTests.test_resummation_criterion_zero_entries_and_boundary
```
