---
title: External tadpole residue source traces
status: derived
created: 2026-11-26
updated: 2026-11-26
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: analytical-derivation
---

# External tadpole residue source traces

## Inputs and assumptions

Use the incoming momentum convention, pole legs and loop-sizing prescription of the [estimate record](../../../couplings-and-decays/finite-mixing/murayama-radial-up-only-estimate-result.md), following Alda:2026nrd and DiLuzio:2026poh, and the compact vertex rules of the [generator-sum record](murayama-loop-generator-trace-products.md). Let $e=(\lambda_e,\xi_e)$ be an external pseudoscalar pole and $\bar e=(\lambda_{\bar e},\xi_{\bar e})$ its charge conjugate, with incoming momenta $p$ and $-p$, so $k_e\cdot k_{\bar e}=-p^2$. Products with the loop momentum are the loop scale squared and do not depend on $p^2$. Source potentials have no derivatives. Units are GeV.

## Derivation and result

The [script](../../../../../../scripts/derive_radial_loop_external_residues.py) differentiates the full vertices with respect to $p^2$. For the closed ALP line only the four-meson kinetic part depends on $p^2$; the splits of the legs into the pairs $(e,a),(\bar e,a)$ contribute, while the split $(e,\bar e),(a,a)$ vanishes because $[\lambda_a,\lambda_a]=0$. For the closed radial line only the ALP kinetic term contributes. The result is

$$
\frac{\partial V(e,\bar e,a,a)}{\partial p^2}=\frac{\langle[\lambda_a,\lambda_e][\lambda_a,\lambda_{\bar e}]\rangle}{12f_\pi^2},\qquad
\frac{\partial V(e,\bar e,\sigma,\sigma)}{\partial p^2}=\frac{2\xi_e\xi_{\bar e}}{f_a^2}.
$$

Multiply either expression by the tree amplitude and the external residue factor $1/2$; the tadpole symmetry factor $1/2$ and the loop measure are excluded. On the ALP leg, $e=\bar e=a$, the first expression vanishes because both commutators vanish; in the diagonal neutral sector of the up-only model it vanishes for $\pi^0$ as well, since diagonal matrices commute. Charged radial angular components vanish. Additional pole-mixing terms allowed for general sources are not calculated here.

## Reproduction and verification

Run `PYTHONPATH=scripts .venv/bin/python scripts/derive_radial_loop_external_residues.py` and `PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_radial_loop_external_residues`. The test checks exactly that both slopes equal the $p^2$ derivatives of the direct expansion of the Lagrangian of the [generator-sum record](murayama-loop-generator-trace-products.md) for general legs $e$, $\bar e$; that they equal the two closed forms displayed above; and that the first slope vanishes for $\lambda_e=\lambda_{\bar e}=\lambda_a$. No loop integral is evaluated.
