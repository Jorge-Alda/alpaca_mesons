---
title: Chiral-scale dependence of the mixed inverse propagator
status: derived
created: 2026-11-10
updated: 2026-11-10
bibtex_keys: [Alda:2026nrd, Kaiser:1998ds, Gasser:1983yg]
provenance: analytical-derivation
---

# Chiral-scale dependence of the mixed inverse propagator

## Inputs and assumptions

Use the leading mixed block and field transformation in [mixed ultraviolet field renormalization](mixed-uv-field-renormalization.md), the local operators in [quadratic Lagrangian](../quadratic-lagrangian.md), the tadpole contraction in [mixed ALP–meson self-energy](mixed-alp-meson-self-energy.md), and the running in [subtraction scheme and scales](../counting/subtraction-scheme-and-scales.md). The conventions are those of Alda:2026nrd and Kaiser:1998ds, with $Q=sK-M^2+\Sigma$, masses and decay constants in GeV and $s$ in GeV$^2$. The source matrices are arbitrary and need not commute. Retain one-loop terms and first order in $1/f_a$, with the loop content and exclusions of [loop classes](../counting/loop-classes.md).

## Derivation

1. Extract each $a\Phi$ operator from the independent quadratic-Lagrangian derivation. For each of the nine generators, its coefficient of $\partial_\mu a\partial^\mu P_i$ multiplies $s$ and its coefficient of $aP_i$ contributes directly. The script prints all nine local entries with independent $L_i^r$, $\Lambda_1$ and $\Lambda_2$.
2. Add the leading mixed block and the finite tadpole contraction of the linked mixed self-energy record. Evaluate the running of the leading term by differentiating with respect to $B_0$, and use $\mu\,dB_0^r/d\mu=-[2m_0^2/(3f_\pi^2)]B_0^r/(16\pi^2)$. Running inside one-loop and NLO terms is beyond the retained order.
3. The tadpole changes by $-2M^2\log(\mu/\mu_0)/(16\pi^2)$ in the generator basis. The local constants change according to their $\Gamma_i$ values. The script evaluates these variations separately, without assuming commuting sources.
4. For the field transformation matrix $R$ in [mixed ultraviolet field renormalization](mixed-uv-field-renormalization.md), the resulting mixed block at scale $\mu$ differs from that at $\mu_0$ by $R^Tb^{(0)}\log(\mu/\mu_0)/(16\pi^2)$. To compare in the same fields, apply the finite transformation $P(\mu)=[1-R\log(\mu/\mu_0)/(16\pi^2)]P(\mu_0)$. It contributes the negative of this difference. This statement is at first order in the loop expansion.

## Verification and scope

Run `.venv/bin/python scripts/derive_mixed_running_block.py`. The unit test checks the cancellation symbolically in all nine entries for arbitrary source components and quark masses. It uses the existing mixed tadpole contraction, whose finite coefficients have not all been independently checked against the vertex expansion; the scale cancellation alone does not establish that independent comparison. The result is conditional on the linked operator basis and subtraction.
