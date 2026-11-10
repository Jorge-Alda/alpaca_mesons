---
title: Perturbative and resummed pole prescription
status: mixed
created: 2026-11-10
updated: 2026-11-10
bibtex_keys: [Alda:2026nrd, Kaiser:1998ds]
provenance: assumption
---

# Perturbative and resummed pole prescription

## Inputs and assumptions

The strict first-order pole and LSZ formulas are given in [derivatives and pole residues](derivatives-and-pole-residues.md), following [the LSZ loop record](../loop-corrections-to-alp-meson-lsz-mixing.md) [@Alda:2026nrd]. The retained Goldstone loops and excluded resonance loops are specified in [loop classes](../counting/loop-classes.md) [@Kaiser:1998ds]. The criterion below is a working prescription, not a probability or an estimate of the full error.

## Comparison and prescription

For the declared truncation estimate ([validity range](../counting/validity-range-and-truncation.md)),

$$
\delta_{\rm est}(s)=\max\!\left({1\over N_c},\ {s\over(4\pi f_\pi)^2},\ {M_K^2\over(4\pi f_\pi)^2}\right),
$$

compare the retained loop self-energy with the LO inverse propagator entry by entry in the chosen field basis:

$$
|\Sigma_{ij}(s)|\leq\delta_{\rm est}(s)\,|Q^{(0)}_{ij}(s)|.
$$

Equality passes. If any entry violates the inequality, use the loop-corrected propagator $[Q^{(0)}(s)+\Sigma(s)]^{-1}$, with the local corrections appropriate to the chosen calculation also included in $Q(s)$, and solve $\det Q(s)=0$ for the poles. For the ALP meson components, the corresponding block expression is $w^{(a)}=-A^{-1}b$. Report these results together with the strict NLO pole and LSZ components. Inverting the loop-corrected matrix retains higher powers of the one-loop correction and is therefore a resummation outside strict NLO.

This entrywise test is basis dependent. A zero LO entry passes only if its loop entry also vanishes; at an LO pole a vanishing diagonal entry can trigger the test even for a small correction. The comparison is consequently a diagnostic rather than a basis-independent measure of convergence. The first-order expansion itself assumes the relevant LO inverse exists, as stated in the [pole derivation](pole-script-derivation.md).

The retained tadpoles have no two-meson cut in $s$ ([thresholds](thresholds-and-analytic-structure.md)). If a pole approaches a threshold or a resonance region associated with excluded effects, this resummation does not supply those effects. In particular, resonance loops remain excluded and a tree-level resonance pole needs the separate width prescription specified in [validity range](../counting/validity-range-and-truncation.md). The estimate $\delta_{\rm est}$ is the factor multiplying the NLO correction in the truncation estimate; the separately estimated excluded resonance-loop uncertainty and its assumptions are given in that record. Neither resummation nor the criterion extends the declared validity range.

## Verification

The illustrative two-meson test in [the pole script derivation](pole-script-derivation.md) gives $0.23375\,\mathrm{GeV}^2$ at first order and $0.23413447\,\mathrm{GeV}^2$ after resummation for the lower pole. These are test values, not fitted results. Run `.venv/bin/python scripts/derive_nlo_loop_self_energies.py --poles`; the focused tests linked in that record check agreement through first order, quadratic scaling of the difference, and equality and zero-entry cases.
