---
title: Counterterm and field-renormalization cancellation in one-loop two-point functions
status: derived
created: 2026-11-10
updated: 2026-11-10
bibtex_keys: [Alda:2026nrd, Kaiser:1998ds, Gasser:1983yg, Gasser:1984gg]
provenance: analytical-derivation
---

# Counterterm cancellation

## Inputs and assumptions

Use the subtraction and $B$ renormalization of [subtraction scheme and scales](../counting/subtraction-scheme-and-scales.md), the $U(3)$ coefficients of [counterterm basis](../counting/counterterm-basis.md), the quadratic operators of [quadratic Lagrangian](../quadratic-lagrangian.md), and the loop coefficients of [meson self-energy](meson-self-energy.md) and [mixed tadpole pole check](mixed-tadpole-pole-check.md). The field and source conventions, counting, units and omitted interactions are those stated in [mixed ALP–meson self-energy](mixed-alp-meson-self-energy.md).

## Cancellation and coefficient comparison

Replace $L_i$ by $L_i^r(\mu)+\Gamma_i\lambda(\mu)$ and $B$ by $B^r(\mu)[1+4\tau\lambda(\mu)/f_\pi^4]$. The source normalizations of Alda:2026nrd and Kaiser:1998ds give $\tau=f_\pi^2m_0^2/6$. Then apply the divergent meson field redefinition given explicitly in [mixed ultraviolet field renormalization](mixed-uv-field-renormalization.md), including its contribution to the full leading mixed block. The quadratic local terms alone leave the nonzero pole shown in [mixed tadpole pole check](mixed-tadpole-pole-check.md); the induced mixed term cancels it and all other calculated mixed poles, entry by entry. The $c_G$ leading mixed term induces no pole because the singlet row of the transformation vanishes. This conclusion is conditional on the listed operator basis and subtraction.

The script uses the nonet $\Gamma_i$ table of [counterterm basis](../counting/counterterm-basis.md). In particular $\Gamma_6=1/16$ and $\Gamma_{25}=0$: the implemented pole cancellation finds no discrepancy with these values after the field redefinition. The $SU(3)$ value $\Gamma_6=11/144$ belongs to a different theory, with the singlet removed. Apparent coefficients inferred by demanding cancellation in the original fields without the divergent transformation are field-coordinate dependent and do not establish a disagreement with the table. The combinations determined from the meson self-energy are listed in [meson self-energy](meson-self-energy.md); they do not independently determine every individual coefficient.

## Scale and verification

Run `.venv/bin/python scripts/derive_nlo_loop_self_energies.py --cancel` and `.venv/bin/python scripts/derive_nlo_loop_self_energies.py --mixed`. The first checks the meson poles and the scale variation at $\mu=0.5,0.77,1.0\,$GeV after the field redefinition; the second checks the calculated mixed poles and compares their local coefficients with the quadratic-operator derivation. The meson test distinguishes changes in the off-shell matrix under a field transformation from pole positions. A separate full mixed-block scale-independence check and independent verification of every finite mixed coefficient have not yet been established; these checks are tracked in Beads issue `alp-gii.2.1`. Interaction vertices must also be transformed when constructing physical amplitudes.
