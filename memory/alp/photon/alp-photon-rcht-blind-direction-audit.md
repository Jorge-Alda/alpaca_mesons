---
title: Real-photon ALP resonance sensitivity to SM blind directions
status: derived
created: 2026-10-13
updated: 2026-10-13
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx]
provenance: analytical-derivation
---

# Real-photon ALP resonance audit

## Inputs and assumptions

Use the [ALP source insertions](../../rcht/odd-parity/two-source-alp-rcht-odd-parity-extension.md),
the [actual two-photon contractions](../../rcht/odd-parity/radiative-tff-amplitude.md),
and the [no-vector photon invariant](alp-photon-invariant-from-mixing.md).
The existing no-vector WZW/contact coefficient has been derived, but the
fully normalized resonance-corrected ALP amplitude has not been assembled.
The SM pion or eta TFF is not that missing calculation.

Assume one ALP at $O(1/f_a)$, the retained single-trace tree action,
electromagnetic-compatible sources and $f_-=0$ at the one-ALP vacuum vertex.
Set both photons real, $s_1=s_2=0$. Preserve
$X=\mathbb M_a^{-1}(\mathcal C)$ and $Y=c_q^{AP}+X$ independently.
For source responses $W_i$, define formal derivative and mass weights

$$
D_i=\operatorname{tr}(\{W_i,Q\}X),\quad
M_i=\operatorname{tr}(\{W_i,Q\}\{Y,B_0m_q\}),
$$

and corresponding pair weights $D_{12},M_{12}$ with $\{W_1,W_2\}$.
All spectral/source factors can be absorbed into these independent weights.
The ALP insertions give stripped mass/derivative relative weight $-2M/D$;
overall physical normalization is deliberately not assigned here.

## Executable projection

Run [the script](../../../scripts/audit_alp_gamma_gamma_blindness.py):

```text
.venv/bin/python scripts/audit_alp_gamma_gamma_blindness.py
```

It substitutes the real virtualities in the actual Lorentz contractions,
retains both photon placements and both VVP assignments, constructs the
independent-spurion polynomials, and differentiates the coefficient rows.
It finds exact zero projection onto $c_7$, $d_4$ and the common
$\delta c_2=\delta c_5=\delta c_6$ shift. Each photon placement separately
annihilates $c_7$ and the common shift; nonideal response weights cannot undo it.

The resulting stripped VJP dependence uses $c_1+c_2-c_5$ and $c_3$;
the stripped VVP dependence uses $d_1-d_3$ and $d_2$. There is no equal-spurion
assumption tying the derivative and mass terms. The WZW/UV photon seed itself
has no dependence on these eleven resonance coefficients.

## Interpretation

The retained real-photon ALP resonance projections do not access a structural
blind direction of the current SM sensitivity calculation. Their complete
normalization and interference with the WZW/contact term still require
[the completion plan](alp-photon-rcht-completion-plan.md). No full predictor,
width, numerical ALP coupling or loop/extra-operator independence is claimed.
