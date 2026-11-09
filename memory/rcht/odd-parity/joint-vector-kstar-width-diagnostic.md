---
title: Joint vector fit K-star width diagnostic
status: numerical
created: 2026-11-09
updated: 2026-11-09
bibtex_keys: [Alda:2026nrd, Ecker:1988te, Ecker:1989yg, ParticleDataGroup:2024cfk]
provenance: numerical-computation
---

# Joint vector fit $K^*$ width diagnostic

## Inputs and conventions

The accepted `basin_probes/even_shared_widths` endpoint and its qualification are described in the [basin-probe record](joint-vector-fit-basin-probes.md). The [joint likelihood baseline](babar-kkpi-global-vector-sector-fit-baseline.md) specifies its physical branch, the shared propagator width, and the data treatment. The [width audit](shared-even-sector-gv-constraint-audit.md) documents the measured $K^{*+}$ and $K^{*0}$ widths, their physical masses and the approximation of zero cross-covariance [@ParticleDataGroup:2024cfk]. Masses and widths are in GeV, in natural units.

The leading interaction, decay-constant convention, flavor factors, and tree-level two-body formula are given in the [conditional $G_V$ record](../foundations/vector-gv-fit.md) [@Ecker:1988te; @Ecker:1989yg; @Alda:2026nrd]. The model assumes $K\pi$ saturation and uses the fitted common tensor mass in the vertex normalization, with physical daughter masses in phase space. The comparison using physical masses replaces the tensor mass by each physical charged or neutral $K^*$ mass only in the external tree-width normalization; the fitted $G_V$ and daughter masses are held fixed. This comparison is not a refit.

## Method

1. Load the qualified endpoint and reconstruct its fitted $G_V$, tensor mass and two $K\pi$ widths. Their arithmetic mean supplies the propagator's common on-shell width. Evaluate the existing ground-state propagator at the squared fitted tensor mass and extract $-\operatorname{Im}(D^{-1})/m_{K^*}$.
2. Recompute each charged and neutral tree-level $K\pi$ width from the two physical daughter channels using its physical parent mass for both phase space and tensor normalization. Take their arithmetic mean for comparison with the common propagator width.
3. For measured widths $\Gamma_{K^{*+}}$ and $\Gamma_{K^{*0}}$ with covariance $C$, minimize the Gaussian sum over an unconstrained common width. The minimum is

$$
\chi^2_{\min}=\frac{(\Gamma_{K^{*+}}-\Gamma_{K^{*0}})^2}{C_{++}+C_{00}-2C_{+0}}.
$$

The two measurements are treated as uncorrelated because their cross-covariance is unavailable. The minimum has one nominal degree of freedom. This common-width result is distinct from the fit of one shared $G_V$ with charge-dependent predicted widths in the [width audit](shared-even-sector-gv-constraint-audit.md).

## Reproduction and verification

```text
.venv/bin/python scripts/diagnose_joint_kstar_width.py
.venv/bin/python -m unittest scripts.test_diagnose_joint_kstar_width
```

The script reports the inputs and width comparisons as JSON. The unit tests compare the common-width formula and fitted width with an independent numerical minimization for a correlated pair, and check equal measurements. The executable is [`diagnose_joint_kstar_width.py`](../../../scripts/diagnose_joint_kstar_width.py).
