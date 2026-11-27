---
title: ALP mass axis in only-up decay plots
status: numerical
created: 2026-11-19
updated: 2026-11-19
bibtex_keys: [DiLuzio:2026poh]
provenance: numerical-computation
---

# ALP mass axis in only-up decay plots

## Inputs and assumptions

The scan uses the leading-order only-up spectrum and conventions in [the input record](../../meson-mixing/finite-mixing/spectrum-and-poles/quadratic-mixing/only-up-finite-mixing-inputs.md). The supplied array of ALP masses is in GeV and corresponds to the increasing $\kappa_u$ scan. The axis assumes that this mass array is increasing.

## Method and output

The plotting routine in `scripts/implementation/murayama.py` interpolates between the supplied masses and $\kappa_u$ in both directions. It places logarithmic major and minor mass ticks only within the supplied mass range. This prevents interpolation outside that range from placing an out-of-range tick at a plot endpoint.

## Verification

Run `.venv/bin/python -m unittest scripts.implementation.test_murayama_alp_axis`. The test checks major and minor ticks against the mass range and checks the interpolation at each major tick.
