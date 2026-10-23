---
title: Conditional omega-pion Dalitz fit and preliminary-null projection
status: numerical
created: 2026-10-23
updated: 2026-10-23
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx, Adlarson:2016hpp, NA60:2016nad]
provenance: numerical-computation
---

# Conditional $ω\to\pi^0\ell^+\ell^-$ fit

## Assumptions

This is a conditional test at the stored central preliminary odd fit and its
single-trace-negative vector-sector point. The even parameters, vector-sector
source residues, and all odd coefficients except the common shift
$\delta c_1=\delta c_5$ are fixed. The finite-$q^2$ vector propagator uses
zero-width tree poles. Physical masses are used in the Dalitz phase-space
weight, while the pseudoscalar mass and source matrix in the odd vertex follow
the tree inputs from Alda:2026nrd. The published bin errors are treated as
diagonal and experiment-specific normalization directions are projected out
in the identifiability diagnostic. These are conditional approximations; the
papers supply no inter-bin covariance matrix.

For Adlarson:2016hpp, the table gives one quoted uncertainty per bin, not a
separate statistical/systematic breakdown. It is stored in the `stat_sigma`
field for compatibility with the data adapter and must be interpreted as the
quoted total error. For NA60:2016nad the statistical and systematic errors
are tabulated separately and added in quadrature for this comparison. The
released spectra are normalized transition-form-factor-squared values, not
raw decay counts. Their normalization at $q^2=0$ is therefore used through
$|G(q^2)/G(0)|^2$; no extra radiative width observation is included.

## Inputs and provenance

The 14 A2 bins and 22 NA60 bins are transcribed into
[`vector_dalitz_measurements.json`](../../../data/odd_parity/vector_dalitz_measurements.json).
Their source extracts are [Adlarson:2016hpp's electron spectrum](../../extracts/Adlarson:2016hpp/omega-pi-electron-spectrum.md)
and [NA60:2016nad's muon spectrum](../../extracts/NA60:2016nad/omega-pi-muon-spectrum.md).
The fixed coupling representative and four-dimensional old-fit null basis
come from [the preliminary-fit record](preliminary-fit-lagrangian-and-null-directions.md).
The finite-virtuality direct and exchange polynomials are the existing
canonical VJP/VVP contractions; their independent symbolic tensor check is
part of the fit script.

## Method

The differential Dalitz spectrum is proportional to the pointlike QED weight
times $|G(q^2)|^2$. For each bin the script compares (i) the published
bin-center value with the prediction and (ii) the pointlike-rate-weighted
average over the bin. It minimizes a one-coordinate diagonal Gaussian
least-squares function for the common $c_1,c_5$ shift, separately for A2,
NA60, and their combined bins.

At the joint best point, finite differences of all binned predictions with
respect to the exact old-fit null basis give a $36\times4$ weighted
Jacobian. Its singular values are also calculated after projecting out two
independent experiment-normalization columns. The null basis is reconstructed
from `summarize_odd_preliminary_fit.main(write_output=False)` and the named
candidate coordinates are checked to span it.

## Result

For bin-center evaluation, A2 alone gives
$\delta c_1=\delta c_5=-0.05553\pm0.01046$ with $\chi^2=1.87$ for 13
nominal degrees of freedom. NA60 alone gives
$-0.09500\pm0.00339$ with $\chi^2=119.17/21$. The joint result is
$-0.09041\pm0.00324$ with $\chi^2=135.30/35$. QED-weighted bin averaging
changes the joint shift only to $-0.09065\pm0.00325$ and gives
$\chi^2=135.70/35$.

The joint old-null Jacobian has singular values
$(238.90,3.32\times10^{-10},2.04\times10^{-16},0)$ in the script's
standardized observable metric. After independent A2 and NA60 normalization
directions are projected out, these are
$(74.47,2.85\times10^{-10},1.50\times10^{-16},0)$. Thus the spectra add one
shape-sensitive direction in the former four-dimensional null subspace,
including after normalization projection. The data do not lift the other
three old null directions in this channel block.

## Interpretation and reproduction

The one new singular direction establishes local shape sensitivity to the
old null subspace, but the poor NA60 and joint fits make the fitted shift
unacceptable as a coupling calibration. The disagreement between the A2 and
NA60 preferred shifts and the large structured NA60 residuals indicate that
the frozen zero-width tree-level model does not describe the combined
spectra. Do not quote the local $0.0032$ error as a physical uncertainty or
as a successful determination of $c_1$ and $c_5$ separately. The fit only
varies their common shift, so their difference remains fixed by the chosen
preliminary representative.

The calculation and its machine-readable output are
[`fit_vector_dalitz.py`](../../../scripts/fit_vector_dalitz.py) and
[`vector_dalitz_conditional_fit.json`](../../../data/vector_sector/vector_dalitz_conditional_fit.json).
Reproduce the result with:

```sh
PYTHONPATH=scripts .venv/bin/python scripts/fit_vector_dalitz.py
```

The input definitions are linked to the source tables above and to the
retained theory convention in [the fit plan](vector-dalitz-form-factor-fit-plan.md).
