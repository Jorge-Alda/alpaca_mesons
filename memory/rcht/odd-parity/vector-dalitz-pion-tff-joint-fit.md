---
title: Conditional joint fit of vector Dalitz and pion single-virtual TFF data
status: numerical
created: 2026-10-24
updated: 2026-10-24
bibtex_keys: [Alda:2026nrd, CLEO:1997fho, BaBar:2009rrj, Adlarson:2016hpp, NA60:2016nad]
provenance: numerical-computation
---

# Scope and inputs

This is a conditional fit of the original seven VJP and four VVP coefficients
to the existing low-energy width block, the A2/NA60
$\omega\to\pi^0\ell^+\ell^-$ spectra, and the spacelike pion transition
form factor. The even sector is fixed to the preliminary
`single_trace_negative` representative. The odd fit refloats its identifiable
coefficient combinations; shared total widths and the published BaBar TFF
normalization nuisance are profiled. No new operators or virtuality cutoff
are introduced.

The A2 and NA60 likelihood inputs are recorded in the
[vector Dalitz fit](vector-dalitz-conditional-fit.md). The 15 CLEO pion points
cover $Q^2=1.64$--$7.90\,\mathrm{GeV}^2$ and the 17 BaBar points cover
$Q^2=4.24$--$34.36\,\mathrm{GeV}^2$. Their tables and diagonal-error
prescriptions are described in the [CLEO extract](../../extracts/CLEO:1997fho/tff-tables-and-errors.md)
and [BaBar extract](../../extracts/BaBar:2009rrj/pi0-tff-table.md). The
source-specific masses and matrices use Alda:2026nrd tree inputs. The complete
hard spacelike ranges are conditional diagnostics because this retained tree
TFF does not enforce the QCD $1/Q^2$ asymptotic falloff.

# Method

The executable
[`fit_joint_dalitz_tff.py`](../../../scripts/fit_joint_dalitz_tff.py)
constructs exact affine amplitude rows in the original odd coefficients,
uses the multi-vector photon response, and converts the canonical coefficients
to the legacy backend convention before evaluating the TFF. It fits the
low-energy widths, Dalitz bin-center values, CLEO and BaBar $Q^2|F|$ points,
and the stated Gaussian nuisance constraints. Experimental covariance is
diagonal where no bin covariance is supplied; experiment blocks are taken
independent. Three deterministic starts are compared, and the profiled
$c_1-c_2+c_5$ coordinate is scanned with the other fit variables reoptimized.

Run:

```text
PYTHONPATH=scripts .venv/bin/python scripts/fit_joint_dalitz_tff.py
```

The data hashes, residuals, covariance, profiles and blockwise fit statistics
are saved in
[`joint_dalitz_tff_conditional_fit.json`](../../../data/vector_sector/joint_dalitz_tff_conditional_fit.json).

# Results

| Included observations | $\chi^2$ | Nominal degrees of freedom |
|---|---:|---:|
| Existing widths | 5.56 | 5 |
| Widths and Dalitz spectra | 136.86 | 40 |
| Widths and pion TFF | 92.12 | 36 |
| Widths, Dalitz spectra and pion TFF | 412.13 | 72 |

The nominal joint-fit probability under the stated diagonal Gaussian model is
$3.66\times10^{-49}$. This is a descriptive diagnostic, not a precision
significance claim: the fit is conditional, the covariance is incomplete,
the hard-range TFF model is not asymptotically matched, and global minimality
is not certified. The best joint solution has a substantially poorer
low-energy block ($\chi^2=212.22$) and NA60 block ($\chi^2=150.44$); the
CLEO and BaBar pion blocks contribute $9.17$ and $22.98$, respectively.

The finite-virtuality rows increase the original odd-coefficient rank from
seven to eight. Dalitz and pion-TFF rows each probe one direction in the old
four-dimensional null space, and their combined projection remains rank one
after the separate experiment-normalization directions are removed. The
other three old null directions remain invisible in these data.

Both data classes constrain the same relevant coefficient combination, but
the conditional fits disagree strongly:

| Fit block | $c_1-c_2+c_5$ | Local conditional error |
|---|---:|---:|
| Widths and Dalitz spectra | $-0.14453$ | $0.00933$ |
| Widths and pion TFF | $2.6364\times10^{-5}$ | $4.9373\times10^{-6}$ |
| Combined | $2.2139\times10^{-5}$ | $4.9666\times10^{-6}$ |

These local errors are not physical calibration uncertainties. The TFF rows
therefore lift the algebraic null but do not improve the accepted description
with the retained model; the strong tension calls for an improved finite- and
high-virtuality description before using the TFF fit to calibrate the ALP
doubly virtual form factor. The ALP source projection itself remains the
separate result in the [ALP null note](../../alp/lepton-pair/alp-isovector-tff-fit-null-projection.md).

# Verification and limits

The two focused fit checks passed, `git diff --check` was clean, and repeated
starting points recovered the reported minimum twice. A third start found a
poorer branch, so the nonlinear global minimum is not proven. The numerical
output reports block residuals, nuisance shifts, coefficient profiles and
the full local covariance for audit. The result is conditional on the fixed
even sector, original-operator truncation, zero-width tree vector response,
independent diagonal data errors, and the selected bin-center Dalitz
approximation.
