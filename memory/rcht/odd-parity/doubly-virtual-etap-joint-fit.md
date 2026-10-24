---
title: Conditional fit including BaBar doubly virtual eta-prime data
status: numerical
created: 2026-10-24
updated: 2026-10-24
bibtex_keys: [Alda:2026nrd, CLEO:1997fho, BaBar:2009rrj, BaBar:2018zpn, Adlarson:2016hpp, NA60:2016nad]
provenance: numerical-computation
---

# Matching and scope

This adds BaBar:2018zpn's double-tag spacelike $\eta'$ transition form
factor to the conditional original-seven-VJP/four-VVP fit already containing
low-energy widths, A2/NA60 $\omega\to\pi^0\ell^+\ell^-$ data, and CLEO/BaBar
single-virtual pion TFF data. The odd coefficients are common RChT parameters;
the meson and ALP source matrices are not interchangeable. For this SM data
fit, the predictor uses Alda:2026nrd's physical tree-level $\eta'$ generator,
mass and $\chi$ matrix. These fitted coefficients can later be inserted into
the separately derived ALP no-mixing source. No meson mixing rotation is
applied to the ALP source.

BaBar:2018zpn defines its scalar with an $e^2$ prefactor and reports
$F_{\eta'}(Q_1^2,Q_2^2)$ in $\mathrm{GeV}^{-1}$. After matching the tensor
normalization, the fit compares absolute values, so an overall epsilon-tensor
sign does not affect the likelihood. The spacelike inputs are
$s_i=-Q_i^2$. The five published points are broad two-dimensional bins whose
form factors were extracted from bin-integrated cross sections and reported
at efficiency-weighted mean virtualities; evaluating the predictor at those
means is an approximation, not an exact bin-folded likelihood. The paper
reports separate asymmetric statistical, systematic and model errors but no
bin covariance. The primary Gaussian result symmetrizes the statistical
error and adds all three error components in quadrature, with a diagonal
covariance. Alternate asymmetric-statistical and upper-limit sensitivities
are also saved. The upper-limit treatment is an assumed Gaussian censoring
surrogate, not the source event-count likelihood.

# Reproduction

The executable
[`fit_joint_dalitz_tff.py`](../../../scripts/fit_joint_dalitz_tff.py)
reads the pinned five-bin input
[`babar_double_virtual_etap.json`](../../../data/odd_parity/babar_double_virtual_etap.json)
and records input hashes, residuals, nuisance profiles, and covariance in
[`joint_dalitz_tff_conditional_fit.json`](../../../data/vector_sector/joint_dalitz_tff_conditional_fit.json).
Run:

```text
PYTHONPATH=scripts .venv/bin/python scripts/fit_joint_dalitz_tff.py
```

# Results

| Fit selection | $\chi^2$ | Nominal degrees of freedom |
|---|---:|---:|
| Widths, Dalitz and singly virtual pion TFF | 412.13 | 72 |
| Widths plus doubly virtual $\eta'$ TFF | 112.72 | 9 |
| All blocks, symmetrized Gaussian errors | 1297.87 | 77 |
| All blocks, asymmetric statistical errors | 1298.91 | 77 |
| All blocks, upper-limit surrogate | 1291.34 | 77 |

Before refitting, the existing width+Dalitz+single-virtual solution predicts
$|F_{\eta'}|=(0.188,0.223,0.189,0.229,0.240)\,\mathrm{GeV}^{-1}$ at the five
BaBar points, versus measured central values
$(0.01432,0.00535,0.00824,0.00607,0.00871)\,\mathrm{GeV}^{-1}$. Its added
double-tag residual is $\chi^2=5.48\times10^4$ under the primary diagonal
Gaussian approximation. Refitting the combined blocks lowers the double-tag
block contribution to $19.34$, but raises the low-energy and NA60 contributions
to $352.75$ and $814.20$. The nominal total-fit probability is extremely
small under these approximations and should be read as a model/data tension
diagnostic, not a precision significance.

The original odd-coefficient rank remains eight, out of eleven. Both
singly- and doubly-virtual rows probe the same one-dimensional component of
the preliminary four-dimensional null space; after normalization projection
the double-virtual rows still add no independent null direction. The other
three null combinations remain unconstrained. The joint fit's
$c_1-c_2+c_5=(3.37\pm0.48)\times10^{-5}$ error is local and conditional; the
poor global description prevents treating it as a physical calibration.

# Interpretation and limits

The data addition does not improve the fit. It demonstrates directly that
the retained tree TFF's double-virtual shape is incompatible with the
measured $\eta'$ form factor when the existing blocks are held in the same
original-operator model. The virtualities are all large, and the model does
not have a validated short-distance completion; this fit must not be used as
a controlled ALP prediction. The mismatch is a reason to complete the
finite- and high-virtuality TFF description and its matching before inferring
an ALP TFF uncertainty from the shared coefficients.

Four focused checks passed, `git diff --check` was clean, and up to ten
deterministic cross-fit starts recovered the reported minima. This is not a
proof of the global nonlinear minimum. Results remain conditional on the
fixed even sector, original operator set, tree spectral response, independent
diagonal data blocks, average-point treatment of broad bins, and nuisance
prescriptions above.
