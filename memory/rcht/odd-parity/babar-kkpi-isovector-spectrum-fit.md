---
title: Conditional fit of the BaBar isovector K star K spectrum
status: numerical
created: 2026-10-30
updated: 2026-10-30
bibtex_keys: [BaBar:2007ceh, ParticleDataGroup:2024cfk]
provenance: numerical-computation
---

# Purpose and inputs

This fit describes the 19 published Table 7 bins for the isovector component
of $e^+e^-\to K^*(892)K$ in [@BaBar:2007ceh]. The measured quantity is the
Dalitz-fit quasi-two-body component. It is not an inclusive $KK\pi$ rate and
is not a direct measurement of the low-energy RChT coefficients. The fitted
interval is $1.40$--$1.78\,\mathrm{GeV}$; the last bin crosses the physical
$\tau$ endpoint, but is retained as published.

The executable [fit](../../../scripts/fit_babar_kkpi_isovector.py) reads the
transcribed [Table 7 bins](../../../data/vector_sector/babar_2007_isovector_kstar_k_tau_region.csv)
and writes the reproducible [fit output](../../../data/vector_sector/babar_kkpi_isovector_fit.json).
The BaBar key was already present in `bibliography.bib`. Masses use the PDG
2024 charged-kaon and pion values; the $K^*(892)$ mass is a stated stable-pole
input of $0.89166\,\mathrm{GeV}$.

# Model and assumptions

The cross section uses the BaBar quasi-two-body resonance form
$$
\sigma(s)=12\pi\,\mathcal P_{K^*K}(s)\left|
\frac{A^{\mathrm{n.r.}}}{s}+
\frac{\sqrt{\Gamma_{ee}^{\rho'}\mathcal B_{K^*K}^{\rho'}}
\sqrt{\Gamma_{\rho'}/\mathcal P_{K^*K}(M_{\rho'}^2)}}
{M_{\rho'}^2-s-i\sqrt{s}\,\Gamma_{\rho'}(s)}\right|^2,
$$
with
$$
\mathcal P_{K^*K}(s)=\left[\frac{\lambda(s,M_{K^*}^2,M_K^2)}{s}\right]^{3/2},
\qquad
\Gamma_{\rho'}(s)=\Gamma_{\rho'}\left[\mathcal B_{4\pi}
\frac{\mathcal P_{4\pi}(s)}{\mathcal P_{4\pi}(M_{\rho'}^2)}
+1-\mathcal B_{4\pi}\right],
$$
and $\mathcal P_{4\pi}(s)=(s-16M_\pi^2)^{3/2}/s$. The baseline follows
Eqs. (11), (12), (15), (18), and (19) of the source HTML for [@BaBar:2007ceh].
Each prediction is averaged uniformly across its energy bin because the
bin-level ISR and acceptance response is unavailable. The $K^*$ is stable in
the phase-space factor. No finite-width $K^*$ convolution is included.

The default fits fix the four-pion fraction to the source's global-fit
central value $\mathcal B_{4\pi}=0.64$ and the resonance phase to zero, while
fitting $M_{\rho'}$, $\Gamma_{\rho'}$, and
$\Gamma_{ee}^{\rho'}\mathcal B_{K^*K}^{\rho'}$. A second variant adds a real
background $A^{\mathrm{n.r.}}/s$, with its signed normalization fitted. A
separate fit frees $\mathcal B_{4\pi}$. For the published statistical errors
the likelihood treats bins as diagonal. The additional 15% Dalitz-fit
uncertainty and 5.1% Table 1 cross-section uncertainty are each shown as
separate common multiplicative Gaussian nuisances; those correlations are
assumptions because BaBar supplies no covariance matrix.

# Results

With the $\rho'$ alone, the conditional fit gives
$$
M_{\rho'}=1.46763\,\mathrm{GeV},\qquad
\Gamma_{\rho'}=0.18797\,\mathrm{GeV},\qquad
\Gamma_{ee}^{\rho'}\mathcal B_{K^*K}^{\rho'}=46.80\,\mathrm{eV},
\qquad \chi^2/\mathrm{dof}=20.37/16.
$$
The conditional local linear errors on the three fitted quantities are
0.01895 GeV, 0.01799 GeV, and 10.65 eV, respectively.
Allowing the real background gives
$$
M_{\rho'}=1.47435\,\mathrm{GeV},\qquad
\Gamma_{\rho'}=0.14682\,\mathrm{GeV},\qquad
\Gamma_{ee}^{\rho'}\mathcal B_{K^*K}^{\rho'}=40.84\,\mathrm{eV},
\qquad \chi^2/\mathrm{dof}=19.96/15,
$$
with signed square root of the background cross section at the reference
energy equal to $-0.2380\,\sqrt{\mathrm{nb}}$. The fit with a free four-pion
fraction reaches its zero lower bound and has a rank-deficient local
Jacobian (rank four for five fitted parameters); the fraction is not
identified by this restricted data sample, and its linearized errors are
suppressed. The background-fit local errors on mass, width, resonance
normalization, and signed background amplitude are 0.01350 GeV, 0.04920 GeV,
12.10 eV, and 0.2943 times the square root of nb. These are conditional on
the model and diagonal statistical errors. Changing the four-pion phase-space
power from $3/2$ to $5/2$ in the background fit changes $\chi^2$ only from
19.9563 to 19.9618. The common systematic nuisances remain essentially zero
at the minimum because the resonance normalization is free; they do not
resolve model or covariance ambiguity. With both common nuisances included,
the background-fit local error on the resonance normalization is 13.72 eV.

Using the source's Table 15 global-fit values as a fixed benchmark, rather
than fitting them to this subset, gives $\chi^2=54.40$ for these 19 bins.
Those Table 15 values came from a global fit to additional datasets,
including the relative isospin phase and other final states. The subset fit
is therefore a conditional description, not a reproduction of the source's
global result. The source gives the background cross section but not its
amplitude sign; the benchmark selects the positive square root. Increasing
the energy-bin quadrature from 12 to 24 points
changes benchmark predictions by less than $1.4\times10^{-15}\,\mathrm{nb}$.

# Interpretation and verification

Invoke the deterministic analysis with
```text
.venv/bin/python scripts/fit_babar_kkpi_isovector.py --output data/vector_sector/babar_kkpi_isovector_fit.json
```
The output stores the bin predictions, residual diagnostics, local Jacobian
ranks, fitted variants, nuisance assumptions, and the fixed Table 15
benchmark. The numerical quadrature refinement above checks the bin
averaging calculation.

This completes a conditional fit of the measured isovector $K^*K$ spectrum.
It does not determine the RChT coefficients. In particular, the $K^*$
component removes the local VPPP terms, so it has no $g_i$ sensitivity; the
full electromagnetic charge map and a coefficient-matched excited-vector
amplitude are also absent. A coefficient extraction needs those physical
ingredients and stability against the allowed resonance and covariance
variants.
