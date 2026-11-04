---
title: Joint BaBar and low-energy vector-sector likelihood baseline
status: numerical
created: 2026-11-04
updated: 2026-11-04
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx, ParticleDataGroup:2024cfk, BaBar:2007ceh, Qin:2024ulb]
provenance: numerical-computation
---

# Joint BaBar and low-energy vector-sector likelihood baseline

This record describes the current stored snapshot of the joint likelihood, not a
converged best fit. The authoritative executable is
[`scripts/fit_babar_kkpi_global_vector_sector.py`](../../../scripts/fit_babar_kkpi_global_vector_sector.py),
and its output is
[`data/vector_sector/babar_kkpi_global_vector_sector_fit.json`](../../../data/vector_sector/babar_kkpi_global_vector_sector_fit.json).
The reproduction command printed by the script is:

```text
.venv/bin/python scripts/fit_babar_kkpi_global_vector_sector.py --output data/vector_sector/babar_kkpi_global_vector_sector_fit.json
```

## Inputs and row inventory

The script assembles one residual vector containing 73 observed summaries and
two normalization constraints. The JSON reports the following six blocks:

| block | count | stored $\chi^2$ |
| --- | ---: | ---: |
| odd low-energy summaries | 12 | 11.154174549835913 |
| masses and leptonic widths | 7 | 10.029707538403345 |
| shared total widths | 7 | 102.45464257297881 |
| $\phi\to K\bar K$ branching summaries | 2 | 16.73735273513131 |
| inclusive BaBar $KK\pi$ bins | 45 | 156.6343800174444 |
| normalization constraints | 2 | 10.974594826121375 |

The first five rows sum to 73; the final block is appended as two pull
coordinates. The script asserts the 12 odd observations and 45 BaBar rows and
constructs the boundaries $(0,12)$, $(12,19)$, $(19,26)$, $(26,28)$,
$(28,73)$, and $(73,75)$.

## Fitted coordinates and nuisances

There are 32 fitted coordinates: six scaled even-sector coordinates
$(m_\rho,m_{K^*},\Delta_V,\kappa_V,F_V,z_V)$, the shared $G_V$ coordinate,
five free total-width pulls for $\rho^0$, $\rho^+$, $\omega$, $\phi$, and
$\eta'$, twelve odd-sector combinations

$$
c_1-c_2+c_5,\quad c_1-c_2-c_5+2c_6,\quad c_1+c_2-c_5,\quad
c_3,\ c_4,\ d_1,\ d_2,\ d_3,\ g_1+2g_2-g_3,\ g_2,\ g_4,\ 2g_4+g_5,
$$

two tensor products $(g_T C_T,g_T D_T)$, four excited-pole residues
$(\beta_{\mathrm{first},KK\pi^0},\beta_{\mathrm{first},K_SK\pi},
\beta_{\mathrm{second},KK\pi^0},\beta_{\mathrm{second},K_SK\pi})$, and two
BaBar normalization pulls. The coordinate count is therefore
$6+1+5+12+2+4+2=32$.

The script maps the twelve odd coordinates to the original odd coefficients by
the representative convention $c_5=c_7=d_4=g_3=0$; this is a coordinate choice,
not a claim that the original sixteen coefficients are all independently
measured.

## Conventions, assumptions, and exclusions

The stored JSON and script specify a negative single-trace source and the
near-ideal inverse-mass branch, with no coefficient or previous-fit
$\chi^2$ priors. The fit uses $K^*$ $K\pi$ saturation and counts each of the
two duplicate experimental $K^*$ widths once. PDG summaries are diagonal except
for the $\phi\to K\bar K$ branching-ratio correlation $-0.72$; unprovided
cross-covariances are set to zero. The ground neutral propagator is fully
rotated, excited multiplets are ideal, vertex masses are averaged, and phase
space uses physical masses. Ground $\rho/K^*$ reference running shapes are
rescaled to shared on-shell widths, with the $K^*$ propagator using the
charged/neutral mean width. Source excited poles are fixed with a $2.5$ GeV
regulator, and there is no $K_2^*$ contribution.

The BaBar part uses the two inclusive below-$1.78$ GeV CSV inputs loaded by the
script, with coherent current integration and no separate BaBar Table 7 isobar
component. That is an intentional duplicate-input exclusion: the isobar
component is not added on top of the inclusive event-rate bins.

## Current stored optimizer state

The JSON snapshot has `success=false`, message “The maximum number of function
evaluations is exceeded.”, `nfev=1500`, and `njev=1453`. It reports

$$
\chi^2=307.98485223991514,\qquad \mathrm{dof}=43,
$$

with nominal $p=9.606791712164976\times10^{-42}$ and
`model_adequacy="inadequate"`. The global minimum is not proven and coefficient
confidence intervals are not certified. The stored scaled optimality and
gradient-infinity norm are both $8967.380245802657$; these values, together
with `success=false`, are why this record must be treated as a nonconverged
baseline snapshot rather than a calibrated fit result.

The JSON also records the six input SHA-256 values under `input_sha256`; those
metadata values identify the exact source snapshot used for the fit. The
numerical statements above were checked against the corresponding JSON fields
and against the script's row-count assertions, parameter declarations, block
boundaries, and output construction.

## Distinction from the older BaBar-only fit

This global likelihood is not a second optimization of the older BaBar-only
inclusive analysis. It combines the 12 odd low-energy summaries, 7
mass/leptonic summaries, 7 shared total-width summaries, 2
$\phi\to K\bar K$ branching summaries, 45 inclusive BaBar bins, and 2
normalization constraints. Its current stored snapshot is therefore

$$
\chi^2=307.98485223991514\ \text{for}\ 43\ \text{dof},
$$

with `success=false` and the evaluation limit reached.

The older record
[`babar-kkpi-unified-vector-coefficient-fit.md`](babar-kkpi-unified-vector-coefficient-fit.md)
describes a separate likelihood containing only the two inclusive BaBar
$KK\pi$ channels (Tables 2 and 3), through the $1.76$--$1.78\,\mathrm{GeV}$
bin. It varies ten $c_i,d_i,g_i$ combinations and normalization nuisances;
it does not include the low-energy, mass, leptonic, shared-width, or
$\phi\to K\bar K$ blocks. Its stored output gives

$$
\chi^2=35.71636391362228\ \text{for}\ 29\ \text{dof},
$$

quoted in the existing record as $\chi^2=35.7164$ for 29 dof with
$p=0.182$. The values are not competing results for one dataset: they refer
to different likelihood scopes and parameterizations. The global value is the
current nonconverged baseline snapshot, whereas the older value is the
converged BaBar-only inclusive fit.

The comparison was verified against both memory records and the outputs: the
global JSON fields `measurement_count`, `constraint_count`, `chi2`, `dof`,
`success`, and `nfev` give $(73,2,307.98485223991514,43,\mathrm{false},1500)$,
while `babar_kkpi_unified_vector_coefficients_fit.json` stores
`fit.chi_squared=35.71636391362228` and
`fit.nominal_degrees_of_freedom=29`.
