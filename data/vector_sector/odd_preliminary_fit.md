# Conditional preliminary fit of the original eleven operators

One fixed experimental-input setup, four distinct source hypotheses, three one-at-a-time input sensitivities and two data-removal checks for the negative single-trace source case, retained as a conditional working example because it describes the primary data; archive extraction without new minimization.

Charged rho and Kstar mass summaries; evaluated phi electronic width; near-ideal inverse mass root; primary twelve odd summaries. These choices define a conditional comparison, not a uniquely preferred source hypothesis.

$\epsilon^{0123}=+1; metric (+---); canonical c_i=-c_i^R and d_i=-d_i^R; old archives remain reference backend inputs.$

Quoted diagonal experimental summaries, fixed central Alda:2026nrd neutral inputs, tree-level real-photon model. Null directions are unmeasured; representatives are not eleven separately measured coefficients. Local prediction errors are sign-branch conditional; no global nonlinear confidence intervals.

| Case | Source hypothesis | Total chi2 | Odd rank | Unmeasured coefficient directions |
|---|---|---:|---:|---:|
| central_universal | universal | 362.14547 | 7 | 4 |
| central_single_trace_positive | single_trace_positive | 346.29371 | 7 | 4 |
| central_single_trace_negative | single_trace_negative | 5.5675233 | 7 | 4 |
| central_double_trace | double_trace | 54.216008 | 7 | 4 |
| change_rho_input | single_trace_negative | 5.6466225 | 7 | 4 |
| change_kstar_input | single_trace_negative | 12.158901 | 7 | 4 |
| change_phi_input | single_trace_negative | 5.6276328 | 7 | 4 |
| without_broad_rho | single_trace_negative | 0.37741791 | 7 | 4 |
| without_constrained_Kstar | single_trace_negative | 3.537427 | 7 | 4 |

## Canonical coefficient representatives

These dimensionless values use a numerical minimum-norm choice for the unmeasured directions. Their zeros or values do not constrain those directions.

| Case | c1 | c2 | c3 | c4 | c5 | c6 | c7 | d1 | d2 | d3 | d4 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| central_universal | -0.001539904 | 0.0003138976 | -0.001790481 | 0.0241735 | 0.001539904 | -0.001853801 | 3.469447e-17 | 0.2515078 | -0.01172403 | 0.2912592 | -0 |
| central_single_trace_positive | -0.002555957 | 0.002333397 | -0.001742148 | 0.02441372 | 0.002555957 | -0.004889354 | -4.857226e-17 | 0.288533 | -0.0104271 | 0.3286617 | -0 |
| central_single_trace_negative | -0.00351401 | -0.02407059 | 0.005638958 | 0.01952114 | 0.00351401 | 0.02055658 | 3.469447e-18 | -0.1296398 | 0.003351725 | 0.06535895 | -0 |
| central_double_trace | 0.01815696 | 0.06570483 | -0.07943627 | -0.02428666 | -0.01815696 | -0.04754787 | -0 | 0.9636437 | -0.4977092 | 0.7723491 | -0 |
| change_rho_input | -0.003581372 | -0.02417257 | 0.005671896 | 0.01949065 | 0.003581372 | 0.0205912 | -3.642919e-17 | -0.1311947 | 0.003442626 | 0.06524004 | -0 |
| change_kstar_input | -0.001073148 | -0.02028599 | 0.004544677 | 0.01995197 | 0.001073148 | 0.01921284 | -1.387779e-17 | -0.07366003 | -0.002292249 | 0.06892076 | -0 |
| change_phi_input | -0.003402205 | -0.02387643 | 0.005569727 | 0.01957747 | 0.003402205 | 0.02047423 | 3.989864e-17 | -0.1271583 | 0.00320972 | 0.06579039 | -0 |
| without_broad_rho | 0.002061877 | -0.01560307 | 0.002856007 | 0.02240751 | -0.002061877 | 0.01766495 | -0 | -0.006206432 | -0.002197154 | 0.07574939 | -0 |
| without_constrained_Kstar | 0.00667794 | 0.03011709 | -0.004183061 | -0.01658478 | -0.00667794 | -0.02343915 | -0 | 0.2065952 | -0.03862109 | -0.03739321 | -0 |

Complete predictions, pulls, available conditional local prediction errors, null bases, nuisance values and input hashes are in `odd_preliminary_fit.json`. No new optimizer search is performed. No individual coefficient uncertainty is assigned.
