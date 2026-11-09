---
title: Sensitivity variations of the joint vector-sector fit
status: numerical
created: 2026-11-09
updated: 2026-11-09
bibtex_keys: [Alda:2026nrd, ParticleDataGroup:2024cfk, BaBar:2007ceh, Qin:2024ulb]
provenance: numerical-computation
---

# Sensitivity variations

## Inputs and provenance

- Variation definitions and their switches: [likelihood options](joint-vector-likelihood-variation-options.md).
- Pole uncertainties: [PDG 2024 excited poles](pdg2024-excited-vector-and-tensor-poles.md).
- Refit procedure and gates: [refit procedure](joint-vector-block-exclusion-and-variation-refits.md).
- Decisions on the pole scheme and on deferring the double-trace source: `openspec/changes/assess-joint-vector-fit-adequacy/design.md`.
- Ledgers:
  - `data/vector_sector/joint_vector_sensitivity_variations.json`, written by `scripts/refit_joint_vector_sensitivity.py --study variations`;
  - the linear pole screen, `data/vector_sector/joint_vector_pole_shift_screen.json`.
- Table printed by `.venv/bin/python scripts/summarize_joint_vector_fit_adequacy.py`.

## Assumptions

- The 5% bin-uncorrelated BaBar systematic is an assumption with no experimental source.
- The pole shifts are $\pm1\sigma$ PDG 2024 uncertainties, many of which are OUR ESTIMATE values. They are applied around the hard-coded baseline values.
- A variation is accepted when its selected run passes the optimizer and quadrature gates. For an unaccepted variation, the $\chi^2$ comes from a refit that did not reach its minimum, so it is only an upper bound on that minimum.

## Coverage

- Run: 47 variations, made of 11 named variations and 36 pole shifts.
- Accepted: 31.
- Not accepted: 16.
  - In 15 of them, all four starts stopped at the 5000-evaluation cap.
  - In `pole_omega_1420_width_minus_one_sigma`, the optimizer converged but the quadrature check failed.
  - `pole_phi_1680_mass_plus_one_sigma` had one of its four starts at the cap, but its selected run is accepted.
- Deferred: `source_single_plus_double_trace`, by user decision on 2026-11-07. It is not part of this study.

## Result: named variations

| variation | accepted | $\chi^2$ | dof | $p$ | worst block: $\chi^2$ (rows) | max $|$shift$|$ in baseline errors |
|---|---|---|---|---|---|---|
| BaBar cut 1.70 GeV | no (cap) | $\le200.67$ | 35 | $3.9\times10^{-25}$ | shared widths 104.0 (7) | 2.7 ($d_3$) |
| BaBar cut 1.60 GeV | no (cap) | $\le156.37$ | 25 | $5.6\times10^{-21}$ | shared widths 99.7 (7) | 9.9 ($d_3$) |
| regulator 2.0 GeV | no (cap) | $\le217.89$ | 43 | $3.1\times10^{-25}$ | shared widths 105.7 (7) | 1.0 ($d_3$) |
| regulator 3.0 GeV | yes | 227.32 | 43 | $6.6\times10^{-27}$ | shared widths 107.0 (7) | 0.2 ($d_3$) |
| tensor products zero | no (cap) | $\le244.50$ | 45 | $3.1\times10^{-29}$ | shared widths 104.7 (7) | 5.1 ($\beta$ second, $K_SK\pi$) |
| universal source | yes | 550.99 | 44 | $8.4\times10^{-89}$ | masses 354.1 (7) | 121 ($z_V$) |
| positive single-trace source | no (cap) | $\le556.25$ | 43 | $2.0\times10^{-90}$ | masses 336.6 (7) | 126 ($z_V$) |
| second inverse-mass branch | no (cap) | $\le1748.6$ | 43 | $<10^{-200}$ | odd 1030.7 (12) | 743 ($\Delta_V$) |
| $\phi\to K\bar K$ correlation zero | yes | 207.00 | 43 | $2.5\times10^{-23}$ | shared widths 76.2 (7) | 4.0 ($G_V$) |
| $\phi$ electronic width, PDG average | yes | 226.91 | 43 | $7.8\times10^{-27}$ | shared widths 108.0 (7) | 1.0 ($z_V$) |
| 5% BaBar systematic (assumption) | yes | 221.79 | 43 | $6.4\times10^{-26}$ | shared widths 106.8 (7) | 0.02 ($d_3$) |

## Result: pole shifts

All 36 pole shifts give $\chi^2$ between 218.3 and 238.2 for 43 dof, with $p<10^{-25}$. The shared total widths stay between 105.3 and 107.3 in every case.

- **Not accepted (cap):** $\rho(1450)$ mass $+$, width $\pm$; $\rho(1700)$ mass $+$, width $-$; $\omega(1420)$ width $+$; $\phi(1680)$ mass $-$, width $\pm$.
- **Not accepted (quadrature):** $\omega(1420)$ width $-$, with $\chi^2=238.2$.
- **Effect of the shifts:** the $K^*(1410)$, $K^*(1680)$, $\phi(2170)$ and $a_2(1320)$ shifts change $\chi^2$ by less than 0.2. The largest changes, up to about 12, come from the $\rho(1450)$, $\omega(1420)$ and $\phi(1680)$ shifts. In every pole shift, the largest parameter change is at most 3.2 baseline errors.

The full per-shift rows, with the worst block and the basin-rerun outcome, are in the ledger and in the summary script output.

## Conclusion

None of the 47 variations gives $p>10^{-20}$. The lowest $\chi^2$ per degree of freedom, about 6.3, comes from the unaccepted 1.60 GeV cut, which still has $p\sim10^{-21}$. Removing the $\phi\to K\bar K$ correlation reduces the shared-width $\chi^2$ to 76, which is the only notable redistribution. The changes to the source hypothesis and to the mass branch make the fit much worse, which supports the baseline choices.
