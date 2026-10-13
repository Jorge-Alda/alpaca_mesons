---
title: Rank of the action-based radiative and neutral TFF channel sets
status: numerical
created: 2026-10-04
updated: 2026-10-13
bibtex_keys: [Ruiz-Femenia:2003jdx, ParticleDataGroup:2024cfk, Aloni:2018vki, Niecknig:2012sj]
provenance: numerical-computation
---

# Rank of the radiative and neutral TFF channel sets

## Inputs and assumptions

The executable rank probe is [`scripts/diagnose_radiative_rank.py`](../../../scripts/diagnose_radiative_rank.py), which uses the action-derived amplitudes in [radiative-sm-amplitude.md](radiative-sm-amplitude.md) and [radiative-tff-amplitude.md](radiative-tff-amplitude.md). It combines the five short-distance OPE rows with normalized coupling-sensitivity rows. No data errors or fit enter.

The illustrative channel set uses canonical ideal vector matrices, exploratory neutral-pseudoscalar matrices and residues, and the leading-order isospin spurion values stated in the SM amplitude record. These assumptions do not fix physical $\eta,\eta'$ inputs. The common internal pole mass is the baseline. Distinct $\rho,\omega,\phi$ pole masses are included only as a phenomenological sensitivity extension; a symmetry-breaking action and its matching have not been supplied.

## Method

For each channel the script evaluates the action-derived amplitude and differentiates its scalar coefficient with respect to $(c_1,\ldots,c_7,d_1,\ldots,d_4)$. Rows are normalized before singular-value decomposition. Rank is determined with the documented numerical tolerance. The TFF rows use the full WZW, VJP and VVP amplitude at the selected virtualities, with the same explicit meson and spurion inputs as the radiative rows.

## Results

- The five OPE constraints have rank 5 of 11.
- Adding the ten radiative channels with a common internal pole gives rank 7 of 11. The additional directions are one correlated $c_3,d_2$ combination and $c_4$ from charged-kaon channels.
- Using distinct internal vector pole masses gives rank 8 of 11. The eighth normalized singular value is $0.00357049$; the remaining null singular values are below $3\times10^{-16}$. This establishes algebraic sensitivity for the chosen extension, not statistical usefulness.
- In both radiative runs, $c_7$, $d_4$, and the common shift $\delta c_2=\delta c_5=\delta c_6$ are unprobed.
- Adding neutral TFF rows with independent exploratory mass-spurion inputs raises the common-pole rank from 7 to 8 by resolving a second $c_3,d_2$ combination; the split-pole rank remains 8. If the TFF inputs instead impose $F_H=2m_P^2F_Q$, the TFF rows lie in the OPE span and add no rank.

## Verification and limits

Run:

```text
.venv/bin/python scripts/derive_radiative_sm_amplitude.py
.venv/bin/python scripts/derive_radiative_tff_amplitude.py
.venv/bin/python scripts/diagnose_radiative_rank.py
```

The scripts assert the canonical $\omega\to\pi^0\gamma$ amplitude, the WZW and resonance normalization, and the stated rank results. Physical $\eta,\eta'$ field and decay-constant inputs, a consistently matched split-pole action, experimental covariance, and channel-composition issues remain prerequisites for a physical fit. These ranks alone imply no saturation or fit conclusion.
