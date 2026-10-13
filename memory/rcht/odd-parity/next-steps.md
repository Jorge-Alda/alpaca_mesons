---
title: Next steps for the odd-parity vector-meson fit
status: mixed
created: 2026-10-04
updated: 2026-10-13
bibtex_keys: [Ruiz-Femenia:2003jdx, ParticleDataGroup:2024cfk, CLEO:1997fho, BaBar:2009rrj, Niecknig:2012sj, Aloni:2018vki]
provenance: assumption
---

# Next Steps

Start-of-session map. Read this, then
[the rank and saturation record](radiative-rank-saturation.md) for the open
question, and the two records it depends on for the machinery.

## Where things stand

Established and verified (all assertions exercised by scripts):

- **SM radiative amplitudes**: action-derived direct and complete neutral-vector
  exchange terms, including linear $f_-$ and $c_4$, reproduce the published
  canonical $\omega\to\pi^0\gamma$ amplitude. The implementation takes
  explicit meson matrices and mass-spurion inputs.
  [Amplitude record](radiative-sm-amplitude.md).
- **Rank and saturation**: the ten radiative channels add only **two** directions
  beyond the five OPE rows under the common internal-pole action baseline
  (rank 7 of 11). The split-pole rank 8 is an explicitly phenomenological
  extension, not a matched SU(3)-breaking action. See
  [the action-derived amplitude](radiative-sm-amplitude.md).
  `scripts/diagnose_radiative_rank.py`, [record](radiative-rank-saturation.md).
- **Space-like TFF amplitude and data connection**: the normalized tree-level
  amplitude includes both photon assignments, WZW and VVP terms and reproduces
  the published pion anomaly factor. CLEO/BaBar Observation factories are
  implemented. Under $F_H=2m_P^2F_Q$, the TFF adds no rank. With independent
  $F_H$ and the action-derived common-pole radiative rows, it raises rank $7\to8$ by
  resolving the other correlated $c_3,d_2$ direction. See
  [the TFF amplitude record](radiative-tff-amplitude.md).

Open, in order of how much they block:

1. **Alda-grounded $\eta,\eta'$ inputs.** The TFF and radiative predictors take
   explicit pseudoscalar matrices, mass-spurion traces, decay constants and
   residues. The indexed Alda extract defines the field rescaling but not the
   component-level $\eta$ mixing and normalization needed for physical fit
   inputs. The current neutral-meson choices remain exploratory.
2. **Single-virtual asymptotics.** The constrained one-resonance TFF approaches
   a nonzero constant at large $Q^2$, so it does not enforce the expected
   $1/Q^2$ falloff. Treat fits over the CLEO/BaBar range as diagnostics until
   this short-distance behavior is addressed.
3. **Split-pole matching.** The observed rank gain with distinct internal
   $\rho,\omega,\phi$ masses is not derived from a flavor-resolved
   SU(3)-breaking action; keep it as a sensitivity extension.

## Space-like transition form factors: amplitude and observations connected

The TFF calculation is complete at tree level in the stated one-vector
approximation. It uses genuinely different off-shell kinematics. Its rank
impact depends on independent mass-spurion inputs and the radiative baseline.

The normalized amplitude in [radiative-tff-amplitude.md](radiative-tff-amplitude.md)
includes both photon placements, Fourier signs and both VVP assignments. Under
the equal-spurion relation $F_H=2m_P^2F_Q$, TFF rows lie in the OPE row span and
add no rank. For independent $F_H$, they still annihilate $c_4,c_7,d_4$ and the
common $\delta c_2=\delta c_5=\delta c_6$ shift, while supplying the other
$c_3,d_2$ combination when combined with common-pole radiative rows.

The full amplitude is implemented in
[`scripts/implementation/pseudoscalar_tff.py`](../../../scripts/implementation/pseudoscalar_tff.py).
CLEO and BaBar data connect through the `Observation` factories in
[`odd_parity_data.py`](../../../scripts/implementation/odd_parity_data.py).
Physical $\eta,\eta'$ fits must supply explicit mass-spurion traces and field
normalizations; none are inferred from meson masses. No fit has been run.

**Data already validated and checksum-pinned, with loaders:**

| file | rows | content |
|---|---|---|
| `data/odd_parity/cleo_tff.csv` | 63 | $\pi^0,\eta,\eta'$, $Q^2\in[1.64,15.3]$ GeV$^2$ [@CLEO:1997fho] |
| `data/odd_parity/babar_pi0_table11.csv` | 17 | $\pi^0$ only, $Q^2\in[4.2,34.4]$ GeV$^2$, stat and syst separate, asymmetric errors [@BaBar:2009rrj] |

Both have checksum-pinned loaders and predictor-backed Observation factories in
`scripts/implementation/odd_parity_data.py`. This completes the TFF data
connection, but does not validate the one-resonance high-$Q^2$ behavior or
establish a fit result.

## Checked: $P\to\gamma\gamma$

The [real-photon check](p-to-gamma-gamma.md) establishes that these processes
are the $s_1=s_2=0$ limit of the complete recorded TFF: WZW contact, both
VJP photon placements and both VVP vector assignments. The width is
$\Gamma=\pi\alpha_{\rm em}^2m_P^3|F_P(0,0)|^2/4$, including the identical-photon
factor, and reproduces the published pion correction. No new operators are
required within the declared single-trace, one-vector tree action.

Equal-spurion inputs add no direction beyond the five OPE constraints;
independent spurions add at most one common $c_3,d_2$ direction (rank 6).
This is weaker than the independent-spurion off-shell TFF sensitivity.
Physical $\eta,\eta'$ predictions still require the Alda-grounded inputs above
and the nonet action remains an assumed extension. No experimental comparison
or fit has been performed.

## Running things

```text
.venv/bin/python scripts/derive_radiative_sm_amplitude.py
.venv/bin/python scripts/derive_radiative_tff_amplitude.py
.venv/bin/python scripts/derive_p_to_gamma_gamma.py
.venv/bin/python scripts/diagnose_radiative_rank.py          # ~30 s
.venv/bin/python scripts/extract_odd_parity_tables.py
.venv/bin/python -m unittest scripts.implementation.test_odd_parity_data
.venv/bin/python -m unittest scripts.implementation.test_pseudoscalar_tff_verification
.venv/bin/python -m unittest scripts.implementation.test_radiative_sm_verification
```

`scripts/implementation/__init__.py` eagerly imports the external `alpaca`
package, so the fit test modules must be loaded against the package directory
without executing that `__init__`; this is a recorded environment limitation, not
a property of the fit code.

## Do not

- Insert the ALP spurion $\mathbb M_a^{-1}(C)$ or $c_q^{AP}$ into a Standard
  Model channel. They belong to the $aV\gamma$ vertex.
- Rescale the nonet states to unit $\mathrm{Tr}(T^2)$; that changes the mixing
  angles and hence every flavor invariant.
- Fit a coefficient to force agreement with Eq. (27); the action-derived
  canonical amplitude reproduces it with the published coupling conventions.
- Treat a PDG average as independent data: $S=2.2$ on $\omega\to\pi^0\gamma$,
  $S=1.7$ on $\rho\to\pi^0\gamma$, $S=1.1$ on $\omega\to\eta\gamma$, and
  $\mathcal B(\eta'\to\rho^0\gamma)$ includes non-resonant $\pi^+\pi^-\gamma$.
