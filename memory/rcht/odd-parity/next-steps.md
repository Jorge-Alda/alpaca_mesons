---
title: Next steps for the odd-parity vector-meson fit
status: mixed
created: 2026-10-04
updated: 2026-10-04
bibtex_keys: [Ruiz-Femenia:2003jdx, ParticleDataGroup:2024cfk, CLEO:1997fho, BaBar:2009rrj, Niecknig:2012sj, Aloni:2018vki]
provenance: assumption
---

# Next Steps

Start-of-session map. Read this, then
[the rank and saturation record](radiative-rank-saturation.md) for the open
question, and the two records it depends on for the machinery.

## Where things stand

Established and verified (all assertions exercised by scripts):

- **SM flavor algebra** for $V\to P\gamma$ and $P\to V\gamma$: per-channel
  $F_Q,F_M$ and exchange invariants, from the meson matrices, $Q$ and
  $\mathcal M$. All fourteen channels distinct; no ALP quantity anywhere.
  `scripts/derive_radiative_sm_flavor_algebra.py`,
  [record](radiative-sm-flavor-algebra.md).
- **Per-operator kinematic factors**, direct and VVP-exchange, from explicit index
  contraction, reproducing the published $\omega\to\pi^0\gamma$ bracket in
  magnitude, plus the assembled per-channel $G$ and width.
  `scripts/derive_radiative_sm_kinematics_probe.py`,
  [record](radiative-width-kinematics.md).
- **Rank and saturation**: the ten radiative channels add only **two** directions
  beyond the five OPE rows (rank 7 of 11). Four null directions, all inherited
  from the OPE block alone. Non-ideal mixing and the two-body $V\to V'P$ decays
  were tested and add **zero** rank.
  `scripts/diagnose_radiative_rank.py`, [record](radiative-rank-saturation.md).
- **Space-like TFF kinematics**: the off-shell factors break the on-shell
  $c_6=c_1-c_2$ identity (giving $c_6=c_2-c_1$), and the three TFF rows raise
  the rank from 7 to **8 of 11**. The TFF route is viable; it lifts one of the
  four null directions.
  `scripts/derive_radiative_tff_kinematics.py`,
  [record](radiative-tff-kinematics.md).

Open, in order of how much they block:

1. **$\mathcal N$ and the absolute $c_3$ coefficient.** The common
   normalization and the $c_3$ functional mismatch block the
   $\sqrt{\chi^2/\text{ndof}}$ half of the saturation probe, i.e. two of its
   three falsification criteria.
2. **The $d_1$ sign discrepancy** in the exchange block: derived
   $(M^2-m_P^2)$, published $m_P^2$ only. Affects $d_1$ alone.
3. **TFF predictor and data connection.** The flavor factors are now derived in
   [the flavor record](radiative-tff-flavor.md). The source does not provide the
   normalized off-shell amplitude; derive the WZW plus resonance contributions,
   including the relevant VVP exchange, before creating CLEO/BaBar observations.
4. **The three remaining null directions** in $(c_2,c_4,c_5,c_6,c_7)$ and $d_4$,
   not lifted by the TFF.

## Space-like transition form factors: decisive test done, route viable

The most credible remaining candidate for lifting the null directions, because
they are the only available observable with genuinely different kinematics in the
$c$-sector: the pseudoscalar is the **initial** state and the photon is
**space-like off-shell** with invariant $Q^2$, so the operator-level kinematic
factors differ from the on-shell case rather than rescaling it.

**Decisive test, 2026-10-05: the rank rises.** The off-shell kinematic factors
were derived by explicit index contraction
([record](radiative-tff-kinematics.md),
`scripts/derive_radiative_tff_kinematics.py`). The off-shell case gives
$c_6=c_2-c_1$, **not** $c_6=c_1-c_2$ as in the on-shell case, so the TFF
sensitivity row is not in the span of the on-shell rows. Adding the three TFF
rows ($\pi^0,\eta,\eta'$ at $Q^2=5$ GeV$^2$) to the Jacobian in
`scripts/diagnose_radiative_rank.py` raises the rank from **7 to 8 of 11**.
Only the $\pi^0$ row adds rank; the $\eta$ and $\eta'$ rows lie in the span of
the OPE block, the ten radiative rows, and the $\pi^0$ TFF row. The TFF route
is **viable**: it lifts one of the four null directions inherited from the OPE
block. The remaining three are not lifted by the TFF.

The flavor factors are derived in [the TFF flavor record](radiative-tff-flavor.md).
The next step is to derive the normalized off-shell amplitude, including the WZW
and relevant VVP contributions, then write the predictor and turn the CLEO and
BaBar rows into `Observation`s. The cited RChT paper does not print this formula.

**Data already validated and checksum-pinned, with loaders:**

| file | rows | content |
|---|---|---|
| `data/odd_parity/cleo_tff.csv` | 63 | $\pi^0,\eta,\eta'$, $Q^2\in[1.64,15.3]$ GeV$^2$ [@CLEO:1997fho] |
| `data/odd_parity/babar_pi0_table11.csv` | 17 | $\pi^0$ only, $Q^2\in[4.2,34.4]$ GeV$^2$, stat and syst separate, asymmetric errors [@BaBar:2009rrj] |

Both have validated loaders in `scripts/implementation/odd_parity_data.py` but
**deliberately create no `Observation`**, because no physical predictor exists.
That is the gap to close.

**What remains to derive.** Steps 1 and 2 of the original plan are done
([record](radiative-tff-kinematics.md)): the off-shell slot assignment was
redone, the alternative contracted pairs all vanish (the uniqueness criterion
holds), and the $c_2,c_5,c_6$ fixed-ratio degeneracy is **broken**
($c_6=c_2-c_1$, not $c_1-c_2$). The remaining work:

3. Derive the fully normalized off-shell amplitude (WZW and resonance terms,
   including relevant VVP exchange); then write a predictor and connect the
   CLEO and BaBar rows as `Observation`s. The flavor trace itself is recorded in
   [radiative-tff-flavor.md](radiative-tff-flavor.md).

## Next step after that: $P\to\gamma\gamma$

$\pi^0\to\gamma\gamma$ and $\eta,\eta'\to\gamma\gamma$ have branching fractions
in the PDG and probe the axial anomaly. They use the same VJP basis with **two**
photons and one pseudoscalar, so the slot content differs and the kinematic
factors are new. This is the place where the four null directions are most likely
to be attacked, since two photon structures give two independent slot structures
rather than one. Note the caveat already recorded: the source's basis is declared
complete for one-pseudoscalar vertices, which $P\to\gamma\gamma$ satisfies, but
the two-photon operator content has not been checked here.

## Running things

```text
.venv/bin/python scripts/derive_radiative_sm_flavor_algebra.py
.venv/bin/python scripts/derive_radiative_sm_kinematics_probe.py
.venv/bin/python scripts/derive_radiative_tff_kinematics.py
.venv/bin/python scripts/diagnose_radiative_rank.py          # ~30 s
.venv/bin/python scripts/extract_odd_parity_tables.py
.venv/bin/python -m unittest scripts.implementation.test_odd_parity_data
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
- Fit a coefficient to force agreement with Eq. (27); the one real mismatch
  ($d_1$) is recorded as found.
- Treat a PDG average as independent data: $S=2.2$ on $\omega\to\pi^0\gamma$,
  $S=1.7$ on $\rho\to\pi^0\gamma$, $S=1.1$ on $\omega\to\eta\gamma$, and
  $\mathcal B(\eta'\to\rho^0\gamma)$ includes non-resonant $\pi^+\pi^-\gamma$.
