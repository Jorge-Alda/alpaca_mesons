---
title: PDG 2024 radiative and C-violating modes of eta prime 958
status: literature
created: 2026-10-04
updated: 2026-10-04
bibtex_keys: [ParticleDataGroup:2024cfk]
provenance: literature
---

# PDG 2024: $\eta'(958)$ Decay Modes Bearing on $P\to V\gamma$

## Why this extract exists

$P\to V\gamma$ is the same operator vertex as $V\to P\gamma$: the seven VJP
operators contain exactly one $V$, one $f_+$ and one $u$, and which of the two
mesonic fields is the decaying external state is a choice of external on-shell
substitution, not a different interaction. The question of whether
$\eta'\to\rho^0\gamma$ and $\eta'\to\omega\gamma$ belong in the radiative
coupling fit is therefore settled by the operator basis, not by a selection
rule. This extract supplies the experimental side.

## Source and extraction method

Cached at `memory/sources/ParticleDataGroup:2024cfk/`:

- `rpp2024-list-eta-prime-958.pdf`
- `rpp2024-list-eta-prime-958.txt`, SHA256
  `575e452ac8c76c107adef1926c64338acbe50cdf5cbdf91533394a24a3df4c38`

retrieved from `https://pdg.lbl.gov/2024/listings/rpp2024-list-eta-prime-958.pdf`
and converted with `pdftotext -layout`, matching the method used for the cached
$\rho(770)$ and $\omega(782)$ listings. Line numbers refer to the `.txt` file.

## The two channels

From the decay-mode table, lines 49--58:

| mode | PDG label | branching fraction | rate [MeV] | line |
|---|---|---|---|---|
| $\rho^0\gamma$ including non-resonant $\pi^+\pi^-\gamma$ | $\Gamma_2$ | $(29.5\pm0.4)\,$% | $0.0554\pm0.0019$ | 49--50 |
| $\rho^0\gamma$, resonant entry | $\Gamma_3$ | value carried inside $\Gamma_2$ | -- | 51 |
| $\omega\gamma$ | $\Gamma_5$ | $(2.52\pm0.07)\,$% | $0.00474\pm0.00020$ | 58 |

$\Gamma_2$ is defined to **include** the non-resonant $\pi^+\pi^-\gamma$
background, and $\Gamma_3$ is the resonant component. The two must not be
treated as independent data, and a genuine $V\to P\gamma$ coupling fit must use
the resonant part, so $\Gamma_3$ alone. The extracted table does not print a
separate fraction for $\Gamma_3$; obtaining it requires the release data files,
not this PDF.

The constrained fit (lines 108--140) reports 52 measurements and one constraint
for 9 parameters, $\chi^2=69.5$ for 44 degrees of freedom, and publishes a
correlation matrix. The rows relevant here are $x_2$ and $x_5$, with
$\delta x_2\delta x_5/\delta x_2\delta x_5=-2$ percent, and $x_2$ correlates at
$-25$ percent with $x_1$ and $-6$ percent with $x_4$.

## $C$- and $P$-violation labels: what the listing does and does not say

The listing carries a separate table headed "Charge conjugation (C), Parity (P),
Lepton family number (LF) violating modes", lines 86--98, containing
$\Gamma_{31}$--$\Gamma_{39}$: $\pi^+\pi^-$, $\pi^0\pi^0$, $\pi^0e^+e^-$,
$\pi^0\rho^0$, $\eta e^+e^-$, $3\gamma$, $\mu^+\mu^-\pi^0$, $\mu^+\mu^-\eta$,
$e\mu$. **Neither $\rho^0\gamma$ nor $\omega\gamma$ appears in that table.**

Two distinct things must not be conflated:

1. Their absence is consistent with $C$-conservation, and for
   $\eta'\to\rho^0\gamma$ an explicit check agrees: the $P$-wave final state has
   $C_\rho C_\gamma(-1)^L = (+1)(-1)(-1)=+1$, equal to $C_{\eta'}=+1$. No
   $C$-based argument excludes this channel.
2. That table is a list of **limits**. Its silence about a mode with a measured
   branching fraction is not by itself a proof of conservation, so the
   $C$-parity bookkeeping for $\eta'\to\omega\gamma$ should be settled against a
   standard reference before any $C$-based selection rule is applied to it. It is
   **not** needed here: the operator-basis argument already admits the channel.

## Consequences for the radiative channel list

- $\eta'\to\rho^0\gamma$ and $\eta'\to\omega\gamma$ are first-class members of
  the radiative class and must be carried in the channel list. The
  vector-initiated orientations $\rho\to\eta'\gamma$ and $\omega\to\eta'\gamma$
  are kinematically **closed** ($m_{\rho^0}=775.3$, $m_\omega=782.7$ MeV, both
  below $m_{\eta'}=957.8$ MeV), so the two orientations are not redundant with
  one another even though they share a flavor invariant.
- These are among the best-determined radiative numbers available, far better
  than the $\rho$ and $\omega$ entries, and they sit in a global PDG fit with a
  published correlation matrix. That is an advantage for a coupled fit and an
  obligation: the correlation matrix is available here, so the
  user-authorized diagonal-error approximation
  ([the dataset plan](../../rcht/odd-parity/odd-parity-global-fit-dataset-plan.md))
  is not needed for them.

## Not claimed

No likelihood or covariance matrix is derived here. Only the printed diagonal
fractions and the printed percentage correlation coefficients are recorded; the
full correlation matrix in the listing is transcribed by line reference, not
copied, and no bin-level or cross-observable covariance is inferred. The
non-resonant/resonant split of $\Gamma_2$ versus $\Gamma_3$ is **not** resolved
by this extract and is recorded as an open data task.