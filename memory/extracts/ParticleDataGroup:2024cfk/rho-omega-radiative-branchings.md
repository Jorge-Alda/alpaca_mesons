---
title: PDG 2024 radiative decay branching fractions for rho(770) and omega(782)
status: literature
created: 2026-10-04
updated: 2026-10-04
bibtex_keys: [ParticleDataGroup:2024cfk]
provenance: literature
---

# PDG 2024 Radiative Branching Fractions: $\rho(770)$ and $\omega(782)$

## Source and extraction method

Cached at `memory/sources/ParticleDataGroup:2024cfk/`:

- `rpp2024-list-rho-770.pdf` / `.txt`
- `rpp2024-list-omega-782.pdf` / `.txt`

retrieved from `https://pdg.lbl.gov/2024/listings/` and converted with
`pdftotext -layout`. The converted omega value reproduces the independently
recorded `data/odd_parity/omega_pi0_gamma_pdg.csv` datum of
$0.0835\pm0.0027$ exactly, which validates the pipeline against a previously
extracted and SHA256-pinned number. Line numbers refer to the `.txt` files.

## $\rho(770)$ radiative modes

| mode | branching fraction | $S$ | line |
|---|---|---|---|
| $\pi^0\gamma$ | $(4.7\pm0.8)\times10^{-4}$ | 1.7 | 645 |
| $\eta\gamma$ | $(3.00\pm0.21)\times10^{-4}$ | -- | 646 |
| $\pi^0\pi^0\gamma$ | $(4.5\pm0.8)\times10^{-5}$ | -- | 647 |
| $\pi^+\pi^-\gamma$ | $(9.9\pm1.6)\times10^{-3}$ | -- | 643 |

Partial widths from the same listing: $\Gamma(\rho\to\pi^0\gamma)=0.070\pm0.012$
keV and $\Gamma(\rho\to\eta\gamma)=0.0447\pm0.0032$ keV, against
$\Gamma_\rho^{\rm total}=147.9\pm1.3$ MeV, reproducing the branching fractions
above.

The $S$ factor of **1.7** on $\pi^0\gamma$ signals that the measurements feeding
it do not agree well. The $\rho\to\pi^0\gamma$ width sub-table states that PDG
does **not** use the listed SND and ND points for averages, fits, or limits.

## $\omega(782)$ radiative modes

| mode | branching fraction | $S$ | line |
|---|---|---|---|
| $\pi^0\gamma$ | $(8.35\pm0.27)\%$ | 2.2 | 168 |
| $\eta\gamma$ | $(4.5\pm0.4)\times10^{-4}$ | 1.1 | 173 |
| $\pi^0\pi^0\gamma$ | $(6.7\pm1.1)\times10^{-5}$ | -- | 187 |
| $\pi^+\pi^-\gamma$ | $<3.6\times10^{-3}$ (95% CL) | -- | 179 |

## Validation against the live PDG website

The 2024 PDF extraction above was cross-checked against the current live
listing, `https://pdg.lbl.gov/2026/listings/rpp2026-list-{rho-770,omega-782}.pdf`,
converted the same way. Three of the four radiative branching fractions are
**identical** across editions, which validates the extraction:

| channel | 2024 (pinned) | 2026 (live) |
|---|---|---|
| $\rho\to\pi^0\gamma$ | $(4.7\pm0.8)\times10^{-4}$, $S=1.7$ | same |
| $\rho\to\eta\gamma$ | $(3.00\pm0.21)\times10^{-4}$ | same |
| $\rho\to\pi^0\pi^0\gamma$ | $(4.5\pm0.8)\times10^{-5}$ | same |
| $\omega\to\eta\gamma$ | $(4.5\pm0.4)\times10^{-4}$, $S=1.1$ | same |
| $\omega\to\pi^0\gamma$ | $(8.35\pm0.27)\%$, $S=2.2$ | **$(8.33\pm0.25)\%$, $S=2.1$** |

The single change is a genuine PDG update, not an extraction artifact. It has a
direct consequence: **the connected $\omega\to\pi^0\gamma$ datum in
[`omega_pi0_gamma_pdg.csv`](../../../data/odd_parity/omega_pi0_gamma_pdg.csv) is
the 2024-edition value** and must stay that way, because `bibliography.bib` pins
`ParticleDataGroup:2024cfk` and the sidecar records a SHA256. The shift is
$0.24$ percentage points against a quoted error of $0.27$, i.e. roughly one
third of the error, which is unremarkable at $S\approx2$ but is a concrete
demonstration that these averages are **not stable at the level of their quoted
errors**.

## Provenance chain to the original experiments

Consulting the live site pays off here: the PDG width sub-tables for
$\rho\to\pi^0\gamma$ and $\omega\to\eta\gamma$ contain **only excluded points**,
each under the banner "We do not use the following data for averages, fits,
limits, etc.". The tabulated branching fractions are therefore not derived from
the listed rows but from earlier determinations cited in their footnotes:

- $\rho\to\pi^0\gamma$ is obtained using $\Gamma_\rho^{\rm total}=147.9\pm1.3$
  MeV together with $\mathcal B(\rho\to\pi\gamma)$ **from Achasov 03**, i.e.
  SND [@Achasov:1999wp]. Held: `bibliography.bib` entry, cached source archive,
  and an existing extract,
  [`radiative-vector-data.md`](../Achasov:1999wp/radiative-vector-data.md).
  That extract is thin and states it must be reconciled with the PDG average
  before use as an input.
- $\omega\to\eta\gamma$ is obtained using $\Gamma_\omega=8.4\pm0.1$ MeV and
  $\mathcal B(\omega\to\eta\gamma)$ **from Dolinsky 89** (ND).
  **Not currently in `bibliography.bib`.** Using this branching fraction rests
  on the PDG citation alone until that measurement is added.

The $S$ factors are the quantitative statement of the same problem: $S=1.7$ on
$\rho\to\pi^0\gamma$, $S=2.1$--$2.2$ on $\omega\to\pi^0\gamma$.

## Usable channels for the coupling fit

Four channels in the $\rho,\omega$ sector have tabulated branching fractions:
$\rho\to\pi^0\gamma$, $\rho\to\eta\gamma$, $\omega\to\eta\gamma$, and the
already-connected $\omega\to\pi^0\gamma$.

Three caveats carry into any use:

1. **These are PDG averages, not independent measurements.** Following
   [the dataset plan](../../rcht/odd-parity/odd-parity-global-fit-dataset-plan.md),
   experiment labels belong in the provenance rather than being treated as
   independent data. Treating four PDG averages as four independent Gaussian
   constraints would overstate the information content.
2. **Consistency factors are not small.** $S=1.7$ on $\rho\to\pi^0\gamma$ and
   $S=2.1$--$2.2$ on $\omega\to\pi^0\gamma$ indicate the underlying measurements
   disagree. The quoted error is an average error and understates the true
   spread; this is a further reason not to read the four channels as four clean
   constraints. Where the width sub-table lists only excluded points, $S$ is the
   only available handle on that spread.
3. **No total-width renormalization is implied.** Each branching fraction is
   compared against the PDG total width of its own resonance, which is an
   independent input to the fit, not a fitted quantity.

## Not claimed

No likelihood, covariance, or correlated treatment is derived here. The
sidecar-based loader pattern in
[`odd_parity_data.py`](../../../scripts/implementation/odd_parity_data.py) must be
extended before any of these rows becomes an `Observation`, and the saturation
diagnostics of
[the probe record](../../rcht/odd-parity/radiative-channel-saturation-probe.md)
must be run before the route is declared viable.