---
title: Why the held KLOE phi three-pion table cannot constrain the VVP couplings
status: derived
created: 2026-10-04
updated: 2026-10-04
bibtex_keys: [KLOE:2003kas, Niecknig:2012sj, GomezDumm:2013sib]
provenance: assumption
---

# Why the Held KLOE $\phi\to3\pi$ Table Cannot Constrain $d_1,d_2$

## Decision

Do **not** build a KLOE $\phi\to3\pi$ predictor for the current first release. The
tabulated material we hold cannot deliver an interpretable constraint on the
VVP couplings $d_1,\dots,d_4$ with the tree-level single-multiplet model fixed by
[the implementation decision record](odd-parity-rcht-implementation-decisions.md).
This is a scoping decision, recorded so the table is not re-litigated as "ready".

## What we actually hold

[`load_kloe_phi_3pi_dalitz_fit`](../../../scripts/implementation/odd_parity_data.py)
returns a **fit-parameter table, not a Dalitz density**. Ten rows, of which four
are blank by isospin/CPT constraint, leaving six numbers [@KLOE:2003kas]:

| row | value | $\pm$ stat | $\pm$ syst | unit | carries $d_i$ information? |
|---|---|---|---|---|---|
| `rho_plus_mass` | 775.8 | 0.5 | 0.3 | MeV | **no** |
| `rho_plus_width` | 143.9 | 1.3 | 1.1 | MeV | **no** |
| `direct_amplitude` | 0.78 | 0.09 | 0.13 | dimensionless | weakly |
| `direct_phase` | $-2.47$ | 0.08 | 0.08 | rad | weakly |
| `omega_pi_amplitude` | $7.1\times10^{-3}$ | 0.6 | 0.8 | dimensionless | no |
| `omega_pi_phase` | $-0.22$ | 0.11 | 0.04 | rad | no |

Three independent reasons this cannot be used as planned:

1. **Two of the six numbers are inputs, not predictions.** $\rho$ mass and width
   enter our model as caller-supplied parameters (the `MV`-type inputs of
   [`OPEParameters`](../../../scripts/implementation/odd_parity_rcht.py)). Fitting
   them would re-measure the $\rho$ resonance, not test any $c_i$ or $d_i$.
   A predictor must never return these two rows.
2. **The high-rank observable is absent.** The Dalitz density — roughly
   $1.8\times10^3$ efficiency-corrected binned KLOE points — is what would carry
   real information about the three-pseudoscalar amplitude. KLOE publishes no
   reusable bin contents in the cached source, and digitizing the published plot
   is forbidden by [the dataset plan](odd-parity-global-fit-dataset-plan.md).
   What remains is a four-number amplitude decomposition in KLOE's own basis.
3. **Our model is documented to misrepresent that basis.** The two numbers with
   any $d_i$ sensitivity are KLOE's `direct` amplitude and phase, i.e. its
   VVP contact term. A tree-level RChT amplitude lacks the crossed-channel
   three-particle cuts that the dispersive analysis shows are essential here
   [@Niecknig:2012sj].

## The crossed-channel problem, in the source's own words

Niecknig et al. solve $\phi\to3\pi$ with a once- or twice-subtracted Omnès
representation whose subtraction constants are the only free parameters, and
find crossed-channel rescattering is essential: dropping it degrades the fit
from $\chi^2/\mathrm{ndof}\simeq1.2$ to $\simeq2.0$ [@Niecknig:2012sj]. They warn
that isobar models, even with finite intermediate widths,

> "do not fully respect analyticity and unitarity constraints, still seem to be
> state of the art"

and that such models are

> "doomed to fail for the decays of $\omega$ or $\phi$ into three pions, where
> the influence of the $\rho$ resonance is already strongly felt ($\omega$) or
> even shows up in the form of resonance bands inside the Dalitz plot ($\phi$)."

Crucially, crossed-channel cuts do not merely rescale the $\rho$ peak: they
**shift and attenuate** it in $\phi\to3\pi$ and **shift and enhance** it in
$\omega\to3\pi$ [@Niecknig:2012sj]. A tree-level amplitude therefore mispredicts
the very shape KLOE measured. Relatedly, Gómez-Dumm et al. warn that replacing
$M^2-s$ by $M^2-s-iM\Gamma(s)$ inside a chiral framework double counts
$\mathrm{Im}A_P(s)$ and violates analyticity at $O(p^6)$ [@GomezDumm:2013sib].

## Consequences

- Connecting the four amplitude rows would produce a **basis-dependent,
  shape-biased** extraction presented with a small statistical error. Under the
  project's evidence discipline that must be labeled `conjectured` at best, and
  it would be indistinguishable from a real $d_i$ determination.
- This does **not** retire the channel. It becomes usable if and only if the
  amplitude is built dispersively — Omnès function plus crossed-channel
  integral, with the $\pi\pi$ P-wave input already available in this repository
  ([the K$\pi$/Omnès records](../../kpi/INDEX.md)). That work is scoped in
  [the dispersive treatment plan](dispersive-three-pion-treatment-plan.md)
  (6--10 weeks for one channel) and is explicitly **not** a first-release item.
  A cheaper staged alternative (c′), which takes the total width rather than the
  Dalitz shape, is specified there.
- The total partial width $\Gamma(\phi\to\pi^+\pi^-\pi^0)$ is a weaker and more
  robust alternative to the density, since a shift-and-attenuation of the $\rho$
  band perturbs an integral far less than a pointwise comparison. It is *not*
  currently tabulated in our CSV; the PDG branching fraction would be the
  natural source. Whether it adds information beyond what the fixed OPE prior
  already fixes in $d_1+8d_2$ and $d_3$ has **not** been checked and must be
  before it is used.

## Verification status

No predictor was written and no fit was run against this table. The loader and
its checksum-validated sidecar are unchanged and remain correct; they simply do
not produce an `Observation`. See
[the data-implementation record](odd-parity-data-implementation.md).