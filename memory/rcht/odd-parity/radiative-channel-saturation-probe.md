---
title: Radiative-channel saturation probe and the rank-is-not-information criterion
status: conjectured
created: 2026-10-04
updated: 2026-10-04
bibtex_keys: [ParticleDataGroup:2024cfk, Ruiz-Femenia:2003jdx, BaBar:2004ytv]
provenance: assumption
---

# Radiative Channels: A Saturation Probe Before Full Implementation

## The problem this record addresses

Adding an observable to the fit adds a **rank direction**, not necessarily
**information**. The distinction has already bitten us once: the single
connected datum $\mathcal B(\omega\to\pi^0\gamma)$ raises the full Jacobian rank
from 5 to 6 of 11, which looks like progress, but its normalized residual is
$\approx0.045$. The datum is reproduced almost exactly by whatever the fixed OPE
block already selects, so it constrains nothing that the OPE prior has not
already fixed.

The same near-degeneracy would defeat the radiative route wholesale, and it is
cheap to test before we build anything.

## Structural advantage of this route

Independently of saturation, the radiative route is the only one that stays
inside the published operator basis. Ruiz--Femenia et al. state their VVP
operators are complete "for constructing vertices with only
**one-pseudoscalar**", and warn that "for a larger number of pseudoscalars
additional operators may emerge" [@Ruiz-Femenia:2003jdx]. $V\to P\gamma$ is a
one-pesudoscalar vertex, so it is in the regime the basis covers. By contrast
their published three-pion $\rho$-exchange amplitude depends only on
$d_1+8d_2-d_3$ and $d_3$, both already fixed individually by the short-distance
matching, so the three-pion route carries no free-$d_i$ information at all — see
[the VVP extract](../../extracts/Ruiz-Femenia:2003jdx/vvp-three-pseudoscalar.md)
and [the treatment plan](dispersive-three-pion-treatment-plan.md). This is a
prior-independent reason to try the radiative channels first.

## Why the risk is real

Every $V\to P\gamma$ width has the same structure: one VMD amplitude plus one
anomalous local term. In the
[coupling prescription record](odd-parity-coupling-prescriptions-and-direct-alp-vgamma.md)
the fitted objects are identifiable *combinations* of $c_i,d_i$, not individual
operator coefficients [@Ruiz-Femenia:2003jdx]. If the radiative widths are all
dominated by the same two combinations, they will be mutually consistent but
mutually uninformative — adding $N$ of them adds at most 1--2 effective rank,
not $N$.

## The proposal

Implement **two or three** radiative channels rather than the whole set, then
measure saturation before committing. Cheapest first, in order:

1. $\mathcal B(\rho\to\pi^0\gamma)$ and $\mathcal B(\omega\to\eta\gamma)$ —
   data is tabulated (see the
   [PDG extract](../../extracts/ParticleDataGroup:2024cfk/rho-omega-radiative-branchings.md)),
   but each channel needs an **operator-level derivation**, not a transcription.
2. $\mathcal B(\phi\to\pi^0\gamma)$ — highest $d_i$ leverage, and the source
   names it as needing flavor-resolved matching, so budget extra work.

**Scope correction, 2026-10-04.** An earlier draft of this record treated the
non-$\phi$ channels as close to mechanical. That was wrong, and the error is
recorded in the
[on-shell radiative extract](../../extracts/Ruiz-Femenia:2003jdx/on-shell-radiative-amplitudes.md):
Ruiz--Femenia et al. treat **exactly one** on-shell vector radiative decay,
$\omega\to\pi^0\gamma$. They state that extending to other channels "requires
that exact SU(3)$_{\mathrm V}$ symmetry is left aside ... keeping distinct masses
for each quark flavour, a rather non trivial task". The two channels they name,
$K^*\to K\gamma$ and $\phi\to\eta\gamma$, are given as *examples*, not as an
exhaustive exclusion list, so $\rho\to\pi^0\gamma$ and $\omega\to\eta\gamma$ are
covered by the same statement.

The requirement splits in two, and the distinction matters for cost:

- The **OPE-eliminated** numerical prediction would need fresh short-distance
  matching per channel. Not needed for the fit, which uses uneliminated
  couplings.
- The **uneliminated** amplitude for a new channel needs a genuine operator-level
  derivation with the right flavor trace and singlet/octet content. Tractable,
  but not transcription. For an ALP source $\mathbb M_a^{-1}(C)$ is a general
  Hermitian matrix, so the SU(3) arguments that would license a quick recipe are
  unavailable.

### Diagnostics to run

For each added channel, on the coupled problem including the fixed OPE block:

- the **normalized residual** $\sqrt{\chi^2_{\text{data}}/\text{ndof}}$ — values
  near zero mean the datum is saturated and carries no constraint;
- the **full-Jacobian singular value spectrum**, to see how many directions
  actually separate from the existing near-null ones;
- the **correlation between the new datum's residual and the existing
  $\omega\to\pi^0\gamma$ residual** — near $+1$ means it is the same
  information twice;
- $\chi^2$ improvement from adding the channel, compared with the $\chi^2$ change
  from *perturbing* the channel's own systematic nuisance.

### Falsification criteria

The radiative route is viable if at least one added channel shows a normalized
residual of order unity or larger, adds a singular value that is well separated
from the near-null cluster, and has residual correlation with the existing
datum below roughly 0.9. If all channels fail these, the radiative route is
reported as saturated and the multi-body dispersive treatment
([the treatment plan](dispersive-three-pion-treatment-plan.md)) becomes the only
remaining route.

## Status and honesty requirements

**Not yet run.** This record is a proposal. Until the probe exists:

- no claim may be made that the radiative route works or fails;
- the $\omega\to\pi^0\gamma$ saturation figure quoted above is a **numerical
  observation from a convergence probe**, not a published result and not a
  fitted $\chi^2$;
- the $O(p^6)$ analyticity caveat of [@GomezDumm:2013sib] concerns
  multi-body propagators and does **not** apply to the two-body
  $V\to P\gamma$ widths used in this probe. Two-body final states are the
  regime where our tree-level baseline is best controlled, which is the main
  reason to try them before (c).

## Consequence either way

A positive probe gives the first release its first unsaturated data, cheaply.
A negative probe is a genuine result: it establishes that no available two-body
dataset discriminates the free directions, which is exactly the evidence needed
to justify scoping (c) as its own project.

## Flavor-level prior, and its limits

[the flavor-algebra record](radiative-flavor-algebra.md) shows that the five
radiative channels share only **three** distinct source patterns, with
$\omega\to\pi^0\gamma\leftrightarrow\rho\to\eta\gamma$ and
$\rho\to\pi^0\gamma\leftrightarrow\omega\to\eta\gamma$ redundant at the flavor
level. All three are reached by $\phi\to\eta\gamma$ and by the crossed $\eta'$
channels as well, so the non-kaon radiative class spans exactly three patterns.
The kaon channels are excluded here because they give $D=M=0$ identically and
constrain nothing. This lowers the prior on the probe succeeding.

Two limits on that prior must be kept in view. First, the off-diagonal entries of
$\mathbb M_a^{-1}(C)$ are reachable only through $|\Delta S|=1$ amplitudes, and
no such data is loaded here — so the radiative class cannot be said to constrain
the full spurion, only its diagonal. Second, the flavor redundancy is not
observable redundancy: channels differ in kinematics and in which operators
contribute, so the Jacobian measurement is still the deciding test and has not
been done.