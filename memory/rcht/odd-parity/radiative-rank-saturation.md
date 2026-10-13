---
title: Rank and saturation measurement of the radiative channel set
status: numerical
created: 2026-10-04
updated: 2026-10-04
bibtex_keys: [Ruiz-Femenia:2003jdx, ParticleDataGroup:2024cfk, Aloni:2018vki, Niecknig:2012sj]
provenance: numerical-computation
---

# Rank and Saturation of the Radiative Channel Set

## Scope

This is the first half of the probe proposed in
[the saturation-probe record](radiative-channel-saturation-probe.md): the
**normalization-free** half. It answers whether the radiative channels carry
information at all, by pure linear algebra, before any fit is built. The other
half, the normalized residual $\sqrt{\chi^2/\text{ndof}}$, needs the common
normalization $\mathcal N$ and the absolute $c_3$ coefficient, which are still
open, and is **not** attempted here.

Reproduce with
[`scripts/diagnose_radiative_rank.py`](../../../scripts/diagnose_radiative_rank.py).

Inputs are local and verified: the per-channel flavor invariants from
[`derive_radiative_sm_flavor_algebra.py`](../../../scripts/derive_radiative_sm_flavor_algebra.py),
the per-operator kinematic factors of
[the width-kinematics record](radiative-width-kinematics.md), the five fixed OPE
rows from `scripts/implementation/odd_parity_rcht.py`, and PDG branching
fractions and total widths from the
[$\rho,\omega$](../../extracts/ParticleDataGroup:2024cfk/rho-omega-radiative-branchings.md)
and [$\eta'$](../../extracts/ParticleDataGroup:2024cfk/eta-prime-radiative-modes.md)
extracts.

## Why this half needs no normalization

Each channel is reduced to the sensitivity vector
$\mathbf s_V=\partial G_{VP\gamma}/\partial(c_1..c_7,d_1..d_4)$ and then
normalized to unit length. Span membership, rank, singular values and residual
cosines are all invariant under row rescaling, so the undetermined common
normalization $\mathcal N$, the open $c_3$ absolute normalization, and any fitted
value are irrelevant to every number below.

## Result

The ten channels are $\omega,\rho,\rho^+,\phi,\eta'\to\pi,\eta,\eta',\rho^0,\omega,K\gamma$.

**The fixed OPE block alone has rank 5 of 11**, deficiency 6, with null
directions

$$
(c_1,c_2,c_3,c_5,c_6),\quad (c_7),\quad (c_4,d_1,d_2),\quad (d_4).
$$

**Adding all ten radiative channels raises the rank only to 7 of 11.** With the
published Eq. (27) sign convention:

| added channel | residual after OPE + previous | verdict |
|---|---|---|
| $\omega\to\pi^0\gamma$ | $4.10\times10^{-1}$ | **informative** |
| $\rho\to\pi^0\gamma$ | $1.83\times10^{-2}$ | **informative, marginal** |
| $\rho\to\eta\gamma$ | $7.7\times10^{-15}$ | saturated |
| $\omega\to\eta\gamma$ | $9.0\times10^{-16}$ | saturated |
| $\eta'\to\rho^0\gamma$ | $6.8\times10^{-16}$ | saturated |
| $\eta'\to\omega\gamma$ | $7.7\times10^{-16}$ | saturated |
| $\phi\to\eta\gamma$ | $9.8\times10^{-16}$ | saturated |
| $K^{*0}\to K^0\gamma$ | $1.0\times10^{-15}$ | saturated |
| $K^{*+}\to K^+\gamma$ | $2.7\times10^{-15}$ | saturated |
| $\rho^+\to\pi^+\gamma$ | $2.7\times10^{-15}$ | saturated |

So **two of ten channels carry information the OPE block does not already
supply**, and the second of those two is marginal. The full-Jacobian singular
spectrum is

$$
8.07,\ 4.14,\ 3.22,\ 1.90,\ 1.19,\ 0.479,\ 0.454,\ 0,\ 0,\ 0,\ 0,
$$

with a clean gap after the seventh value and a **null cluster of four**.

### What stays unprobed

The right null vectors of the full Jacobian are, to numerical precision,

$$
(c_2,c_4,c_5,c_6,c_7)\ \text{threefold degenerate},\qquad (d_4).
$$

The degeneracy is **inherited from the fixed OPE block alone**, and this is a hard
statement, verified in the script rather than argued:

- the five OPE rows restricted to $(c_2,c_5,c_6)$ have rank 2, so one combination
  of those three is free before any channel is added;
- $c_4$ and $c_7$ appear in **no** OPE row and in **no** channel row;
- $d_4$ appears in **no** OPE row and in **no** channel row.

That accounts for exactly four free directions, matching the four observed null
vectors. So the radiative channels add **no** constraint on these combinations
beyond what the OPE relations already leave free. $c_4,c_7,d_4$ are structurally
unreachable because $\mathcal O^4_{VJP}=0$ at the vacuum,
$\mathcal O^7_{VJP}$ vanishes on shell and $\mathcal O^4_{VVP}=0$; they cannot be
reached by any radiative channel in this basis, at any mixing angle, mass or sign
convention.

Separately, the $c_6$ entry equals $c_1-c_2$ in every row, identically, because
$-2M^2=-(M^2-m^2)-(M^2+m^2)$. Verified to $1.1\times10^{-16}$ on
unit-normalized rows.

### Two channels with no exchange, and two with only $d_2$

From the flavor algebra: $F_Q^{\rm exch}=F_M^{\rm exch}=0$ for
$\rho^+\to\pi^+\gamma$, so that channel has no exchange at all, and
$F_Q^{\rm exch}=0$ with $F_M^{\rm exch}\propto(m_u-m_d)$ for
$\rho\to\pi^0\gamma$, so only $d_2$ can enter there, and only through isospin
breaking. Conversely $\omega\to\eta\gamma$ and $\eta'\to\omega\gamma$ have
$F_Q^{\rm exch}=0$, so they probe the exchange only through $d_2$.

### Robustness, and a real convention dependence

Row rescaling cannot affect the above, but the **sign convention** is a
component-wise flip on $c_2,c_5$ and therefore not a rescaling. Both conventions
were run:

| convention | total rank | informative channels |
|---|---|---|
| published Eq. (27) signs | 7 of 11 | $\omega\to\pi^0\gamma$, $\rho\to\pi^0\gamma$ |
| raw derivation signs | 8 of 11 | the same two, plus $\rho\to\eta\gamma$ |

So the exact count is convention dependent and must not be quoted without
saying which. The convention-independent statement is the one that matters: **at
most eight of eleven directions are ever reached, at least three are always
unprobed, and the majority of channels are saturated in both conventions.** This
is exactly the degeneracy the
[optimization note](odd-parity-coupling-prescriptions-and-direct-alp-vgamma.md)
already flagged as a property of the parameterization, now quantified for the
radiative class.

## The two-body $V\to V'P$ decays: relevant, but they cannot lift any degeneracy

$\phi\to\omega\pi^0$ and $\phi\to\rho\pi$ were tested and added to the Jacobian.
They are genuinely relevant, and the answer on degeneracies is negative for a
structural reason.

**They are inside the stated validity of the basis.** The VVP operators contain
two vectors and one pseudoscalar, which is exactly the field content of
$V\to V'P$ with **one** pseudoscalar, the regime the source declares complete
[@Ruiz-Femenia:2003jdx, source lines 431--433]. So these are direct tree-level
probes of $\mathcal O^1,\mathcal O^2,\mathcal O^3_{VVP}$, i.e. of $d_1,d_2,d_3$.

**There is no charge matrix.** The second *vector* occupies the slot that $Q$
occupies in $f_+=-2eQF$, so the invariant is
$\mathrm{tr}(\{T_V,T_{V'}\})T_P)$ for $d_1,d_3$ and its mass-weighted version
for $d_2$. The kinematic factor is $(q'\cdot k)=(M_{V'}^2-m_P^2)/2$, i.e. the
*final-state* vector mass, not the photon-momentum combination of the radiative
case.

**What they probe, numerically:**

| channel | $\mathrm{tr}(\{T_V,T_{V'}\})T_P)$ | probes |
|---|---|---|
| $\phi\to\omega\pi^0$ | $0$ exactly | $d_2$ alone |
| $\phi\to\eta\pi^0$ | $0$ exactly | $d_2$ alone |
| $\phi\to\rho^0\pi^0$ | $-0.195262$ | $d_1,d_2,d_3$ |
| $\phi\to\rho^+\pi^-$ | $-0.195262$ | $d_1,d_2,d_3$ |

$\phi$ and $\omega$ (and $\eta$) all lie in the span of the singlet and
$\lambda^8$, and products of that span are orthogonal to the isovector $\pi^0$.
So the charge-type invariant of $\phi\to\omega\pi^0$ vanishes **exactly**, and
that channel isolates $d_2$ with a single coupling. The neutral and charged
$\rho\pi$ modes give the *same* invariant, because the trace is cyclic and both
pick out the same flavor component of $T_\phi$.

**They cannot lift any of the four null directions.** Every $V\to V'P$ row has
zero $c$-sector content, and $d_4=0$ on shell because $\mathcal O^4_{VVP}$ forces
$q$ twice into the $\epsilon$. All four null vectors live in $(c_2,c_4,c_5,c_6,c_7)$
or in $d_4$, so such rows are automatically orthogonal to every one of them.
Measured: the dot products with the four null vectors are $1.6\times10^{-17}$ and
$1.9\times10^{-16}$, and

$$
\operatorname{rank}\bigl[\text{OPE} + 10\,\text{radiative} + 4\,(V\to V'P)\bigr]=7
=\operatorname{rank}\bigl[\text{OPE} + 10\,\text{radiative}\bigr] \quad \text{of 11}.
$$

**Zero added rank.** This is the same conclusion the mixing scan reached, and for
the same kind of reason: the degeneracy lives in the $c$-sector, and only
*c*-sector observables can touch it.

**But they do add something, in the $d$-sector.** The radiative exchange measures
only two combinations, $(d_1+d_3)$ and $2d_2$, and the $d_1$ coefficient carries
the unresolved sign discrepancy. $\phi\to\omega\pi^0$ isolates $d_2$ alone,
with a different kinematic function, no photon, no internal propagator and no
$c_i$ interference. Since the OPE matching fixes only $d_1+8d_2$ and $d_3$
individually, $d_2$ is the least-constrained coefficient in the whole
parameterization, and this is the cheapest clean handle on it. That is a real
gain, just not a gain in rank.

## Verdict against the falsification criteria

The criteria in [the probe record](radiative-channel-saturation-probe.md) require
a channel with normalized residual of order unity or larger, a separated
singular value, and residual correlation with the $\omega\to\pi^0\gamma$ residual
below roughly $0.9$. The first and third of those need $\mathcal N$ and are
**not yet testable**. What is established is the middle one, and negatively for
eight of ten channels: they add no separated singular value at all, because they
add no rank.

The radiative route is therefore **not** declared saturated and **not** declared
viable. It is declared to be carrying at most two additional directions here,
which is fewer than the six the OPE block leaves open, and which would have to be
weighed against the multi-body dispersive alternative
([the treatment plan](dispersive-three-pion-treatment-plan.md)).

## What would and would not break the result

Three of the four candidates were conjectures on my part. Two are now measured.

**Tested and it does NOT help: non-ideal mixing.** The $\omega$–$\phi$ angle
$\xi$ was scanned from $10^\circ$ to $45^\circ$ (ideal is $30^\circ$) and the
$\eta$–$\eta'$ angle from $0$ to $-35^\circ$:

| $\xi$ | $\theta$ | rank | informative |
|---|---|---|---|
| $30^\circ$ | $-19.47^\circ$ | 7 | 2 |
| $30^\circ$ | $0^\circ$ | 7 | 2 |
| $30^\circ$ | $-35^\circ$ | 7 | 2 |
| $10^\circ,20^\circ,40^\circ,45^\circ$ | $-19.47^\circ$ | 7 | 2 |

**Nothing changes.** The reason is structural: mixing enters only through the
per-channel scalar flavor invariants $F_Q,F_M,F_Q^{\rm exch},F_M^{\rm exch}$,
which rescale a channel's direct and exchange blocks without changing the internal
ratios that generate the degeneracy. So the earlier suggestion that non-ideal
mixing might help is **wrong and is withdrawn**. The same reasoning rules out a
different $M_V$ per state: it changes only the scalar kinematic prefactors per
channel.

**Still conjectured, not tested:**

- **singlet-specific or nonet-symmetric operators.** These would add structure the
  present basis has no slot for, and are the only candidate within the operator
  approach that can reach $c_4$ or $c_7$. Nothing here tests it.
- **dispersive treatment.** For the two-body radiative widths there is no
  intermediate propagator to modify, so what changes would be the amplitude and
  phase structure and the field-normalization residue, not a propagator. The
  propagator question belongs to the multi-body route, where the internal
  $\rho\pi$ line is genuinely replaced by a dispersive kernel
  ([propagator conventions](../../extracts/Ruiz-Femenia:2003jdx/timelike-propagator-mixing.md)).
- **transition form factors** $F_{P\gamma^*\gamma}(Q^2)$ and $P\to\gamma\gamma$,
  which carry the same operator basis but with a *space-like* off-shell photon
  and the pseudoscalar as initial state, hence different kinematic factors and
  genuinely different sensitivity rows. The TFF kinematics has now been derived
  ([the TFF-kinematics record](radiative-tff-kinematics.md)): the off-shell case
  gives $c_6=c_2-c_1$, breaking the on-shell $c_6=c_1-c_2$ identity, and the
  three TFF rows raise the rank from 7 to **8 of 11**. The TFF route is viable
  and lifts one of the four null directions. $P\to\gamma\gamma$ remains untested.
- **two-body $V\to V'P$ decays such as $\phi\to\omega\pi^0$ and
  $\phi\to\rho\pi$**: tested above, and they add **zero** rank. They cannot
  help, because they couple only to the $d_i$ while all four null directions
  live in the $c$-sector and in $d_4$.

## Not claimed

No $\chi^2$, no fit, no fitted coupling and no branching-fraction prediction is
reported here. The PDG rows used for the informative/saturated labels are the six
in the two local extracts; $\phi\to\eta\gamma$ and the two kaon channels have no
verified local branching fraction and are included in the rank analysis only
because rank does not need data. No covariance is used or implied: the
$S=2.2$ on $\omega\to\pi^0\gamma$, $S=1.7$ on $\rho\to\pi^0\gamma$ and $S=1.1$
on $\omega\to\eta\gamma$ mean the quoted errors understate the true spread, and
$\mathcal B(\eta'\to\rho^0\gamma)$ includes non-resonant $\pi^+\pi^-\gamma$, so it
is not a pure $V\to P\gamma$ signal.