---
title: Standard Model flavor algebra of the vector radiative channels
status: derived
created: 2026-10-04
updated: 2026-10-04
bibtex_keys: [Ruiz-Femenia:2003jdx, Aloni:2018vki]
provenance: analytical-derivation
---

# Flavor Algebra of $V\to P\gamma$ in the Standard Model

## Scope

This record establishes the flavor invariant reached by each radiative decay in
the class, for **both** orientations $V\to P\gamma$ and $P\to V\gamma$. This is the flavor half. The companion kinematic half, the factor each operator
contributes to the scalar $G$ of $\mathcal M=i\epsilon_{\alpha\beta\rho\sigma}
\epsilon_V^\alpha\epsilon_\gamma^\beta q^\rho k^\sigma G$, is in
[the width-kinematics record](radiative-width-kinematics.md); the two together
give the per-channel width.

Reproduce with
[`scripts/derive_radiative_sm_flavor_algebra.py`](../../../scripts/derive_radiative_sm_flavor_algebra.py),
which prints every table below and exercises every assertion.

## Why no ALP quantity appears

$V\to P\gamma$ is a purely Standard Model process. The final-state
pseudoscalar is a real meson, so the flavor algebra is built only from

- the meson state flavor matrices,
- the charge matrix $Q=\mathrm{diag}(2/3,-1/3,-1/3)$, which is what
  $f_+^{\mu\nu}=-2eQF^{\mu\nu}$ supplies at the Goldstone vacuum,
- the quark-mass matrix $\mathcal M=\mathrm{diag}(m_u,m_d,m_s)$, which is what
  the $\chi_-$ meson weight supplies.

$\mathbb M_a^{-1}(C)$ and $c_q^{AP}$ belong to the $aV\gamma$ vertex, i.e. to
the ALP in the final state, and must not be inserted here.

## Conventions

Octet fields enter as $\Phi=\lambda^a\phi^a/\sqrt2$ and
$V_{\mu\nu}=\lambda^aV^a_{\mu\nu}/\sqrt2$, and the singlet is added as
$\mathrm{diag}(1,1,1)/\sqrt3$. The source's own $\Phi$ matrix gives
$\Phi_{12}=\pi^+$, $\Phi_{13}=K^+$, $\Phi_{23}=K^0$, so the unit field of octet
adjoint $a$ is $\lambda^a/\sqrt2$ and of the singlet
$\mathrm{diag}(1,1,1)/\sqrt3$ [@Ruiz-Femenia:2003jdx, source lines 198--232].

**States are not rescaled to unit $\mathrm{Tr}(T^2)$.** Doing so silently
changes the $\eta$-$\eta'$ and $\omega$-$\phi$ mixing angles; the size of the
error is quantified below. Ideal mixing is imposed as a statement about
*states* built from those unit fields [@Ruiz-Femenia:2003jdx, source lines
784--797]:

$$
T_\omega=\sqrt{\tfrac23}\,T_{\eta_1}+\tfrac{1}{\sqrt3}T_{\eta_8}
=\tfrac16\mathrm{diag}(2\sqrt2+1,\,2\sqrt2+1,\,2\sqrt2-2),
$$

$$
T_\phi=-\tfrac{1}{\sqrt3}\,T_{\eta_1}+\sqrt{\tfrac23}\,T_{\eta_8}
=\tfrac16\mathrm{diag}(\sqrt2-2,\,\sqrt2-2,\,-2-2\sqrt2).
$$

The $\eta$-$\eta'$ rotation is one orthogonal rotation of the same unit fields,
$T_\eta=\cos\theta\,T_{\eta_8}+\sin\theta\,T_{\eta_1}$ and
$T_{\eta'}=-\sin\theta\,T_{\eta_8}+\cos\theta\,T_{\eta_1}$, evaluated at the
$|\sin\theta|=1/3$, $\cos\theta=2\sqrt2/3$ value of [@Aloni:2018vki]. The
rotated matrices are orthogonal only if $\sin^2\theta+\cos^2\theta=1$; the
script asserts this.

$\mathrm{Tr}(T_\omega^2)=5/6$ and $\mathrm{Tr}(T_\phi^2)=2/3$, while
$\mathrm{Tr}(T_{\pi^0}^2)=1$ and $\mathrm{Tr}(T_{\eta_8}^2)=1/2$. For
$\pi^\pm,K^0,K^+$ the unit field is a single off-diagonal element, so
$\mathrm{Tr}(T^2)=0$ and no unit-trace requirement can be imposed on them.

For the pseudoscalar and the vector of the same quark content the crossing
partner matrix is $T^\star=T^T$. This applies to the kaons and charged pions
only. It does **not** apply to nonet channels, where $T_V$ and $T_P$ are
different states that are not each other's crossing partners; there both
matrices are carried independently.

### $P\to V\gamma$ is the same vertex

The seven VJP operators contain exactly one $V$, one $f_+$ and one $u$. Which of
the two mesonic fields is the decaying external state is a choice of external
on-shell substitution, not a different interaction, so **$P\to V\gamma$ is
first-class in this basis**. It is therefore admitted here, and it is admitted for
a reason that does not rest on any symmetry argument: there is no selection rule
against it. The ordering of the trace becomes $\mathrm{tr}(\{T_P,Q\}T_V)$,
which equals $\mathrm{tr}(\{T_V,Q\}T_P)$ for the diagonal nonet states, so the
two orientations of a given meson pair share a flavor invariant. They are **not**
redundant as observables, because their kinematics differ and one orientation is
often closed where the other is open; see the openness table.

## Method

Every one of the seven VJP operators and the four VVP operators is a single
trace of a triple matrix product, with the fields in the order written in the
operator. The meson enters either as $u_\mu\Rightarrow-\partial_\mu\Phi$, whose
one-meson content is the state matrix $T_P$ with one momentum, or through
$\chi_-$, whose one-meson content is mass weighted. Expanding
$u=\exp(i\Phi/(\sqrt2F))$ gives

$$
\chi_-=4iB_0p-\frac{2i\sqrt2B_0}{F}\{\mathcal M,\Phi\},
$$

so the mass-weighted meson structure is the matrix $\{T_P,\mathcal M\}$. The two
flavor invariants are therefore

$$
F_Q^{VP}=\mathrm{tr}\bigl(\{T_V,Q\}\,T_P\bigr),\qquad
F_M^{VP}=\mathrm{tr}\bigl(\{T_V,Q\}\,\{T_P,\mathcal M\}\bigr),
$$

with the charge matrix between the vector and the meson, following the order
$\{V,f_+\}$ in the operators. $F_Q$ multiplies $c_1,c_2,c_5,c_6,c_7$ and
$d_1,d_3,d_4$; $F_M$ multiplies $c_3$ and $d_2$. For the exchange operators $Q$
is replaced by the internal octet vector matrix.

The $F_M$ value is a quark-mass combination, **not** a meson mass squared. Any
expression that displays $m_P^2$ in the $c_3$ or $d_2$ slot is in exact $SU(3)$,
where $F_M=2m F_Q$; see the limit below.

## Results

Evaluated at $\sin\theta=-1/3$, $\cos\theta=2\sqrt2/3$ and at
$m_u,m_d,m_s=2.2,4.2,20.4$ MeV for orientation:

| channel | $F_Q$ | $F_M$ [MeV] |
|---|---|---|
| $\omega\to\pi^0\gamma$ | $0.902369$ | $0.238414$ |
| $\rho\to\pi^0\gamma$ | $0.333333$ | $0.006144$ |
| $\rho^+\to\pi^+\gamma$ | $0.333333$ | $0.098310$ |
| $\omega\to\eta\gamma$ | $0.101728$ | $0.128137$ |
| $\rho\to\eta\gamma$ | $0.112735$ | $0.029786$ |
| $\phi\to\eta\gamma$ | $-0.400466$ | $-0.743289$ |
| $\omega\to\eta'\gamma$ | $0.240091$ | $-0.055876$ |
| $\rho\to\eta'\gamma$ | $0.905883$ | $0.239342$ |
| $\phi\to\eta'\gamma$ | $0.147089$ | $0.354175$ |
| $K^{*0}\to K^0\gamma$ | $-0.666667$ | $-0.755760$ |
| $K^{*+}\to K^+\gamma$ | $0.333333$ | $0.347158$ |
| $\eta'\to\rho^0\gamma$ | $0.905883$ | $0.239342$ |
| $\eta'\to\omega\gamma$ | $0.240091$ | $-0.055876$ |
| $\eta'\to\phi\gamma$ | $0.147089$ | $0.354175$ |

Exact closed forms for two:

$$
F_Q^{\omega\pi^0}=\frac{\sqrt2}{6}+\frac23,\qquad
F_M^{\omega\pi^0}=\frac{\sqrt2\,(1+2\sqrt2)}{9}\,(m_d+2m_u).
$$

**Kinematic openness.** The initial state must be heavier than the final-state
meson. With $m_{\eta'}=0.9578$ GeV:

- $\omega\to\eta'\gamma$ and $\rho\to\eta'\gamma$ are **closed**
  ($0.7827$ and $0.7753 < 0.9578$ GeV);
- $\eta'\to\rho^0\gamma$, $\eta'\to\omega\gamma$ are **open**;
- $\eta'\to\phi\gamma$ is **closed** ($1.0195 > 0.9578$ GeV), the only
  closed channel among the three pseudoscalar-initiated ones.

So of the eleven $V\to P\gamma$ channels nine are open, and all three
$\eta'\to V\gamma$ channels considered are handled, two open and one closed.
The script asserts the paired statements: $\eta'\to\rho^0\gamma$ and
$\eta'\to\omega\gamma$ are open exactly where their vector-initiated
counterparts are closed.

### Redundancy

Of the fourteen channels, **eleven have pairwise distinct $(F_Q,F_M)$**; the
three duplicates are precisely the orientation pairs
$\omega\leftrightarrow\eta'$, $\rho^0\leftrightarrow\eta'$ and
$\phi\leftrightarrow\eta'$, where the shared invariant follows from the
symmetry of the trace for diagonal states. **Flavor degeneracy here is not
observable degeneracy**: each such pair differs in external masses, in the
kinematic factors of [the width-kinematics record](radiative-width-kinematics.md),
and in phase space, and exactly one member of each pair is open.

Two further near-degeneracies are worth recording. $\rho\to\pi^0\gamma$ and
$\rho^+\to\pi^+\gamma$ share $F_Q=1/3$ but differ in $F_M$ at the $10\%$ level,
and $K^{*+}\to K^+\gamma$ shares $F_Q=1/3$ with both. Only in exact $SU(3)$ do
all channels collapse, because then $F_M=2mF_Q$ for **every** channel
identically, leaving a single flavor invariant. That collapse is precisely the
assumption under which the source publishes one channel, and it is why its result
is not predictive for the rest.

### The unit-trace trap

Rescaling each state matrix independently to unit $\mathrm{Tr}(T^2)$ is not a
harmless normalization; it changes the mixing angles. Singlet-to-octet weight
ratios $|a/b|$ in $T=aT_{\eta_1}+bT_{\eta_8}$:

| state | unit fields (used here) | unit-trace rescaling |
|---|---|---|
| $\omega$ | $\sqrt2$ | $1$ |
| $\phi$ | $\sqrt2/2$ | $1/2$ |
| $\eta$ | $\sqrt2/4$ | $1/4$ |
| $\eta'$ | $2\sqrt2$ | $2$ |

Every ratio is wrong by a factor $\sqrt2$ under that rescaling. Both $\eta$ and
$\eta'$ remain octet-dominant and singlet-dominant respectively, as the
observed mass ordering requires, so the qualitative statement survives while the
quantitative one does not, which is what makes the error easy to miss.

## Consequences for the fit

- No flavor-level degeneracy should be assumed among these channels. The
  [saturation probe](radiative-channel-saturation-probe.md) should therefore be
  re-scoped over all fourteen channels, both orientations, with the three
  orientation pairs counted once each in flavor but twice in kinematics.
- The $P\to V\gamma$ channels are the best-measured data in this class:
  $\mathcal B(\eta'\to\rho^0\gamma)=(29.5\pm0.4)\,\%$ and
  $\mathcal B(\eta'\to\omega\gamma)=(2.52\pm0.07)\,\%$, inside a global PDG
  fit with a published correlation matrix
  ([the eta' extract](../../extracts/ParticleDataGroup:2024cfk/eta-prime-radiative-modes.md)).
  For those rows the user-authorized diagonal-error approximation is not needed,
  and the $\Gamma_2$/$\Gamma_3$ resonant/non-resonant split must be resolved
  before use.
- The kaon channels are nonzero and are the only ones whose mass-weighted
  structure $\{T_P,\mathcal M\}$ is non-diagonal, so they probe $m_d+m_s$ and
  $m_u+m_s$ rather than the $u,d$ combinations. Nothing excludes them.

## Assumptions and limits

- Tree level, one resonance multiplet, Goldstone vacuum, $f_-^{\mu\nu}=0$ at the
  vacuum, electromagnetic-compatible sources.
- Ideal $\omega$-$\phi$ mixing and the $\eta$-$\eta'$ angle of
  [@Aloni:2018vki] are declared approximations; a single common parametrization
  of both nonets should replace them before a fit.
- Absolute normalizations of $F_Q$ and $F_M$ are convention dependent, fixed
  here by the source's field normalization. Only their ratios and their
  inequalities carry physical content.
- No $C$- or $P$-based selection rule has been applied to any channel. None is
  needed: every channel listed is in the operator basis on kinematic grounds
  alone.
- No isospin limit is imposed. $Q$ and $\mathcal M$ are kept fully general.
- Results are at the flavor level only. Equal flavor invariants would still not
  imply equal observables, since channels differ in kinematics and in which
  operators contribute; the three orientation pairs are the cases where this
  bites. The kinematic factors are in
  [the width-kinematics record](radiative-width-kinematics.md).