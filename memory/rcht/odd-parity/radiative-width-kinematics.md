---
title: Operator-level kinematic factors of V to P gamma in RChT
status: derived
created: 2026-10-04
updated: 2026-10-04
bibtex_keys: [Ruiz-Femenia:2003jdx, Ecker:1989yg, ParticleDataGroup:2024cfk]
provenance: analytical-derivation
---

# Operator-Level Kinematic Factors of $V\to P\gamma$

## Scope

This record fixes, operator by operator, the kinematic factor multiplying the
scalar $G$ in the radiative amplitude. Combined with the per-channel flavor
invariants of
[the SM flavor-algebra record](radiative-sm-flavor-algebra.md) it is what a
per-channel width expression needs. The published calibration channel is
$\omega\to\pi^0\gamma$ [@Ruiz-Femenia:2003jdx, Eq. (27)].

Reproduce with
[`scripts/derive_radiative_sm_kinematics_probe.py`](../../../scripts/derive_radiative_sm_kinematics_probe.py);
every assertion in this record is exercised there.

## Conventions

Natural units, metric $(+,-,-,-)$, $\epsilon^{0123}=+1$. External vector in the
Ecker normalization $V^{\mu\nu}=(q^\mu\epsilon^\nu-q^\nu\epsilon^\mu)/M_V^{\rm ext}$
with $q\cdot\epsilon=0$; photon massless,
$F^{\mu\nu}=r^\mu\epsilon_\gamma^\nu-r^\nu\epsilon_\gamma^\mu$, with $r^2=0$ and
$r\cdot\epsilon_\gamma=0$ [@Ecker:1989yg; @Ruiz-Femenia:2003jdx, source lines
823--824]. Meson momentum $k$, $k^2=m_P^2$, $q=k+r$.

$M_V$ is the **octet** internal mass, appearing in the propagator pole $1/M_V^2$
and in the $c_i/M_V$ prefactor; $M_V^{\rm ext}$ is the **external** on-shell mass.
They are not identified. Only the magnitude of a kinematic factor is physical
here, since the amplitude enters the width squared, so one overall sign is a
convention and is fixed against the published sign of $i\mathcal M^{\rm direct}$.

Contracting the lowered tensors with $\epsilon_{\mu\nu\rho\sigma}$ is equivalent
to the all-upper form by exactly one factor, $\prod_a g_{aa}=\det\eta=-1$,
because $\epsilon$ is fully antisymmetric and every slot index appears once. No
kinematic factor depends on that sign.

## The single invariant

Writing $S_{\mathrm{pub}}=\epsilon_{\alpha\beta\rho\sigma}
\epsilon^\alpha_V\epsilon^\beta_\gamma q^\rho k^\sigma$, **all** admissible
structures of this vertex are proportional to it. Reordering the $\epsilon$
indices changes only the sign, and the six orderings of $(q,r,k)$ in the three
momentum slots give
$S_{\mathrm{pub}}$, $-S_{\mathrm{pub}}$,
$\epsilon(q,\epsilon_V,r,k)=0$, and their negatives. In particular

$$
\epsilon(q,\epsilon_V,\epsilon_\gamma,k)=-S_{\mathrm{pub}},
$$

so a single reference structure suffices. The invariant space is
**one-dimensional**, and $S_{\mathrm{pub}}$ itself carries an explicit factor
$(M^{\rm ext\,2}-m_P^2)$:

$$
S_{\mathrm{pub}}\propto (M^{\rm ext\,2}-m_P^2)\,
\bigl(g_1v_y-g_2v_x\bigr).
$$

So the $m_P$ dependence lives partly in the structure and partly in $G$.

## The three kinematic relations

$$
q\cdot k=\frac{M^{\rm ext\,2}+m_P^2}{2},\qquad
r\cdot k=q\cdot r=\frac{M^{\rm ext\,2}-m_P^2}{2},\qquad
M^{\rm ext\,2}=(q\cdot k)+(r\cdot k),\qquad
m_P^2=(q\cdot k)-(r\cdot k).
$$

## Per-operator factors

Each operator has four free $\epsilon$ slots and one contracted index pair. The
$\{V,f_+\}$ anticommutator contributes a factor 2, because the two orderings of
the same triple matrix product are equal by cyclicity of the trace. The
surviving factors, in units of $S_{\mathrm{pub}}$ and up to one overall sign:

| operator | surviving contracted pair | factor | factor re-expressed |
|---|---|---|---|
| $\mathcal O^1_{VJP}$ | $(r\cdot k)$ | $-(M^{\rm ext\,2}-m_P^2)$ | $-2\,r\cdot k$ |
| $\mathcal O^2_{VJP}$ | $(q\cdot k)$ | $-(M^{\rm ext\,2}+m_P^2)$ | $-2\,q\cdot k$ |
| $\mathcal O^3_{VJP}$ | none | $2$ | pure number |
| $\mathcal O^5_{VJP}$ | $(r\cdot q)$ | $-(M^{\rm ext\,2}-m_P^2)$ | $-2\,r\cdot q$ |
| $\mathcal O^6_{VJP}$ | $(q^2)$ | $-2M^{\rm ext\,2}$ | $-2\,q^2$ |
| $\mathcal O^4_{VJP}$ | -- | $0$ | $f_-^{\mu\nu}=0$ |
| $\mathcal O^7_{VJP}$ | -- | $0$ | vanishes for the external on-shell $V^{\mu\nu}$ |

## Why the assignment is unique

Each operator offers a second way to form its contracted pair, and in every case
that alternative vanishes identically, so nothing is ambiguous:

- $\epsilon(q,\epsilon_V,r,k)=0$, because $q=r+k$ puts two identical momenta in
  the $\epsilon$;
- $\epsilon_V\cdot q=0$ by transversality of the external vector;
- keeping $\epsilon_\gamma\cdot k$ as a free-slot factor forces two identical
  momenta into the $\epsilon$.

This is a genuine selection rule rather than a preference, and it is what makes
the derivation reproduce the published combination.

## The on-shell step that makes $\chi_-$ work

$\chi_-$ is a meson wavefunction with **no Lorentz index**, so $\mathcal O^3$ and
$\mathcal O^2_{VVP}$ appear to have no final-state momentum and therefore no
phase-space suppression. The photon is on shell, so $r=q-k$ exactly, and

$$
\epsilon(q,\epsilon_V,r,\epsilon_\gamma)
=-\epsilon(q,\epsilon_V,k,\epsilon_\gamma)=+S_{\mathrm{pub}},
$$

the $q$-part of the substitution vanishing because $q$ would appear twice. The
structure is therefore the reference one and the explicit
$(M^{\rm ext\,2}-m_P^2)$ dependence is available to these operators.

The remaining $m_P^2$ in the $c_3$ and $d_2$ terms is **not** kinematic: it comes
from the quark-mass flavor factor $F_M=\mathrm{tr}(\{T_V,Q\}\{T_P,\mathcal M\})$,
which in the exact-$SU(3)$ limit satisfies $F_M=2mF_Q$ and $m\propto m_P^2$.
That the derived factor for $\mathcal O_3$ is a pure number is exactly what the
published $8m_\pi^2$ requires.

## Calibration against Eq. (27)

The published direct bracket is
$(c_2-c_1+c_5-2c_6)M_\omega^2+(c_1+c_2+8c_3-c_5)m_\pi^2$. Derived and published
magnitudes, in units of $S_{\mathrm{pub}}$:

| coupling | derived | published | agreement |
|---|---|---|---|
| $c_1$ | $-(M^2-m^2)$ | $m^2-M^2$ | exact |
| $c_2$ | $-(M^2+m^2)$ | $M^2+m^2$ | exact |
| $c_5$ | $-(M^2-m^2)$ | $M^2-m^2$ | exact |
| $c_6$ | $-2M^2$ | $-2M^2$ | exact |
| $c_3$ | $2$, times $F_M$ | $8m_\pi^2$, times $F_M$ | consistent |

All four reproduce up to one overall sign. The signs split as
$c_1,c_2,c_3,c_6$ against $c_5$, i.e. the single remaining freedom is the
**relative** sign of $\mathcal O^5$ with respect to $\mathcal O^1,\mathcal O^2,
\mathcal O^3$. That relative sign is physical, since the amplitude enters the
width squared, and it is fixed by the published result rather than by the
kinematics. It is convention dependent here only through the position of the
meson Lorentz index between $\nabla_\alpha u^\sigma$, $u^\nu$ and
$(\nabla_\alpha V^{\mu\nu})$.

## The four VVP exchange operators

The photon reaches the internal octet line through the even-parity mixing
$\mathcal L_2(V)=\tfrac{F_V}{2\sqrt2}\langle Vf_+\rangle$, so the internal vector
plays the role $f_+$ played in the direct diagrams. The VVP vertex has
$q=p+k$ with $k$ the meson, so the internal line carries the **photon** momentum
$r$ and the pole is $1/M_V^2$ with the octet mass. Slots: $\mu,\nu$ from the
external vector, $\rho$ from the internal line, $\sigma$ from the meson,
$\alpha$ the contracted pair. Applying the identical procedure:

| operator | surviving pair | factor in units of $S_{\rm pub}$ |
|---|---|---|
| $\mathcal O^1_{VVP}$ ($\to d_1$) | $(r\cdot k)$ | $-(M^2-m_P^2)$ |
| $\mathcal O^3_{VVP}$ ($\to d_3$) | $(q\cdot r)$ | $-(M^2-m^2)$ |
| $\mathcal O^2_{VVP}$ ($\to d_2$) | none | $2$ |
| $\mathcal O^4_{VVP}$ ($\to d_4$) | -- | $0$ |

Every alternative assignment vanishes identically, as for the direct operators.
$\mathcal O^4_{VVP}=0$ because its derivative index is the $\sigma$ slot, which
already carries $q$, forcing two identical momenta into the $\epsilon$; hence
$d_4$ never enters, exactly as the published bracket shows. Since
$r\cdot k=q\cdot r=(M^2-m_P^2)/2$, $\mathcal O^1_{VVP}$ and $\mathcal O^3_{VVP}$
carry the **same** factor.

So the derived exchange combination, before the common factor
$eF_V/(M_V^{\rm ext}M_V^2F)$, is

$$
(d_1+d_3)\bigl(M^{\rm ext\,2}-m_P^2\bigr)+2d_2 ,
$$

with $d_2$'s $m_P^2$ entering through $F_M$ as for $\mathcal O^3_{VJP}$.

**Comparison with the published $d_3M^2+(d_1+8d_2-d_3)m^2$.** $d_3$ agrees
exactly: derived $(M^2-m^2)$, published $M^2-m^2$. $d_2$ agrees in structure: a
pure number times the mass-weighted flavor factor, versus $8m_P^2$, the $8$
being the octet contraction $\dim 8$. **$d_1$ is a genuine residual
discrepancy**: the derivation gives it the same $(M^2-m_P^2)$ factor as $d_3$,
whereas the published bracket places $d_1$ only in the $m_P^2$ part. Since
$\mathcal O^1_{VVP}$ has exactly the slot structure of $\mathcal O^1_{VJP}$, and
$\mathcal O^1_{VJP}$ reproduces the published $c_1$, this is not explained by the
bookkeeping. It is reported, not tuned away, and it affects $d_1$ alone.

### Which channels exchange at all

The exchange flavor invariants use the internal **octet** matrix
$T_\rho=\lambda^3/\sqrt2$ in place of $Q$:

$$
F_Q^{\rm exch}=\mathrm{tr}(\{T_V,T_\rho\}T_P),\qquad
F_M^{\rm exh}=\mathrm{tr}(\{T_V,T_\rho\}\{T_P,\mathcal M\}).
$$

Two structural results follow, and both are predictions of the flavor algebra:

- For the **isovector** final states $\rho\to\pi^0\gamma$ and
  $\rho^+\to\pi^+\gamma$, $F_Q^{\rm exch}=0$, because
  $\{T_\rho,T_\rho\}=2T_\rho^2\propto\mathrm{diag}(1,1,0)$ is orthogonal to the
  isovector. So $d_1,d_3$ do not enter those channels at all. For
  $\rho^+\to\pi^+\gamma$ also $F_M^{\rm exch}=0$, i.e. **no exchange
  whatsoever**; for $\rho\to\pi^0\gamma$ $F_M^{\rm exch}\propto(m_u-m_d)$, pure
  isospin breaking, vanishing in exact $SU(3)$.
- The **kaon channels do exchange**: $F_Q^{\rm exch}=-\sqrt2/2$ and $+\sqrt2/2$
  for $K^{*0}$ and $K^{*+}$, because $\{E_{32},E_{23}\}$ contracts onto $E_{22}$.
  So the $\rho$-exchange is present for the kaons.

For $\omega,\phi,\eta'$ the exchange is likewise nonzero.

## Assembled width for any channel

With $\kappa_i$ the table above, $F_Q^{VP},F_M^{VP}$ the direct invariants of
[the SM flavor-algebra record](radiative-sm-flavor-algebra.md) and
$F_Q^{\rm exch},F_M^{\rm exch}$ the exchange ones,

$$
G_{VP\gamma}=\frac{\mathcal N}{M_V^{\rm ext}}
\Bigl[\kappa_{\rm dir}\,F_Q^{VP}\,(c_1+c_2+c_5+c_6)
+2c_3F_M^{VP}\Bigr]
+\frac{eF_V}{M_V^{\rm ext}M_V^2F}\Bigl[(d_1+d_3)F_Q^{\rm exh}
+2d_2F_M^{\rm exh}\Bigr],
$$

$$
\Gamma(V\to P\gamma)=\frac{|G_{VP\gamma}|^2\bigl(M_V^{{\rm ext}\,2}-m_P^2\bigr)^3}
{96\pi M_V^{{\rm ext}\3}} .
$$

$\mathcal N$ is an undetermined common normalization, fixed once by calibrating
on $\omega\to\pi^0\gamma$ against the PDG datum
[@ParticleDataGroup:2024cfk]; it multiplies the direct block only and does not
affect the ratio of two direct-only channels. The same formula holds for
$P\to V\gamma$ with the roles of $M_V^{\rm ext}$ and $m_P$ interchanged.

## Width

With $i\mathcal M=i\epsilon_{\alpha\beta\rho\sigma}\epsilon^\alpha_V
\epsilon^\beta_\gamma q^\rho k^\sigma G$,

$$
\Gamma(V\to P\gamma)=\frac{|G|^2\bigl(M_V^{{\rm ext}\,2}-m_P^2\bigr)^3}
{96\pi M_V^{{\rm ext}\,3}} .
$$

The exchange contribution carries the internal octet propagator. In this class
the internal vector carries the **photon** momentum, so $(q-k)^2=r^2=0$ exactly
and the pole is $1/M_V^2$ with the octet mass, not $1/(M_V^2-2m_P^2)$, matching
the published $1/M_V^2$ explicitly.

## Scope: same vertex, both orientations

$P\to V\gamma$ is the same operator vertex with the two mesonic fields exchanged
as external states, so all of the above holds verbatim once the roles of
$M_V^{\rm ext}$ and $m_P$ are interchanged. One thing does change: the external
on-shell substitution for the vector, since the vector is then a final-state
particle. The two orientations share a flavor invariant but differ in kinematics
and phase space, and one is usually closed where the other is open.

## What remains open

- The relative sign of $\mathcal O^5$ against the other three direct operators,
  and the absolute sign overall. Both are fixed by the published channel and are
  convention dependent in this record.
- The residual $d_1$ discrepancy in the exchange block, stated above. It is the
  only place where the derived combination does not reproduce the published
  grouping.
- The common normalization $\mathcal N$, and the assumption that the internal
  line is represented in the Ecker form with the octet mass. The latter is
  corroborated by $\mathcal O^3_{VVP}$ reproducing the published $d_3$ exactly, but
  it is an assumption about the off-shell internal vector, whose momentum is
  massless while its pole sits at $M_V$.
- Tree level, narrow on-shell external states, canonical $Z_V=1$, single
  resonance multiplet, one common $M_V$, ideal mixing, photon-vector mixing
  through $\mathcal L_2(V)$ only. Each is an assumption. No width prescription,
  resonance self-energy or fitted coupling is inferred.