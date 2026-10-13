---
title: Operator-level kinematic factors of the space-like TFF P to gamma* gamma
status: mixed
created: 2026-10-05
updated: 2026-10-05
bibtex_keys: [Ruiz-Femenia:2003jdx, CLEO:1997fho, BaBar:2009rrj]
provenance: analytical-derivation
---

# Operator-Level Kinematic Factors of the Space-like TFF $P\to\gamma^*\gamma$

## Scope

This record fixes, operator by operator, the kinematic factor multiplying the
scalar $G$ in the space-like transition form factor $F_{P\gamma^*\gamma}(Q^2)$,
where the pseudoscalar is the initial state and one photon is space-like
off-shell with invariant $Q^2$. It is the off-shell companion of
[the on-shell width-kinematics record](radiative-width-kinematics.md), and it
answers the decisive question posed in
[the next-steps record](next-steps.md): does the off-shell kinematics break the
on-shell fixed-ratio degeneracy in $(c_2,c_5,c_6)$?

Reproduce with
[`scripts/derive_radiative_tff_kinematics.py`](../../../scripts/derive_radiative_tff_kinematics.py);
every assertion in this record is exercised there. The rank consequence is
exercised in section 6 of
[`scripts/diagnose_radiative_rank.py`](../../../scripts/diagnose_radiative_rank.py).

## Conventions

Natural units, metric $(+,-,-,-)$, $\epsilon^{0123}=+1$. The "vector" in the
VJP operators is the **off-shell** photon (momentum $\ell$, polarization
$\epsilon_\ell$, $\ell^2=-Q^2$); the "photon" $f_+$ is the **real** photon
(momentum $k$, polarization $\epsilon_k$, $k^2=0$); the "pseudoscalar" $u$ is
$P$ (momentum $p$, $p^2=m_P^2$). The off-shell photon field tensor is
$V^{\mu\nu}=\ell^\mu\epsilon_\ell^\nu-\ell^\nu\epsilon_\ell^\mu$ with **no**
$1/M$ factor, because the off-shell photon is not a resonance with a pole mass.
The real photon field tensor is $F^{\mu\nu}=k^\mu\epsilon_k^\nu-k^\nu\epsilon_k^\mu$.

The on-shell relation $r=q-k$ does **not** hold: the off-shell photon has
$\ell^2\neq0$. Instead $p=k+\ell$.

### Kinematics

Pseudoscalar rest frame:

$$
p=(m_P,0,0,0),\qquad
k=\Bigl(\tfrac{m_P^2+Q^2}{2m_P},0,0,\tfrac{m_P^2+Q^2}{2m_P}\Bigr),
$$

$$
\ell=p-k=\Bigl(\tfrac{Q^2-m_P^2}{2m_P},0,0,-\tfrac{m_P^2+Q^2}{2m_P}\Bigr),
\qquad \ell^2=-Q^2.
$$

The reference structure is $S_{\rm TFF}=\epsilon(p,k,\epsilon_\ell,\epsilon_k)$,
which equals $\epsilon(\ell,k,\epsilon_\ell,\epsilon_k)$ because
$\ell=p-k$ and $\epsilon(k,k,\cdot,\cdot)=0$. All factors are reported in units
of $S_{\rm TFF}$.

## Per-operator factors

Each operator has four free $\epsilon$ slots and one contracted index pair. The
$\{V,f_+\}$ anticommutator contributes a factor 2. The slot mapping is
$q\to\ell$ (off-shell photon), $\epsilon_V\to\epsilon_\ell$, $r\to k$ (real
photon), $\epsilon_\gamma\to\epsilon_k$, $k\to p$ (pseudoscalar). The surviving
factors, in units of $S_{\rm TFF}$:

| operator | surviving contracted pair | factor |
|---|---|---|
| $\mathcal O^1_{VJP}$ | $(k\cdot p)$ | $Q^2+m_P^2$ |
| $\mathcal O^2_{VJP}$ | $(\ell\cdot p)$ | $m_P^2-Q^2$ |
| $\mathcal O^3_{VJP}$ | none | $-2$ |
| $\mathcal O^5_{VJP}$ | $(k\cdot\ell)$ | $Q^2+m_P^2$ |
| $\mathcal O^6_{VJP}$ | $(\ell^2)$ | $-2Q^2$ |
| $\mathcal O^4_{VJP}$ | -- | $0$ |
| $\mathcal O^7_{VJP}$ | -- | $0$ |

Every alternative contracted pair vanishes identically, exactly as in the
on-shell case, so the assignment is unique.

## The decisive output: the on-shell identity is broken

The on-shell rows all satisfy the identity $c_6=c_1-c_2$ (in the $c$-sector),
because $-2M^2=-(M^2-m^2)-(M^2+m^2)$. The off-shell factors give instead

$$
c_1-c_2=(Q^2+m_P^2)-(m_P^2-Q^2)=2Q^2,\qquad c_6=-2Q^2,
$$

so $c_6=c_2-c_1$, **not** $c_6=c_1-c_2$. The TFF sensitivity row is therefore
**not** in the span of the on-shell rows. This is the decisive output, and it
is the reason the TFF route was expected to lift the null directions.

## Rank consequence

Adding the three TFF sensitivity rows ($\pi^0,\eta,\eta'$ at $Q^2=5$ GeV$^2$,
flavor factors set to 1 as a placeholder) to the Jacobian in
[`scripts/diagnose_radiative_rank.py`](../../../scripts/diagnose_radiative_rank.py)
raises the rank from **7 of 11** to **8 of 11**. Only the $\pi^0$ TFF row adds
rank; the $\eta$ and $\eta'$ TFF rows lie in the span of the OPE block, the ten
radiative rows, and the $\pi^0$ TFF row. The singular value spectrum develops a
new separated value $0.373$ above the null cluster.

The TFF route is therefore **viable**: the off-shell kinematics lifts one of the
four null directions inherited from the OPE block. The remaining three null
directions are not lifted by the TFF, because the three TFF rows span only one
new direction in the $(c_2,c_5,c_6)$ subspace.

## Caveats and open items

- **Flavor factors are placeholders.** The TFF flavor structure differs from the
  on-shell case: the off-shell photon is a flavor singlet, so the invariant is
  not $\mathrm{tr}(\{T_V,Q\}T_P)$ with $T_V$ a meson matrix. The $t_2$ operator
  in the chiral Lagrangian [@Ruiz-Femenia:2003jdx, source lines 440--444] has
  the flavor structure $\mathrm{tr}(Q\{Q,T_P\})$, which is nonzero for $\pi^0$
  and must be derived for $\eta,\eta'$ before a fit. For the rank test the
  flavor factors are overall per-row factors and do not affect the result.
- **VVP contribution not included.** The $t_2$ operator also carries a $d_3$
  term [@Ruiz-Femenia:2003jdx, source lines 467--468]. The $d$-sector is zero in
  the TFF rows here; adding it would not change the $c$-sector rank.
- **Off-shell vector normalization.** The off-shell photon field tensor carries
  no $1/M$ factor, unlike the on-shell external vector. This is a convention
  that does not affect the rank but must be fixed before a fit.
- **Data not yet connected.** The CLEO and BaBar TFF data are checksum-pinned
  with loaders in `scripts/implementation/odd_parity_data.py` but create no
  `Observation`, because no physical predictor exists yet. That is the next
  step after the flavor structure is derived.

## Not claimed

No $\chi^2$, no fit, no fitted coupling and no TFF prediction is reported here.
The rank result is normalization-free and does not depend on the undetermined
common normalization $\mathcal N$, the open $c_3$ absolute normalization, or any
fitted value.
