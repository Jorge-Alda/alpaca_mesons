---
title: Symmetry selection rules and widths for neutral ALP decays to a pseudoscalar and photon, omega, or phi
status: mixed
created: 2026-09-24
updated: 2026-09-24
bibtex_keys: [Ecker:1988te, Ecker:1989yg, Alda:2026bar, ParticleDataGroup:2024cfk]
provenance: analytical-derivation
---

# Neutral $a\to PV$ Selection Rules and Widths

## Scope, conventions, and inputs

The executable derivation is [`scripts/derive_neutral_alp_pv_selection_rule_width.py`](../../scripts/derive_neutral_alp_pv_selection_rule_width.py). It uses natural units, the mostly-minus metric, and the convention $U=\exp(i\sqrt{2}\Phi/F)$ from [the RChT summary](../rcht/foundations/rcht-basics.md). At leading order, $F=f_\pi$ is set explicitly to

$$
f_\pi=\frac{93}{1000}\,\mathrm{GeV}.
$$

The derivative $WPa$ and mixing-induced $\gamma Pa$ operators are inputs from [the generic-generator expansion](../rcht/foundations/rcht-vpp-vpa-photon-pa-expansion.md). The meson and photon quantum numbers and physical masses are taken from [the PDG extract](../extracts/ParticleDataGroup:2024cfk/light-neutral-meson-quantum-numbers.md) [@ParticleDataGroup:2024cfk].

Assumptions are explicit: the ALP is a self-conjugate $J^{PC}=0^{-+}$ state and its strong/electromagnetic interactions conserve $C$ and $P$; $P$, the neutral vector generator $V$, and the photon-mixing generator $\mathcal T_\gamma$ are diagonal; the $3\times3$ ALP spurion $C_A$ is otherwise arbitrary. Generic pseudoscalar and vector mixing angles are retained. The massive-vector width assumes the narrow-state approximation, $m_a>m_P+m_X$, and $m_X>0$.

## State-level symmetry analysis

The PDG assignments are

$$
J^{PC}(\pi^0,\eta,\eta')=0^{-+},
\qquad
J^{PC}(\omega,\phi)=1^{--},
\qquad
C_\gamma=-1.
$$

Angular momentum conservation in $0\to0+1$ requires orbital angular momentum $L=1$. Parity is therefore conserved:

$$
P_f=P_PP_V(-1)^L=(-1)(-1)(-1)=-1=P_a.
$$

For two distinguishable neutral $C$ eigenstates there is no exchange factor, and

$$
C_f=C_PC_V=(+1)(-1)=-1\ne C_a=+1.
$$

Thus all six $a\to P\omega$ and $a\to P\phi$ channels are forbidden by $C$ conservation. Equivalently, $(CP)_a=-1$ while $(CP)_f=+1$, so they are also forbidden when $CP$ is conserved. This state-level conclusion is independent of ideal mixing and of the RChT truncation.

The real-photon channels obey a stronger statement. Lorentz covariance and the Ward identity leave no nonzero on-shell amplitude for a spin-zero particle decaying to another spin-zero particle and one real photon. Therefore $a\to P\gamma$ vanishes even if $C$ is explicitly violated. This does not apply to $a\to P\gamma^*$ with an off-shell photon.

## Generator conventions

The normalized diagonal basis is

$$
T_{\pi^0}=\frac{1}{\sqrt2}\operatorname{diag}(1,-1,0),\qquad
T_{\eta_8}=\frac{1}{\sqrt6}\operatorname{diag}(1,1,-2),\qquad
T_{\eta_0}=\frac{1}{\sqrt3}\mathbf1,
$$

with $\operatorname{Tr}(T_iT_j)=\delta_{ij}$. The mixing conventions are

$$
T_\eta=\cos\theta_P T_{\eta_8}-\sin\theta_P T_{\eta_0},
\qquad
T_{\eta'}=\sin\theta_P T_{\eta_8}+\cos\theta_P T_{\eta_0},
$$

$$
T_\omega=\cos\theta_V T_{\eta_0}+\sin\theta_V T_{\eta_8},
\qquad
T_\phi=\sin\theta_V T_{\eta_0}-\cos\theta_V T_{\eta_8}.
$$

At ideal mixing, $\cos\theta_V=\sqrt{2/3}$ and $\sin\theta_V=1/\sqrt3$, giving $T_\omega=\operatorname{diag}(1,1,0)/\sqrt2$ and the chosen phase $T_\phi=\operatorname{diag}(0,0,1)$.

## Flavor derivation

For arbitrary diagonal matrices $P$ and $D$ and an arbitrary spurion $C_A$, the script evaluates the full component matrices and obtains

$$
\operatorname{Tr}\!\left(D[P,C_A]\right)=0.
$$

Consequently, both tensor-origin structures relevant to diagonal $a\to P\omega/\phi$ and mixing-induced $a\to P\gamma$ vanish:

$$
\operatorname{Tr}\!\left(V[P,C_A]\right)=0,
\qquad
\operatorname{Tr}\!\left(\mathcal T_\gamma[P,C_A]\right)=0.
$$

Because $Q$ and $P$ are diagonal, the direct Goldstone-current structure also vanishes without assuming $[Q,C_A]=0$:

$$
\operatorname{Tr}\!\left([Q,C_A]P\right)=0.
$$

These are derived selection rules only for the stated diagonal external generators; they are not claims about charged or other nondiagonal channels.

In particular, at leading electromagnetic order the tensor-origin coefficient for a physical massive vector is

$$
C_{PV}^{(G_V)}
=-\frac{2iG_V}{f_\pi f_a}
\left\langle\frac{T_V}{M_V}[T_P,C_A]\right\rangle=0,
$$

while the mixing-induced photon coefficient is

$$
C_{P\gamma}^{(G_V)}
=-\frac{2iG_V}{f_\pi f_a}
\left\langle\mathcal T_\gamma[T_P,C_A]\right\rangle=0.
$$

The direct Goldstone-current coefficient is also zero. Hence the explicit $U(3)$ calculation reproduces the state-level selection rule for arbitrary $\eta$--$\eta'$ and diagonal $\omega$--$\phi$ mixing, even when $C_A$ itself is not diagonal.

## Lorentz amplitude and width before the flavor zero

Introduce a formal coefficient $C_X$ through

$$
\mathcal L_{PXa}=C_X X_{\mu\nu}\partial^\mu P\,\partial^\nu a.
$$

For $a(p)\to P(k)X(q,\epsilon)$ with $p=k+q$, direct contraction of the field strength gives, up to an irrelevant overall S-matrix phase,

$$
\mathcal M=C_X\left[(q\cdot k)(\epsilon^*\cdot p)
-(q\cdot p)(\epsilon^*\cdot k)\right].
$$

The replacement $\epsilon^*\to q$ gives zero identically. For an on-shell massive vector, $q^2=m_X^2$ and $q\cdot\epsilon=0$, so

$$
\mathcal M(a\to PX)=-C_Xm_X^2(\epsilon^*\cdot k).
$$

With

$$
\lambda=\lambda(m_a^2,m_P^2,m_X^2),
$$

the script derives

$$
\sum_\lambda|\epsilon^{(\lambda)}\cdot k|^2
=\frac{\lambda}{4m_X^2},
\qquad
\sum_\lambda|\mathcal M|^2
=\frac{|C_X|^2m_X^2\lambda}{4},
$$

and therefore

$$
\boxed{
\Gamma(a\to PX)
=\frac{|C_X|^2m_X^2}{64\pi m_a^3}
\lambda^{3/2}(m_a^2,m_P^2,m_X^2)
}.
$$

For a real photon, $q^2=0$, $q\cdot\epsilon=0$, $\epsilon\cdot p=\epsilon\cdot k$, and $q\cdot p=q\cdot k$. Hence

$$
\mathcal M(a\to P\gamma)=0
$$

even for a formal nonzero $C_\gamma$. This kinematic and gauge-invariant zero is logically independent of the diagonal flavor-trace zero.

## Phenomenological channel table

Using the PDG central masses [@ParticleDataGroup:2024cfk], the on-shell thresholds and widths in the $C$-conserving theory are

| channel | threshold | result |
|---|---:|---:|
| $a\to\pi^0\gamma$ | $134.9768\,\mathrm{MeV}$ | $\Gamma=0$ |
| $a\to\eta\gamma$ | $547.862\,\mathrm{MeV}$ | $\Gamma=0$ |
| $a\to\eta'\gamma$ | $957.78\,\mathrm{MeV}$ | $\Gamma=0$ |
| $a\to\pi^0\omega$ | $917.6368\,\mathrm{MeV}$ | $\Gamma=0$ |
| $a\to\eta\omega$ | $1330.522\,\mathrm{MeV}$ | $\Gamma=0$ |
| $a\to\eta'\omega$ | $1740.44\,\mathrm{MeV}$ | $\Gamma=0$ |
| $a\to\pi^0\phi$ | $1154.4378\,\mathrm{MeV}$ | $\Gamma=0$ |
| $a\to\eta\phi$ | $1567.323\,\mathrm{MeV}$ | $\Gamma=0$ |
| $a\to\eta'\phi$ | $1977.241\,\mathrm{MeV}$ | $\Gamma=0$ |

The thresholds are only kinematic conditions and do not remove the symmetry zero. Below an $\omega$ or $\phi$ threshold, an off-shell vector can contribute to a higher-multiplicity final state, which is a different process and requires its propagator and decay amplitude.

## Explicitly $C$-violating extension

If an additional $C$- and $CP$-violating operator gives a nonzero massive-vector coefficient $C_X$, the formal width derived above applies:

$$
\Gamma_{\not C}(a\to PX)
=\frac{|C_X|^2m_X^2}{64\pi m_a^3}
\lambda^{3/2}(m_a^2,m_P^2,m_X^2).
$$

This is not a prediction of the $C$-conserving RChT Lagrangian. Its channel dependence requires a specified $C$-violating flavor spurion. No analogous real-photon width is possible: electromagnetic gauge invariance still forces $\Gamma(a\to P\gamma)=0$.

## Evidence status and verification

- **Literature inputs:** the RChT conventions and operators linked above, and the PDG quantum numbers and masses [@ParticleDataGroup:2024cfk].
- **Derived:** the $C$, $P$, and $CP$ comparison; generator normalization; all diagonal flavor traces; the Ward identity; on-shell amplitudes; polarization sum; and massive two-body width.
- **Assumptions:** a self-conjugate $0^{-+}$ ALP, $C$- and $P$-conserving strong/electromagnetic interactions, diagonal neutral external generators, and the stated narrow-state/on-shell conditions.
- **Approximations:** none in the algebra. The input operator itself belongs to the previously recorded leading-order RChT truncation.

Run

```text
.venv/bin/python scripts/derive_neutral_alp_pv_selection_rule_width.py
```

The script prints each SymPy derivation step in LaTeX and ends with diagnostic assertions.
