---
title: Intrinsic-parity selection rule for the cubic vertices in the model with a radial mode
status: derived
created: 2026-11-19
updated: 2026-11-19
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: analytical-derivation
---

# Cubic vertices of pseudoscalars and the ALP: intrinsic-parity selection rule

**Result.** In the polar coordinates of the model of DiLuzio:2026poh [@DiLuzio:2026poh] at $\sigma=0$, every vertex with an odd number of pseudoscalar and ALP legs vanishes if and only if the dressed source obeys
$$
\widehat M(-a,\sigma)=\widehat M(a,\sigma)^\dagger
\quad\Longleftrightarrow\quad
\widehat S_{rs}=0\ \text{for odd } s,\qquad \widehat P_{rs}=0\ \text{for even } s .
$$
The rule holds in the original coordinates and in the basis without ALP–meson mixing of the draft (`eq:murayama-chiral-shift`, `eq:murayama-T`). Then $\mathcal V^{(3)}\equiv0$ and the term $\sum\mathcal V^{(3)}G\,\mathcal V^{(3)}$ of the [radial-exchange record](finite-mixing/alp-three-pseudoscalar-radial-exchange.md) is absent; only the contact $\mathcal V^{(4)}$ and the $\sigma$ exchange remain. The condition is intrinsic parity of the ALP couplings, not flavour conservation.

## Conventions and units

Conventions of Alda:2026nrd [@Alda:2026nrd] and of the [radial-exchange record](finite-mixing/alp-three-pseudoscalar-radial-exchange.md): $\widehat M(a,\sigma)=\widehat S+i\widehat P=\sum_{r,s}(\widehat S_{rs}+i\widehat P_{rs})(\sigma/f_a)^r(a/f_a)^s$, $\Phi=\sum_kP_k\lambda_k$, $U=e^{i\Phi/f_\pi}$, $\chi=2B_0\widehat M$. Vertices omit the factor $i$. Masses, fields, $B_0$, $f_\pi$, $f_a$ and source coefficients in GeV.

## Assumptions

1. Tree level, LO $U(3)$, polar connection-free representative (assumptions 1–4 of the [radial-exchange record](finite-mixing/alp-three-pseudoscalar-radial-exchange.md)).
2. The anomalous singlet term $-\tfrac{m_1^2}{12}\langle\Phi\rangle^2$ does not depend on $a$ or $\sigma$ (assumption 3 there). **Not checked:** with an $a$-dependent topological source $\Theta(a)$ the rule additionally needs $\Theta$ odd in $a$.

## Inputs

- Polar Lagrangian and source expansion: [radial-exchange record](finite-mixing/alp-three-pseudoscalar-radial-exchange.md), [general-source radial record](../radial-mode/couplings-and-decays/finite-mixing/radial-decay-polar-cartesian-general-sources.md).
- Change of fields of the chiral shift (no quadratic term): [quartic-vertex record](finite-mixing/murayama-unmixed-basis-quartic-vertex.md), section "Check".
- Only-up coefficients: [only-up radial-decay record](../radial-mode/couplings-and-decays/finite-mixing/murayama-sigma-two-pseudoscalar-decays.md).
- Draft chiral shift and linear change of fields: `latex/murayama_matching.tex`, `eq:murayama-chiral-shift`, `eq:murayama-T`.

## Step 1: the $\mathbb Z_2$ transformation

At $\sigma=0$ the polar Lagrangian is
$$
\mathcal L=\frac{f_\pi^2}{4}\big\langle\partial_\mu U\partial^\mu U^\dagger+\chi(a)U^\dagger+\chi(a)^\dagger U\big\rangle-\frac{m_1^2}{12}\langle\Phi\rangle^2+\frac12(\partial a)^2 .
$$
Under $(\Phi,a)\to(-\Phi,-a)$, $U\to U^\dagger$. The kinetic trace is invariant by cyclicity, and $\langle\Phi\rangle^2$ and $(\partial a)^2$ are even. The source term becomes $\langle\chi(-a)U+\chi(-a)^\dagger U^\dagger\rangle$, equal to the original for every $U$ if and only if $\chi(-a)=\chi(a)^\dagger$. With $\sigma$ unchanged, the same holds at every $\sigma$, so the condition above makes the complete polar Lagrangian invariant under $(\Phi,a,\sigma)\to(-\Phi,-a,\sigma)$. Vertices with an odd number of $\Phi$ and $a$ legs vanish; $\sigma PP$ vertices are allowed.

## Step 2: symbolic check of the cubic terms (SymPy)

With noncommutative Hermitian symbols $\Phi$, $D=\partial_\mu\Phi$, $S_n$, $P_n$ (derivatives in $x=a/f_a$) and trace words identified up to cyclic order, the cubic part of $\langle\chi U^\dagger+\chi^\dagger U\rangle$ (with $f_\pi=1$, $\chi=\sum_nx^n(S_n+iP_n)$) is
$$
2x^3\langle S_3\rangle+2x^2\langle P_2\Phi\rangle-x\langle S_1\Phi^2\rangle-\tfrac13\langle P_0\Phi^3\rangle ,
$$
and it vanishes identically for $S_1=S_3=P_0=P_2=0$. The cubic part of $\langle\partial U\partial U^\dagger\rangle$ vanishes for any $\Phi$.

## Check: the rule in the basis without ALP–meson mixing

- **Basis without ALP–meson mixing.** The chiral shift of `eq:murayama-chiral-shift`, $U=g_RU'g_R$ with $g_R=\exp\big(\tfrac{iaB_0}{f_a}\mathbb M^{-1}_{m_a^2}(\widehat P_{01})\big)$, obeys $g_R(-a)=g_R(a)^\dagger$. Hence $(\Phi',a)\to(-\Phi',-a)$ gives $U\to g_R^\dagger U'^\dagger g_R^\dagger=U^\dagger$, and $\langle\Phi\rangle=\langle\Phi'\rangle+a\tfrac{2B_0f_\pi}{f_a}\langle\mathbb M^{-1}_{m_a^2}(\widehat P_{01})\rangle$ changes sign: the $\mathbb Z_2$ of Step 1 is the same transformation in the new fields. The change of fields of `eq:murayama-T` and the meson-only rotation are linear. The rule therefore holds in the basis without ALP–meson mixing, for external and internal legs. Equivalently, the shift has no quadratic term, so its cubic vertices are those of the original fields with the linear part of the change of fields inserted ([quartic-vertex record](finite-mixing/murayama-unmixed-basis-quartic-vertex.md), section "Check"). No coefficients of the transformed source are needed.

## Which sources obey the rule

| source | condition | $\mathcal V^{(3)}$ |
|---|---|---|
| $\operatorname{diag}[m_q(1+\sigma/f_a)^{n_q}e^{-in_qa/f_a}]$, any real $n_q$, $\theta=0$ (only-up: $n=(1,0,0)$) | holds (SymPy: all Taylor coefficients) | $0$ |
| $(1+\sigma/f_a)\,e^{iac/f_a}\,m\,e^{iac/f_a}$ with Hermitian, flavour off-diagonal $c$ (axial-type, flavour violating) | holds | $0$ |
| only-up plus $\widehat S_{01}=\epsilon\operatorname{diag}(1,-1,0)$ (flavour diagonal, allowed by $\langle\widehat S_{01}\rangle=0$) | violated | $\neq0$ |
| only-up plus $\widehat P_{02}=\epsilon\operatorname{diag}(0,1,-1)$ (flavour diagonal) | violated | $\neq0$ |

Flavour conservation is therefore neither necessary nor sufficient. In the linear parametrization of Alda:2026nrd the first two orders are $\widehat S_1=-i[c_q^{VS},m_q]$ and $\widehat P_2=-i[c_q^V,\{2c_q^P-c_q^A,m_q\}]-i\{c_q^A,[2c_q^S-c_q^V,m_q]\}$ ([finite-mixing record](finite-mixing/alp-finite-mixing-three-pseudoscalar-amplitudes.md)). Both vanish for flavour-diagonal couplings, so there flavour-diagonal couplings satisfy the rule to this order. Higher orders of that parametrization were not checked here.

## Verification

[Script](../../../scripts/derive_murayama_only_up_three_pseudoscalar_amplitude.py) and [unit tests](../../../scripts/test_derive_murayama_only_up_three_pseudoscalar_amplitude.py) (`SymbolicTest`, `CubicVertexTest`); `DraftParityTest` checks `eq:murayama-parity` and the cubic source terms `eq:murayama-cubic-parity-trace` (Fourier extraction of the cubic part of $\tfrac{f_\pi^2}{4}\langle\chi U^\dagger+\chi^\dagger U\rangle$ for random sources with $\widehat P_{00}\neq0$, relative residual below $10^{-12}$). Numerical checks use the exact polar Lagrangian of the [radial-exchange script](../../../scripts/derive_alp_three_pseudoscalar_radial_exchange.py) (eleven fields, all derivative terms). $m_u,m_d,m_s=0.0022,0.0047,0.095$ GeV, $B_0=2.6$ GeV, $f_\pi=0.0924$ GeV, $f_a=0.4$ GeV, $m_1^2=0.75$ GeV$^2$, $\lambda=0.8$, $\epsilon=0.002$ GeV. "max $|\mathcal V^{(3)}|$" runs over all 220 cubic vertices of $(a,P_0,\dots,P_8)$ at random momenta.

| check | result |
|---|---|
| SymPy: parity of the PQ-charge Taylor coefficients ($r\le2$, $s\le4$) | no violation |
| SymPy: cubic trace words; with $S_1=S_3=P_0=P_2=0$; kinetic term | as displayed; none; none |
| SymPy: $\log(e^Xe^Ye^X)$ has no even orders | exact |
| $g_R(-a)=g_R(a)^\dagger$ and $U(-\Phi',-a)=U(\Phi',a)^\dagger$, random sources (seed 2026), `test_chiral_shift_commutes_with_the_parity` | $<10^{-13}$ |
| max $|\mathcal V^{(3)}|$ / max $|\mathcal V_{\sigma PP}|$: $n=(1,0,0)$, $(0,1,0)$, $(2,0,1)$, axial $c=0.5\lambda_6+0.3\lambda_4$ | $4.4\times10^{-15}$, $4.3\times10^{-15}$, $2.9\times10^{-16}$, $2.9\times10^{-16}$ |
| same ratio with $\widehat P_{02}$ / $\widehat S_{01}$ added | $0.19$ / $0.38$ |
| max $|\mathcal V^{(3)}|$ in the basis without ALP–meson mixing, ten random triples of $(a',P')$ directions: the four parity-conserving sources / the two parity-violating sources (`test_basis_without_mixing_keeps_the_rule`) | $\le1.7\times10^{-14}$ / $\ge0.055$ |
| direct $\sum\mathcal V^{(3)}G\mathcal V^{(3)}$, only-up, eight final states incl. $\pi^0\pi^+\pi^-$, $\pi^0K^+K^-$, $\eta K^+K^-$, $\pi^0K^0\bar K^0$ | $\le10^{-60}$ (exactly zero vertices) |
| same with $\widehat P_{02}$ / $\widehat S_{01}$, neutral legs | $1.3\times10^{-5}$ / $-1.3\times10^{-5}$ against total $-2.75\times10^{-3}$ |

```sh
.venv/bin/python scripts/derive_murayama_only_up_three_pseudoscalar_amplitude.py
.venv/bin/python -m unittest scripts.test_derive_murayama_only_up_three_pseudoscalar_amplitude
```
