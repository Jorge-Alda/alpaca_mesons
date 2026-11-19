---
title: Only-up a_phys -> P_i P_j P_k amplitude with radial-mode exchange, channel by channel
status: derived
created: 2026-11-19
updated: 2026-11-19
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: analytical-derivation
---

# Only-up model: $a_{\rm phys}\to P_iP_jP_k$ with $\sigma$ exchange

**Result.** In the model of DiLuzio:2026poh [@DiLuzio:2026poh] with $(n_u,n_d,n_s)=(1,0,0)$, the exact tree-level LO amplitude in polar coordinates is
$$
\boxed{\mathcal M(a_{\rm phys}\to P_iP_jP_k)=\mathcal V^{(4)}_{aijk}+\sum_{(r|uv)}\frac{g_{ar}(s_{uv})\,g_{uv}(s_{uv})}{m_\sigma^2-s_{uv}}}
$$
with no exchange of mesons or of the ALP. The couplings, the contact and the surviving channels are given below for every final state.

## Conventions and units

Conventions of Alda:2026nrd [@Alda:2026nrd] and of the draft (`latex/murayama_matching.tex`): $\mathbb M_s$ (`eq:murayama-operator`), pole legs with meson matrix $\tfrac{2B_0f_\pi}{f_a}\mathbb M^{-1}_{m_r^2}(\widehat P_{01})/\sqrt{N(m_r^2)}$ and ALP component $1/\sqrt{N(m_r^2)}$ (`eq:murayama-pole-direction`), $N(m_r^2)$ (`eq:murayama-pole-norm`), $g_{ij}$ (`eq:murayama-polar-vertex`), $m_\sigma^2=B_0m_uf_\pi^2/f_a^2+2\lambda f_a^2$ (`eq:murayama-only-up-vacuum`). $(r|uv)\in\{(i|jk),(j|ik),(k|ij)\}$, $s_{uv}=(p_u+p_v)^2$. $g_{xy}(q^2)$ is the coupling $g_{xy}$ of the draft with $m_\sigma^2\to q^2$ (radial-exchange record, not a new quantity). Vertices omit the factor $i$; $S=1+iT$. GeV, GeV$^2$.

## Assumptions

Those of the [radial-exchange record](alp-three-pseudoscalar-radial-exchange.md): tree level, LO $U(3)$, no expansion in $f_\pi/f_a$, polar connection-free coordinates, stationary radial potential, stable poles, **no widths** (the formula needs $s_{uv}\neq m_\sigma^2$).

## Inputs

- General amplitude $\mathcal V^{(4)}+\sum\mathcal V^{(3)}G\mathcal V^{(3)}+\sum g\,g/(m_\sigma^2-s)$ and off-shell $g_{xy}(q^2)$: [radial-exchange record](alp-three-pseudoscalar-radial-exchange.md).
- $\mathcal V^{(3)}\equiv0$ for this model: [parity record](murayama-cubic-vertex-parity-rule.md).
- $\mathcal V^{(4)}$ in invariant form (original fields, original source coefficients) and its forms (4a)–(4c) for diagonal sources: [quartic-vertex record](murayama-unmixed-basis-quartic-vertex.md).
- Zero contact for four distinct neutral poles: [cancellation record](murayama-distinct-neutral-contact-cancellation.md).
- Only-up coefficients $\widehat S_{10}=m_ue_1$, $\widehat P_{01}=\widehat P_{11}=-m_ue_1$, $\widehat S_{02}=\widehat S_{12}=-\tfrac{m_u}{2}e_1$, $\widehat S_{20}=0$, pole equation `eq:murayama-only-up-pole`: [only-up radial-decay record](../shared/murayama-sigma-two-pseudoscalar-decays.md).

## Step 1: no meson or ALP exchange

The only-up source satisfies $\widehat M(-a,\sigma)=\widehat M(a,\sigma)^\dagger$, so all cubic vertices of pseudoscalars and the ALP vanish ([parity record](murayama-cubic-vertex-parity-rule.md)). The radial-exchange amplitude therefore loses its $\sum\mathcal V^{(3)}G\mathcal V^{(3)}$ term.

## Step 2: radial couplings

**Neutral legs with an ALP component** ($a_{\rm phys}$, $\pi^0$, $\eta$, $\eta'$). From the trace-form vertex of the radial-exchange record, $\big(\mathbb M^{-1}_{m_x^2}(\widehat P_{01})\big)_{uu}=-m_u\langle e_1\mathbb M^{-1}_{m_x^2}(e_1)\rangle$ with $e_1=\operatorname{diag}(1,0,0)$, and the pole equation, SymPy gives
$$
g_{xy}(q^2)=\frac{1}{f_a\sqrt{N(m_x^2)\,N(m_y^2)}}\Big[q^2-m_x^2-m_y^2+\frac{m_x^2\,m_y^2\,f_a^2}{B_0m_uf_\pi^2}\Big].
$$
It uses only the pole masses and normalizations. At $q^2=m_\sigma^2$ it is the on-shell $g_{xy}$ of the draft. As a check, $g_{xy}(q^2)$ is unchanged by the chiral shift to the basis without ALP–meson mixing, which does not touch $\sigma$ and has no quadratic term ([quartic-vertex record](murayama-unmixed-basis-quartic-vertex.md), section "Check").

**Legs without ALP component.** $\widehat S_{10}$, $\widehat P_{11}$, $\widehat S_{12}$ are proportional to $e_1$ and the neutral columns are diagonal, so
$$
g_{\pi^+\pi^-}=g_{K^+K^-}=\frac{B_0m_u}{f_a}\ \ (\text{independent of }q^2),\qquad g_{K^0\bar K^0}=0,
$$
and $g$ vanishes between a neutral leg and a charged or $K^0$-type leg.

## Step 3: channel by channel

| final state | $\mathcal V^{(4)}$ | $\sigma$ channels |
|---|---|---|
| three distinct neutral daughters ($\pi^0\eta\eta'$) | $0$ ([cancellation record](murayama-distinct-neutral-contact-cancellation.md)) | all three: $(i|jk)$, $(j|ik)$, $(k|ij)$ |
| neutral daughters with repetition ($3\pi^0$, $\pi^0\pi^0\eta$, ...) | form (4a) or (4b) | all three |
| $P^0\pi^+\pi^-$, $P^0K^+K^-$ ($P^0=\pi^0,\eta,\eta'$; $i=P^0$) | form (4c), pair $(u,d)$ or $(u,s)$ | only $(i|jk)$: $\dfrac{g_{ai}(s_{jk})}{m_\sigma^2-s_{jk}}\dfrac{B_0m_u}{f_a}$ |
| $P^0K^0\bar K^0$ | form (4c), pair $(d,s)$ | none |

Hence for three distinct neutral daughters
$$
\mathcal M(a_{\rm phys}\to\pi^0\eta\eta')=\sum_{(r|uv)}\frac{g_{ar}(s_{uv})\,g_{uv}(s_{uv})}{m_\sigma^2-s_{uv}},
$$
with all four couplings from the closed form of Step 2. For $P^0\pi^+\pi^-$,
$$
\mathcal M=\mathcal V^{(4)}_{aijk}+\frac{B_0m_u}{f_a}\,\frac{1}{f_a\sqrt{N(m_a^2)N(m_{P^0}^2)}}\,\frac{s_{+-}-m_a^2-m_{P^0}^2+m_a^2m_{P^0}^2f_a^2/(B_0m_uf_\pi^2)}{m_\sigma^2-s_{+-}} .
$$

Sample point ($f_a=0.4$ GeV, other inputs as in the [parity record](murayama-cubic-vertex-parity-rule.md), $s_{12}=0.3$, $s_{13}=0.5$ GeV$^2$, outside the physical region; algebraic check only):

| final state | $\mathcal M$ | $\mathcal V^{(4)}$ | $\sigma$ channels $(1|23)$, $(2|13)$, $(3|12)$ |
|---|---|---|---|
| $\pi^0\eta\eta'$ | $-2.764\times10^{-3}$ | $2\times10^{-17}$ | $-2.65\times10^{-3}$, $-1.43\times10^{-4}$, $2.78\times10^{-5}$ |
| $\pi^0\pi^+\pi^-$ | $6.796$ | $6.798$ | $-2.04\times10^{-3}$, $0$, $0$ |
| $\pi^0K^0\bar K^0$ | $0.4357$ | $0.4357$ | $0$, $0$, $0$ |

## Implementation

`scripts/implementation/only_up_decays.py`: `three_pseudoscalar(inputs, mesons, radial_mass_squared=None)` and `three_pseudoscalar_width(inputs, mesons, radial_mass_squared=None, *, epsrel=1e-5)`. With the default `None` the radial mode is omitted and the results are unchanged. With `radial_mass_squared` $=m_\sigma^2$ (GeV$^2$) the class `RadialExchange` adds the $\sigma$ term coherently. It uses the general trace-form vertex of the radial-exchange record on the same leg columns as the amplitude without $\sigma$ (`prepare_amplitude(...).external`, added to `calculate_alp_three_pseudoscalar_finite.py`), so the relative sign is fixed. The width raises `ValueError` if $m_\sigma^2$ (or a meson pole with nonzero residue) lies in the Dalitz range of a channel with nonzero couplings. For $\lambda>0$ the input must satisfy $m_\sigma^2>B_0m_uf_\pi^2/f_a^2$; this is not checked.

Sample widths ($\kappa_u=0.3$, $f_a=0.0104$ GeV, $m_a=0.711$ GeV): $\Gamma(\pi^0\pi^+\pi^-)=3.113\times10^{-4}$ GeV without $\sigma$ and $1.950\times10^{-4}$ GeV with $m_\sigma^2=2.5$ GeV$^2$. $\Gamma(3\pi^0)=4.80\times10^{-10}$ GeV without and $2.81\times10^{-6}$ GeV with $\sigma$. Here $f_a<f_\pi$, so the weak-mixing counting of the radial-exchange record (the $\sigma$ term suppressed by $(f_\pi/f_a)^2$) does not apply.

## Verification

[Script](../../../scripts/derive_murayama_only_up_three_pseudoscalar_amplitude.py), [unit tests](../../../scripts/test_derive_murayama_only_up_three_pseudoscalar_amplitude.py) (`OnlyUpAmplitudeTest`), [implementation tests](../../../scripts/implementation/test_only_up_decays.py) (`RadialExchangeTests`). "Direct" is the tree amplitude of the exact eleven-field polar Lagrangian with charged legs entered as complex directions. `DraftOnlyUpAmplitudeTest` evaluates `eq:murayama-only-up-offshell`, the relation $2B_0\big(\mathbb M^{-1}_{m_x^2}(\widehat P_{01})\big)_{uu}+1=m_x^2f_a^2/(B_0m_uf_\pi^2)$, `eq:murayama-only-up-vacuum`, `eq:murayama-only-up-cartesian-vertex` at $q^2=m_\sigma^2$ and `eq:murayama-only-up-three` with the channel table and the forms (4a), (4c) of the quartic-vertex record literally. The result agrees with the direct amplitude for all eight final states to the $10^{-9}$ precision of $N(m^2)$ near meson-like poles.

| check | residual |
|---|---|
| direct $\sigma$ exchange per channel vs trace-form $g\,g/(m_\sigma^2-s)$, eight final states (relative to $|\mathcal M|$) | $\le6.3\times10^{-16}$ |
| forbidden channels (charged pairs, $K^0\bar K^0$) | $\le5\times10^{-31}$ |
| closed-form $g_{xy}(q^2)$ vs trace form, all pairs of neutral poles, $q^2=0,0.77$ GeV$^2$ (relative) | $\le7\times10^{-10}$ (precision of $\mathbb M_{m^2}^{-1}$ near meson-like poles, $N\sim10^6$) |
| SymPy: closed form of $g_{xy}(q^2)$ | exact |
| $g_{\pi^+\pi^-}$, $g_{K^+K^-}$ vs $B_0m_u/f_a$ | $<10^{-14}$ |
| assembled formula (invariant $\mathcal V^{(4)}$ + trace-form $\sigma$ exchange) vs direct, eight final states (relative) | $\le3.9\times10^{-12}$ |
| implementation amplitude (with and without $\sigma$) vs direct, $\kappa_u=0.3$, $m_\sigma^2=2.5$ GeV$^2$, five final states incl. $\pi^0\pi^+\pi^-$, $\pi^0K^0_1K^0_1$ (relative) | $<10^{-9}$ |
| implementation width with $\sigma$ vs `scipy.dblquad` of the same $|\mathcal M|^2$ ($\pi^0\pi^+\pi^-$) | $<10^{-6}$ relative |
| $m_\sigma^2=10^8$ GeV$^2$ reproduces the width without $\sigma$ | $<10^{-6}$ relative |

```sh
.venv/bin/python scripts/derive_murayama_only_up_three_pseudoscalar_amplitude.py
.venv/bin/python -m unittest scripts.test_derive_murayama_only_up_three_pseudoscalar_amplitude
.venv/bin/python -m unittest scripts.implementation.test_only_up_decays
```
