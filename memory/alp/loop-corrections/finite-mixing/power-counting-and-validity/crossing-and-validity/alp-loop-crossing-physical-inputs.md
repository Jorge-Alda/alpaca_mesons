---
title: Loop-to-tree crossing values of f_pi/f_a for ALP lines with physical inputs
status: numerical
created: 2026-11-21
updated: 2026-11-22
bibtex_keys: [Alda:2026nrd, Gasser:1983yg]
provenance: numerical-computation
---

# Loop-to-tree crossing values of $f_\pi/f_a$ for ALP lines with physical inputs

**Superseded crossing values.** The numerical values below describe the fixed coupling representative and unmixed approximation used here. The current finite-mixing reliability criterion compares the [projected one-loop size](../../three-meson-amplitudes/alp-three-meson-one-loop-estimate.md) with the exact on-shell tree amplitude through the [implicit crossing equation](alp-loop-crossing-exact-tree-condition.md#accepted-crossing-equation). Its [validity conditions](alp-loop-crossing-validity.md) apply to that comparison.

**Status note (2026-11-21):** this record uses a representative of the ALP couplings ($c_q^P=0$, fixed $c_q^A$) and off-shell flavour-basis entries in the unmixed approximation. The frame-invariant definition (pole-mass shift, dressed invariants only) is in the [invariant crossing record](alp-invariant-loop-crossing-meson-poles.md); class (i) here corresponds to its class A with $\widehat S_2=-4c^2m_q$ and $Y=0$.

Recomputation of the regime map (section L4) of the [one-loop power-counting record](../loop-counting/alp-loop-power-counting.md) with physical inputs, with and without the $\overline{\rm MS}$ logarithms. The earlier L4 numbers were illustrative.

## Conventions and units

Units GeV; $r=f_\pi/f_a$; notation of [@Alda:2026nrd]. Dressed frame, flavour basis and **unmixed approximation** (ALP and meson propagators diagonal in $(P_u,P_d,P_s,a)$, $m_0$ neglected in internal meson lines) exactly as in the [power-counting record](../loop-counting/alp-loop-power-counting.md), sections L2-L3. $B_0$ always denotes the low-energy constant [@Gasser:1983yg]; the dimensionless $\overline{\rm MS}$ bubble function is written $\bar B(s;m_1,m_2)$ here and called "the bubble function" in words, to avoid confusion with $B_0$.

- Tadpole function: $A_0(m^2)=\frac{m^2}{16\pi^2}\big[1-\ln\frac{m^2}{\mu^2}\big]$. Size estimate: $|A_0|\to m^2/(16\pi^2)$.
- Bubble function: $\bar B(s;m_1,m_2)=-\int_0^1dx\,\ln\frac{x m_1^2+(1-x)m_2^2-x(1-x)s-i0}{\mu^2}$ (same normalization as L3, $\operatorname{Im}\bar B=\pi\sqrt{1-4m^2/s}$ for equal masses). **Evaluation point (choice):** threshold of the intermediate state, $s=(m_1+m_2)^2$, where $\operatorname{Im}\bar B=0$ and
$$
\bar B_{\rm thr}=2-\frac{2\,[m_1\ln(m_1/\mu)+m_2\ln(m_2/\mu)]}{m_1+m_2}.
$$
Size estimate: $|\bar B|\to1$.
- Source convention: $c_q^P=0$; primary $c_q^A=c\,\mathbf 1$ ("universal"), second $c_u^A=c$, $c_d^A=c_s^A=0$ ("up only"); $c\in\{1,3\}$. $\kappa_q=c^A_q(c^A_q-2c^P_q)$, $\kappa^{(4)}_q=(c^A_q)^3(c^A_q-4c^P_q)$, $c^{AP}_q=c^A_q-c^P_q$.

## Inputs

| input | value | source |
|---|---|---|
| $m_{\pi^0}$ | $0.13498$ GeV | stated input |
| $f_\pi$ | $0.09207$ GeV | stated input, as in the [power-counting record](../loop-counting/alp-loop-power-counting.md) |
| $(m_u,m_d,m_s)$ primary | $(2.2,4.7,95)$ MeV | default of `charge_model_parameters` in `scripts/derive_murayama_only_up_three_pseudoscalar_amplitude.py` |
| $(m_u,m_d,m_s)$ alternative | $(2.16,4.70,95)$ MeV | $m_u,m_d$ of the [partial-wave unitarity record](../../../../decay-amplitudes/vertices-and-kinematics/kinematics-and-validity/finite-mixing/alp-lo-partial-wave-unitarity.md); $m_s$ added (assumption) |
| $B_0$ | $m_{\pi^0}^2/(m_u+m_d)=2.6405$ (alternative $2.6559$) GeV | LO relation [@Gasser:1983yg] |
| meson-line masses | $m_{qq}^2=2B_0m_q$, $m_{qr}^2=B_0(m_q+m_r)$ | LO, quark masses (not kaon masses) |
| $\mu$ | $0.77$ GeV | choice |

Primary set: $m_{uu}^2=0.01162$, $m_{dd}^2=0.02482$, $m_{ss}^2=0.50170$ GeV$^2$. For comparison $2m_{K^\pm}^2-m_{\pi^0}^2=0.46922$ GeV$^2$ ($m_{K^\pm}=0.49368$ GeV): $m_{ss}^2$ from $2B_0m_s$ is $7\%$ larger, which lowers the $s$-flavour crossings of (ii) by $3.3\%$ and of (iii) by $1.7\%$ relative to kaon-based inputs.

## Method

Each loop/tree ratio is a monomial in $r$; set $|{\rm ratio}|=1$ (closed forms; ratios from L2-L3 of the power-counting record):

- (i) ALP line in the meson entry $qq$, fixed (tuned) $m_a$: ratio $=2\kappa_q r^2|A_0(m_a^2)|/f_\pi^2$, so $r_*=f_\pi/\sqrt{2\kappa_q|A_0(m_a^2)|}$, size estimate $r_*=4\pi f_\pi/(\sqrt{2\kappa_q}\,m_a)$. Independent of $B_0$ and of the quark masses (unit-tested).
- (ii) $aP$ bubble in $P_qP_q\to P_qP_q$: ratio $=4(c^{AP}_q)^2r^2m_{qq}^2|\bar B|/(16\pi^2f_\pi^2)$ with $\bar B$ at $s=(m_a+m_{qq})^2$; $r_*=4\pi f_\pi/(2|c^{AP}_q|m_{qq}\sqrt{|\bar B|})$.
- (iii) $aa$ bubble: ratio $=4\kappa_q^2r^4m_{qq}^2|\bar B|/(16\pi^2f_\pi^2)$ with $\bar B$ at $s=4m_a^2$, $\bar B_{\rm thr}=2-\ln(m_a^2/\mu^2)$; $r_*=[16\pi^2f_\pi^2/(4\kappa_q^2m_{qq}^2|\bar B|)]^{1/4}$.
- (iv) ALP line in $M^2_{aa}$ relative to $m_a^2$: ratio $=V^L_{aaaa}|A_0(m_a^2)|/(2m_a^2)=4r^4\sum_q\kappa^{(4)}_qm_{qq}^2\,|1-\ln(m_a^2/\mu^2)|/(16\pi^2f_\pi^2)$; with the size estimate it is $V^L_{aaaa}/(32\pi^2)$, independent of $m_a$.
- (v) meson-only loops: ratio $m_P^2/(16\pi^2f_\pi^2)$ (times $|1-\ln(m_P^2/\mu^2)|$ for the tadpole, $|\bar B_{\rm thr}|$ for the $PP$ bubble); independent of $r$, no crossing.
- Untuned (i): $m_a^2=M^2_{aa}|_{\rm quark}=2r^2\sum_q\kappa_qm_{qq}^2$ (no bare or topological ALP mass). Size estimate $r_*=[16\pi^2f_\pi^2/(4\kappa_q\sum_r\kappa_rm_{rr}^2)]^{1/4}$; with the logarithm the smallest root of ratio $=1$ is found numerically (grid plus `scipy.optimize.brentq`).

$A_0(m_a^2)=0$ at $m_a=\sqrt e\,\mu=1.270$ GeV and $\bar B(4m_a^2)=0$ at $m_a=e\mu=2.093$ GeV; near these points the "full" crossings are large only because the subtracted loop accidentally vanishes. The remainder is a polynomial absorbed into counterterms of the next order, so these zeros are not physical suppressions.

## Results: crossing $r_*=f_\pi/f_a$ (primary quark masses, universal convention)

Columns: $m_a$ in GeV. Entries of (ii), (iii), (iv) with the size estimate are independent of $m_a$ (repeated). "S" = size estimate, "L" = $\overline{\rm MS}$ logarithms kept, $\mu=0.77$ GeV.

$c=1$:

| class | | 0.1 | 0.3 | 0.5 | 1.0 | 2.0 |
|---|---|---|---|---|---|---|
| (i) ALP line, meson entry ($uu$ = $ss$) | S | 8.18 | 2.73 | 1.64 | 0.818 | 0.409 |
| | L | 3.63 | 1.61 | 1.20 | 1.18 | 0.429 |
| (ii) $aP$ bubble, $u$ | S | 5.37 (all $m_a$) | | | | |
| | L | 2.19 | 2.55 | 2.91 | 3.88 | 8.60 |
| (ii) $aP$ bubble, $s$ | S | 0.817 (all $m_a$) | | | | |
| | L | 0.502 | 0.499 | 0.521 | 0.615 | 1.03 |
| (iii) $aa$ bubble, $u$ | S | 2.32 (all $m_a$) | | | | |
| | L | 1.48 | 1.65 | 1.78 | 2.10 | 4.22 |
| (iii) $aa$ bubble, $s$ | S | 0.904 (all $m_a$) | | | | |
| | L | 0.575 | 0.644 | 0.695 | 0.820 | 1.65 |
| (iv) ALP line in $M^2_{aa}$ | S | 0.888 (all $m_a$) | | | | |
| | L | 0.591 | 0.681 | 0.760 | 1.07 | 0.909 |

$c=3$: every crossing above is divided by exactly 3 (all ratios are homogeneous in $c$: $c^2r^2$ or $c^4r^4$). Explicitly, S: (i) $2.73, 0.909, 0.545, 0.273, 0.136$; (ii) $u$ $1.79$, $s$ $0.272$; (iii) $u$ $0.772$, $s$ $0.301$; (iv) $0.296$. L: (i) $1.21, 0.535, 0.400, 0.395, 0.143$; (ii) $u$ $0.730, 0.850, 0.969, 1.29, 2.87$; (ii) $s$ $0.167, 0.166, 0.174, 0.205, 0.342$; (iii) $u$ $0.492, 0.550, 0.594, 0.700, 1.41$; (iii) $s$ $0.192, 0.215, 0.232, 0.273, 0.549$; (iv) $0.197, 0.227, 0.253, 0.356, 0.303$.

**Up-only convention** ($c_u^A=c$): (i) $uu$, (ii) $u$, (iii) $u$ unchanged; (i) $ss$, (ii) $s$, (iii) $s$ vanish (no crossing); (iv) only $m_{uu}^2$ contributes and coincides with (iii) $u$ in the size estimate: $2.32$ ($c=1$), $0.772$ ($c=3$); with logarithms $1.54, 1.78, 1.98, 2.79, 2.37$ ($c=1$). Off-diagonal entries of (i) acquire the factor $2m_u/(m_u+m_q)$ and become quark-mass dependent.

**Untuned case (i)** ($u$ entry; $m_a$ at the crossing in brackets): universal $c=1$: S $0.888$ ($m_a=0.921$ GeV), L $1.345$; $c=3$: S $0.296$, L $0.448$. Up only $c=1$: S $2.32$ ($m_a=0.353$ GeV), L $1.74$; $c=3$: S $0.772$, L $0.581$. The illustrative value of the power-counting record was $0.90$ (kaon-based $m_{ss}^2$). With logarithms the universal crossing lies beyond the zero of $A_0$ ($m_a^2>e\mu^2$).

**Meson-only loops (v)** (no crossing; primary set), ratio S / tadpole L / $PP$ bubble L: $uu$ $0.0087/0.043/0.052$; $dd$ $0.019/0.077/0.096$; $ss$ $0.375/0.437/0.812$; $\pi^0$ $0.014/0.061/0.075$. The $s\bar s$ line is close to tree size by itself.

**Alternative quark masses** $(2.16,4.70,95)$ MeV: largest relative change of any crossing $0.8\%$ ((ii) $u$), $\le0.3\%$ for the other classes, $0$ for (i).

## Regime statement

With physical inputs and the universal convention, ALP-line loops reach tree size at $r_*\approx0.8/c$ (size estimate) or $0.5/c$-$0.6/c$ (logarithms kept, $m_a\lesssim1$ GeV). The first class to cross is:

- $m_a\lesssim1$ GeV: the $aP$ bubble with an $s\bar s$ line, $r_*=0.82/c$ (S) or $0.50/c$-$0.62/c$ (L). The size-estimate boundary with class (i) is $m_a=\sqrt2\,m_{ss}=1.00$ GeV, independent of $c$. Class (iv) and the $s$-flavour $aa$ bubble follow closely ($0.89/c$, $0.90/c$ in S).
- $m_a\gtrsim1$ GeV: the ALP line in a meson entry (i), $r_*=4\pi f_\pi/(\sqrt2\,c\,m_a)$, e.g. $0.41/c$ at $2$ GeV.

In the up-only convention the strange-flavour classes are absent: below $m_a=0.353$ GeV (S, independent of $c$) the $aa$ bubble with $u\bar u$ and class (iv) cross first ($2.32/c$); above it class (i) crosses first. With logarithms class (i) is first already from $m_a\approx0.3$ GeV ($1.61$ vs $1.65$ at $c=1$). For $f_\pi/f_a\gg1$ the expansion parameter is $c\,f_\pi/f_a$, as stated in the power-counting record; the crossings scale as $1/c$.

## Convention dependence (flagged)

- **Quark-mass set:** (ii), (iii), (iv) depend on it (through $m_{qq}^2$), at the sub-percent level here; (i) does not (exactly independent of $B_0$). The untuned (i) depends on it through $m_a^2$. $2B_0m_s$ vs $2m_K^2-m_\pi^2$ is a $7\%$ effect on $m_{ss}^2$.
- **$c_q^A$ vs $c_q^P$ representative:** (i), (iii), (iv) depend on $\kappa_q$ or $\kappa^{(4)}_q$, not on $c^{AP}_q$ alone, and vanish in the representative $c_q^A=0$. Only (ii) depends on $c^{AP}_q$ alone.
- **Flavour structure (universal vs up only):** removes the $s$-flavour classes and changes (iv) by a factor $2.6$.
- **$\mu$:** all "L" entries of (i), (iii), (iv) and of (ii) depend on $\mu$; the size-estimate entries do not. The full (i) and (iv) crossings are artificially large near $m_a=1.27$ GeV and (iii) near $m_a=2.09$ GeV.
- **Bubble-function evaluation point:** (ii) and (iii) "L" entries depend on the choice $s=(m_1+m_2)^2$; another $s$ changes $\operatorname{Re}\bar B$ and adds $\operatorname{Im}\bar B$.
- **Unmixed approximation** (all classes): valid when $|m_a^2-m_P^2|$ is large compared with the mixing entries, see the power-counting record.

## Verification

[`test_derive_alp_loop_crossing_table.py`](../../../../../../scripts/test_derive_alp_loop_crossing_table.py) (14 tests): threshold bubble function against numerical integration (`scipy.integrate.quad`) and the equal-mass formula; crossings of (i), (ii), (iii), untuned (i) against `regime_map` of the power-counting script; ratios of (i) and (iv) against `unmixed_closed_forms` with $A_0$ inserted, ratios of (ii), (iii) against `bubble_to_tree_ratios`; ratio $=1$ at every crossing (S and L, $c=1,3$) and smallest-root property for the untuned case; independence of (i) from $B_0$ and quark masses; up-only switch-off of the strange classes; (iv) equal to untuned (i) for universal couplings; LO relations for $B_0$ and $m_{qr}^2$.

Remark: `illustrative_regime_numbers` of the power-counting script substitutes $\sum\kappa m^2$ with $c^2$ also for the class "ALP line in $M^2_{aa}$", whose coefficient is $c^4$; harmless at its default $c=1$, wrong for $c\neq1$. The new script uses $\kappa^{(4)}_q=c^4$.

## Reproducible calculation

[`derive_alp_loop_crossing_table.py`](../../../../../../scripts/derive_alp_loop_crossing_table.py) (imports `loop_function_A0` from `derive_alp_loop_power_counting.py`).

```sh
.venv/bin/python scripts/derive_alp_loop_crossing_table.py
.venv/bin/python -m unittest discover -s scripts -p test_derive_alp_loop_crossing_table.py
```
