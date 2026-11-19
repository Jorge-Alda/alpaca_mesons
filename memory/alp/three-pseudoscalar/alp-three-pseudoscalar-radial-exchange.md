---
title: Radial-mode exchange in the ALP decay to three pseudoscalars
status: mixed
created: 2026-11-18
updated: 2026-11-19
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: analytical-derivation
---

# Radial-mode exchange in $a_{\rm phys}\to P_iP_jP_k$

This record gives the exact tree-level LO amplitude of $a_{\rm phys}\to P_iP_jP_k$ in the model of DiLuzio:2026poh [@DiLuzio:2026poh] with general Hermitian source coefficients, including the exchange of the radial mode $\sigma$. The result is the amplitude of the [finite-mixing three-pseudoscalar record](alp-finite-mixing-three-pseudoscalar-amplitudes.md) plus three $\sigma$ exchanges written in invariant traces. The derivation is analytical (SymPy for the symbolic steps). The amplitude is also checked against a direct numerical tree calculation from the exact Lagrangian.

## Conventions and units

Natural units, mostly-minus metric, $S=1+iT$, conventions of Alda:2026nrd [@Alda:2026nrd]. Masses, fields, $B_0$, $f_\pi$, $f_a$, $\widehat S_{rs}$ and $\widehat P_{rs}$ are in GeV; $s_{ij}$ and squared masses in GeV$^2$. The source expansion, meson field, polar and Cartesian fields, radial potential and $\mathbb M_s$ are those of [the general-source radial-decay record](../shared/radial-decay-polar-cartesian-general-sources.md):

- $\widehat M(a,\sigma)=\widehat S+i\widehat P$, $\widehat S=\sum_{r,s}\widehat S_{rs}(\sigma/f_a)^r(a/f_a)^s$ and likewise for $\widehat P$; the first index counts $\sigma$, the second $a$.
- $\Phi=\sum_kP_k\lambda_k$, $\langle\lambda_k\lambda_l\rangle=2\delta_{kl}$, $U=e^{i\Phi/f_\pi}$, $\chi=2B_0\widehat M$, meson Lagrangian $\tfrac{f_\pi^2}{4}\langle u_\mu u^\mu+\chi_+\rangle-\tfrac{m_1^2}{12}\langle\Phi\rangle^2$.
- Polar kinetic term $\tfrac12(\partial\sigma)^2+\tfrac12(1+\sigma/f_a)^2(\partial a)^2$, $V_\varphi=\tfrac{m_\varphi^2}{2}(f_a+\sigma)^2+\tfrac\lambda4(f_a+\sigma)^4$; Cartesian fields from $f_a+s+ip=(f_a+\sigma)e^{ia/f_a}$.
- $\mathbb M_s(X)=B_0\{\widehat S_{00},X\}+\tfrac{m_1^2}{3}\langle X\rangle\mathbf 1-sX$; $N(m^2)=1+\tfrac{2B_0^2f_\pi^2}{f_a^2}\langle\mathbb M_{m^2}^{-1}(\widehat P_{01})^2\rangle$.
- Vertices omit the factor $i$. External momenta are all outgoing; the parent carries $-p$. With $s_{ij}=(p_i+p_j)^2$, $s_{ij}+s_{ik}+s_{jk}=m_a^2+m_i^2+m_j^2+m_k^2$.
- $\mathcal L\supset-\tfrac12\sigma\sum_{xy}g_{xy}P_xP_y$, as in [the only-up radial-decay record](../shared/murayama-sigma-two-pseudoscalar-decays.md).

## Assumptions

1. **No tadpoles.** $\widehat P_{00}=0$ and $\langle\widehat S_{01}\rangle=0$. The radial tadpole is removed by the stationarity condition $m_\varphi^2=B_0f_\pi^2\langle\widehat S_{10}\rangle/f_a^2-\lambda f_a^2$.
2. **No quadratic radial mixing.** $\widehat P_{10}=0$ and $\langle\widehat S_{11}\rangle=0$.
3. **Anomalous singlet term.** It does not depend on $\sigma$ or $a$, so $c_{GA}=0$ and $\Theta_n=0$ in the language of Alda:2026nrd.
4. **Order.** Tree level, leading $U(3)$ chiral order, no expansion in $f_\pi/f_a$. Use the polar representative without a radial connection.
5. **Poles.** External states are stable, simple poles. The trace forms require a nonzero angular (ALP-field) component of the leg. Exact meson poles with no ALP component are treated through their meson matrix, as in [the $K^0$ record](../shared/alp-lo-k0-doublet-exact-pole.md).
6. **No width (approximation).** The tree propagators are used without widths. This requires $s_{xy}\neq m_\sigma^2$ and $s_{xy}$ away from the poles of $G$; see the limits.

## Inputs

- Source expansion, stationarity, $m_\sigma^2$, polar cubic Lagrangian, pole columns: [general-source radial-decay record](../shared/radial-decay-polar-cartesian-general-sources.md).
- On-shell $g_{ij}$, only-up coefficients, $N(m^2)$, pole-aligned change of fields: [only-up radial-decay record](../shared/murayama-sigma-two-pseudoscalar-decays.md).
- Amplitude without $\sigma$: $\mathcal V^{(4)}$, $\mathcal V^{(3)}$, $G(s)=(M^{2\prime}-sK')^{-1}$, step-2 columns: [finite-mixing three-pseudoscalar record](alp-finite-mixing-three-pseudoscalar-amplitudes.md).
- Equality of on-shell amplitudes in polar and Cartesian coordinates: [polar and Cartesian record](../shared/alp-polar-cartesian-amplitude-equivalence.md).
- Model: DiLuzio:2026poh, [model extract](../../extracts/DiLuzio:2026poh/model-lagrangian-and-field-normalization.md).

## Choice of polar coordinates

The complete on-shell amplitude is the same in polar and Cartesian coordinates ([polar and Cartesian record](../shared/alp-polar-cartesian-amplitude-equivalence.md); checked numerically below). The polar coordinates are used for two reasons.

**(1) The part without $\sigma$ is an existing result.** At $\sigma=0$ the polar Lagrangian is the LO chiral Lagrangian with the $a$-dependent source $\widehat M(a,0)$, no derivative couplings of $a$, the kinetic term $\tfrac12(\partial a)^2$ and no bare ALP mass. This is the Lagrangian of the finite-mixing three-pseudoscalar record in its derivative-source-free coordinates, with the substitutions of Step 2. In Cartesian coordinates $a=p-sp/f_a-p^3/(3f_a^2)+\dots$ and $\sigma=s+p^2/(2f_a)+\dots$. SymPy gives the Cartesian minus polar quartic contacts at $s=0$ (with $a\to p$):
$$
\Delta\mathcal L=\Big[\frac{B_0f_\pi^2}{24f_a^4}\big(-16\langle\widehat S_{02}\rangle-3\langle\widehat S_{10}\rangle+12\langle\widehat S_{12}\rangle+6\langle\widehat S_{20}\rangle\big)-\frac\lambda4\Big]p^4
+\frac{B_0f_\pi}{6f_a^3}\langle(3\widehat P_{11}-2\widehat P_{01})\Phi\rangle p^3-\frac{B_0}{4f_a^2}\langle\widehat S_{10}\Phi^2\rangle p^2 .
$$
The Cartesian cubic vertices also differ ([general-source record](../shared/radial-decay-polar-cartesian-general-sources.md), Step 3). Hence in Cartesian coordinates the contact and the exchanges without $s$ are not those of the existing record, and the separation into "no $\sigma$" and "$\sigma$" parts is different.

**(2) The polar $\sigma$ term decouples as $1/m_\sigma^2$.** The polar $\sigma$ vertices do not contain $\lambda$. Since $m_\sigma^2=2\lambda f_a^2+\dots$, the polar $\sigma$ exchanges fall as $1/m_\sigma^2$. The Cartesian $sp^2$ vertex contains $-\lambda f_a$ and the contact contains $-\tfrac\lambda4p^4$. Each of these Cartesian pieces grows linearly with $\lambda$, and only their sum stays finite. Seed 2026 of the test, with the three daughters with the largest ALP component, gives:

| $\lambda$ | $m_\sigma^2$ (GeV$^2$) | polar contact | polar $\sigma$ term | $m_\sigma^2\times$ polar $\sigma$ term | Cartesian contact | Cartesian $s$ exchange | total (both) |
|---|---|---|---|---|---|---|---|
| 0.8 | 0.263 | 2.1143 | $7.35\times10^{-2}$ | 0.0193 | 2.1897 | $-1.89\times10^{-3}$ | 2.19140 |
| 8 | 2.567 | 2.1143 | $6.64\times10^{-2}$ | 0.170 | 1.9366 | 0.2441 | 2.18427 |
| 80 | 25.61 | 2.1143 | $7.45\times10^{-3}$ | 0.1908 | $-0.594$ | 2.716 | 2.12536 |
| 800 | 256.0 | 2.1143 | $7.72\times10^{-4}$ | 0.1975 | $-25.90$ | 28.02 | 2.11868 |
| 8000 | 2560 | 2.1143 | $7.74\times10^{-5}$ | 0.1983 | $-279.0$ | 281.1 | 2.11799 |

The Cartesian exchange of the angular and meson fields is $3.63\times10^{-3}$ at every $\lambda$. The polar and Cartesian totals agree to $\le4\times10^{-14}$ (relative) at every $\lambda$. At $\lambda=8000$ the two Cartesian pieces are each about 130 times the amplitude.

## Diagrams

At tree level a four-point amplitude has one contact or one propagator.

| diagram | status | reason |
|---|---|---|
| contact $\mathcal V^{(4)}_{aijk}$ | included | no $\sigma$ leg |
| exchange of the nine mesons and the ALP in the three channels, full $G(s)$ | included | as in the existing record |
| $\sigma$ exchange in the three channels | included | Step 3 |
| contacts or vertices with external $\sigma$ | not present | $\sigma$ appears only as an internal line |
| $\sigma$ tadpole | removed | stationarity (assumption 1) |
| $a$, $\Phi$ tadpoles | absent | assumption 1 |
| $\sigma$–$a$, $\sigma$–$\Phi$ mixing on a propagator | absent | SymPy: polar quadratic terms are $B_0f_\pi^2\langle\widehat S_{11}\rangle\sigma a/f_a^2$ and $B_0f_\pi\sigma\langle\widehat P_{10}\Phi\rangle/f_a$, zero by assumption 2; the polar kinetic term has no $\sigma a$ quadratic part |
| loops, widths, higher chiral orders | not included | assumptions 4 and 6 |

Hence the propagator matrix is block diagonal: $\sigma$ with $i/(q^2-m_\sigma^2)$, and $G$ for the mesons and the ALP.

## Step 1: separation of the amplitude

With no quadratic radial mixing,
$$
\mathcal M=\mathcal M\big|_{\text{no }\sigma\text{ line}}+\sum_{(r|uv)}\frac{V_{\sigma ar}\,V_{\sigma uv}}{m_\sigma^2-s_{uv}},
$$
with $(r|uv)\in\{(i|jk),(j|ik),(k|ij)\}$. The sign follows from $(iV)\,\tfrac{i}{q^2-m_\sigma^2}\,(iV)=iVV/(m_\sigma^2-q^2)$, the same convention as the plus sign of the existing record. Since $V_{\sigma xy}=-g_{xy}$, the product of two vertices equals the product of two couplings $g$.

## Step 2: amplitude without $\sigma$ line; mapping onto the existing record

The first term is the boxed step-4 amplitude of the [finite-mixing three-pseudoscalar record](alp-finite-mixing-three-pseudoscalar-amplitudes.md) with the following substitutions. All were checked (symbolically and against the existing script's `component_quadratic_blocks`):

| existing record | here | check |
|---|---|---|
| $\widehat M_n=\widehat S_n+i\widehat P_n$ ($x=a/f_a$ derivatives) | $n!\,(\widehat S_{0n}+i\widehat P_{0n})$ | Taylor coefficients |
| $m_q$ | $\widehat S_{00}$ | same $\mathbb M_s$ |
| $\mathcal C=-B_0\{m_q,c_q^{AP}\}-\frac{m_0^2}{N_f}c_{GA}\mathbf 1$ | $B_0\widehat P_{01}$ | $M^2_{Ia}=-f_\pi\langle\mathcal C\lambda_I\rangle/f_a=-B_0f_\pi\langle\widehat P_{01}\lambda_I\rangle/f_a$ |
| $c_{GA}$, $\Theta_n$ | $0$ | assumption 3 |
| $m_0^2/N_f$ | $m_1^2/3$ | singlet term |
| bare ALP mass $m_{a,0}^2$ | $0$ | no bare mass in the model |
| $M^2_{aa}=m_{a,0}^2-\frac{f_\pi^2}{4}\langle\partial_a^2\chi_+|_0\rangle+\frac{2f_\pi^2m_0^2}{N_ff_a^2}c_{GA}^2$ | $-\frac{2B_0f_\pi^2}{f_a^2}\langle\widehat S_{02}\rangle$ | $\partial_a^2\chi_+|_0=8B_0\widehat S_{02}/f_a^2$ |
| $[c_q^{VS},m_q]$ in $\mathcal V^{(3)}$ | $i\widehat S_{01}$ | $\widehat S_1=-i[c_q^{VS},m_q]$ |

The second-derivative correction to $M^2_{aa}$ in the existing record therefore reduces to $\langle\widehat S_{02}\rangle$ alone. The normalization agrees: $2(f_\pi/f_a)\mathbb M^{-1}_{s}(\mathcal C)=(2B_0f_\pi/f_a)\mathbb M^{-1}_{s}(\widehat P_{01})$, so the factor $2B_0$ of the radial records and the factor $2$ of the existing record are the same. The existing amplitude then contains $\widehat c_q^A=-B_0\mathbb M^{-1}_{m_a^2}(\widehat P_{01})$. Its pole-aligned coefficients follow from the Taylor law of that record:
$$
S_n'+iP_n'=n!\sum_{r+s+t=n}\frac{1}{r!\,s!}\,(i\widehat c_q^A)^r\big(\widehat S_{0t}+i\widehat P_{0t}\big)(i\widehat c_q^A)^s ,
$$
with $\widehat P=P_1'=\widehat P_{01}+\{\widehat c_q^A,\widehat S_{00}\}$ in $T^P$. The anomalous terms with $\Theta_1'=\langle\mathbb M^{-1}_{m_a^2}(\mathcal C)\rangle\neq0$ multiply $\Theta_2'=\Theta_3'=0$ and vanish.

## Step 3: off-shell $\sigma$ vertex

The polar cubic Lagrangian ([general-source record](../shared/radial-decay-polar-cartesian-general-sources.md), Step 3) is
$$
\mathcal L_{\sigma PP}=\frac{\sigma}{f_a}\Big[-\frac{B_0}{2}\langle\widehat S_{10}\Phi^2\rangle+\frac{B_0f_\pi}{f_a}a\langle\widehat P_{11}\Phi\rangle+\frac{B_0f_\pi^2}{f_a^2}a^2\langle\widehat S_{12}\rangle+(\partial a)^2\Big].
$$
For two legs $x,y$ with meson column $\lambda_x$ (the notation of step 2 of the existing record; not a generator unless stated) and ALP-field coordinate $\xi_x$ in the original polar coordinates, the derivative term gives $-2\xi_x\xi_yk_x\cdot k_y/f_a$ in $V$. For outgoing daughters, $2k_x\cdot k_y=q^2-m_x^2-m_y^2$. For the parent pair, $k_a=-p$, $k_r=p_r$, so $2k_a\cdot k_r=-2p\cdot p_r=q^2-m_a^2-m_r^2$ with $q^2=(p-p_r)^2=s_{uv}$. Hence in both vertices of a channel
$$
f_ag_{xy}(q^2)=\frac{B_0}{2}\langle\widehat S_{10}\{\lambda_x,\lambda_y\}\rangle-\frac{B_0f_\pi}{f_a}\langle\widehat P_{11}(\xi_x\lambda_y+\xi_y\lambda_x)\rangle-\frac{2B_0f_\pi^2}{f_a^2}\xi_x\xi_y\langle\widehat S_{12}\rangle+\xi_x\xi_y(q^2-m_x^2-m_y^2),
$$
with $q^2=s_{uv}$. This is the on-shell $g_{xy}$ of the radial-decay records with $m_\sigma^2$ replaced by $q^2$; SymPy confirms both statements.

The pole columns in the original polar coordinates are $\lambda_x=\tfrac{2B_0f_\pi}{f_a}\mathbb M_{m_x^2}^{-1}(\widehat P_{01})/\sqrt{N(m_x^2)}$ and $\xi_x=1/\sqrt{N(m_x^2)}$ for every leg with an ALP component, including the parent. A leg without ALP component, such as a charged meson or an exact $K^0$-type pole, has as $\lambda_x$ its LO meson matrix and $\xi_x=0$.

## Step 4: the vertex is unchanged by the pole alignment

The existing record works in pole-aligned coordinates, while Step 3 uses the original polar columns. The two are consistent:

- The pole-aligned change of fields does not involve $\sigma$. Its chiral part $e^{i\Phi/f_\pi}=e^{iaY/f_a}e^{i\Phi'/f_\pi}e^{iaY/f_a}$ has no quadratic term: SymPy gives $\log(e^{tA}e^{tB}e^{tA})=t(2A+B)+O(t^3)$ for non-commuting $A,B$. The rest of the change is linear.
- A quadratic term would also leave the $\sigma xy$ vertex unchanged, because $\sigma$ has no quadratic mixing with $a$ or $\Phi$.

Hence the projected $\sigma xy$ vertex in the aligned coordinates equals the original vertex on the transformed columns, which are the same physical pole vectors.

## Result

$$
\boxed{\begin{aligned}
\mathcal M(a_{\rm phys}\to P_iP_jP_k)={}&\mathcal V^{(4)}_{aijk}+\sum_{(r|uv)}\sum_{\ell,m}\mathcal V^{(3)}_{ar\ell}\,G_{\ell m}(s_{uv})\,\mathcal V^{(3)}_{muv}
+\sum_{(r|uv)}\frac{g_{ar}(s_{uv})\,g_{uv}(s_{uv})}{m_\sigma^2-s_{uv}},\\[1mm]
f_a\sqrt{N(m_x^2)N(m_y^2)}\;g_{xy}(q^2)={}&\frac{2B_0^3f_\pi^2}{f_a^2}\Big\langle\widehat S_{10}\big\{\mathbb M_{m_x^2}^{-1}(\widehat P_{01}),\mathbb M_{m_y^2}^{-1}(\widehat P_{01})\big\}\Big\rangle
-\frac{2B_0^2f_\pi^2}{f_a^2}\Big\langle\widehat P_{11}\big(\mathbb M_{m_x^2}^{-1}(\widehat P_{01})+\mathbb M_{m_y^2}^{-1}(\widehat P_{01})\big)\Big\rangle\\
&-\frac{2B_0f_\pi^2}{f_a^2}\langle\widehat S_{12}\rangle+q^2-m_x^2-m_y^2,\\[1mm]
m_\sigma^2={}&2\lambda f_a^2+\frac{B_0f_\pi^2}{f_a^2}\big(\langle\widehat S_{10}\rangle-2\langle\widehat S_{20}\rangle\big).
\end{aligned}}
$$

Here $(r|uv)$ runs over $(i|jk)$, $(j|ik)$, $(k|ij)$, and $x,y\in\{a,i,j,k\}$ with $m_x^2$ the pole masses. The first line has $\mathcal V^{(4)}$, $\mathcal V^{(3)}$, $G$ and the step-2 columns of the [finite-mixing record](alp-finite-mixing-three-pseudoscalar-amplitudes.md), with the substitutions of Step 2. For a leg without ALP component, replace $\tfrac{2B_0f_\pi}{f_a}\mathbb M^{-1}_{m_x^2}(\widehat P_{01})/\sqrt{N(m_x^2)}\to\lambda_x$ and $1/\sqrt{N(m_x^2)}\to0$ in the general form of Step 3. The residue at $s_{uv}=m_\sigma^2$ is $g_{ar}\,g_{uv}$ with the on-shell couplings of [the radial-decay records](../shared/murayama-sigma-two-pseudoscalar-decays.md).

Note on notation: $g_{xy}(q^2)$ is the existing coupling $g_{ij}$ of the radial-decay records with $m_\sigma^2\to q^2$, not a new quantity.

## Power counting and LO truncation

Hold $\widehat S_{rs}$, $\widehat P_{rs}$ and $m_\sigma^2$ fixed as $f_a$ grows. For meson-like daughters, $\lambda_r=O(1)$ and $\xi_r=O(f_\pi/f_a)$. For the ALP-like parent, $\lambda_a=O(f_\pi/f_a)$ and $\xi_a=O(1)$. Then:

- $g_{uv}$ for two meson-like daughters is $O(1/f_a)$. It is dominated by $\tfrac{B_0}{2f_a}\langle\widehat S_{10}\{\lambda_u,\lambda_v\}\rangle$. The $\widehat P_{11}$, $\widehat S_{12}$ and $q^2$ terms are suppressed by $(f_\pi/f_a)^2$.
- $g_{ar}$ is $O(f_\pi/f_a^2)$. Three terms contribute at this order: the $\widehat S_{10}$ term with $\lambda_a$, the $\widehat P_{11}$ term with $\xi_a\lambda_r$, and the kinetic term $\xi_a\xi_r(q^2-m_a^2-m_r^2)$. The $\widehat S_{12}$ term is suppressed by a further $(f_\pi/f_a)^2$.
- The $\sigma$ term is therefore $O(f_\pi/f_a^3)$. Relative to the $O(1/f_a)$ leading amplitude of the existing record it is $O\big((f_\pi/f_a)^2\big)$, times a ratio of the order of $m^2/(m_\sigma^2-s)$. It does not change the $O(1/f_a)$ amplitude. At fixed $\lambda$ instead of fixed $m_\sigma$, $m_\sigma^2\simeq2\lambda f_a^2$ and the $\sigma$ term falls as $f_a^{-5}$.

The formula itself is not truncated. The counting fails near degenerate poles, for $O(1)$ ALP–meson mixing, and near $s_{uv}=m_\sigma^2$.

Numerical check (seed 2026, $m_\sigma^2=1$ GeV$^2$, three meson-like daughters):

| $f_a$ (GeV) | $f_a\,\mathcal M_{\text{no }\sigma}$ | $f_a^3\times$ ($\sigma$ term) (GeV$^2$) | $\sigma$ term / $\mathcal M_{\text{no }\sigma}$ |
|---|---|---|---|
| 4 | 0.110955 | $-6.86\times10^{-4}$ | $-3.9\times10^{-4}$ |
| 8 | 0.111401 | $-7.18\times10^{-4}$ | $-1.0\times10^{-4}$ |
| 16 | 0.111645 | $-7.50\times10^{-4}$ | $-2.6\times10^{-5}$ |
| 32 | 0.111719 | $-7.61\times10^{-4}$ | $-6.6\times10^{-6}$ |

## Only-up model

With the coefficients of [the only-up record](../shared/murayama-sigma-two-pseudoscalar-decays.md) ($\widehat S_{10}=m_ue_1$, $\widehat P_{11}=\widehat P_{01}=-m_ue_1$, $\widehat S_{12}=\widehat S_{02}=-\tfrac{m_u}{2}e_1$, $\widehat S_{20}=0$), the vertex for two neutral (diagonal) columns is
$$
f_ag_{xy}(q^2)=B_0m_u\Big[(\lambda_x)_{11}+\frac{f_\pi}{f_a}\xi_x\Big]\Big[(\lambda_y)_{11}+\frac{f_\pi}{f_a}\xi_y\Big]+\xi_x\xi_y(q^2-m_x^2-m_y^2),
$$
and $g_{\pi^+\pi^-}=g_{K^+K^-}=B_0m_u/f_a$ for charged pairs (no ALP component).

Numerical observation (only-up model): for the parent and three distinct neutral mixed daughters ($\pi^0$, $\eta$, $\eta'$-like), the LO amplitude without $\sigma$ vanishes to $10^{-16}$ in both the direct calculation and the existing formula. A toy potential $\sum_qm_q(1-\cos\theta_q)+\tfrac A2\langle\Phi\rangle^2$ with the same structure (only $\theta_u$ contains $a$) also gives zero, while random linear forms do not. The radial exchange is then the whole tree amplitude for such final states. **Update (2026-11-19):** this is now derived. The cubic vertices vanish by intrinsic parity ([parity record](murayama-cubic-vertex-parity-rule.md)). The remaining contact is a third divided difference of a linear function ([cancellation record](murayama-distinct-neutral-contact-cancellation.md)). The assembled only-up amplitude, channel by channel, is in [the only-up amplitude record](murayama-only-up-three-pseudoscalar-amplitude.md); the quartic vertex in the draft's basis without ALP–meson mixing is in [the quartic-vertex record](murayama-unmixed-basis-quartic-vertex.md).

Until 2026-11-18, `scripts/implementation/only_up_decays.py` (`three_pseudoscalar_width`, through `calculate_alp_three_pseudoscalar_finite.decay_width`) did **not** include $\sigma$ exchange: its docstring stated that the radial scalar is omitted, and it received no $\lambda$, $m_\sigma$ or $\widehat S_{10}$ input. Since 2026-11-19 `three_pseudoscalar` and `three_pseudoscalar_width` take an explicit `radial_mass_squared` ($m_\sigma^2$); the default `None` keeps the former result without $\sigma$. See [the only-up amplitude record](murayama-only-up-three-pseudoscalar-amplitude.md).

## Limits

- If $s_{uv}$ reaches $m_\sigma^2$ inside the Dalitz region ($m_\sigma<m_a-m_r$ for some channel), the tree propagator diverges. A width of $\sigma$ (from $\sigma\to PP$, $\sigma\to aa$, ...) is then needed; this is outside the present calculation. The same holds for meson or ALP poles of $G$.
- Exceptional poles: if $\mathbb M_{m_x^2}$ is singular, as for the state of a degenerate pair of $B_0\{\widehat S_{00},\cdot\}$ with no ALP component, use its LO meson matrix as $\lambda_x$ and $\xi_x=0$. Degenerate eigenspaces need an orthonormal basis. Near a degeneracy the weak-mixing counting fails.
- Without stationarity of $V_\varphi$ the polar and Cartesian tree amplitudes differ, because the $\sigma$ tadpole is not removed (check (c)).
- Tree level and LO only: no loops, rescattering, scalar resonances other than $\sigma$, or higher chiral orders.

## Verification

[Script](../../../scripts/derive_alp_three_pseudoscalar_radial_exchange.py) and [unit tests](../../../scripts/test_derive_alp_three_pseudoscalar_radial_exchange.py). The direct calculation builds the polar Lagrangian, or the Cartesian one through the exact map, as a polynomial in four formal variables with vanishing squares. This is exact to quartic order and includes the nonlinear meson kinetic term $\tfrac{f_\pi^2}{4}\langle\partial U\partial U^\dagger\rangle$ and all derivative vertices. The script takes the poles from the Hessian and contracts the contact and the three exchanges with the full eleven-field propagator ($\sigma$, $a$, nine mesons). The sources are random non-commuting Hermitian matrices for all $r+s\le4$ except those fixed by assumptions 1–2 (seeds 2026, 7, 11). Each seed gives 7 poles with an ALP component and 3 without. The analytic side calls `compact_amplitude` of the existing script and adds the trace-form $\sigma$ term.

| check | residual |
|---|---|
| symbolic: second-order term of the symmetric map; polar radial mixing $\propto\langle\widehat S_{11}\rangle,\widehat P_{10}$; off-shell vertex = conjectured form = on-shell $g$ with $m_\sigma^2\to q^2$; parent-pair momentum | exact (0) |
| symbolic: mapping $\mathcal C=B_0\widehat P_{01}$, $M^2_{aa}$; Cartesian quartic contacts | exact |
| polar Hessian vs existing `component_quadratic_blocks` with the mapping | $1\times10^{-16}$ GeV$^2$ |
| (a) direct polar vs analytic, 3 seeds × 3 daughter sets × 2 Dalitz points (relative) | $\le1.8\times10^{-11}$ |
| (a) direct $\sigma$ exchange vs trace form (relative to amplitude) | $\le1.0\times10^{-12}$ |
| (a) size of the $\sigma$ term in these points (relative) | $8\times10^{-6}$ to $1.4$ |
| (b) direct polar vs direct Cartesian (relative) | $\le7.6\times10^{-15}$ |
| (c) $m_\varphi^2$ shifted by $0.05$ GeV$^2$: polar vs Cartesian (relative) | $0.51$ (breaks, as expected) |
| (d) residue at $s_{12}=m_\sigma^2(1\pm10^{-6})$ vs $g_{a3}g_{12}$ of the $\sigma\to PP$ record (relative) | $3.9\times10^{-10}$ |
| (e) large $\lambda$: table above; polar = Cartesian (relative) | $\le4\times10^{-14}$ |
| (f) only-up: direct polar vs analytic; polar vs Cartesian; factorized vertex (relative) | $9\times10^{-11}$; $5\times10^{-15}$; $<10^{-12}$ |
| (g) weak mixing: table above | $f_a^3\sigma$ term and $f_a\mathcal M_{\text{no }\sigma}$ converge |

The precision of (a) is limited by the numerical pole vectors and the resolvent; it is not a missing contribution, since the $\sigma$ term exceeds the residual by at least $10^3$ in every point.

```sh
.venv/bin/python scripts/derive_alp_three_pseudoscalar_radial_exchange.py
.venv/bin/python -m unittest scripts.test_derive_alp_three_pseudoscalar_radial_exchange
```
