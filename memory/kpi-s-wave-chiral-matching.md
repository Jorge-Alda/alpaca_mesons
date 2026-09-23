---
title: Chiral matching of the K pi S-wave production kernel
status: derived
created: 2026-09-23
updated: 2026-09-23
bibtex_keys:
  - Alda:2026bar
  - Aloni:2018vki
  - BaBar:2015kii
  - Omnes:1958hv
  - Watson:1952ji
  - Buettiker:2003pp
  - Pelaez:2016tgi
  - Pelaez:2020gnd
  - Aston:1987ir
  - Jamin:2000wn
  - Jamin:2001zq
  - Jamin:2006tj
  - VonDetten:2021rax
provenance: analytical-derivation
---

# Chiral Matching of the K Pi S-Wave Production Kernel

## Scope

This record derives the production kernel that replaces the unexplained scalar-sector normalization in the $a\to S_{K\pi}K$ ansatz of Aloni et al. The phenomenological $\mathcal F_{SPP}$ factor is omitted. The result follows from the chiral-invariant leading-order ALP Lagrangian and no-mixing prescription of Alda and Rigolin. [@Aloni:2018vki; @Alda:2026bar]

The executable derivation is [`scripts/derive_kpi_s_wave_matching.py`](../scripts/derive_kpi_s_wave_matching.py).

## Inputs and assumptions

The matrix-level interaction is the result derived in [the leading-order $a\Phi^3$ record](alp-phi3-lo-chipt-derivation.md), whose literature inputs are recorded in [the Alda Appendix A extract](extracts/Alda:2026bar/appendix-a-expansion.md) [@Alda:2026bar]:

$$
\mathcal L_{a\Phi^3}
=\frac{\partial_\mu a}{6f_af_\pi}
\left\langle\widehat c_q^A[\Phi,[\Phi,\partial^\mu\Phi]]\right\rangle
-\frac{B_0a}{6f_af_\pi}
\left\langle\Phi^3\{\widehat c_q^P,m_q\}\right\rangle.
$$

The following assumptions are used:

- The calculation is leading order in the chiral expansion and linear in $1/f_a$.
- The hatted couplings are in the Alda no-ALP-meson-mixing basis.
- The isospin limit $m_u=m_d=\widehat m$ and $m_{K^+}=m_{K^-}=m_K$ is used.
- The coupling matrices are diagonal, with entries $A_i=(\widehat c_q^A)_{ii}$ and $P_i=(\widehat c_q^P)_{ii}$.
- The exact isospin-limit meson generators are used; physical $\pi^0$--$\eta$--$\eta'$ mixing can instead be included by evaluating the same trace formulas with the fitted generators.
- Natural units, the mostly-minus metric, and $S=1+iT$ are used.

For

$$
a(p)\to K^+(p_1)\pi^0(p_2)K^-(p_3),
$$

define

$$
s=(p_1+p_2)^2,
\qquad t=(p_2+p_3)^2,
\qquad u=(p_1+p_3)^2,
$$

so that

$$
s+t+u=m_a^2+2m_K^2+m_\pi^2.
$$

The generators are

$$
\lambda_{K^+}=\sqrt2E_{13},
\qquad
\lambda_{K^-}=\sqrt2E_{31},
\qquad
\lambda_{\pi^0}=\operatorname{diag}(1,-1,0).
$$

## Flavor traces

Using the ordered-generator definitions in [the generic-generator record](alp-three-meson-generic-generators.md), SymPy gives

$$
T^A_{K^+\pi^0;K^-}=2(A_s-A_u),
$$

$$
T^A_{K^+K^-;\pi^0}=4(A_u-A_s),
$$

and

$$
T^A_{\pi^0K^-;K^+}=2(A_s-A_u).
$$

Their sum vanishes, as required by the Jacobi identity. For the symmetry-breaking term,

$$
\sum_{\sigma\in S_3}
\lambda_{\sigma(K^+)}\lambda_{\sigma(\pi^0)}
\lambda_{\sigma(K^-)}
=\operatorname{diag}(4,0,2),
$$

and hence

$$
T^P_{K^+\pi^0K^-}
=8P_u\widehat m+4P_sm_s.
$$

## Full leading-order amplitude

Define

$$
\Delta_A=A_u-A_s.
$$

Substitution into the three-distinct-meson amplitude gives

$$
\boxed{
\mathcal M_{\mathrm{LO}}(s,t,u)
=\frac{1}{6f_af_\pi}
\left[
\Delta_A\left(s+t-2u+2m_\pi^2-2m_K^2\right)
-4B_0\left(2\widehat mP_u+m_sP_s\right)
\right].
}
$$

This is the chiral-invariant contact amplitude. It is not proportional to the momentum polynomial used by Aloni et al.

## Two-band decomposition

Eliminating $u$ gives

$$
\mathcal M_{\mathrm{LO}}(s,t)
=\frac{1}{6f_af_\pi}
\left[
3\Delta_A(s+t)-2\Delta_A m_a^2-6\Delta_A m_K^2
-4B_0(2\widehat mP_u+m_sP_s)
\right].
$$

The symmetric decomposition into the two $K\pi$ bands is

$$
\boxed{
\mathcal M_{\mathrm{LO}}(s,t)=P_{\mathrm{LO}}(s)+P_{\mathrm{LO}}(t),
}
$$

where

$$
\boxed{
P_{\mathrm{LO}}(s)=
\frac{
\Delta_A(3s-m_a^2-3m_K^2)
-2B_0(2\widehat mP_u+m_sP_s)
}{6f_af_\pi}.
}
$$

This is the Alda-consistent per-band scalar production kernel. Under the requirements that both charge-conjugate bands use the same one-variable kernel and that no additional contact remainder is retained, the decomposition is fixed by the full leading-order amplitude.

If the overall factor is written in the convention of Aloni et al.,

$$
\mathcal M=\frac{f_\pi}{f_a}
\left[\widetilde P_{\mathrm{LO}}(s)+\widetilde P_{\mathrm{LO}}(t)\right],
$$

then

$$
\boxed{
\widetilde P_{\mathrm{LO}}(s)=
\frac{
\Delta_A(3s-m_a^2-3m_K^2)
-2B_0(2\widehat mP_u+m_sP_s)
}{6f_\pi^2}.
}
$$

Thus the chiral replacement for the Aloni normalization is a coupling- and $s$-dependent production kernel, not a universal constant such as $(8.2/\mathrm{GeV})^2$.

## Why an Aloni-like constant cannot be calculated

For one band, the Aloni momentum factor is

$$
(p\cdot p_3)(p_1\cdot p_2)
=\frac14(m_a^2+m_K^2-s)(s-m_K^2-m_\pi^2).
$$

This is quadratic in $s$, whereas $P_{\mathrm{LO}}(s)$ is linear in $s$ and contains two independent coupling structures, $\Delta_A$ and $2\widehat mP_u+m_sP_s$. Consequently no constant prefactor can make the Aloni band amplitude reproduce the chiral result for all $s$.

There is a second independent obstruction: the BaBar $S_{K\pi}$ production amplitude has an arbitrary overall magnitude and phase. BaBar fixes these by a reference-bin convention rather than by a chiral or scattering normalization. Therefore a coefficient multiplying the raw BaBar function is convention dependent. [@BaBar:2015kii]

One can force equality at a selected point $s_\star$, but the resulting coefficient

$$
C(s_\star)=
\frac{P_{\mathrm{LO}}(s_\star)}
{S_{K\pi}^{\mathrm{BaBar}}(s_\star)}
$$

depends on $m_a$, the hatted ALP couplings, the arbitrary BaBar normalization and phase, and the chosen matching point. It is not a universal hadronic prefactor and does not reproduce the full leading-order $s$ dependence.

## Consistent final-state-interaction prescription

Introduce a dimensionless scalar final-state-interaction function normalized at a low-energy subtraction point,

$$
F_0^{K\pi}(s_0)=1.
$$

The minimal matched amplitude is

$$
\boxed{
\mathcal M_{S\text{-wave}}(s,t)
=P_{\mathrm{LO}}(s)F_0^{K\pi}(s)
+P_{\mathrm{LO}}(t)F_0^{K\pi}(t).
}
$$

Setting $F_0^{K\pi}=1$ reproduces the complete Alda-consistent leading-order contact amplitude exactly. A normalized elastic Omnes function is one possible realization [see the Omnes extract](extracts/Omnes:1958hv/phase-representation.md) [@Omnes:1958hv]:

$$
\Omega(s;s_0)=
\exp\left\{
\frac{s-s_0}{\pi}
\int_{s_{\mathrm{th}}}^{\infty}
\frac{\delta_0^{1/2}(s')\,ds'}
{(s'-s_0)(s'-s-i0)}
\right\},
\qquad \Omega(s_0;s_0)=1.
$$

Then $F_0^{K\pi}=\Omega$ at the elastic level. More generally the dispersive representation is $\Omega(s)$ times a subtraction polynomial whose low-energy expansion is matched to $P_{\mathrm{LO}}(s)$.

## Where BaBar enters

The Omnes construction and direct use of the BaBar amplitude are two alternative prescriptions. They should not be applied simultaneously to the same band.

### Baseline dispersive prescription: BaBar is not required

The phase $\delta_0^{1/2}(s)$ in the Omnes integral is the phase of the $I=1/2$ $K\pi\to K\pi$ scattering partial wave. It should therefore be taken from a scattering analysis, such as LASS or a modern dispersive coupled-channel fit, rather than from the BaBar magnitude. The calculation proceeds as follows:

1. Compute $P_{\mathrm{LO}}(s)$ from the hatted Alda couplings.
2. Obtain the scattering phase $\delta_0^{1/2}(s)$.
3. Construct $\Omega(s;s_0)$ with $\Omega(s_0;s_0)=1$.
4. Use $P_{\mathrm{LO}}(s)\Omega(s;s_0)+P_{\mathrm{LO}}(t)\Omega(t;s_0)$.

In this baseline prescription, the BaBar result is used only as an external validation of the predicted production-wave phase and line shape. It does not determine the normalization or prefactor.

Watson's theorem explains the possible phase comparison: in an elastic region, and in the absence of additional process-dependent phases, a production amplitude with given isospin and angular momentum has the phase of the corresponding scattering channel, up to the production convention [see the Watson extract](extracts/Watson:1952ji/final-state-phase.md) [@Watson:1952ji]. BaBar observes phase behavior similar to LASS below about $1.5\,\mathrm{GeV}$ after an arbitrary constant phase shift, while the measured magnitudes differ substantially. This supports using scattering-phase information more strongly than treating the BaBar magnitude as universal. [@BaBar:2015kii]

### Empirical alternative: use BaBar instead of Omnes

If a data-driven isobar model is preferred over a dispersive construction, define a convention-independent relative BaBar shape at a physical matching point $s_\star$,

$$
\widehat S_{\mathrm{BaBar}}(s;s_\star)
=\frac{S_{K\pi}^{\mathrm{BaBar}}(s)}
{S_{K\pi}^{\mathrm{BaBar}}(s_\star)},
\qquad
\widehat S_{\mathrm{BaBar}}(s_\star;s_\star)=1.
$$

Then a phenomenological amplitude is

$$
\mathcal M_{\mathrm{emp}}(s,t)
=P_{\mathrm{LO}}(s)\widehat S_{\mathrm{BaBar}}(s;s_\star)
+P_{\mathrm{LO}}(t)\widehat S_{\mathrm{BaBar}}(t;s_\star).
$$

The ratio removes BaBar's arbitrary overall complex normalization. This prescription uses both the BaBar magnitude and phase, but it matches the undressed chiral kernel only at $s=s_\star$. Since BaBar measures only the physical decay region and fixes its convention near $1.45\,\mathrm{GeV}$, it cannot by itself supply the subthreshold normalization $F_0^{K\pi}(0)=1$. The empirical prescription is therefore a model choice, not a consequence of chiral symmetry or unitarity.

An improved empirical matching would fit a subtraction polynomial multiplying an interpolation of the BaBar wave so that its low-energy behavior agrees with the chiral kernel. That requires additional assumptions about analytic continuation and about the process dependence of the BaBar production amplitude.

### Approximate use of only the BaBar phase

The BaBar phase may be used as an approximate input to $\Omega$ only if one assumes Watson's theorem applies to the $\eta_c$ production wave in the chosen range and that crossed-channel and production phases are negligible. This is weaker than importing the full BaBar amplitude, because its nonuniversal magnitude is discarded. It is nevertheless less controlled than using a genuine scattering phase, especially after inelastic channels open.

Therefore the recommended hierarchy is:

1. Use a coupled-channel scattering solution for a controlled calculation.
2. Use an elastic scattering phase in a normalized Omnes function below important inelasticities.
3. Use the normalized BaBar complex shape as a phenomenological alternative when reproducing the observed isobar line shape is more important than controlled chiral matching.

## Phase input and the meaning of effective elasticity

Write the elastic-channel partial-wave $S$ matrix as

$$
S_{11}(s)=\eta(s)e^{2i\delta(s)},
\qquad 0\leq\eta(s)\leq1.
$$

Here $\delta$ is the $K\pi$ phase shift and $\eta$ is the elasticity. In flux-normalized conventions, $1-\eta^2$ measures the probability removed from the $K\pi$ channel into all other channels. A scalar one-channel Omnes function uses $\delta(s)$ and assumes $\eta=1$ on the relevant right-hand cut. If $\eta<1$, unitarity couples the $K\pi$ form factor to form factors for the other open channels, and the scalar Omnes function must be replaced by a Muskhelishvili--Omnes matrix. [@Watson:1952ji; @Omnes:1958hv; @Jamin:2001zq]

The phrase ``below important inelasticities'' is therefore not identical to ``below the first inelastic threshold.'' A channel can be open but weakly coupled, leaving $\eta$ close to one. Conversely, a nearby strongly coupled threshold can invalidate a one-channel treatment rapidly. There is no universal numerical boundary.

For numerical work, an explicit tolerance should be stated. One possible conservative operational assumption is

$$
1-\eta^2<0.05,
$$

which corresponds to $\eta>\sqrt{0.95}\simeq0.975$. The five-percent choice is an analysis assumption, not a literature result. It should be varied as a systematic uncertainty.

For the $I=1/2$ scalar $K\pi$ wave, the relevant two-body thresholds are approximately

$$
\sqrt{s_{K\eta}}\simeq1.04\,\mathrm{GeV},
\qquad
\sqrt{s_{K\eta'}}\simeq1.45\,\mathrm{GeV}.
$$

Thus strict two-body elasticity ends at the $K\eta$ opening. Phenomenologically, however, the $K\eta$ coupling is small in the cited coupled-channel analyses, while $K\eta'$ is the important coupled channel around $K_0^*(1430)$. The broad $K_0^*(1430)$ overlaps the $K\eta'$ threshold, so a single-channel Omnes function should not be used as a controlled description of its complete line shape. [@Jamin:2000wn; @Jamin:2001zq; @VonDetten:2021rax]

### Recommended phase sources

For threshold through $\sqrt{s}=0.967\,\mathrm{GeV}$, use the $I=1/2$ Roy--Steiner phase of Buttiker, Descotes-Genon, and Moussallam. This is the canonical low-energy solution and imposes elastic unitarity in its domain. [See the source extract](extracts/Buettiker:2003pp/kpi-roy-steiner-phase.md) [@Buettiker:2003pp]

For an implementable phase from threshold into the resonance region, use the dispersively constrained parametrization of Pelaez and Rodas, checking its stated isospin and partial-wave conventions. It is constrained through $1.6\,\mathrm{GeV}$ and is a more convenient primary input than interpolating old LASS points directly. [See the source extracts](extracts/Pelaez:2016tgi/kpi-dispersive-parametrization.md) and [the later review](extracts/Pelaez:2020gnd/kpi-dispersive-review.md) [@Pelaez:2016tgi; @Pelaez:2020gnd]

Use LASS as high-energy scattering information and a cross-check, not as an unqualified table of pure-$I=1/2$ phase and inelasticity. The measured charge channel requires an $I=3/2$ subtraction, elasticity below $1.29\,\mathrm{GeV}$ was imposed, and ambiguities occur at high energy. [See the LASS extract](extracts/Aston:1987ir/lass-kpi-scattering.md) [@Aston:1987ir]

If the calculation includes the full $K_0^*(1430)$ region, use a coupled-channel solution rather than choosing a phase for a scalar Omnes integral. The canonical option is the $K\pi$--$K\eta$--$K\eta'$ solution of Jamin, Oller, and Pich. A modern phenomenological option is the construction of von Detten et al., which combines a dispersive low-energy phase with coupled-channel resonance dynamics and LASS-region information. [See the coupled-channel extracts](extracts/Jamin:2000wn/kpi-coupled-scattering.md), [the scalar-form-factor solution](extracts/Jamin:2001zq/coupled-scalar-form-factors.md), and [the modern inelastic construction](extracts/VonDetten:2021rax/kpi-form-factor-inelastic.md) [@Jamin:2000wn; @Jamin:2001zq; @VonDetten:2021rax]

For the present ALP problem, the practical choice is therefore:

1. If every physical $K\pi$ invariant mass remains below about $1.3\,\mathrm{GeV}$, use the Pelaez--Rodas phase, with the Roy--Steiner result fixing or validating the low-energy part, and assign an uncertainty for treating the weak $K\eta$ channel as elastic.
2. If the integration reaches the $K_0^*(1430)$ and $K\eta'$ region, replace the scalar Omnes factor by a coupled-channel matrix. This also requires deriving the Alda-consistent production kernels for $K\eta$ and $K\eta'$, since a coupled-channel solution acts on a vector of source amplitudes.

Even when the external value of $s$ lies in the elastic region, the Omnes integral extends to infinity. Above the last trustworthy phase input, a high-energy phase continuation is still required. Additional subtractions suppress this sensitivity, and varying the continuation and matching scale provides the associated systematic uncertainty. [@Omnes:1958hv]

## Limitations

- Above the $K\eta$ and $K\eta'$ thresholds, a one-channel Omnes factor is an approximation; the systematic replacement is a coupled-channel Omnes matrix acting on a vector of production kernels.
- Exact three-body unitarity includes crossed-channel rescattering. A Khuri--Treiman-type equation is required beyond the independent-band approximation.
- Dressing the complete $s$-channel projection and then adding a separately dressed crossed band would double count the crossed contribution. The object to dress in the two-isobar prescription is the single-band kernel $P_{\mathrm{LO}}$, not the projection of the full contact amplitude.
- If $A_u\ne A_d$ or $P_u\ne P_d$, the source has isospin-breaking components. The $I=1/2$ and $I=3/2$ $K\pi$ production kernels must then be separated before applying the BaBar $I=1/2$ phase.

A charge-basis implementation with split thresholds, $I=1/2$--$I=3/2$ mixing, strong $m_d-m_u$ effects, and electromagnetic caveats is developed in [the isospin-breaking record](kpi-isospin-breaking.md).

## Verification

Run

```bash
.venv/bin/python scripts/derive_kpi_s_wave_matching.py
```

The script prints the flavor tensors, complete invariant amplitude, symmetric band decomposition, full-amplitude $K\pi$ $S$-wave projection, and the nonconstant coefficient that would be required to force the Aloni momentum polynomial onto the chiral kernel. Assertions verify the flavor traces, Jacobi identity, Mandelstam reduction, exact reconstruction $P_{\mathrm{LO}}(s)+P_{\mathrm{LO}}(t)=\mathcal M_{\mathrm{LO}}$, and angular projection.
