---
title: Role of the K pi S wave in ALP three-meson decays
status: mixed
created: 2026-09-23
updated: 2026-09-23
bibtex_keys:
  - Aloni:2018vki
  - BaBar:2015kii
  - Fariborz:1999gr
  - BaBar:2014asx
provenance: literature
---

# Role of the K pi S Wave in ALP Three-Meson Decays

## Inputs and provenance

The ALP decay model and its definition of $S_{K\pi}$ are recorded in [the Aloni et al. extract](../extracts/Aloni:2018vki/kpi-s-wave.md) [@Aloni:2018vki]. The experimental $K\pi$ production amplitude is recorded in [the BaBar MIPWA extract](../extracts/BaBar:2015kii/kpi-mipwa.md) [@BaBar:2015kii]. The scalar-nonet couplings are recorded in [the Fariborz--Schechter extract](../extracts/Fariborz:1999gr/kappa-couplings.md) [@Fariborz:1999gr]. Direct evidence concerning the scalar $K\eta$ channel is recorded in [the BaBar $K\eta$ extract](../extracts/BaBar:2014asx/keta-scalar-channel.md) [@BaBar:2014asx].

## Identification of the amplitude

In the ALP model, $S_{K\pi}(m_{K\pi})$ is the empirical mass-dependent shape and phase of the $I=1/2$, $J=0$ $K\pi$ production amplitude measured in $\eta_c\to K\bar K\pi$. It is not a single $K_0^*(1430)$ Breit-Wigner and is not the normalization of a $K\pi$ scattering partial wave. BaBar fixes an arbitrary reference magnitude and phase in its binned production amplitude; Aloni et al. separately set the magnitude of their decay term using $\kappa$-sector parameters and the form factor $\mathcal F$. [@Aloni:2018vki; @BaBar:2015kii; @Fariborz:1999gr]

The term ``$S$ wave'' means that the $K\pi$ pair has orbital angular momentum $J=L=0$. For a spin-zero parent, the corresponding $Y_0^0$ angular factor is constant. Therefore this scalar contribution has no independent $H_\pm$ helicity components and no nontrivial longitudinal, parallel, or perpendicular transversity decomposition. This is a consequence of the angular structure recorded by BaBar, rather than a separately stated ALP-model amplitude. [@BaBar:2015kii]

## Effect on the ALP decay amplitude

The model writes

$$
\mathcal A(a\to K\bar K\pi)=\frac{f_\pi}{f_a}
\left[\mathcal A(a\to S_{K\pi}(K\pi)K)
+\mathcal A(a\to a_0(K\bar K)\pi)\right].
$$

For $a\to K^+K^-\pi^0$, the scalar contribution contains

$$
\begin{aligned}
\mathcal A_{S_{K\pi}K}\propto{}&
(p_a\cdot p_{K^-})(p_\pi\cdot p_{K^+})S_{K\pi}(m_{K^+\pi})\\
&+(p_a\cdot p_{K^+})(p_\pi\cdot p_{K^-})S_{K\pi}(m_{K^-\pi}).
\end{aligned}
$$

It consequently controls the scalar line shape and strong phase in both charged $K\pi$ bands. Because the two terms occur in one amplitude, they interfere with each other. Because $\mathcal A_{S_{K\pi}K}$ and $\mathcal A_{a_0\pi}$ are added before squaring, the $K\pi$ scalar contribution also interferes with the $a_0(K\bar K)\pi$ contribution. It therefore affects the $K\bar K\pi$ Dalitz distribution and the phase-space-integrated width. These interference statements follow directly from expanding the coherent squared amplitude $|\mathcal A_{S_{K\pi}K}+\mathcal A_{a_0\pi}|^2$. [@Aloni:2018vki]

No $S_{K\pi}$ factor occurs in the paper's displayed amplitudes for $a\to\eta^{(\prime)}\pi\pi$, vector-vector, or vector-pseudoscalar channels. Its role in that model is confined to the scalar quasi-two-body part of $a\to K\bar K\pi$. [@Aloni:2018vki]

## Charge channels

The explicit formula for $K^+K^-\pi^0$ is a representative charge channel, not the only channel to which $S_{K\pi}$ applies. Aloni et al. state that the corresponding amplitude has the same value for the six physical final states

$$
K^+K^-\pi^0,
\quad K^+K_S\pi^-,
\quad K^+K_L\pi^-,
\quad K^-K_S\pi^+,
\quad K^-K_L\pi^+,
\quad K_SK_L\pi^0,
$$

with different $U(3)$ generator expressions for the different charge assignments. They consequently use

$$
\Gamma(a\to K\bar K\pi)
=6\,\Gamma(a\to K^+K^-\pi^0).
$$

This equality is the flavor-symmetric prescription adopted in the paper. It neglects small charged-neutral kaon and pion mass differences. [@Aloni:2018vki]

In the flavor basis, these physical channels arise from exactly the four modes

$$
K^+K^-\pi^0,
\qquad K^0\bar K^0\pi^0,
\qquad \pi^+K^-K^0,
\qquad \pi^-K^+\bar K^0.
$$

The mixed charged-neutral mode $\pi^-K^+\bar K^0$ appears experimentally as the two channels $K^+K_S\pi^-$ and $K^+K_L\pi^-$, while its charge conjugate $\pi^+K^-K^0$ gives $K^-K_S\pi^+$ and $K^-K_L\pi^+$. The neutral flavor mode $K^0\bar K^0\pi^0$ supplies the listed $K_SK_L\pi^0$ channel. Thus all three additional modes in the question are included in the paper's $K\bar K\pi$ rate. The same empirical $S_{K\pi}$ shape and phase are used after the appropriate charge substitutions, but the charged formula's flavor trace must be replaced by the corresponding $U(3)$ generator expression. [@Aloni:2018vki]

## Relation to K eta

There is a physical scalar $K\eta$ channel: BaBar observed $K_0^*(1430)\to K\eta$ and measured

$$
\frac{\mathcal B(K_0^*(1430)\to K\eta)}
{\mathcal B(K_0^*(1430)\to K\pi)}
=0.092\pm0.025^{+0.010}_{-0.025}.
$$

The scalar-nonet interaction used as an input by Aloni et al. also contains a distinct $\kappa K\eta$ coupling. Thus a scalar $K\eta$ amplitude is physically and model-theoretically meaningful. [@Fariborz:1999gr; @BaBar:2014asx]

However, arXiv:1811.03474 does not define an $S_{K\eta}(m_{K\eta})$ function or an $a\to K\bar K\eta$ amplitude. The available BaBar $K\eta$ analysis uses a $K_0^*(1430)$ Breit-Wigner and a constant nonresonant term, rather than the binned model-independent construction used for $K\pi$. Therefore it is not a directly interchangeable $K\eta$ analogue of the imported $S_{K\pi}$ function. [@Aloni:2018vki; @BaBar:2014asx; @BaBar:2015kii]

The measured $K\pi$ mass interval crosses the $K\eta$ and $K\eta'$ thresholds. It is reasonable to interpret the empirical $K\pi$ shape as already containing whatever inelastic effects those openings induce in the observed production process. This is an interpretation, not a coupled-channel result: the binned $K\pi$ amplitude does not by itself determine a transition amplitude into $K\eta$. Constructing $S_{K\eta}$ consistently would require an additional coupled-channel production or scattering model. [@BaBar:2015kii]

## Chiral-consistent normalization

The derivation of the production kernel that replaces the unexplained Aloni scalar prefactor is recorded in [the chiral-matching record](kpi-s-wave-chiral-matching.md). The result is an $s$-dependent kernel fixed by the hatted Alda couplings, not a universal constant. A normalized final-state-interaction function may multiply this kernel without changing its leading-order chiral normalization. [@Alda:2026nrd; @Aloni:2018vki; @BaBar:2015kii]
