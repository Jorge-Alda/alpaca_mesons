---
title: Comparison of vector-meson extensions of chiral perturbation theory for ALPs
status: mixed
created: 2026-09-23
updated: 2026-09-23
bibtex_keys: [Alda:2026bar, Ecker:1988te, Ecker:1989yg, Bando:1987br, Fujiwara:1984mp, Harada:2003jx, Ruiz-Femenia:2003jdx, Ovchynnikov:2025gpx, Bai:2025fvl, Balkin:2025enj]
provenance: literature
---

# Vector-Meson Extensions for the Alda--Rigolin Programme

## Executive comparison

The three options are not mutually exclusive descriptions of the same scope.

| Framework | What is added | Main theoretical input | Main strength | Main limitation | Best role here |
|---|---|---|---|---|---|
| Resonance Chiral Theory (RChT) | Explicit vector nonet, preferably antisymmetric tensors, coupled to CCWZ blocks | Chiral symmetry, large-$N_c$, resonance saturation, QCD short-distance constraints | Systematic operator basis in even and odd intrinsic parity; clean matching back to chiral low-energy constants | Many couplings; finite resonance truncation and tree-level widths are model assumptions | Primary dynamical EFT for $\rho,\omega,K^*,\phi$ |
| Hidden Local Symmetry (HLS) | Vectors as gauge bosons of a redundant local factorization of the chiral field | Hidden gauge organization plus derivative expansion | Compact $VPP$, vector--photon mixing, and anomalous interactions; efficient phenomenology | Gauge structure imposes correlations not required by QCD; $a=2$, ideal mixing, universality, and strict VMD are optional assumptions | Economical implementation and cross-check; useful if a global HLS fit is adopted |
| Background fields in gauged WZW | ALP, SM gauge fields, and vector/axial mesons inserted as external chiral sources in the anomaly functional | Anomaly matching and Bardeen-type counterterms | Most transparent route to basis-independent odd-parity ALP vertices | Does not by itself supply vector dynamics, $VPP$, masses, propagators, widths, or all invariant homogeneous terms | Anomaly module and non-negotiable consistency check, not a standalone resonance model |

This classification follows the foundational RChT, HLS, and anomalous-background-field constructions. [@Ecker:1988te; @Ecker:1989yg; @Bando:1987br; @Fujiwara:1984mp; @Bai:2025fvl]

## 1. Resonance Chiral Theory

### Construction and assumptions

In the representation best adapted to standard short-distance constraints, vectors are antisymmetric tensors $V_{\mu\nu}$ transforming as $V_{\mu\nu}\to hV_{\mu\nu}h^\dagger$. The leading even-parity terms are

$$
\mathcal L_V^{\rm even}=\mathcal L_{\rm kin}(V)
+\frac{F_V}{2\sqrt2}\langle V_{\mu\nu}f_+^{\mu\nu}\rangle
+\frac{iG_V}{\sqrt2}\langle V_{\mu\nu}u^\mu u^\nu\rangle+\cdots.
$$

Odd parity requires the WZW action plus independent $VJP$, $VVP$, and, for multiparticle applications, further resonance operators. [@Ruiz-Femenia:2003jdx]

The predictive assumptions are explicit:

- leading order in $1/N_c$, so resonances are initially narrow and amplitudes are tree level;
- truncation to the lightest vector nonet (“single-resonance” or lowest-meson dominance);
- resonance saturation of selected chiral low-energy constants;
- QCD short-distance constraints imposed on chosen Green functions and form factors;
- additional prescriptions for finite widths, loops, $SU(3)$ breaking, and $\omega$--$\phi$ mixing.

Different spin-1 representations are physically equivalent only after all local terms required by field redefinitions are retained. At a finite operator truncation, Proca, antisymmetric-tensor, massive-Yang--Mills, and HLS implementations can therefore produce different ultraviolet behavior. This is a truncation effect, not evidence for physically distinct vector particles. [@Ecker:1989yg; @Ruiz-Femenia:2003jdx]

### Contact with data

The short-distance-constrained antisymmetric theory gives a useful but not precision-level description without subleading corrections. A historical example is $\omega\to\pi\gamma$, predicted within $5$--$30\%$ of the measured width quoted in the source. Conversely, $\rho$ exchange alone gives only about $1.4\,\mathrm{MeV}$ for $\omega\to3\pi$ versus $7.52\,\mathrm{MeV}$, demonstrating that strict VMD fails and a direct term is needed. [@Ruiz-Femenia:2003jdx]

Thus “agreement with data” should be assessed channel by channel after a correlated fit. Chiral symmetry and short-distance constraints reduce parameter freedom, but they do not determine all couplings, $SU(3)$ breaking, widths, or final-state phases.

### Extension of `Alda:2026bar`

1. Promote the vector octet to a $U(3)$ nonet and let $V_{\mu\nu}$ transform with the same CCWZ compensator $h$ used in `Alda:2026bar`.
2. Insert the Alda covariant blocks $u_\mu$, $\Gamma_\mu$, $\chi_\pm$, $f_\pm^{\mu\nu}$ and the topological singlet $X_{TC}$ into the complete one-vector and two-vector operator bases.
3. For the even-parity sector, impose the Alda no-mixing parametrization
   $$
   \widehat c_q^A=-\mathbb M_a^{-1}(C),\qquad
   \widehat c_q^P=-c_q^{AP}-\mathbb M_a^{-1}(C),\qquad
   \widehat c_q^V=c_q^{VS},
   $$
   and derive direct $aVP$, $aVPP$, and $aVV$ vertices from invariant sources. This is a derived programme based on the transformation law of the CCWZ blocks, not yet a completed calculation. [@Alda:2026bar]
4. Treat the odd-parity sector separately: include WZW, $VJP$, $VVP$, and the anomalous counterterm before applying any field redefinition. The direct-only no-mixing rule established outside WZW cannot simply be assumed to hold in the anomalous sector.
5. Fit $F_V,G_V$ and odd-parity combinations to $\rho\to\pi\pi$, $V\to P\gamma$, $e^+e^-\to PP$, $e^+e^-\to VP$, $\omega/\phi\to3P$, and transition-form-factor data, while enforcing low-energy anomaly normalization and selected QCD asymptotics.
6. Add finite widths through self-energies or dispersive propagators in a way that preserves the low-energy expansion; do not insert widths into tree propagators before checking Ward identities.

**Assessment:** this is the strongest default for a reusable, systematically improvable Alda extension, especially if both even- and odd-parity amplitudes are required.

## 2. Hidden Local Symmetry

### Construction and assumptions

HLS factorizes the chiral field and introduces a local redundancy whose gauge bosons are identified with $\rho$, $\omega$, $K^*$, and $\phi$. At leading order the central parameters are $F_\pi$, $g$, and $a=F_\sigma^2/F_\pi^2$. The familiar $a=2$ choice yields tree-level KSRF relations, coupling universality, and pion-form-factor vector dominance, but $a=2$ is not fixed by chiral symmetry. [@Bando:1987br; @Harada:2003jx]

The anomalous action contains homogeneous terms with free coefficients in addition to WZW. Therefore HLS does not imply strict VMD unless those coefficients are specialized. [@Fujiwara:1984mp]

At loop level HLS is an EFT with $g=O(p)$ and $m_V=O(p)$. A calculation that uses tree-level HLS with measured nondegenerate vector masses, constant widths, ideal mixing, and $a=2$ is consequently a phenomenological model, not a controlled fixed-order HLS calculation. [@Harada:2003jx]

### Contact with data

The success of $a\simeq2$ and KSRF-type relations explains why HLS efficiently reproduces leading $\rho$ phenomenology and the pion electromagnetic form factor. However, anomalous channels require non-VMD contact terms. The covariant ALP implementation of Balkin et al. explicitly retains a $\gamma PPP$ term to improve $\omega\to3\pi$, while fitting high-mass form factors to $e^+e^-$ data. [@Balkin:2025enj]

That continuation is useful but strongly model dependent: it assumes form factors are real, flavor independent within a Lorentz class, and functions only of the hardest scale. Exact basis-dependent cancellations also rely on the universal HLS vector mass but are imposed after measured masses and widths are substituted. [@Balkin:2025enj]

### Extension of `Alda:2026bar`

1. Factorize the same $U(3)$ field $U=u^2$ used by Alda into HLS variables $\xi_{L,R}$ and define their covariant derivatives with the general ALP left/right sources, including flavor-changing entries.
2. Introduce a full vector nonet gauge field and retain general $a$, $g$, kinetic terms, $SU(3)$-breaking operators, and singlet mixing. Do not hard-wire $a=2$ or ideal $\omega$--$\phi$ mixing before fitting.
3. Express all ALP source dependence through $c_q^{AP}$, $c_q^{VS}$, $c_{GA}$, and $X_{TC}$, or prove explicitly that amplitudes are invariant under the Alda quark-field redefinitions.
4. Include the complete anomalous HLS functional, not only a selected $PVV$ term. Fit its homogeneous coefficients to radiative and three-pseudoscalar vector decays.
5. Verify that low-energy expansion of HLS exchange plus contact terms reproduces the pure Alda amplitudes and does not double count their local $O(p^4)$ or higher terms.
6. Compare amplitudes numerically against an RChT implementation using the same masses, widths, and fitted channels. Differences after matching quantify truncation uncertainty.

**Assessment:** HLS is the shortest route to a broad event-generator implementation, but it is less transparent than RChT for separating symmetry consequences from VMD-inspired parameter choices.

## 3. Background Fields in the WZW Functional

### Construction and assumptions

The background-field method treats the ALP derivative source and vector or axial-vector mesons as left/right chiral gauge sources inside the gauged WZW functional. The associated counterterm is essential for gauge invariance, anomaly matching, and cancellation of auxiliary chiral rotations. [@Bai:2025fvl]

This method fixes anomaly-related structures, not complete resonance dynamics. Promoting a background vector to a propagating meson additionally requires kinetic and mass terms, non-anomalous interactions such as $VPP$, symmetry-breaking terms, mixing, and a width prescription. Gauge-invariant homogeneous odd-parity terms remain possible and generally carry free coefficients. Consequently, “background fields for WZW” is best viewed as an anomalous-sector construction that can be combined with RChT or HLS, not as a third complete alternative.

The Bai et al. implementation includes three flavors, the anomalous $U(1)_A$ sector, vectors and axial vectors, and demonstrates auxiliary-phase cancellation. It uses leading-order and idealized inputs, so its benchmark widths establish consistency rather than precision agreement with a global body of vector-meson data. [@Bai:2025fvl]

### Extension of `Alda:2026bar`

1. Promote the general Alda left/right ALP sources, including non-diagonal $c_q^A,c_q^V,c_q^P,c_q^S$, to background connections in the full $U(3)$ gauged WZW functional.
2. Retain the anomalous $U(1)_A$ topological source rather than requiring a rotation with $\operatorname{Tr}\kappa=c_G$.
3. Derive and include the Bardeen-type counterterm before choosing the no-mixing variables. Then transform the WZW functional, counterterm, and explicit partonic anomalous couplings together.
4. Prove field-redefinition invariance for a minimal test set: $a\to\gamma\gamma$, $a\to V\gamma$, $a\to VV$, $V\to a\gamma$, and the five-field WZW amplitudes $a\to4P$.
5. Couple this anomaly module to either RChT or HLS dynamics for propagating vectors. Subtract local terms already generated by resonance exchange to prevent double counting.

**Assessment:** this is the preferred way to complete the WZW part of the Alda programme and to derive $a\to4P$, but it cannot replace a dynamical vector EFT.

## 4. Recommended staged programme

The following recommendation is an inference from the literature rather than an established result.

1. **Use RChT as the primary dynamical layer.** Its CCWZ transformation law interfaces directly with `Alda:2026bar`, and its operator basis makes assumptions and short-distance constraints auditable.
2. **Use the full background-field WZW functional as the anomaly layer.** This is needed regardless of whether the propagating vectors are represented by RChT or HLS fields.
3. **Implement HLS as a matched alternative.** Fit it to the same calibration set and use the RChT--HLS spread, after matching low-energy constants and pole inputs, as one estimate of resonance-model uncertainty.
4. **Validate in ascending complexity:** SM $V\to PP$ and $V\to P\gamma$; pion/kaon electromagnetic form factors; $\omega,\phi\to3P$; pseudoscalar transition form factors; then ALP amplitudes.
5. **Require four checks for every ALP amplitude:** independence of quark-field redefinitions, correct soft/chiral limit, anomaly normalization where applicable, and acceptable short-distance scaling.

## 5. Evidence status and open uncertainties

- **Literature established:** the structures and transformation properties of RChT, HLS, and the gauged WZW functional; the need for anomalous counterterms; failures of strict VMD in selected anomalous channels. [@Ecker:1988te; @Bando:1987br; @Fujiwara:1984mp; @Ruiz-Femenia:2003jdx; @Bai:2025fvl]
- **Recent literature implementation:** invariant or covariant ALP resonance amplitudes exist, but they use restricted flavor assumptions or phenomenological form-factor continuations. [@Ovchynnikov:2025gpx; @Balkin:2025enj; @Bai:2025fvl]
- **Derived recommendation:** RChT plus a background-field WZW anomaly module is the most transparent baseline for the general, flavor-changing, no-mixing framework of `Alda:2026bar`.
- **Unresolved:** a controlled error budget near and above $1\,\mathrm{GeV}$; consistent loop/width treatment; complete $U(3)$ odd-parity operator basis with the Alda topological source; coupled-channel final-state interactions; and matching of the exclusive sum to inclusive QCD.

Recent invariant phenomenological models do not yet close the last issue: one implementation finds that its exclusive invariant width above roughly $1\,\mathrm{GeV}$ does not approach the perturbative gluonic result, signaling missing channels or interactions. [@Ovchynnikov:2025gpx]

## 6. Does vector--photon mixing arise automatically?

No. One must distinguish a bilinear vector--photon transition from the stronger assumption of strict vector-meson dominance (VMD), in which photon--hadron interactions are required to proceed only through vector exchange.

- **RChT:** after electromagnetism is included as an external source, the allowed operator $F_V\langle V_{\mu\nu}f_+^{\mu\nu}\rangle/(2\sqrt2)$ generates neutral-vector--photon transitions. This is a standard leading operator, but it is not forced merely by introducing a vector field: it must be retained and its coefficient must be fixed or fitted. Field redefinitions can move contributions between an explicit bilinear transition and local photon vertices, provided the complete accompanying operator set is kept. [@Ecker:1988te; @Ecker:1989yg]
- **HLS:** electromagnetic gauging of the leading HLS Lagrangian naturally produces a bilinear transition through the combination of the hidden gauge field and the electromagnetic external field. Thus mixing is built into the usual gauged leading-order construction, as in $f_\pi^2\langle(gV_\mu-e_\mu)^2\rangle$. Nevertheless, strict VMD follows only for special parameter choices, such as the familiar tree-level $a=2$ limit in the even-parity pion form factor and corresponding choices of anomalous homogeneous coefficients. [@Bando:1987br; @Fujiwara:1984mp; @Harada:2003jx; @Balkin:2025enj]
- **Background-field WZW:** this construction naturally supplies anomaly-related vertices involving photons and vectors, such as $PV\gamma$ and related odd-parity interactions, but it does not by itself supply the ordinary even-parity vector--photon two-point transition. A separate dynamical vector Lagrangian, for example RChT or HLS, is required. [@Bai:2025fvl]
- **Generic Proca or massive-Yang--Mills realizations:** symmetry permits vector--photon transition operators, but their explicit appearance and form are representation- and basis-dependent. They are not guaranteed by the statement that vector mesons have been added. Equivalence to other spin-1 representations requires the local terms induced by field redefinitions. [@Ecker:1989yg]

Therefore, among the methods compared here, conventional gauged HLS includes vector--photon mixing most automatically; RChT includes it naturally through a standard allowed operator; and the background-field WZW method alone does not generate the required non-anomalous mixing. None of these facts alone implies strict VMD.
