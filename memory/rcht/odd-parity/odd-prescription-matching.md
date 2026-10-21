---
title: Aloni and Bai vector-meson prescriptions in the original odd basis
status: derived
created: 2026-10-21
updated: 2026-10-21
bibtex_keys: [Aloni:2018vki, Bai:2025fvl, Ruiz-Femenia:2003jdx, Fujiwara:1984mp]
provenance: analytical-derivation
---

# Aloni:2018vki and Bai:2025fvl matching

## Scope and conventions

This record matches the specified low-energy neutral two-vector and one-vector/one-photon amplitudes to the repository's original eleven odd coefficients $(c_1,\ldots,c_7,d_1,\ldots,d_4)$. It uses the canonical convention $\epsilon^{0123}=+1$, metric $(+---)$, $c_i=-c_i^R$, and $d_i=-d_i^R$. The Aloni phenomenological suppression function is set to $\mathcal F(m_a)=1$, as appropriate to the requested pure chiral regime. There is no numerical fit.

The matching is at amplitude level with common vector mass $M_V$, common $F_V$, commuting neutral flavor generators, independent pseudoscalar-mass and quark-spurion insertions, and the displayed source field normalizations. It does not identify an on-shell pole-amplitude match with equality of the complete off-shell covariant actions.

## Literature inputs

Aloni:2018vki, Supplemental Eq. $\mathcal L_{\rm VMD}$, uses

$$
\mathcal L_{\rm VMD}\supset \frac{g_{VVP}}{4}\langle \boldsymbol P V^{\mu\nu}\widetilde V_{\mu\nu}\rangle,
\qquad g_{VVP}=\frac{3g^2}{8\pi^2 f_\pi},\qquad g\simeq\sqrt{12\pi},
$$

with $\boldsymbol P=\boldsymbol\Phi/\sqrt2$ and $\boldsymbol V=\boldsymbol V_{\rm ideal}/\sqrt2$ in its displayed matrix convention. Its ALP dependence in this term is through the ALP--pseudoscalar overlap matrix. The same Lagrangian supplies vector--photon mass mixing. [@Aloni:2018vki, Supplemental Material, Eq. $\mathcal L_{\rm VMD}$]

Bai:2025fvl uses the complete gauged WZW functional, its Bardeen counterterm, vector/axial background connections, the ALP shifts of those connections, and pseudoscalar mixing. For the neutral vector sector, $\mathbb B_V=\mathbb B_L+\mathbb B_R$ and $v=\mathbb B_V/2=g\boldsymbol V_{\rm ideal}/\sqrt2$. The convention conversion is $\boldsymbol\Phi_{\rm Bai}=2\boldsymbol P_{\rm Aloni}$ and $f_{\pi,\rm Bai}=\sqrt2 f_{\pi,\rm Aloni}$. [@Bai:2025fvl, Eqs. (2.8), (2.12), (2.5), and (3.1)]

## Neutral VVP projection

For $F_Q=\operatorname{tr}(\{T_i,T_j\}T_P)$ and independent mass-spurion factor $F_H$, the repository's neutral VVP numerator is

$$
\frac{H(s,t)}{F_Q}=(d_1-d_3)m_P^2+d_3(s+t)+4d_2\frac{F_H}{F_Q}.
$$

Matching a mass-independent $PV\widetilde V$ pole vertex for independent $m_P^2$ and $F_H$ sets $d_1=d_3$, $d_2=0$. Matching the Lorentz coefficient $G_{PVV}/F_Q$ then gives, for a common pole mass,

$$
d_1=d_3=-\frac{F G_{PVV}}{4\sqrt2},\qquad d_2=0,
$$

where $F=f_{\pi,\rm Aloni}$ and $G_{PVV}$ multiplies the canonical tensor $\epsilon(\varepsilon_i,\varepsilon_j,q_i,q_j)$. The independent coefficients $d_4$ and the four other VVP-insensitive directions remain free in the combined projection below.

The source normalizations obtained by expanding the displayed Lagrangians and converting their field-strength tensors to the canonical Lorentz tensor are

$$
\begin{array}{c|cc}
 & G_{PVV}/F_Q & G_{PV\gamma}/F_Q\\
\hline
\text{Aloni:2018vki} & \dfrac{3\sqrt2 g^2}{64\pi^2F} & \dfrac{3eg}{64\pi^2F}\\[1ex]
\text{Bai:2025fvl, pseudoscalar pole} & -\dfrac{\sqrt2 N_Cg^2}{16\pi^2F} & -\dfrac{N_Ceg}{8\pi^2F}
\end{array}
$$

The $PV\gamma$ coefficient in the RChT matching is the direct VJP piece after subtracting the VVP exchange contribution. Solving both neutral pole projections, with $M_i=M_j=M_V$, yields these affine families:

$$
\begin{array}{c|ccccc}
 & c_1 & c_2 & c_3 & d_1=d_3 & d_2\\
\hline
\text{Aloni:2018vki} & c_5-c_2 & c_6+\dfrac{3\sqrt2g}{512\pi^2}\left(1-\dfrac{F_Vg}{M_V}\right) & 0 & -\dfrac{3g^2}{256\pi^2} & 0\\[2ex]
\text{Bai:2025fvl pole} & c_5-c_2 & c_6+\dfrac{\sqrt2N_Cg}{128\pi^2}\left(\dfrac{F_Vg}{M_V}-2\right) & 0 & \dfrac{N_Cg^2}{64\pi^2} & 0
\end{array}
$$

In each row $c_4,c_5,c_6,c_7,d_4$ remain free. The exact SymPy system has rank six and five free coefficient directions. These constraints match the specified neutral pole amplitudes; the free directions are not assigned values or uncertainties.

## Bai complete ALP source and the independent-source limitation

The Bai functional also contains a direct ALP axial-source contact and an anomalous rotation term. With $\mathcal C=N_C/(48\pi^2)$, the neutral cubic expansion gives a pseudoscalar term proportional to $-3\mathcal C\sqrt2g^2\boldsymbol\Phi/F$, a direct source term proportional to $N_Cg^2(k_L-k_R)/(24\pi^2f_a)$, and the anomalous rotation term proportional to $N_Cg^2c_{gg}(\theta_L-\theta_R)/(48\pi^2f_a)$. In the displayed matrix normalization, the combined coefficient matrix multiplying $a\,dV\wedge dV$ is

$$
K_{\rm Bai}^{\rm eff}=\frac{N_Cg^2}{48\pi^2Ff_a}\left[2F(k_L'-k_R')+Fc_{gg}(\theta_L-\theta_R)-3\sqrt2f_a\sum_P\theta_{Pa}T_P\right].
$$

The flavor trace with the two vector generators is taken after assembling this matrix. The anomalous term cancels the basis-rotation dependence of the direct source and pseudoscalar-pole pieces. The executable calculation verifies the rotation residual is exactly zero.

The direct axial-source projection of the original VVP operators, with $\mathcal U_a=-(k_L-k_R)/f_a$, requires

$$
d_3^{\rm direct\ source}=-\frac{N_Cg^2}{96\pi^2},
$$

whereas the neutral pseudoscalar pole match requires $d_3^{\rm pole}=N_Cg^2/(64\pi^2)$. Their difference is $5N_Cg^2/(192\pi^2)$. For Aloni:2018vki, the displayed Lagrangian has no independent direct axial-source contact, while its pole match has $d_3=-3g^2/(256\pi^2)$. Thus the pole coefficient families do not extend to a single complete matching of all independent source operators in either framework without additional local source terms. They still give explicit, useful matches for the requested physical neutral decay amplitudes.

## Numerical coefficient representatives

The companion fit-style JSON stores numerical coefficient dictionaries, not
the affine formulae above. Its displayed representative sets the five
undetermined coordinates $c_4,c_5,c_6,c_7,d_4$ to zero; this convention does
not determine those coefficients physically. Keep the paper decay-constant
conventions distinct: $f_{\pi,\mathrm{Aloni}}=0.092\,\mathrm{GeV}$ and
$f_{\pi,\mathrm{Bai}}=0.130\,\mathrm{GeV}$. They differ by $\sqrt2$ in the
definitions of $U$, while both give the neutral-pion entry
$(\log U)_{11}=i\pi^0/(0.092\,\mathrm{GeV})+\cdots$. Thus the common RChT
normalization in the coefficient equations is
$F=f_{\pi,\mathrm{Aloni}}=f_{\pi,\mathrm{Bai}}/\sqrt2$.

Use $g=\sqrt{12\pi}$ and $N_C=3$. For Aloni:2018vki, set the common vector
mass from the cited Fujiwara relation $M_V^2=a g^2f_{\pi,\mathrm{Aloni}}^2$
with $a=2$. This gives $M_V=0.7988557850\,\mathrm{GeV}$. For Bai:2025fvl,
use $M_V=0.770\,\mathrm{GeV}$ as the repository operator normalization
because the paper does not specify a common dynamical vector pole mass for
this matching.

## Matching the even-parity photon--vector mixing

The RChT photon source contains both $F_V$ and the spurion correction $z_V$.
For a neutral vector generator $T$, define its charge overlap
$q=\langle TQ\rangle$ and spurion overlap
$h=\langle T\{\chi_+,Q\}\rangle$. The tensor--photon bilinear from the
quadratic RChT source is

$$
\mathcal L_{V\gamma}^{\mathrm{RChT}}
=-\frac{e}{\sqrt2}(F_Vq+z_Vh)V_{\mu\nu}F^{\mu\nu}.
$$

Eliminating the tensor in the first-order quadratic Lagrangian gives

$$
V_{\mu\nu}=\frac{P_{\mu\nu}+\sqrt2e(F_Vq+z_Vh)F_{\mu\nu}}{M_V^2},
\qquad
\mathcal L_{\mathrm{red}}
=-\frac{[P_{\mu\nu}+\sqrt2e(F_Vq+z_Vh)F_{\mu\nu}]^2}{4M_V^2}.
$$

Writing $P_\mu=M_VR_\mu$ and shifting to the canonically normalized vector
field yields the Proca mass bilinear

$$
\mathcal L_{R\gamma}^{\mathrm{RChT}}
=-\sqrt2M_Ve(F_Vq+z_Vh)R_\mu A^\mu.
$$

Aloni:2018vki instead writes
$-M_V^2(e/g)\langle\boldsymbol V_\mu Q\rangle A^\mu$, with
$\boldsymbol V=V_{\rm ideal}/\sqrt2$. Its coefficient for the same generator
is $-M_V^2eq/(\sqrt2g)$. The matching equation is therefore

$$
\sqrt2M_VeF_Vq=\frac{M_V^2eq}{\sqrt2g},
\qquad
F_V=\frac{M_V}{2g}=\frac{f_{\pi,\mathrm{Aloni}}}{\sqrt2}
=0.06505382387\,\mathrm{GeV},
\qquad z_V=0.
$$

This determines $F_V$ from the displayed Aloni mass mixing; it does not use
the provisional-fit $F_V$. Bai:2025fvl's displayed WZW/background functional
has no dynamical even-parity quadratic photon--vector mixing term. For the
requested prescription representative, the missing term is assigned
$F_V=0$ and $z_V=0$; this is the absence of these terms in that prescription,
not a measurement of a physical vector-current residue. Bai's odd-parity
WZW $aV\gamma$ vertex remains nonzero and is matched separately.

With the five null coordinates chosen as zero, the resulting numerical
representatives have $(c_1,c_2,d_1,d_3)=(-0.002577520392,0.002577520392,
-0.044762327745,-0.044762327745)$ for Aloni and
$(0.041240326273,-0.041240326273,0.179049310978,0.179049310978)$ for Bai.
The drop-in JSON stores the five-dimensional coefficient null basis and the
matching inputs. Its even-parameter slots use the common vector mass, zero
mass-breaking and spurion-source parameters, the prescription-matched $F_V$,
and $z_V=0$. These are conditional representatives, not fits or uncertainty
estimates.

## Full double-vector-dominance form factor

As a separate check of the Aloni vector--photon mechanism, matching the complete neutral $P\to\gamma^*\gamma^*$ tensor-RChT amplitude (WZW plus VJP plus VVP) to the anomaly-normalized common-mass double-pole form

$$
F_{P\gamma^*\gamma^*}(s,t)=F_{P\gamma\gamma}(0,0)\frac{M_V^4}{(M_V^2-s)(M_V^2-t)}
$$

gives the rank-seven family

$$
c_1=c_3=d_2=0,\qquad
c_2=c_5=c_6-\frac{\sqrt2M_VN_C}{128\pi^2F_V},\qquad
d_1=d_3=\frac{M_V^2N_C}{64\pi^2F_V^2},
$$

with $c_4,c_6,c_7,d_4$ free. This is the full double-pole form-factor realization; it is not the same assertion as matching only the $PV\widetilde V$ vertex or matching an independent ALP axial source.

## Reproducibility

[`scripts/derive_odd_prescription_matching.py`](../../../scripts/derive_odd_prescription_matching.py) carries out the quadratic photon--vector normalization, Lorentz-component conversion, neutral flavor projections, source expansion, affine solves and symbolic residual checks. Inputs include the [Aloni:2018vki VMD extract](../../extracts/Aloni:2018vki/vector-meson-dominance.md), the [Bai:2025fvl matching extract](../../extracts/Bai:2025fvl/no-axial-vector-and-rcht-matching.md), and the [canonical VVP operators](../../extracts/Ruiz-Femenia:2003jdx/vvp-three-pseudoscalar.md). Run:

```text
.venv/bin/python scripts/derive_odd_prescription_matching.py --derive-only
```
