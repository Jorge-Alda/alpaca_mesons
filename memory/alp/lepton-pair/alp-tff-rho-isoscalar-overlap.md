---
title: Two-pion and mixed rho-isoscalar overlap in the ALP TFF
status: derived
created: 2026-11-06
updated: 2026-11-06
bibtex_keys: [Alda:2026nrd, Messerli:2025rnv, Ecker:1988te, Ecker:1989yg, Ruiz-Femenia:2003jdx]
provenance: analytical-derivation
---

# Scope and inputs

This record compares the two-pion discontinuity with the mixed
$\rho I$ VVP contribution, $I=\omega,\phi$, in the low-energy ALP transition
form factor (TFF). It uses Alda:2026nrd's no-mixing ALP sources and the
source-resolved VVP numerator and photon ordering of the
[tree-level isoscalar decomposition](alp-isoscalar-tff-source-decomposition.md),
the charged two-pion cut normalization and production-wave convention of
the [pion current and production record](alp-pion-vector-and-production-form-factors.md),
and the [isoscalar pole-residue record](alp-pion-isoscalar-pole.md)
[@Alda:2026nrd; @Messerli:2025rnv; @Ecker:1988te; @Ecker:1989yg;
@Ruiz-Femenia:2003jdx]. Use natural units, the real-epsilon crossed-current
convention, and $\mathcal F_{a\pi\pi\gamma^*}=f_1/e$.
No numerical width, production phase, or subtraction function is selected.

# Sources and vector generators

Use the Alda:2026nrd no-mixing sources
$$
U_a=\frac{2\hat c_q^A}{f_a}=-\frac{2}{f_a}\mathbb M_a^{-1}(\mathcal C),\qquad
Y_a=-\frac{4 B_0}{f_a}\{\hat c_q^P,m_q\}
=\frac{4 B_0}{f_a}\{c_q^{AP}+\mathbb M_a^{-1}(\mathcal C),m_q\}.
$$

The TFF records use vector matrices $T_i$ normalized by
$\operatorname{tr}(T_iT_j)=\delta_{ij}$. Elsewhere in the repository,
$\lambda_i$ usually denotes pseudoscalar flavor matrices. For this note's
requested vector notation, explicitly define $\lambda_i=\sqrt2\,T_i$, so
$\operatorname{tr}(\lambda_i\lambda_j)=2\delta_{ij}$. In the ideal
neutral basis,
$$
\lambda_\rho=\operatorname{diag}(1,-1,0),\qquad
\lambda_{\omega_I}=\operatorname{diag}(1,1,0),\qquad
\lambda_{\phi_I}=\operatorname{diag}(0,0,\sqrt2),
$$
and the physical $\omega,\phi$ generators use the same mixing-angle
rotation as their $T_i$ counterparts. The signed photon residues keep the
project's original normalization:
$$
R_i=\operatorname{tr}(T_iS_\gamma)
=\frac{1}{\sqrt2}\operatorname{tr}(\lambda_iS_\gamma).
$$
This relation is needed because the response residues were defined with
the normalized $T_i$ generators.

# Mixed term and its two photon placements

The project's ordered VVP numerator is $P_{ij}$. In the requested
$\lambda_i$ normalization, its mixed entry is
$$
P_{\rho I}(s_1,s_2)=2\left\{
\operatorname{tr}(\lambda_\rho\lambda_I U_a)
\left[(d_1-d_3)m_a^2+d_3(s_1+s_2)\right]
-4d_2\operatorname{tr}(\lambda_\rho\lambda_I Y_a)\right\}.
$$
The established notation for the complete mixed contribution, including
both photon placements, is
$$
\begin{aligned}
F_{a\gamma^*\gamma^*}^{(\rho I),\mathrm{dress}}(s_1,s_2)
=-2\sum_{I=\omega,\phi}\Big[&
W_\rho^{\mathrm{dress}}(s_1)W_I^{\mathrm{pole}}(s_2)
P_{\rho I}(s_1,s_2)\\
&+W_I^{\mathrm{pole}}(s_1)W_\rho^{\mathrm{dress}}(s_2)
P_{I\rho}(s_1,s_2)\Big].
\end{aligned}
$$
The first and second terms are the two ordered VVP photon placements. They
are not matrix components or separate physical form factors. In the first
version of this note, $M_{12}$ and $M_{21}$ were temporary names for these
two terms, respectively; they were introduced here, are not used elsewhere
in the project, and have been replaced by the established
$F_{a\gamma^*\gamma^*}^{(\rho I),\mathrm{dress}}$ notation.

The broad rho uses the dressed response
$$
W_\rho^{\mathrm{dress}}(s)=
\frac{R_\rho}{M_{0\rho}^2-s-\Delta_\rho^{\mathrm{self}}(s+i0)},
\qquad
\operatorname{Im}\Delta_\rho^{\mathrm{self}}(s+i0)
=+\sqrt{s}\,\Gamma_{\rho\to\pi\pi}(s).
$$
On the upper rim its denominator is
$M_{0\rho}^2-s-\operatorname{Re}\Delta_\rho^{\mathrm{self}}(s)
-i\operatorname{Im}\Delta_\rho^{\mathrm{self}}(s+i0)$. It does not use a
zero-width denominator $M_\rho^2-s-i0$. The real part, mass and residue
conditions follow the [rho self-energy record](alp-rho-pipi-self-energy.md).
By contrast, the first-approximation $I=\omega,\phi$ leg is the selected
narrow stable pole,
$$
W_I^{\mathrm{pole}}(s)=\frac{R_I}{M_I^2-s-i0}.
$$
Thus the zero-width denominator belongs to the narrow isoscalar leg in this
approximation, not to the broad rho. A finite-width isoscalar response
requires the pole/continuum replacement discussed in the
[mixed-vector strategy](alp-tff-mixed-vector-strategy.md).

For diagonal sources, $P_{I\rho}$ is obtained by exchanging the photon
invariants in $P_{\rho I}$. The affine numerator obeys
$$
\frac{P_{\rho I}(s_1,s_2)}{M_I^2-s_2}
=\frac{P_{\rho I}(s_1,M_I^2)}{M_I^2-s_2}
-2d_3\operatorname{tr}(\lambda_\rho\lambda_I U_a).
$$
Therefore the first photon placement has a rho cut in $s_1$ multiplying
both a narrow $I$ pole in $s_2$ and a term regular in $s_2$; the other
placement has the exchanged structure. The regular term contributes
$+4d_3 R_I W_\rho^{\mathrm{dress}}(s_1)
\operatorname{tr}(\lambda_\rho\lambda_I U_a)$ to the mixed TFF.

For the rho response,
$$
\operatorname{Disc}W_\rho^{\mathrm{dress}}(x)=
\frac{2iR_\rho\,\operatorname{Im}\Delta_\rho^{\mathrm{self}}(x+i0)}
{[M_{0\rho}^2-x-\operatorname{Re}\Delta_\rho^{\mathrm{self}}(x)]^2+
 [\operatorname{Im}\Delta_\rho^{\mathrm{self}}(x+i0)]^2},
\qquad
\operatorname{Disc}W_I^{\mathrm{pole}}(x)=2\pi iR_I\delta(M_I^2-x).
$$
For analytic $P_{\rho I}$, the double discontinuity of the complete mixed
term is
$$
\begin{aligned}
\operatorname{Disc}_{s_1}\operatorname{Disc}_{s_2}
F_{a\gamma^*\gamma^*}^{(\rho I),\mathrm{dress}}
=-2\sum_I\Big[&
P_{\rho I}(s_1,s_2)
\operatorname{Disc}W_\rho^{\mathrm{dress}}(s_1)
\operatorname{Disc}W_I^{\mathrm{pole}}(s_2)\\
&+P_{I\rho}(s_1,s_2)
\operatorname{Disc}W_I^{\mathrm{pole}}(s_1)
\operatorname{Disc}W_\rho^{\mathrm{dress}}(s_2)\Big].
\end{aligned}
$$
Each placement has a rho two-pion single discontinuity and rho-cut times
narrow-pole double support. It does not have two continuum cuts in the
narrow-$I$ variable. The regular term contributes to the rho single cut
but not to that double discontinuity.

# Two-pion cut and matching conditions

With $K_{\pi\pi}(x)=x\sigma_\pi^3(x)/(96\pi)$, the two-pion unitarity
relation in the stated convention is
$$
\operatorname{Im}_x F^{2\pi}_{a\gamma^*\gamma^*}(x,y)
=K_{\pi\pi}(x)[F_\pi^V(x)]^*
\mathcal F_{a\pi\pi\gamma^*}(x,y).
$$
The per-ordering double-spectral density uses
$K_{\pi\pi}/2=x\sigma_\pi^3/(192\pi)$ and includes the exchanged photon
placement separately. When Schwarz reflection applies, the corresponding
discontinuity is $2i$ times this imaginary part. A double density does not
specify the full single-cut function: production subtraction functions,
crossed-channel left-hand cuts, and terms without a discontinuity in the
second photon virtuality must also be matched.

In [the isoscalar pole-residue record](alp-pion-isoscalar-pole.md),
$C_I(x)$ denotes the coefficient of the isoscalar pole in the full
$f_1(x,y)$, before converting to $\mathcal F=f_1/e$. The mixed VVP
subamplitude has its own projected coefficient $C_I^{(\rho I)}(x)$; it is
only one part of $C_I$. The pole-supported single-cut relation for that
component is
$$
2iK_{\pi\pi}(x)[F_\pi^V(x)]^*\frac{C_I^{(\rho I)}(x)}{e}
=-2R_I\operatorname{Disc}W_\rho^{\mathrm{dress}}(x)
P_{\rho I}(x,M_I^2).
$$
This is the two-pion unitarity relation for the factorized $\rho I$
subamplitude when $F_\pi^V$ and $\mathcal F^{(\rho I)}$ use the same
$\rho\to\pi\pi$ vertex, dressed rho response, and VVP source numerator. It
does not identify the full $C_I$ with the mixed VVP term. The corresponding
double discontinuity is this relation multiplied by the narrow-pole jump
$2\pi iR_I\delta(M_I^2-y)$, with the exchanged photon placement included.
The affine numerator also gives a single-rho-cut term
$+4d_3R_I\operatorname{Disc}W_\rho^{\mathrm{dress}}(x)
\operatorname{tr}(\lambda_\rho\lambda_I U_a)$, which is not in the double
density and must be included in the single-cut reconstruction.

# Cut-replacement prescription

The analytic supports are not disjoint: the $\pi\pi$ self-energy cut of the
dressed rho is the same two-pion intermediate state used in the physical
two-pion unitarity contribution. Retain the physical two-pion cut once and
remove the rho-cut projection from the separately added mixed term. This
prescription does not require identifying the full physical production wave
with the mixed model amplitude.

Choose $N$ subtractions at $s=0$, with $N$ large enough for convergence,
and the same upper spectral limit $\Lambda^2$ as the two-pion construction:
$$
W_\rho^{\pi\pi,[N]}(s;\Lambda)=
\frac{s^N}{2\pi i}\int_{4M_\pi^2}^{\Lambda^2}
\frac{dx\,\operatorname{Disc}W_\rho^{\mathrm{dress}}(x)}
{x^N(x-s-i0)}.
$$
Define the mixed term's two-pion cut projection by replacing each dressed
rho response in both photon placements:
$$
\begin{aligned}
F_{a\gamma^*\gamma^*}^{(\rho I),\pi\pi,[N]}(s_1,s_2;\Lambda)
=-2\sum_{I=\omega,\phi}\Big[&
W_\rho^{\pi\pi,[N]}(s_1;\Lambda)W_I^{\mathrm{pole}}(s_2)
P_{\rho I}(s_1,s_2)\\
&+W_I^{\mathrm{pole}}(s_1)W_\rho^{\pi\pi,[N]}(s_2;\Lambda)
P_{I\rho}(s_1,s_2)\Big].
\end{aligned}
$$
In each line the VVP numerator multiplies both response factors: the first
product is $W_\rho^{\pi\pi,[N]}(s_1;\Lambda)
W_I^{\mathrm{pole}}(s_2)P_{\rho I}(s_1,s_2)$, and the second is its
exchanged placement.

The matched representation is
$$
F_{a\gamma^*\gamma^*}^{\mathrm{matched}}
=F_{a\gamma^*\gamma^*}^{2\pi,\mathrm{cut-only},[N]}
+F_{a\gamma^*\gamma^*}^{(\rho I),\mathrm{dress}}
-F_{a\gamma^*\gamma^*}^{(\rho I),\pi\pi,[N]}.
$$
$F^{2\pi,\mathrm{cut-only},[N]}$ contains the physical two-pion single and
double dispersive contributions, including the $\rho I$ production
subamplitude. Its $\rho I$ subtraction functions are assigned to the
explicit mixed remainder below, not included a second time. The last two
terms have zero net rho $\pi\pi$
discontinuity when the full rho cut is included and there are no additional
first-sheet poles or cuts. Thus the physical two-pion discontinuity is kept
once. If a supplied $F^{2\pi}$ already includes the $\rho I$ subtraction
functions, remove those same functions before adding the remainder below.

For full support ($\Lambda=\infty$) and a rho response with no other
first-sheet singularities,
$$
W_\rho^{\mathrm{dress}}(s)-W_\rho^{\pi\pi,[N]}(s;\infty)
=p_{N-1}(s),\qquad
p_{N-1}(s)=\sum_{k=0}^{N-1}
\frac{(W_\rho^{\mathrm{dress}})^{(k)}(0)}{k!}s^k.
$$
The explicit rho-cut-free mixed remainder is therefore
$$
-2\sum_I\Big[
p_{N-1}(s_1)W_I^{\mathrm{pole}}(s_2)P_{\rho I}(s_1,s_2)
+W_I^{\mathrm{pole}}(s_1)p_{N-1}(s_2)P_{I\rho}(s_1,s_2)\Big].
$$
It is an isoscalar-pole/subtraction function, not an additional $\pi\pi$
cut. For finite $\Lambda$, the difference also contains the $x>\Lambda^2$
rho-cut tail; retain it explicitly as a higher-energy term or omit it with
the cutoff limitation stated. Additional rho cuts or isolated poles must
also be kept separately.
For clarity, the tail response used in this replacement is
$$
W_\rho^{\mathrm{tail},[N]}(s;\Lambda)=
\frac{s^N}{2\pi i}\int_{\Lambda^2}^{\infty}
\frac{dx\,\operatorname{Disc}W_\rho^{\mathrm{dress}}(x)}
{x^N(x-s-i0)}.
$$
$F^{(\rho I),\mathrm{tail},[N]}$ is obtained by replacing both rho
responses in the ordered mixed term with this tail response.

The affine decomposition shows that the cut projection contains both the
pole-supported double-spectral term and the regular single-rho-cut term in
each photon placement. The regular remainder is independent of the
isoscalar virtuality, so it cancels from the following double-subtracted
combination; its rho single cut is retained by the single-spectral part of
the cut projection:
$$
F^{(\rho I),\mathrm{dress}}(s_1,s_2)
-F^{(\rho I),\mathrm{dress}}(s_1,0)
-F^{(\rho I),\mathrm{dress}}(0,s_2)
+F^{(\rho I),\mathrm{dress}}(0,0).
$$
# Result and limitations

The disjointness question has a definite answer: the supports are not
disjoint. The dressed rho term contains the $\pi\pi$ intermediate-state cut
also present in the physical two-pion unitarity construction. Its
narrow-pole double support is the same $\pi\pi$ cut paired with the
narrow-$I$ pole in the other photon variable. This statement concerns
analytic support; it does not say that the full physical production
amplitude equals the mixed model term.

The cut-replacement prescription is the matching/subtraction rule: retain
the physical two-pion cut once, subtract the mixed model's complete
$\pi\pi$ projection from its explicit addition, and include the resulting
rho-cut-free isoscalar-pole/subtraction remainder once. The projection
includes both photon placements, the double density, and the regular
single-rho-cut terms. For a finite cutoff the omitted rho tail is not
removed by this rule and must be supplied or explicitly left outside the
approximation.

The analytic prescription is complete, but its numerical input is not:
source-resolved complex $F_\pi^V$ and $a\gamma^*\to\pi\pi$ amplitudes, the
rho self-energy/residue convention, the Omnès subtraction and left-hand-cut
functions, and the cutoff are still needed to evaluate the physical
two-pion contribution. An already assembled two-pion TFF must also expose
its subtraction functions so the mixed isoscalar-pole remainder is assigned
only once.

# Reproducibility

Run the following commands:

    .venv/bin/python scripts/derive_alp_tff_rho_isoscalar_overlap.py
    .venv/bin/python scripts/derive_alp_tff_rho_isoscalar_overlap.py --test

The script verifies the pole-plus-regular decomposition, rho jump sign,
both ordered double discontinuities, the subtracted cut projection and its
finite-cutoff remainder, the single-rho-cut regular term, and the
generator/residue normalization. It uses symbolic quantities only and
assigns no physical production residues, phases, or subtraction functions.

