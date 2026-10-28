---
title: Background-field construction of ALP WZW interactions
status: literature
created: 2026-09-23
updated: 2026-10-28
bibtex_keys: [Bai:2025fvl]
provenance: literature
---

# Background Fields for ALP WZW Interactions

Bai et al. embed the derivative ALP--quark couplings by shifting left- and right-handed background fields,

$$
\mathbb B_{L/R}^\mu\to\mathbb B_{L/R}^\mu+k_{L/R}\frac{\partial^\mu a}{f_a},
$$

and use the full gauged Wess--Zumino--Witten functional together with the required counterterm. The counterterm is needed for fundamental gauge invariance, anomaly matching, and cancellation of auxiliary chiral-rotation parameters. [@Bai:2025fvl, Secs. 2--3, Eqs. (2.4)--(2.5) and (3.1)--(3.7)]

In this construction vector and axial-vector mesons are fictitious background gauge fields above confinement and are treated as dynamical below the QCD scale. This step supplies their anomalous couplings, but a complete dynamical model still requires kinetic, mass, even-parity, symmetry-breaking, and width terms. [@Bai:2025fvl, discussion following Eq. (2.8)]

The three-flavor construction includes the anomalous $U(1)_A$ sector and demonstrates auxiliary-phase cancellation without imposing $\operatorname{Tr}\kappa_q=1$. It agrees with the HLS pion-shift result for the vector-like vector-meson interactions considered there, while extending the background-field treatment to axial-vector mesons and the associated vector-like rotation phases. [@Bai:2025fvl, Sec. 3 after Eq. (3.7)]

The phenomenological study adopts approximations including idealized flavor assignments, $g\simeq\sqrt{12\pi}$, a fixed approximate $\eta$--$\eta'$ mixing pattern, and leading-order large-$N_c$ treatment of the topological susceptibility. Its benchmark decay predictions are therefore a consistency demonstration, not a global fit establishing precision agreement with all vector-meson data. [@Bai:2025fvl, Sec. 2 and Table 1]

## Full gauged functional for meson-containing vertices

Appendix A, Eq. (A11), gives the complete gauged functional as the ungauged
five-dimensional term $\Gamma_0(U)$, the four-dimensional trace below, and
the source counterterm $\Gamma_c$ [@Bai:2025fvl]:
$$
\begin{aligned}
\Gamma_{\rm WZW}^{\rm full}
=\Gamma_0(U)+C\int\operatorname{Tr}\Big\{\,&
 A_L\alpha^3+A_R\beta^3
-\frac{i}{2}\big[(A_L\alpha)^2-(A_R\beta)^2\big]\\
&+i\big(A_LUA_RU^\dagger\alpha^2
-A_RU^\dagger A_LU\beta^2\big)\\
&+i\big(dA_R\,dU^\dagger A_LU-dA_L\,dU\,A_RU^\dagger\big)\\
&+i\big[(dA_LA_L+A_LdA_L)\alpha
+(dA_RA_R+A_RdA_R)\beta\big]\\
&+A_L^3\alpha+A_R^3\beta
-(dA_LA_L+A_LdA_L)UA_RU^\dagger\\
&+(dA_RA_R+A_RdA_R)U^\dagger A_LU
+A_LUA_RU^\dagger A_L\alpha
+A_RU^\dagger A_LUA_R\beta\\
&+i\Big[A_L^3UA_RU^\dagger-A_R^3U^\dagger A_LU
-\frac12(UA_RU^\dagger A_L)^2\Big]\Big\}
+\Gamma_c,
\end{aligned}
$$
Products of differential forms are wedge products,
$C=N_c/(48\pi^2)$, $A_{L,R}=\mathbb A_{L,R}+\mathbb B_{L,R}$,
$\alpha=dU\,U^\dagger$, and $\beta=U^\dagger dU$. The counterterm
$\Gamma_c(\mathbb A_{L/R},\mathbb B_{L/R})$ is given separately in Eq. (A4);
it contains external sources but no $U$. It therefore has no Goldstone-field
vertex. For an ALP insertion, the background currents are shifted by
$\mathbb B_{L/R}\to\mathbb B_{L/R}+\mathbf k_{L/R}(a)\,da/f$
[@Bai:2025fvl, Eq. (16)]. The ALP--two-Goldstone--photon vertex follows from
expanding this full source functional with one photon source, one ALP
background source, and two Goldstone fields.

## Convention comparison with Alda:2026nrd

Equation (A10) uses a pion kinetic prefactor $f_{\pi,\mathrm{Bai}}^2/8$;
Alda:2026nrd uses $f_\pi^2/4$ for the same trace-normalized kinetic term.
Matching canonical pion normalization gives
$f_{\pi,\mathrm{Bai}}=\sqrt2 f_\pi$. With this conversion, the charged-pion
exponent coefficients agree and there is no residual $\sqrt2$ in the
$a\gamma\pi^+\pi^-$ coefficient [@Bai:2025fvl, Eq. (A10);
@Alda:2026nrd, Sec. III and App. A].

The external-current mapping in the calculation uses the Alda vielbein at
zero pion field to identify the left/right source difference with
$2\widehat c_q^A\,da/f_a$. Expansion of Eq. (A11) in the chosen exterior
orientation gives the raw source scalar for $\pi^0\to\gamma\gamma$. Comparing
it with the signed canonical TFF requires both $E=-R$ and the indexed
external-pion phase; this fixes the common A11 multiplier to $-1$ in the
project convention. The project $a\to\gamma\gamma$ contact-plus-WZW/mixing
coefficient is invariant under Alda:2026nrd's chiral source redefinition
with this relative sign. Bai:2025fvl does not state a numerical
$\epsilon^{0123}$ convention and footnote 20 notes a sign difference between
two cited WZW conventions, so the signed pion match is needed to translate
Eq. (A11) [@Bai:2025fvl; @Alda:2026nrd].
