---
title: Strong tree amplitude for an ALP to a narrow vector and two pseudoscalars
status: derived
created: 2026-10-21
updated: 2026-10-21
bibtex_keys: [Alda:2026nrd, Ecker:1988te, Ecker:1989yg, Ruiz-Femenia:2003jdx]
provenance: analytical-derivation
---

# Scope and inputs

This record derives the strong tree amplitude for $a(p_a)\to V(k,\varepsilon_V)
P_1(p_1)P_2(p_2)$ with $V=\omega$ or $\phi$ treated as a narrow external
state with the fit-derived near-ideal neutral-vector generator; it uses its
real on-shell LSZ state, not a Breit--Wigner factor. Natural units and the
metric $(+---)$ are used, with
$\epsilon^{0123}=+1$. All meson masses, vector masses, vector source residues,
and odd coefficients are those of the central provisional fit where relevant;
no additional $\kappa_V$ interaction is included. The ALP sources are already
in Alda:2026nrd's no-mixing basis [@Alda:2026nrd]. Define
$$
\chi=2B_0m_q,\qquad
\mathcal U_a=\frac{2c_q^A}{f_a},\qquad
\mathcal Y_a=-\frac{2}{f_a}\{c_q^P,\chi\},\qquad
\mathcal U_i=-\frac{\sqrt2}{F}\lambda_i^\dagger,\qquad
\mathcal Y_i=-\frac{\sqrt2}{F}\{\lambda_i^\dagger,\chi\}.
$$
Here $F=f_\pi$ at the retained leading chiral order, $c_q^A,c_q^P$ are the
no-mixing ALP matrices, and $\lambda_i^\dagger$ is the outgoing-state
pseudoscalar generator. The no-mixing vector-source matrix is
$c_q^{VS}=c_q^V-c_q^S$, with the no-mixing scalar source set to zero. The
even coupling $G_V$ is not among the parameters of the provisional odd-sector
fit, so it remains an explicit symbolic input and the amplitude is linear in
$G_V$. The odd action contains
the original seven $VJP$ and four $VVP$ operators only
[@Ruiz-Femenia:2003jdx].

The executable calculation is
[`scripts/derive_alp_vector_two_pseudoscalar.py`](../../../scripts/derive_alp_vector_two_pseudoscalar.py):

```text
.venv/bin/python scripts/derive_alp_vector_two_pseudoscalar.py
```

It derives the two-vector odd numerator, inverts the six-component
antisymmetric-tensor Hessian on the $VPP$ current, reduces all three odd
exchange contractions to one Lorentz tensor, derives the even scalar-source
exchange vertices, and checks the fit-null projections. Inputs are the
[no-mixing $aVV$ record](../../rcht/odd-parity/alp-vv-amplitude.md),
[the Alda:2026nrd-normalized $VPP$ and $VPa$ vertices](../../rcht/foundations/rcht-vpp-vpa-photon-pa-expansion.md),
[the original odd-operator fit](../../rcht/odd-parity/odd-preliminary-fit.md),
and the [fit-null audit](../../rcht/odd-parity/alp-target-amplitude-audit.md).

# Odd-parity exchange amplitude

Set $s=(p_1+p_2)^2$, $t=(k+p_2)^2$, $u=(k+p_1)^2$, and
$p_a=k+p_1+p_2$. For each Hermitian normalized internal vector generator
$T_R$ with mass $M_R$, define the following numerator for a pseudoscalar
source $X$ with mass $m_X$ and complementary vector virtuality $r^2$:
$$
\mathcal N_R(X;r^2)=
4\left[
d_1m_X^2\,\frac{\operatorname{tr}(\{T_V,T_R\}\mathcal U_X)}{2}
-d_3(m_X^2-M_V^2-r^2)\,\frac{\operatorname{tr}(\{T_V,T_R\}\mathcal U_X)}{2}
-4d_2\,\frac{\operatorname{tr}(\{T_V,T_R\}\mathcal Y_X)}{2}
\right].
$$
Here $M_V$ is the external vector mass; for the ALP source take
$X=a$, $m_X=m_a$, $r^2=s$, $\mathcal U_X=\mathcal U_a$,
$\mathcal Y_X=\mathcal Y_a$. For an outgoing pseudoscalar $P_i$, take
$X=P_i$, $m_X=m_i$, $r^2$ as shown, $\mathcal U_X=\mathcal U_i$,
$\mathcal Y_X=\mathcal Y_i$. This expression includes both vector assignments
at the $VVP$ vertex. The $d_4$ contraction is zero on these routed tensor
currents.

With the canonical tensor convention and the script's $S=1+iT$ phases, the
odd exchange contribution is
$$
\begin{aligned}
\mathcal M_{\rm odd}={}&
\epsilon_{\mu\nu\rho\sigma}
\varepsilon_V^\mu p_1^\nu p_2^\rho k^\sigma
\sum_{T_R}\Bigg[
-\frac{\sqrt2G_V\,\mathcal N_R(a;s)
\operatorname{tr}(T_R[\lambda_1^\dagger,\lambda_2^\dagger])}
{F^2M_V(s-M_R^2)}\\
&\quad+\frac{2G_V\,\mathcal N_R(P_2;t)
\operatorname{tr}(T_R[\lambda_1^\dagger,c_q^A])}
{Ff_aM_V(t-M_R^2)}
-\frac{2G_V\,\mathcal N_R(P_1;u)
\operatorname{tr}(T_R[\lambda_2^\dagger,c_q^A])}
{Ff_aM_V(u-M_R^2)}\Bigg].
\end{aligned}
$$
The first term is $a\to VT_R^*\to V P_1P_2$. The other two are
$a\to P_1T_R^*\to P_1VP_2$ and
$a\to P_2T_R^*\to P_2VP_1$. Their determinant contractions reduce to the
displayed common epsilon tensor. The full tensor propagator was used; replacing
it by a Proca pole alone would omit its local part and give an incorrect extra
virtuality factor [@Ecker:1989yg].

# Even-parity scalar-source contribution

Write $c_q^{VS}$ for the no-mixing vector-source matrix. The allowed flavor
rotation that removes its pure-gradient source shifts the scalar mass matrix
by $-ia[c_q^{VS},\chi]/f_a$. With no additional vector-mass spurion operator,
the resulting $aP_iP_R$ scalar vertex followed by the even $VPP$ vertex gives
the complete $c_q^{VS}$ contribution in this rotated representation:
$$
\begin{aligned}
\mathcal M_{VS}={}&-
\frac{iG_VM_V}{\sqrt2F^2f_a}
\sum_{\lambda_R}\Bigg[
\frac{
\operatorname{tr}(\{\lambda_2^\dagger,\lambda_R\}
[c_q^{VS},\chi])
\operatorname{tr}(T_V[\lambda_R,\lambda_1^\dagger])
}{u-m_R^2}\, (\varepsilon_V\cdot p_1)\\
&\qquad\qquad+
\frac{
\operatorname{tr}(\{\lambda_1^\dagger,\lambda_R\}
[c_q^{VS},\chi])
\operatorname{tr}(T_V[\lambda_R,\lambda_2^\dagger])
}{t-m_R^2}\, (\varepsilon_V\cdot p_2)
\Bigg].
\end{aligned}
$$
The sum runs over canonically normalized pseudoscalar mass eigenstates
$\lambda_R$ with masses $m_R$. This is a parity-even Lorentz structure and
must be kept separate from $\mathcal M_{\rm odd}$. In the original source
variables the same result is distributed among the direct $G_VaVPP$ contact
and pseudoscalar exchanges; do not add that contact to the rotated expression.
For diagonal $c_q^{VS}$, diagonal $\chi$, and $V=\omega,\phi$, both
$[c_q^{VS},\chi]$ and the covariant-vector kinetic $aVV$ trace vanish, and
the direct contact cancels the derivative-source pseudoscalar exchanges.

No $\kappa_V$ vector-mass interaction is included, as requested. If that
operator is restored in a later calculation, its ALP variation and the
resulting vector exchange must be derived and added consistently.

# Fit null directions and vector poles

The exact four-dimensional null space of the archived original eleven-operator
fit is spanned by
$$
\Delta c_1=\Delta c_5,\qquad
\Delta c_2=\Delta c_5=\Delta c_6,\qquad
\Delta c_7,\qquad \Delta d_4.
$$
The script compares this span with the stored central single-trace fit basis;
the largest projected discrepancy is $5.6\times10^{-15}$. The complete
strong odd exchange row annihilates all four directions exactly. The $c_i$
directions are absent at zero electromagnetic order. Electromagnetic VJP
conversion into a massive vector would depend on the first direction, so it
is excluded from this strong amplitude; no fit-selected representative for
that null is needed here.

Sum every internal vector in the retained vector nonet whose flavor trace is
nonzero. The internal pole set is selected by the flavor vertices and
kinematics, not by whether a resonance is broad or narrow. It therefore
includes narrow $\omega$ or $\phi$ poles when their $VPP$ trace and threshold
permit them, as well as broad $\rho$ or $K^*$ poles when allowed. For a
strict tree expression the denominators are $r^2-M_R^2+i0$; near an internal
resonance replace $M_R^2-i0$ by the chosen dressed pole denominator and state
the width prescription. The external narrow $\omega$ or $\phi$ is on shell
and receives no propagator denominator.

This result is tree level and retains the original odd eleven-operator basis,
the fitted even inputs, the leading $G_V$ interaction, and no additional
$\kappa_V$ term. It does not include electromagnetic conversions, loops,
rescattering, or vector multiplets beyond the fitted nonet. A numerical width
also requires the chosen internal line shapes and the physical channel's
phase space.
