---
title: Flavour-violating ALP decays into two distinct pseudoscalar mesons
status: derived
created: 2026-10-15
updated: 2026-10-15
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Scope and inputs

These are tree-level amplitudes at leading chiral order and first order in
$1/f_a$, from the vector and scalar ALP sources in the even-intrinsic-parity
meson Lagrangian. Units are natural. No resonance exchange, rescattering,
electromagnetic correction or weak neutral-kaon mixing is included.

Inputs are the [two-pseudoscalar expansion](alp-phi2-vector-scalar-lo-chipt-derivation.md),
the [source-verified Appendix-A definitions](../../extracts/Alda:2026nrd/appendix-a-expansion.md),
and the [neutral-meson generators](../../extracts/Alda:2026nrd/neutral-meson-generators.md).
The archived source `main_arXiv.tex`, lines 264--274 and 290, defines
$\bm{c_q^{VS}}=\bm{c_q^V}-\bm{c_q^S}$; lines 2042--2046 fix
$\mathbb C_{\bm\varphi}[\bm A]=[\bm\varphi,\bm A]$; lines 2142--2146
define the meson mass operator and its eigenvectors. [@Alda:2026nrd]

The executable calculation is
[derive_alp_flavour_violating_two_pseudoscalars.py](../../../scripts/derive_alp_flavour_violating_two_pseudoscalars.py).
It derives the interactions by series expansion, extracts both meson
placements, crosses the all-incoming vertex, and evaluates the outgoing
flavour traces. Its important steps and results are printed with SymPy's
LaTeX renderer. The corrected preceding expansion script remains executable.

# Interaction and outgoing-state convention

The source-defined commutator ordering gives

$$
\mathcal L\big|_{a\Phi^2}=
-\frac{i\partial_\mu a}{2f_a}
\left\langle\partial^\mu\bm\Phi[\bm{c_q^V},\bm\Phi]\right\rangle
-\frac{iB_0a}{2f_a}
\left\langle\bm\Phi^2[\bm{c_q^S},\bm{m_q}]\right\rangle.
$$

In particular, the vector sign follows the explicit vielbein in the source
at line 2109. The reversed commutator definition formerly used in the input
record and script has been corrected there. This calculation audits this
two-pseudoscalar sector; it does not establish that unrelated derivations
using the old definition have been rechecked.

Use $\bm\Phi=\sum_P\bm{\lambda_P}P$ with the source's generator and field
normalizations. The Fourier convention is $e^{-ik\cdot x}$. The vertex is
$i\mathcal M$, with $S_{fi}=i(2\pi)^4\delta^{(4)}(p_a-p_i-p_j)\mathcal M$.
For a decay, $k_a=p_a$, $k_i=-p_i$, $k_j=-p_j$.
The matrix for an outgoing meson is $\bm{\lambda_i}^\dagger$: an outgoing
$K^+$, for example, is created by the $K^-$ field in $\bm\Phi$.
A flavour-field monomial cannot be identified with the outgoing channel
without this conjugation.

For distinct final mesons the executed result is

$$
i\mathcal M(a\to P_iP_j)=\frac{1}{2f_a}\left[
-(m_i^2-m_j^2)
\left\langle\bm{c_q^V}[\bm{\lambda_i}^\dagger,\bm{\lambda_j}^\dagger]\right\rangle
+B_0\left\langle
\{\bm{\lambda_i}^\dagger,\bm{\lambda_j}^\dagger\}
[\bm{c_q^S},\bm{m_q}]\right\rangle\right].
$$

This includes both meson placements with no identical-particle factor.
The unequal-mass contraction is
$p_a\cdot(p_j-p_i)=m_j^2-m_i^2$. No assumption $m_i=m_j$ is made.

# Field-redefinition invariance and meson masses

Alda's vector field redefinition shifts both $\bm{c_q^V}$ and
$\bm{c_q^S}$ by the same $\bm\delta$. Choosing
$\bm\delta=-\bm{c_q^S}$ removes the scalar source and leaves
$\bm{\widehat c_q^V}=\bm{c_q^{VS}}$.
The SymPy calculation retains the singlet term in

$$
\mathbb M(\bm A)=B_0\{\bm{m_q},\bm A\}
+\frac{m_0^2}{3}\langle\bm A\rangle\bm 1,
\qquad \mathbb M(\bm{\lambda_i})=m_i^2\bm{\lambda_i}.
$$

For eigenstates of this same LO operator it derives

$$
(m_i^2-m_j^2)
\left\langle\bm{c_q^S}[\bm{\lambda_i}^\dagger,\bm{\lambda_j}^\dagger]\right\rangle
=B_0\left\langle
\{\bm{\lambda_i}^\dagger,\bm{\lambda_j}^\dagger\}
[\bm{c_q^S},\bm{m_q}]\right\rangle.
$$

The singlet contribution cancels in this identity, without setting $m_0^2$
to zero. The complete vector-plus-scalar amplitude is therefore

$$
\boxed{
i\mathcal M(a\to P_iP_j)=
-\frac{m_i^2-m_j^2}{2f_a}
\left\langle\bm{c_q^{VS}}
[\bm{\lambda_i}^\dagger,\bm{\lambda_j}^\dagger]\right\rangle.}
$$

This equality requires LO masses and generators from the same mass operator.
Arbitrary physical masses and an independently chosen neutral mixing angle
do not establish it. Inserting measured phase-space masses or measured mass
differences is an additional approximation beyond this exact LO statement.

# Electrically neutral flavour-changing channels

Take the quark order $(u,d,s)$ and Hermitian couplings. The complex element
$\bm{c_q^{VS}}_{ds}$ is arbitrary, with
$\bm{c_q^{VS}}_{sd}=(\bm{c_q^{VS}}_{ds})^*$.
The amplitudes for strangeness $+1$ select the $ds$ entry, not the $sd$
entry that appears when extracting the corresponding unconjugated field
monomials.

For any real diagonal neutral generator, including the full isospin-breaking
LO $\pi^0,\eta,\eta'$ generators, the directly evaluated expression is

$$
i\mathcal M(a\to K^0P)=\frac{1}{\sqrt2f_a}\left[
-\bm{c_q^V}_{ds}
\big((\bm{\lambda_P})_{dd}-(\bm{\lambda_P})_{ss}\big)(m_{K^0}^2-m_P^2)
+B_0\bm{c_q^S}_{ds}(m_s-m_d)
\big((\bm{\lambda_P})_{dd}+(\bm{\lambda_P})_{ss}\big)\right].
$$

For LO eigenstates this becomes

$$
i\mathcal M(a\to K^0P)=
-\frac{\bm{c_q^{VS}}_{ds}}{\sqrt2f_a}
\big((\bm{\lambda_P})_{dd}-(\bm{\lambda_P})_{ss}\big)
(m_{K^0}^2-m_P^2).
$$

This is the formula for exact LO isospin breaking. Rounded generator entries
must not be treated as exact eigenvectors in the field-redefinition test.

For the familiar one-angle form, impose $m_u=m_d$ and use
$\bm{\lambda_\eta}=\cos\theta_P\bm{\lambda_8}-\sin\theta_P\bm{\lambda_0}$,
$\bm{\lambda_{\eta'}}=\sin\theta_P\bm{\lambda_8}+\cos\theta_P\bm{\lambda_0}$,
where $\bm{\lambda_8}=\operatorname{diag}(1,1,-2)/\sqrt3$ and
$\bm{\lambda_0}=\sqrt{2/3}\bm1$. The angle and neutral masses must diagonalize
the same LO mass operator. In the scalar-free field basis the four results are

| Outgoing channel | $i\mathcal M$ |
|---|---|
| $K^+\pi^-$ | $-\bm{c_q^{VS}}_{ds}(m_{K^+}^2-m_{\pi^+}^2)/f_a$ |
| $K^0\pi^0$ | $+\bm{c_q^{VS}}_{ds}(m_{K^0}^2-m_{\pi^0}^2)/(\sqrt2f_a)$ |
| $K^0\eta$ | $-\sqrt{3/2}\,\bm{c_q^{VS}}_{ds}\cos\theta_P(m_{K^0}^2-m_\eta^2)/f_a$ |
| $K^0\eta'$ | $-\sqrt{3/2}\,\bm{c_q^{VS}}_{ds}\sin\theta_P(m_{K^0}^2-m_{\eta'}^2)/f_a$ |

The charged result also holds without the isospin approximation when its
charged LO masses are used. All four conjugate channels $K^-\pi^+$,
$\bar K^0\pi^0$, $\bar K^0\eta$, $\bar K^0\eta'$ obey

$$
i\mathcal M(a\to\bar P_i\bar P_j)
=-\big[i\mathcal M(a\to P_iP_j)\big]^*,
\qquad
\mathcal M(a\to\bar P_i\bar P_j)
=\mathcal M(a\to P_iP_j)^*.
$$

Thus their $i\mathcal M$ values follow from the table by changing the overall
sign and conjugating $\bm{c_q^{VS}}_{ds}$. $K^0$ and $\bar K^0$ here are
strong-interaction flavour states. Predictions for $K_S,K_L$ require a separate
weak-mixing projection, which is not performed here.

# Widths, checks and next step

For any one distinct-state channel above threshold,

$$
\Gamma(a\to P_iP_j)=
\frac{\sqrt{[m_a^2-(m_i+m_j)^2][m_a^2-(m_i-m_j)^2]}}
{16\pi m_a^3}\,|\mathcal M(a\to P_iP_j)|^2,
\qquad m_a>m_i+m_j.
$$

There is no additional $1/2!$ factor. The conjugate channels have equal widths
at this order. Diagonal-source particle--antiparticle channels vanish on shell
when particle and antiparticle masses agree.

Reproduce the derivation and checks from the repository root:

```text
.venv/bin/python scripts/derive_alp_flavour_violating_two_pseudoscalars.py
.venv/bin/python scripts/derive_alp_flavour_violating_two_pseudoscalars.py --test
.venv/bin/python scripts/derive_alp_phi2_lo_chipt.py
```

The eight checks cover raw trace reductions, all eight flavour-changing
channels, conjugate amplitudes, the exact mass-operator identity, charged
field-redefinition invariance, nondegenerate neutral eigenvectors with nonzero
singlet mass, diagonal-source zeros and the all-incoming vertex.
The LO amplitude calculation is complete. Further precision requires higher
chiral orders and final-state interactions; measured-mass substitutions and
weak neutral-kaon projections must be stated explicitly when introduced.
No manuscript change is made by this calculation.
