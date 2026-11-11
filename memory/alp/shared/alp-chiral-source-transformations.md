---
title: Finite chiral transformations of meson sources and topological coefficient
status: derived
created: 2026-11-11
updated: 2026-11-11
bibtex_keys: [Alda:2026nrd, Bonnefoy:2022rik]
provenance: analytical-derivation
---

# Finite chiral transformations of meson sources and topological coefficient

## Inputs, conventions and assumptions

Use natural units, $x=a/f_a$, and the restricted partonic theory, source
transformation laws and transport equations in the [general chiral-redefinition
record](alp-general-chiral-redefinition-invariants.md). In particular, the
flavor matrices need not commute, $g_{R,L}(0)=\mathbf1$, and
$W_\eta^\dagger W_\eta=\mathbf1$. The theory omits dipoles,
higher-derivative quark bilinears and interactions with two or more ALP
derivatives. The source convention and leading meson theory follow
`Alda:2026nrd` [@Alda:2026nrd]; the restriction of the partonic operator set
is also discussed in `Bonnefoy:2022rik` [@Bonnefoy:2022rik]. The anomalous
Jacobian is an input from the [QCD topological record](alp-qcd-topological-cp.md),
not independently derived here. The [dressed topological record](alp-dressed-topological-invariant.md)
provides the existing invariant coefficient and related fourth-order checks.

Write $r_\mu=(\partial_\mu x)C_R$, $\ell_\mu=(\partial_\mu x)C_L$ and
$\chi=2B_0M$ before dressing. The matrices $C_{R,L}$ are dimensionless,
$M$ has mass dimension one and $\chi$ mass dimension two.

## Source transformations and transport

1. Apply the partonic connection law to both sources and use
$\partial_\mu g_\eta=(\partial_\mu x)\partial_xg_\eta$:
$$
\begin{aligned}
r_\mu'&=g_Rr_\mu g_R^\dagger+ig_R\partial_\mu g_R^\dagger,\\
\ell_\mu'&=g_L\ell_\mu g_L^\dagger+ig_L\partial_\mu g_L^\dagger,\\
\chi'&=g_R\chi g_L^\dagger,\qquad U'=g_RUg_L^\dagger.
\end{aligned}
$$
The $U$ transformation is the conventional chiral transformation of the
Goldstone matrix. Two successive redefinitions compose as
$g_\eta=g_\eta^{(2)}g_\eta^{(1)}$.

2. The transport equations $\partial_x W_\eta=iC_\eta W_\eta$ give
$W_\eta'=g_\eta W_\eta$. Thus $\widehat U=W_R^\dagger UW_L$ and
$\widehat M=W_R^\dagger MW_L$ are unchanged. The dressed spurion is
$2B_0\widehat M$. Differentiating the three factors gives
$$
\partial_\mu\widehat U
=W_R^\dagger(\partial_\mu U-ir_\mu U+iU\ell_\mu)W_L
=W_R^\dagger(D_\mu U)W_L.
$$
This identity also permits independent spacetime dependence of $U$; only
$C_{R,L}$ and $g_{R,L}$ are restricted to depend on $x$.

## Singlet anomaly and the meson source

Use $\Theta$ in the dimensionless normalization of `Alda:2026nrd`, so
$\Theta=(4\pi/\alpha_s)\Theta_{\rm raw}$ for the coefficient multiplying
$G\widetilde G$ without the prefactor $\alpha_s/(4\pi)$. The Jacobian law,
with the logarithm continuous from $x=0$, is
$$
\Theta'=\Theta-\frac{1}{2i}\log\frac{\det g_R}{\det g_L}.
$$
For $g_R=g_L^\dagger=e^{i\theta(x)}\mathbf1$, with real arbitrary
$\theta(0)=0$, the determinant gives
$$
\Theta'=\Theta-N_f\theta,\qquad
\Theta_{\rm raw}'=\Theta_{\rm raw}-\frac{\alpha_s}{4\pi}N_f\theta.
$$
Locally on the same continuous branch, $U'=e^{2i\theta}U$ gives
$\operatorname{tr}\log u'=\operatorname{tr}\log u+iN_f\theta$.
Assuming the extension of the `Alda:2026nrd` linear source to a general
function, as stated in the [dressed topological record](alp-dressed-topological-invariant.md),
$$
X_{TC}=2i\bigl(\operatorname{tr}\log u+i\Theta(x)\bigr)
$$
is therefore invariant. This calculation treats the QCD topological
combination; it does not determine the full Wess-Zumino term.

## Reproducible calculation and checks

[`derive_alp_chiral_source_transformations.py`](../../../scripts/derive_alp_chiral_source_transformations.py)
uses the product rule for the source transformations and the derivative of
$\widehat U$, then substitutes two distinct noncommuting $U(2)$ redefinitions
with nonlinear generators. It computes exact polynomial coefficients through
$x^3$ for connection composition, transport covariance, invariant
$\widehat U$ and $\widehat M$, and the derivative identity. The examples use
fixed matrices; the all-orders identities use the product rule and are not
inferred from the examples. The script also computes the singlet determinant
shift, $X_{TC}$ cancellation and normalization conversion from the stated
Jacobian input. The branch choice is local and must be continued consistently.

```sh
.venv/bin/python scripts/derive_alp_chiral_source_transformations.py
.venv/bin/python -m unittest discover -s scripts -p test_derive_alp_chiral_source_transformations.py
```
