---
title: Partonic quark-bilinear operator generators at arbitrary ALP order
status: derived
created: 2026-11-03
updated: 2026-11-03
bibtex_keys: [Alda:2026nrd, Bonnefoy:2022rik]
provenance: analytical-derivation
---

# Partonic quark-bilinear operator generators at arbitrary ALP order

## Assumptions and scope

Take a neutral pseudoscalar $a$ and partonic quark fields in specified
representations of the unbroken gauge group. Impose the ALP shift symmetry on
the invariant part of the EFT; explicit breaking is carried by spurions with
specified shift transformations. The field is odd under parity, but no
parity- or CP-invariance condition is imposed on its interactions. This note
stays at the partonic operator level.
The distinction between shift-symmetric derivative couplings and
spurion-dependent fermion couplings is also made in
`Bonnefoy:2022rik` [@Bonnefoy:2022rik]. Its dimension-five operator count
belongs to its electroweak theory and is not imported here
[@Bonnefoy:2022rik, Secs. 2--3].

There is a constructive recursion for every ALP insertion number $n$ and
normalized operator dimension $d$. It gives an overcomplete generating set.
An independent finite basis additionally requires fixing the unbroken gauge
representations, the explicit-breaking spurions and their shift charges, and
a maximum $d$, followed by quotienting integration-by-parts, equations-of-
motion, and algebraic relations. Canonical dimension is the only dimension
grading used here.

## Generating set at fixed ALP insertion number

At any insertion number $n$, include all local Lorentz and gauge contractions
of the following quark-bilinear classes:

- Same-chirality bilinears $\bar q_L\Gamma q_L$ and
  $\bar q_R\Gamma q_R$. Their elementary Dirac structures are vector/axial
  type. They may contain any allowed QCD/QED gauge-covariant derivatives, field
  strengths, mass spurions, and $n$ ALP insertions.
- Opposite-chirality bilinears $\bar q_R\Gamma q_L$ and their Hermitian
  conjugates. Their elementary Dirac structures are scalar, pseudoscalar,
  or tensor type. Their coefficient matrices carry the corresponding
  left/right flavor indices and may contain the same covariant insertions.

In a manifestly shift-symmetric basis, each ALP insertion occurs through
derivatives of $a$, such as $\partial_{\mu_1}\cdots
\partial_{\mu_k}a/f_a$. Undifferentiated $a$ can also occur in an
explicit-breaking spurion; its allowed dependence is restricted by that
spurion's shift transformation. An unspecified breaking function leaves
infinitely many independent Taylor coefficients. An ALP-dependent field
redefinition can generate undifferentiated powers $a/f_a$ even from the
manifestly shift-symmetric terms. Those coefficients are fixed by the
field redefinition and are not new explicit-breaking couplings. The general
transformation laws and invariant combinations are in the [chiral
redefinition record](alp-general-chiral-redefinition-invariants.md).

In particular, the leading current and mass-source terms used by
`Alda:2026nrd` are only a small subset of this generating set. Higher
derivative bilinears, tensor bilinears, and bilinears with gauge-field
insertions carry independent Wilson coefficients until a symmetry, matching
condition, or chosen truncation relates them.

All derivatives in this construction are the gauge-covariant derivatives
acting on partonic fields in the chosen gauge theory.

## Induction in ALP insertion number

Let $\mathcal B_{n,d}$ denote the overcomplete manifest-basis set of
open-Lorentz-index quark-bilinear monomials with $n$ ALP insertions and
normalized mass dimension $d$, including the factors of $1/f_a$. The seed $\mathcal B_{0,d}$
contains all ALP-free quark bilinears of dimension $d$: quark fields,
gauge-covariant derivatives, field strengths, allowed mass or Higgs
spurions, and Dirac matrices. Keep Lorentz indices open until all insertions
have been made. At each step form
$$
\mathcal B_{n,d}=\operatorname{span}\!\left[
\bigcup_{k\geq1}\left(
\frac{\partial_{\mu_1}\cdots\partial_{\mu_k}a}{f_a}
\otimes\mathcal B_{n-1,d-k}\right)
\;\cup\;
\bigcup_{r=1}^{n}\left(
\left(\frac{a}{f_a}\right)^r S_r
\otimes\mathcal B_{n-r,d-d(S_r)}\right)
\right].
$$
Here the first union inserts one shift-symmetric ALP jet, which has grade one
and dimension $k$ after division by $f_a$. In the second union, $S_r$ runs
over every allowed explicit-breaking spurion Taylor coefficient of grade
$r$ and dimension $d(S_r)$; $a/f_a$ is dimensionless. A term is omitted if
its residual dimension is negative. Each tensor product includes every
allowed insertion position and retains open Lorentz and gauge indices. Project
onto gauge singlets and contract Lorentz indices only after tensor
generation. For exact shift symmetry, omit the second union. If the breaking
spurions and their shift transformations are not specified, that union is
an unrestricted family, so no finite list follows. If $\Delta$ is the
canonical dimension of the numerator before extracting $1/f_a^n$, then
$d=\Delta-n$.

In a manifestly shift-symmetric basis, every ALP insertion contributes at
least one unit to $d$. Since a quark bilinear starts at dimension three,
fixed $d$ allows only $n\leq d-3$. Explicit-breaking blocks $(a/f_a)^rS_r$
have no dimension cost from $a/f_a$, so arbitrary $r$ can occur at the same
$d$; their coefficients are related only after the spurion shift
representation and its invariant functions have been specified. The
ALP-dependent field-redefinition orbit is another source of $a/f_a$ powers
at fixed $d$; its coefficients are correlated and derived below.

Why keep open indices? Removing one ALP jet from an order-$n$ Lorentz scalar
can expose an uncontracted vector or higher-rank tensor. A scalar-only
$\mathcal B_{n-1,d}$ would have discarded it and would not seed every scalar
at order $n$. After the recursion, impose gauge singlet conditions, contract
Lorentz indices, and quotient by integration by parts, the equations of
motion at the chosen EFT order, Clifford and Schouten identities, Bianchi
identities, and commutators of covariant derivatives. Reduce at each $(n,d)$
with the relations expanded consistently in the ALP grading; do not assume
that multiplying a reduced scalar basis at $n-1$ automatically gives an
independent basis at $n$.

Here $\Lambda$ denotes the suppression scale for operators above dimension
four; each flavor matrix and dimensionless Wilson coefficient is independent
for every allowed operator.

For orientation, in a broken-phase QCD+QED realization where the indicated
bilinears are gauge invariant:

- $n=0$: $\bar q_L i\!\not\!Dq_L$ and
  $\bar q_R M q_L+\mathrm{h.c.}$ are dimension-four seeds; a dipole
  $\bar q_R\sigma^{\mu\nu}T^Aq_LG^A_{\mu\nu}+\mathrm{h.c.}$ is a
  dimension-five seed.
- $n=1$: $(\partial_\mu a/f_a)\bar q_L\gamma^\mu C_Lq_L$ is dimension four;
  $(a/f_a)\bar q_R S_1q_L+\mathrm{h.c.}$ is also dimension four when
  $S_1$ has mass dimension one.
- $n=2$: $(\partial_\mu a\,\partial^\mu a/f_a^2)\bar q_R C q_L$
  has dimension five and coefficient proportional to $1/\Lambda$;
  $(a/f_a)^2\bar q_R S_2q_L+\mathrm{h.c.}$ has dimension four for a
  dimension-one spurion $S_2$.

These are representatives in the generating set, not claims of independence
after equations of motion or integration by parts. In an electroweak-gauge-
invariant theory, a chirality-changing bilinear must also contain fields such
as the Higgs doublet needed to make a gauge singlet. The same recursion
applies with those gauge representations and spurions in the seed.

As a closed but restricted example, omit quark and gauge-field derivatives,
field strengths, and explicit-breaking spurions, and allow only first ALP
gradients. Then Lorentz invariance gives, for even $n\geq2$,
$$
\Lambda^{1-n}
\left(\frac{\partial_\mu a\,\partial^\mu a}{f_a^2}\right)^{n/2}
\bar q(C_S+i\gamma_5C_P)q,
$$
and, for odd $n\geq1$,
$$
\Lambda^{1-n}
\left(\frac{\partial_\mu a\,\partial^\mu a}{f_a^2}\right)^{(n-1)/2}
\frac{\partial_\mu a}{f_a}
\bar q\gamma^\mu(C_V+C_A\gamma_5)q.
$$
The $C$ coefficients are arbitrary flavor matrices compatible with the gauge
representation. These representatives recur from $n$ to $n+2$ by adding
$(\partial a)^2/(f_a^2\Lambda^2)$. This simple recursion does not exhaust
the full basis, because higher ALP derivatives and derivatives or field
strengths on the quark bilinear generate other Lorentz tensors.

## Reproducible recursion

Run the symbolic ALP-grade/canonical-dimension generator through $n=2$,
$d=5$ with:

```sh
.venv/bin/python scripts/derive_alp_partonic_operator_induction.py --order 2 --dimension 5
```

The script performs coefficient extraction for the recursive grading. Lorentz
contractions, gauge representations, insertion positions, and reductions by
equations of motion or integration by parts remain symbolic; its path counts
are not counts of independent operators.

