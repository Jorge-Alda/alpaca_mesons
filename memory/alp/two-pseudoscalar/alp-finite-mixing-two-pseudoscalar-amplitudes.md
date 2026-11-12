---
title: Pole projection for one and two pseudoscalars
status: derived
created: 2026-11-12
updated: 2026-11-12
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Pole projection for one and two pseudoscalars

## Inputs and assumptions

Use the CCWZ Lagrangian, source functions, Taylor derivatives, field normalization and aligned stationary vacuum of the [joint-vertex calculation](../shared/alp-lo-joint-vertex-expansion.md), with the conventions of `Alda:2026nrd` [@Alda:2026nrd]. Natural units and tree-level LO chiral perturbation theory are used. All momenta are incoming, including crossed outgoing momenta. The vertex $\mathcal V_{IJK}$ is the symmetrized derivative of the Lagrangian without its factor $i$, evaluated on the corresponding momenta. Its complete source and momentum dependence is provided by the linked CCWZ calculation. The normalized pole vectors, their assumptions and crossing conventions are given in the [pole-projected construction](../shared/alp-lo-pole-projected-tree-amplitudes.md). The calculation does not discard ALP components of neutral daughters or expand the parent pole.

## One-pseudoscalar mixing

The two-point function is the inverse of $D(s)=sK-M^2$ and its pole residue is $v^{(n)}v^{(n)\mathsf T}$. The off-diagonal two-point function describes mixing, not an independent on-shell decay between distinct nondegenerate mass eigenstates. For a pole with nonzero ALP component, its meson-to-ALP overlap is given analytically in the [operator-inverse record](../shared/alp-lo-finite-mixing-operator-inverse.md), with the general-source modification to the pole equation in the [joint-vertex record](../shared/alp-lo-joint-vertex-expansion.md). For a pole with vanishing ALP component, use its full pole vector rather than dividing by that component.

## Two-pseudoscalar amplitude

For parent $a_{\rm phys}$ and daughters $P_i,P_j$, the exact LO tree amplitude is

$$
\mathcal M(a_{\rm phys}\to P_i P_j)
=\sum_{I,J,K}v_I^{(a)}v_J^{(i)}v_K^{(j)}\mathcal V_{IJK}.
$$

The sums run over all meson and ALP directions. They include the $aPP$, $aaP$, and $aaa$ vertices and, where present, pure-meson cubic vertices. The derivative with respect to the fields includes both placements for identical or distinct daughters. Under simultaneous exchange of daughter directions and momenta the expression is unchanged. No factor $1/2$ enters the amplitude for identical daughters; the phase-space integration carries the usual $1/2!$.

For charged outgoing particles the field directions are conjugate to their outgoing flavour generators. Electric charge prevents charged mesons from mixing with the neutral ALP. An electrically neutral vector source may contain complex Hermitian $d$–$s$ entries; its $|\Delta S|=1$ contributions and conjugate channels are handled by the same contraction. For neutral kaons, $K^0$ and $\bar K^0$ are flavour directions, not automatically separate mass eigenstates when the ALP mixes with the $P_6,P_7$ system. Use the actual neutral pole vectors; weak $K_S,K_L$ projections require additional weak dynamics. The exact zero-ALP neutral-kaon pole and its assumptions are described in the [neutral-block record](../shared/alp-lo-neutral-block-reduction.md).

### Derivation for the neutral $K\pi$ channel

The following explicit $K^0\pi^0$ formula is a **first-order weak-mixing
result**, evaluated with unmixed external flavour states and isospin-symmetric
LO masses. Its explicit $1/f_a$ is not a resummation of neutral ALP--meson
mixing. At finite mixing the physical daughter poles need not be $K^0$ and
$\pi^0$, and the full pole contraction below is required.

Start from the CCWZ building blocks $u_\mu$ and $\chi_+$ and retain one
ALP field and two meson fields. The [stepwise CCWZ expansion](alp-phi2-vector-scalar-lo-chipt-derivation.md) gives

$$
\left.\mathcal L_{\rm LO}\right|_{a\Phi^2}=
-\frac{i\partial_\mu a}{2f_a}
\langle\partial^\mu\bm\Phi[\bm{c_q^V},\bm\Phi]\rangle
-\frac{iB_0a}{2f_a}\langle\bm\Phi^2[\bm{c_q^S},\bm{m_q}]\rangle.
$$

Differentiate with respect to both meson fields, cross the outgoing fields,
and use $p_a=p_i+p_j$ and the two daughter on-shell conditions. For
generators and quark-mass insertions from the same LO meson mass operator,
the scalar and vector terms combine as derived in the
[flavour-changing amplitude record](alp-flavour-violating-two-pseudoscalar-amplitudes.md).

For $P_i=K^0$ and $P_j=\pi^0$, use quark order $(u,d,s)$,
$\bm{\lambda_{K^0}}^\dagger=\sqrt2\begin{pmatrix}0&0&0\\0&0&0\\0&1&0\end{pmatrix}$ and
$\bm{\lambda_{\pi^0}}=\operatorname{diag}(1,-1,0)$. Crossing the outgoing kaon
therefore selects $\bm{\lambda_{K^0}}^\dagger$, and direct multiplication
gives

$$
\langle\bm{c_q^{VS}}[\bm{\lambda_{K^0}}^\dagger,\bm{\lambda_{\pi^0}}]\rangle
=-\sqrt2\,(c_q^{VS})_{ds}.
$$

Evaluating the quark-mass terms directly gives the flavour-state check
$$
i\mathcal M(a\to K^0\pi^0)=\frac{B_0(m_s-m_d)(c_q^{VS})_{ds}}{\sqrt2 f_a},\qquad
i\mathcal M(a\to\bar K^0\pi^0)=-\left[i\mathcal M(a\to K^0\pi^0)\right]^*.
$$
The conjugate relation assumes a Hermitian source. The check retains the
complex $ds$ coupling and its Hermitian conjugate, and does not include
weak $K^0$--$\bar K^0$ mixing.

For finite mixing, let $i$ and $j$ label the actual neutral daughter poles
that approach the chosen kaon and pion states in a specified unmixed limit.
Their complete tree-level LO result is the earlier contraction,

$$
\mathcal M(a_{\rm phys}\to P_iP_j)=
\sum_{I,J,K}v_I^{(a)}v_J^{(i)}v_K^{(j)}
\mathcal V_{IJK}(p_a,-p_i,-p_j),
$$

with the pole masses in the on-shell momenta and with every allowed neutral
meson and ALP direction included in each sum. The pole vectors and the cubic
vertex are computed from the same CCWZ source functions. Thus terms from
daughter ALP components and from parent meson components are retained; no
truncation in $1/f_a$ is made. For a degenerate kaon subspace, choose and
specify a basis of physical poles before assigning a flavour label. The
parent normalization is the $K$-normalized residue vector; for a pole with
nonzero ALP component it is given by $v_a^{(a)}=1/\sqrt{K_{aa}+2K_{ia}w_i^{(a)}+K_{ij}w_i^{(a)}w_j^{(a)}}$, with repeated meson indices summed, as in the [operator-inverse record](../shared/alp-lo-finite-mixing-operator-inverse.md). The
expression displayed for $K^0\pi^0$ above follows by taking the unmixed
external directions and keeping the term linear in $1/f_a$ under the
weak-mixing assumptions below; it is not a general finite-mixing channel
formula. A numerical value for a finite-mixing channel requires the source
functions and the chosen pole identification. The [explicit pole-aligned
calculation](alp-explicit-finite-mixing-neutral-amplitude.md) evaluates the
mixed neutral-pole amplitude in flavour components with vector, axial,
scalar and pseudoscalar sources, as well as the exact decoupled-kaon branch.
The [pole-aligned matrix result](alp-pole-aligned-neutral-matrix-amplitude.md)
uses the established operator inverse and gives the cubic contact with the
vector, axial, scalar and pseudoscalar sources explicitly transformed.

## Weak-mixing and charged limit

For nondegenerate poles, hold the source Taylor derivatives with respect to $x=a/f_a$ fixed as $f_a$ varies. Each ALP insertion then supplies one power of $1/f_a$. Assume also that the vacuum has no pure-meson cubic vertex. The meson admixture of the ALP-like parent and the ALP admixtures of meson-like daughters are then first order in $1/f_a$. The first nonzero term is

$$
\mathcal M(a_{\rm phys}\to P_iP_j)=\sum_{r,s}\mathcal V_{ars}\,v_a^{(a)}v_r^{(i)}v_s^{(j)}+O(f_a^{-3}),
$$

with the leading unmixed meson directions. The remainder assumes the stated power counting; it is not uniform near a degeneracy. Additional ALP derivatives of source functions can contribute at finite mixing even if absent in a linear-spurion special case.

For $K^+\pi^-$, with the same LO quark-mass relation used in the [charged vector-redefinition record](../shared/alp-charged-u3-vector-redefinition.md), the corresponding differentiated Lagrangian gives

$$
i\mathcal M(a_{\rm phys}\to K^+\pi^-)
=-\frac{B_0(m_s-m_d)(c_q^{VS})_{ds}}{f_a}\,v_a^{(a)}.
$$

This equality holds under that record's parity-even vacuum and cubic-vertex assumptions. Its exact parent normalization is $v_a^{(a)}=1/\sqrt{K_{aa}+2K_{ia}w_i^{(a)}+K_{ij}w_i^{(a)}w_j^{(a)}}$ on the applicable branch. The charged vector redefinition changes the off-shell contact but its contribution vanishes when both daughters are on shell. The general flavour-changing source and outgoing-state conventions are given in the [two-pseudoscalar record](alp-flavour-violating-two-pseudoscalar-amplitudes.md).

## Reproducibility and limits

[`scripts/derive_alp_finite_mixing_two_pseudoscalars.py`](../../../scripts/derive_alp_finite_mixing_two_pseudoscalars.py) differentiates a generic cubic polynomial to display all field-direction components, and exposes the existing CCWZ vertex evaluator for arbitrary supplied flavour directions and source derivatives. It derives the nondegenerate weak-mixing term under the stated counting and runs the established charged-redefinition derivation. The tests check daughter exchange symmetry, ALP components on daughters, and the $K^0\pi^0$ source-conjugation and LO mass relation. A direct charged CCWZ check retains both real and imaginary parts of the Hermitian $d$–$s$ vector entry and reproduces the charged result with an arbitrary incoming ALP overlap.

One restricted meson-plus-ALP example derives both the quadratic matrices and cubic vertex from the same CCWZ sources. Its parent and daughter have nonzero admixtures. Its exact finite-pole amplitude approaches the weak-mixing expression with a residual proportional to the cube of the inverse ALP scale in the tested range. This illustrative example does not establish decoupling of omitted flavour directions or real physical phase space. It uses a nonlinear-source-capable setup rather than deriving all source-consistent weak limits for the full neutral block. No loop, resonance, width, or weak neutral-kaon calculation is included.

```sh
.venv/bin/python scripts/derive_alp_finite_mixing_two_pseudoscalars.py
.venv/bin/python -m unittest discover -s scripts -p test_derive_alp_finite_mixing_two_pseudoscalars.py
```
