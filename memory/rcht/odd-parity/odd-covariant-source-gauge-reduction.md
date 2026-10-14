---
title: Independent-source frame and joint covariant trace-jet reduction
status: derived
created: 2026-10-14
updated: 2026-10-15
bibtex_keys: [Alda:2026nrd, Ecker:1988te, Ruiz-Femenia:2003jdx]
provenance: analytical-derivation
---

Canonical tensor and odd-coupling conventions follow [the Levi-Civita audit](levi-civita-convention.md). The linked original scripts and saved artifacts use the explicitly identified reference backend; apply that record's signed coefficient and tensor adapter when interpreting them in the canonical convention.

# Inputs and declared quotient

Use [the covariant identities](odd-covariant-closure.md),
[the source-preserving LO Goldstone EOM](odd-goldstone-eom-multiplier-reduction.md),
[the derivative ceiling](odd-breaking-lorentz-reduction.md),
[the leading quotient](odd-leading-covariant-operator-count.md), and
[Alda's topological action](../../extracts/Alda:2026nrd/lo-topological-action.md).
Use the specified VJP/VVP seed classes, at most one quark spurion and two
flavor traces. Preserve explicit vector fields. No vector EOM, on-shell
conditions, physical-mass substitution or loop calculation is used.
This is not an enumeration of arbitrary independent multi-Goldstone seeds.

# Local source frame

After deriving the covariant Goldstone EOM with independent sources, choose
the local chiral frame $U=I$. Write $v=(r+\ell)/2$ and
$a_{\rm ext}=(r-\ell)/2$. The script's primitive matrices `a0` through `a3`
are components of this external axial source, **not the ALP field**.
Source profiles, their derivatives and the explicit tensor vector are arbitrary.

The building blocks in this frame are

$$
\Gamma_\mu=-iv_\mu,\qquad u_\mu=2a_{{\rm ext},\mu},
\qquad \nabla_\mu K=\partial_\mu K-i[v_\mu,K],
$$

$$
f_{+\mu\nu}=2\left(\partial_\mu v_\nu-\partial_\nu v_\mu
-i[v_\mu,v_\nu]-i[a_{{\rm ext},\mu},a_{{\rm ext},\nu}]\right),
$$

$$
f_{-\mu\nu}=-2\left(\partial_\mu a_{{\rm ext},\nu}
-\partial_\nu a_{{\rm ext},\mu}
-i[v_\mu,a_{{\rm ext},\nu}]-i[a_{{\rm ext},\mu},v_\nu]\right).
$$

Keep $X=\chi_+$ as an arbitrary matrix-valued scalar source. With
$Y=i\chi_-$ and $Z_{TC}=m_0^2X_{TC}/3$, the LO EOM quotient is represented by

$$
Y=2\nabla\cdot u+2Z_{TC}I.
$$

This algebraically solves the source-dependent EOM for the independent
pseudoscalar source; it does not set that source to zero. $Z_{TC}$ and all
its derivatives are retained. In particular, contacts induced by this
replacement remain tied to the parent operator coefficient. Independent
new topological coefficient functions are not added.

The computation is in natural units with metric $(+---)$ and
$\epsilon^{0123}=+1, \epsilon_{0123}=-1$. Derivative symbols are formal derivative eigenvalues;
there are no implicit Fourier factors of $i$ in the polynomial IBP step.

# Ordered derivatives and actual matrix algebra

[The preparation script](../../../scripts/prepare_odd_covariant_blocks.py)
polarizes the two derivative slots independently. The generic Lorentz ranks
before chiral identities are nine for each one-$X$ derivative or ordered
two-derivative block. VVP parity subdivisions are $(4,5)$ for one-$X$
derivatives and $(5,4)$ for ordered two-derivative scalar blocks. These are
preparation ranks, not independent action counts.

[The joint script](../../../scripts/reduce_odd_covariant_higher.py) evaluates
the complete source-frame operators, including ordered covariant derivatives.
Their commutators and both nonlinear Bianchi identities are therefore retained
automatically. Matrix-valued Lorentz components remain abstract three-by-three
matrices. SymPy supplies the noncommutative polynomial algebra, differentiation,
trace expressions and exact linear algebra. A small adapter uses standard
SymPy noncommutative symbols for differentiation because general coordinate
and Frechet derivatives of `MatrixSymbol`, `MatPow` and `Trace` are unavailable.
It is not a replacement flavor-component calculation.

[The verification script](../../../scripts/verify_odd_source_gauge.py)
checks the product rule, commuting primitive partials, curvature, the $u$
curl, both full Bianchi identities and the EOM. It also proves equivalence
of the faster source-polynomial path to the abstract matrix definitions.

# Integration by parts and a common projection

Primitive partial jets commute. A trace monomial is represented by its
cyclic word, derivative labels and scalar-source factors. The script averages
the finite cyclic/permutation orbit, including stabilizers of repeated fields,
and eliminates one derivative momentum using total momentum conservation.
Exact `Poly` coefficients then define its source-functional normal form.
Products of traces and derivatives of $Z_{TC}$ are included. Explicit tests
verify that total derivatives vanish in this normal form.

The first calculation gives 34 candidate relations in the 197 higher
covariant seeds. These are candidates because that initial projection keeps
the lowest source grade of each seed separately. It is not yet a common
projection for the whole action.

In particular, a higher-derivative $f_-X$ relation can generate the very same
four-field source grade as a $Vf_+uX\nabla$ operator. The
[completion script](../../../scripts/complete_odd_higher_reduction.py)
matches all four such candidate remainders to the 80-dimensional lowest-grade
one-$X$ VJP span and subtracts the corresponding **full covariant operators**.
Two candidates then vanish exactly and two retain higher-field remainders.
Ignoring this subtraction would overcount the joint space.

After this common projection, the source projection has rank 163. Its
remaining higher-field map has an abstract trace-word upper rank of 25.
The [exact three-flavor Gram closure](odd-three-flavor-trace-gram.md) establishes
that all 25 directions survive at three flavors. The result is recorded in
the [complete count](odd-higher-covariant-operator-count.md).

# Reproduction

```text
.venv/bin/python scripts/prepare_odd_covariant_blocks.py
.venv/bin/python scripts/verify_odd_source_gauge.py
.venv/bin/python scripts/reduce_odd_covariant_higher.py
.venv/bin/python scripts/complete_odd_higher_reduction.py
.venv/bin/python scripts/certify_odd_higher_rank.py
.venv/bin/python scripts/close_odd_three_flavor_quotient.py
```

The deterministic cache under `/private/tmp/alpaca-odd-covariant-cache` is
keyed by the calculation script hash. It permits resuming interrupted work;
the scripts can regenerate it from scratch. Results and exact relation vectors
are saved under `data/vector_sector/`, including
[the common projection](../../../data/vector_sector/odd_covariant_higher_reduction.json).

No source gradients, curvature contacts, higher-field remainders or two-spurion
EOM images are discarded to enforce a desired count. The higher-derivative
Goldstone images are retained computationally with their fixed coefficients;
they do not become freely adjustable additional seeds.
