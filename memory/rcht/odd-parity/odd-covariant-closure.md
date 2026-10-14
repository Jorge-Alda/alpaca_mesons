---
title: Covariant closure and explicit-vector counting convention
status: derived
created: 2026-10-14
updated: 2026-10-14
bibtex_keys: [Alda:2026nrd, Ecker:1988te, Ecker:1989yg]
provenance: analytical-derivation
---

# Counting convention

User decision on 2026-10-14: preserve explicit vector fields and the fixed
derivative ceiling; reduce chiral identities and source-preserving Goldstone
EOM. Resonance EOM are not an equivalence relation in this count. Independent
topological coefficient functions remain excluded, but fixed contacts induced
by the Goldstone EOM remain attached to their parent coefficients.

Inputs are the [bounded catalogue](odd-breaking-vertex-catalogue.md),
[source EOM audit](odd-trace-source-eom-audit.md), and the building blocks and
tensor kinetic action in [RChT basics](../foundations/rcht-basics.md).
Source and field conventions remain those of Alda:2026nrd.

# Executed identities

Run `.venv/bin/python scripts/derive_odd_covariant_closure.py` from the repository
root. The deterministic SymPy script derives the identities with noncommuting
three-by-three matrix symbols, verifies them exactly, and saves the
[transcript](../../../data/vector_sector/odd_covariant_closure.md) and
[hashed results](../../../data/vector_sector/odd_covariant_closure.json).
These are algebraic derivations from the cited definitions, not new literature
claims. Write $X=\chi_+$ and $Y=i\chi_-$.

The dressed connections $\Omega_R=\Gamma-iu/2$ and
$\Omega_L=\Gamma+iu/2$ have curvatures $-i u^\dagger F_Ru$ and
$-i uF_Lu^\dagger$, respectively. Their sum and difference give

$$
\Gamma_{\mu\nu}=\frac14[u_\mu,u_\nu]-\frac{i}{2}f_{+\mu\nu},
\qquad \nabla_\mu u_\nu-\nabla_\nu u_\mu=-f_{-\mu\nu}.
$$

The corresponding dressed Bianchi identities give

$$
\sum_{\mathrm{cyc}}\nabla_\mu f_{\pm\nu\rho}
=-\frac{i}{2}\sum_{\mathrm{cyc}}[u_\mu,f_{\mp\nu\rho}],
\qquad
[\nabla_\mu,\nabla_\nu]K=[\Gamma_{\mu\nu},K].
$$

In particular, a linearized homogeneous Bianchi relation carries a covariant
$u f_-$ remainder. Commuting derivatives carries both $f_+K$ and $[u,u]K$
remainders. A vanishing one-pseudoscalar vertex projection does not make these
remainders zero as operators.

Define the independently differentiated source blocks

$$
\chi_\mu^\pm=u^\dagger(D_\mu\chi)u^\dagger
\pm u(D_\mu\chi^\dagger)u.
$$

Differentiating the dressed spurions then yields

$$
\nabla_\mu X=\chi_\mu^++\frac12\{u_\mu,Y\},
\qquad
\nabla_\mu Y=i\chi_\mu^--\frac12\{u_\mu,X\},
$$

$$
\nabla_\mu\nabla_\nu Y
=i\nabla_\mu\chi_\nu^-
-\frac12\{\nabla_\mu u_\nu,X\}
-\frac12\{u_\nu,\chi_\mu^+\}
-\frac14\{u_\nu,\{u_\mu,Y\}\}.
$$

Thus replacing the scalar spurion by a constant vacuum matrix before source
differentiation loses source gradients and a three-odd-block interaction.
Those terms must enter the remainder columns of a covariant reduction matrix.

# Why vector EOM are excluded

The script derives the six-component antisymmetric-tensor quadratic Euler
operator in metric $(+---)$, with Fourier derivatives $\partial_\mu\to iq_\mu$.
Its canonical wedge mode has eigenvalue $M_V^2-q^2$. Writing the action as
$V^T H(M_V^2+T)V/2+V^T HJ+\epsilon V^THO$, a field shift
$V\mapsto V-\epsilon O/M_V^2$ removes the last term but generates

$$
-\frac{\epsilon}{M_V^2}V^THTO
-\frac{\epsilon}{M_V^2}J^THO.
$$

Here $H$ is the tensor-component metric, $T$ contains two derivatives, $J$ is
the independent even source, and $O$ is the odd tensor current. The induced
terms are a vector interaction with two extra derivatives and a source-only
contact. This exact quadratic illustration demonstrates nonclosure of the
fixed-ceiling space under resonance-EOM elimination; it does not enumerate
the full interacting resonance EOM.

# Count status and remaining calculation

The saved **205** is a source-resolved vertex-tensor count. The **15** is the
additional formal scalar electromagnetic amplitude rank. Neither is the
number of new independent covariant operators. The present script establishes
the identity dictionary and selected quotient convention; it does not yet
construct the full reduction matrix. Consequently no subtraction from 205
is justified here.

The remaining calculation is to lift the catalogue to explicit covariant
representatives, append their Bianchi/curvature/spurion and Goldstone-EOM
contact columns, and compute the exact quotient rank including the embedded
original eleven operators. Induced remainders share the parent coefficient;
they are not automatically new freely adjustable couplings. Report a final
new-operator count only after that augmented rank calculation.
