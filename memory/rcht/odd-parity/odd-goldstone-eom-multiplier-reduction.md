---
title: Goldstone EOM multiplier images and the one-spurion ceiling
status: derived
created: 2026-10-14
updated: 2026-10-14
bibtex_keys: [Alda:2026nrd, Ecker:1988te]
provenance: analytical-derivation
---

# Inputs and scope

Use [Alda's LO topological action](../../extracts/Alda:2026nrd/lo-topological-action.md),
[the neutral source derivation](odd-trace-source-eom-audit.md),
[covariant closure conventions](odd-covariant-closure.md), and
[the flavor enumeration](odd-breaking-flavor-enumeration.md).
The reduction uses the LO Goldstone action, independent sources, three flavors,
and natural units. It does not include variations of even resonance interactions
in the Goldstone EOM. Vector fields are preserved as the user requested.

# Executable calculation

Run [the script](../../../scripts/reduce_odd_goldstone_eom.py):

```text
.venv/bin/python scripts/reduce_odd_goldstone_eom.py
```

The [transcript](../../../data/vector_sector/odd_goldstone_eom_reduction.md)
and [hashed matrices](../../../data/vector_sector/odd_goldstone_eom_reduction.json)
are deterministic exact SymPy outputs. Flavor matrices are abstract; only
finite cyclic trace labels reuse the earlier catalogue canonicalization.

For a covariant axial tangent $\delta U=i u\eta u$, use the induced variations
$\delta u_\mu=-\nabla_\mu\eta$ (up to compensator commutators),
$\delta\chi_+=-\{\eta,Y\}/2$ and $\delta X_{TC}=-\operatorname{tr}\eta$,
where $Y=i\chi_-$. Compensator commutators vanish in the action traces.
Integrating the kinetic variation by parts gives

$$
\delta S_2=\int d^4x\left[
\frac{F^2}{2}\operatorname{tr}(\eta\nabla\cdot u)
-\frac{F^2}{4}\operatorname{tr}(\eta Y)
+\frac{F^2m_0^2X_{TC}}{6}\operatorname{tr}\eta\right].
$$

The script extracts the tangent coefficients and divides by the kinetic
coefficient, obtaining

$$
\nabla\cdot u=\frac{Y}{2}-\frac{m_0^2X_{TC}}{3}I.
$$

This extends the previous neutral-linear identity to the specified covariant
LO Goldstone action. External sources have not been set to zero.

# Exact multiplier matrices

Let $F_j(A,B,Z)$ be the catalogue's charge-even, at-most-double-trace
flavor functions. Replace the odd slot $Z$ by the divergence above and keep
the common epsilon/tensor factor. The two images are $F_j(A,B,Y)/2$ and
$-m_0^2X_{TC}F_j(A,B,I)/3$. Repeat with one extra $X=\chi_+$ slot.
For VVP keep the Bose-symmetric flavor sector appropriate to this wedge
contraction. This calculation concerns these divergence multipliers, not all
Lorentz structures of the catalogue.

| Sector | Multipliers | Mass-image rank | Topological-image rank | Contact-free dimension |
| --- | ---: | ---: | ---: | ---: |
| Leading VJP | 4 | 4 | 2 | 2 |
| Leading VVP | 3 | 3 | 2 | 1 |
| One-$X$ VJP | 10 | 10 | 4 | 6 |
| One-$X$ VVP | 7 | 7 | 3 | 4 |

All ranks follow from exact row reduction. The augmented mass/contact map is
injective in each sector. The three leading contact-free combinations recover
the earlier neutral-linear relations, now as flavor identities multiplying
the covariant LO EOM. Ordering and normalization of the explicit null vectors
are pinned in the JSON rather than transferred between catalogues silently.

# A further closure constraint

For a one-$X$ divergence operator, the mass image contains both $X=\chi_+$
and $Y=i\chi_-$: it has **two quark-spurion insertions**. The saved action
space permits at most one. The mass-image ranks 10 and 7 are full, so no
nonzero combination of the respective multipliers cancels that entire
two-spurion image. The six/four contact-free combinations cancel only the
topological image; they still generate independent two-spurion terms.

Consequently these seventeen replacements are not reductions wholly inside
the saved one-spurion space. One may retain their two-spurion images with
fixed parent coefficients in an enlarged representative space, or define a
separate spurion-order quotient that discards them. Such a discard has not
been adopted here. Under the current exact retained-remainder convention,
they cannot be counted as seventeen redundant one-spurion coefficients.

This result does not require another user decision to continue: preserve the
remainders in the augmented matrix, just as for Bianchi and source contacts.
It also demonstrates why simply subtracting all contact-free dimensions
from 205 would be incorrect.

# Remaining count calculation

The global count still requires the covariant lift of the Lorentz/IBP basis
and its remainder columns, including the operators with two derivatives on
$Y$ or $f_-$. The original eleven must be embedded into that same quotient.
Neither seven leading multipliers nor seventeen one-$X$ multipliers can be
subtracted globally without this embedding. No final new-operator count or
low-energy fit is asserted by this subsector calculation.
