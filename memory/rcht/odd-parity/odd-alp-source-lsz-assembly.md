---
title: LO diagonal ALP source profile and parentwise LSZ assembly
status: mixed
created: 2026-10-14
updated: 2026-10-15
bibtex_keys: [Alda:2026nrd, Bai:2025fvl, Bauer:2020jbp]
provenance: analytical-derivation
---

# Inputs and assumptions

Use Alda:2026nrd's [source expansion and no-mixing prescription](../../extracts/Alda:2026nrd/appendix-a-expansion.md),
its [LO topological action](../../extracts/Alda:2026nrd/lo-topological-action.md),
the [normalized neutral kernels](odd-power-counted-channel-maps.md), and the
explicit assumptions of the [photon invariant record](../../alp/photon/alp-photon-invariant-from-mixing.md).
Neutral UV coefficients and the quark spurion are diagonal here. Retain one
ALP at tree level and order $1/f_a$. Set $w=m_a^2$, $g=m_0^2/3$,
$c^{AP}=c^A-c^P$ and $c_{GA}=c_G+\operatorname{tr}c^A$.

# Derivation

Run

```text
.venv/bin/python scripts/derive_odd_alp_lsz_assembly.py
```

The executable solves the LO mass-operator equations rather than introducing
independent source weights:

$$
(\chi_i-w)S_i+gT_S=-\chi_i c_i^{AP}-gc_{GA},\qquad
T_S=\sum_i S_i,\qquad S=\mathbb M_a^{-1}(C).
$$

Solving the individual equations and then the trace equation gives

$$
T_S=-\frac{\sum_i\chi_i c_i^{AP}/(\chi_i-w)
+gc_{GA}\sum_i1/(\chi_i-w)}{1+g\sum_i1/(\chi_i-w)},
\qquad
S_i=\frac{-\chi_i c_i^{AP}-g(c_{GA}+T_S)}{\chi_i-w}.
$$

The apparent individual denominators are removable when the full neutral
resolvent determinant is nonzero. The executable verifies this by multiplying
each reduced component by $\det[\operatorname{diag}\chi+g\mathbf1\mathbf1^T-wI]$.
The actual domain excludes neutral mass poles and requires perturbative
ALP–meson overlaps; finite inversion alone does not establish LSZ validity.
The LO masses must come from this same operator.

With Alda:2026nrd's field conversion, the meson replacement is
$\Phi_R/a=\sqrt2F(c^A+S)/f_a$. Direct insertion and this replacement give

$$
U_{\rm dir}=2c^A/f_a,\quad U_{\rm mix}=-2(c^A+S)/f_a,
\quad U=-2S/f_a,
$$

$$
Y_{\rm dir}=-2\{c^P,\chi\}/f_a,\quad
Y_{\rm mix}=2\{c^A+S,\chi\}/f_a,\quad
Y=2\{c^{AP}+S,\chi\}/f_a.
$$

The topological contact is

$$
Z=Z_{\rm dir}+Z_{\rm mix}
=-\frac{2g}{f_a}[c_{GA}+\operatorname{tr}S],
\qquad Y+2wU-2ZI=0.
$$

Every retained cubic resonance kernel is linear in $U,Y$. The script verifies
direct-plus-meson assembly for all 241 raw parents in four neutral external
projections (964 identities), as well as the exact source EOM. Applying the
no-mixing profile already includes the meson LSZ contribution; adding it again
would double count.

# Photon contact and remaining matching

Using the **Bauer photon-operator normalization**, together with the
separately declared quark-rotation and leading-nonet WZW inputs,

$$
C_{a\gamma\gamma}=c_\gamma
+2N_c\operatorname{tr}[Q^2(c^A+S)]
$$

is invariant under the declared coefficient rotation and equals the transformed
contact in the no-mixing basis. These algebraic identities pass. The positive
local photon coefficient convention is supplied by [@Bauer:2020jbp, Eqs. (39)/(54)],
as recorded in the [source digest](../../extracts/Bauer:2020jbp/photon-operator-normalization.md).
Alda:2026nrd excludes the WZW sector; the meson contribution is a separate repository
derivation, not an equation supplied by that paper. The
[shared photon/action phase audit](../../alp/photon/alp-photon-odd-action-phase-matching.md)
uses the single canonical tensor $\epsilon^{0123}=+1$, with
$\epsilon_{0123}=-1$ for $(+---)$. In the common cubic tensor $E$ the scalar
contact is $-C_{a\gamma\gamma}/(4\pi^2f_a)$, to be added to the canonical
resonance gradients. The [convention audit](levi-civita-convention.md)
records the action-preserving conversion of saved reference-backend outputs.
The SM [TFF convention](radiative-tff-amplitude.md)
uses $R=-E$ and the opposite external-meson phase; its LSZ overlap must
follow that phase as well. The convention conversion does not establish
independent UV/local anomaly matching for the enlarged RChT action. No
numerical ALP photon predictor, independent hadronic anomaly coefficient,
or measurement rank is inferred here.

# Outputs

[Exact profile artifact](../../../data/vector_sector/odd_alp_lsz_assembly.json)
and [step-by-step transcript](../../../data/vector_sector/odd_alp_lsz_assembly.md)
include input hashes. Neutral resonance source assembly and charged cubic
kernels and signed photon-phase matching are now available. The next concrete
gate is fixed-scenario low-energy measurement Jacobians.
