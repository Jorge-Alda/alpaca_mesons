---
title: Flavor selection rules for an ALP to a neutral vector and two pseudoscalars
status: derived
created: 2026-10-21
updated: 2026-10-21
bibtex_keys: [Alda:2026nrd, Ecker:1988te, Ecker:1989yg, Ruiz-Femenia:2003jdx]
provenance: analytical-derivation
---

# Scope and inputs

This record derives flavor support for the strong tree amplitudes in the
[$a\to VPP$ amplitude record](alp-vpp-amplitude.md), with external
$V=\omega,\phi$. The external vector generator is diagonal and isospin
symmetric, $T_V=\operatorname{diag}(v_n,v_n,v_s)$; ideal $\omega$--$\phi$ mixing is
considered separately. Pseudoscalar generators use the conventions of the
[generic-generator record](alp-three-meson-generic-generators.md), with generic
$\eta$--$\eta'$ mixing. Charge is conserved, and only neutral ALP spurions are
considered. In the even sector the vector-scalar source is restricted to its
neutral $|\Delta S|=1$ entries $c_{ds}$ and $c_{sd}$; charged $us,su$ source
entries cannot contribute to a neutral ALP and neutral final state.

The executable trace derivation is
[`scripts/derive_alp_vpp_selection_rules.py`](../../../scripts/derive_alp_vpp_selection_rules.py):

```text
.venv/bin/python scripts/derive_alp_vpp_selection_rules.py
```

# Even scalar-source amplitude

The two flavor factors are
$$
\operatorname{tr}\!\left(\{\lambda_j,\lambda_R\}
[c_q^{VS},\chi]\right),\qquad
\operatorname{tr}\!\left(T_V[\lambda_R,\lambda_i]\right).
$$
Since $T_V$ and $\chi$ are diagonal, the second trace requires an off-diagonal
intermediate pseudoscalar conjugate in flavor flow to $P_i$. The first trace
then requires the total final pseudoscalar strangeness to match the neutral
off-diagonal source. Exact symbolic enumeration of all three exchange
routings gives, for either external vector,

| Source entry | Allowed $P_1P_2$ channels |
|---|---|
| $c_{ds}$ | $K^0\pi^0$, $K^0\eta$, $K^0\eta'$, $K^+\pi^-$ |
| $c_{sd}$ | $\bar K^0\pi^0$, $\bar K^0\eta$, $\bar K^0\eta'$, $K^-\pi^+$ |

The nonzero coefficients are proportional to
$(\chi_d-\chi_s)(v_n-v_s)$, up to channel and exchange factors. Thus the even
amplitude also vanishes in the $d$--$s$ degenerate-mass limit or for a vector
generator proportional to the identity. These are additional zeros beyond
the charge and strangeness rules.

# Odd amplitude with diagonal ALP sources

For diagonal $c_q^A$ and $c_q^P$, the ALP-to-vector-pair trace in the first
exchange routing selects diagonal internal vectors. The $VPP$ commutator trace
then requires the two pseudoscalars to carry conjugate flavor transitions.
The other two routings obey the same overall flavor-flow rule after summing
their allowed intermediate-vector components. The allowed pairs are
$$
\pi^+\pi^-,\qquad K^+K^-,\qquad K^0\bar K^0.
$$
Pairs made only from diagonal pseudoscalars, such as $\pi^0\eta$, vanish
because their commutator is zero. For ideal mixing,
$T_\phi=\operatorname{diag}(0,0,1)$ gives zero for $\phi\pi^+\pi^-$, while
$\omega\pi^+\pi^-$ and both neutral-vector $K\bar K$ channels are allowed.
With nonideal isoscalar mixing, $\phi\pi^+\pi^-$ is allowed if the physical
$T_\phi$ has a nonzero nonstrange component.

# Odd amplitude with neutral $|\Delta S|=1$ ALP sources

For off-diagonal neutral entries $c_{ds}$ and/or $c_{sd}$ in $c_q^A$ or
$c_q^P$, the anticommutator trace in $\mathcal N_R$ selects neutral
flavor-changing internal-vector directions, and the $VPP$ or $VPa$ traces
enforce the same net strangeness in all three exchange routings. The support is
the same for either external $\omega$ or $\phi$:

| ALP source entry | Allowed $P_1P_2$ channels |
|---|---|
| $c_{ds}$ | $K^0\pi^0$, $K^0\eta$, $K^0\eta'$, $K^+\pi^-$ |
| $c_{sd}$ | $\bar K^0\pi^0$, $\bar K^0\eta$, $\bar K^0\eta'$, $K^-\pi^+$ |

The amplitude is linear in the ALP source entries. The executable derivation
checks $c_q^A$-only and $c_q^P$-only cases separately; each gives the same
channel support and source orientation shown in the table. Imposing exact
strong isospin, $\chi_u=\chi_d$, does not remove any full channel after summing
the three routings. Individual trace terms can still vanish, and special
pseudoscalar mixing angles or tuned coefficient cancellations can remove a
channel from the full amplitude. If both Hermitian-conjugate entries are
present, both table rows are possible.

# Verification and scope

The script forms the trace factors using explicit $3\times3$ generators,
enumerates unordered pseudoscalar pairs, and sums the real components of
charged internal fields before identifying structural zeros. It verifies the
support for generic diagonal isoscalar vectors and separately substitutes
ideal $\omega$ and $\phi$ generators. It separately evaluates axial-only and
pseudoscalar-only off-diagonal ALP sources and the exact-isospin limit. These
rules apply to the tree-level
amplitudes in the parent record and do not include electromagnetic
conversions, loops, rescattering, or additional flavor-breaking operators.
