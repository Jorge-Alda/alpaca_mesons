---
title: Joint power-counted neutral cubic resonance kernels and source EOM audit
status: derived
created: 2026-10-14
updated: 2026-10-15
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx, Ecker:1988te, Ecker:1989yg]
provenance: analytical-derivation
---

Canonical tensor and odd-coupling conventions follow [the Levi-Civita audit](levi-civita-convention.md). The linked original scripts and saved artifacts use the explicitly identified reference backend; apply that record's signed coefficient and tensor adapter when interpreting them in the canonical convention.

# Inputs and scope

Use [the 228-dimensional bounded action](odd-higher-covariant-operator-count.md),
[the parent filtration](odd-hilbert-power-counting.md),
[the derivative prescription](odd-resonance-derivative-counting.md),
[the higher stripped projection](odd-higher-two-boson-combinations.md),
[the source frame and EOM](odd-covariant-source-gauge-reduction.md),
[the normalized leading scalar kernels](odd-trace-scalar-kernels.md), and
[the stage-4 photon response](vector-stage4-amplitudes.md).
The on-shell source-family rank specialization also uses the checksum-pinned
[Alda:2026nrd central quark inputs](alda-neutral-amplitude-inputs.md).

The executable
[derive_odd_power_counted_channels.py](../../../scripts/derive_odd_power_counted_channels.py)
maps all 44 leading plus 197 higher **raw parents**, before selecting a
quotient basis or imposing the power-counting filter. It retains Alda:2026nrd's
field/source conventions, canonical wedge tensors, $(+---)$,
$\epsilon^{0123}=+1, \epsilon_{0123}=-1$, all-incoming momenta and $X=2\chi=4B_0m_q$.
Quark masses remain symbolic and the neutral rho/omega/positive-strange-phi
generator span is unchanged. Charged source matrices are not covered.

These are normalized local cubic resonance-sector gradients. They are not
complete physical ALP predictors: anomaly/local matching, direct-plus-meson-
pole LSZ equivalence and chosen ALP UV profiles still require assembly.
Ranks across an arbitrary source family are not measurement Jacobian ranks.

# Fourier/source normalization

For one scalar field $a$ define coefficient matrices by

$$
u_\mu|_a=i p_\mu U a,\qquad (i\chi_-)|_a=Y a,
\qquad Z_{TC}|_a=Z_a a.
$$

Here $U$ is a coefficient matrix, not the nonlinear Goldstone field.
It has mass dimension minus one; $Y$ and $Z_a$ have dimension one.
The script uses arbitrary diagonal coefficients until the appropriate
specialization. For a canonical SM scalar pseudoscalar and for the physical
ALP source pair, respectively,

$$
U_{\rm SM}=-\frac{\sqrt2}{F}T_P,
\qquad Y_{\rm SM}=\frac{\sqrt2}{F}\{T_P,\chi\},
$$

$$
U_a=-\frac{2}{f_a}\mathbb M_a^{-1}(\mathcal C),
\qquad Y_a=\frac{2}{f_a}
\{c_q^{AP}+\mathbb M_a^{-1}(\mathcal C),\chi\}.
$$

No identification of these two matrices is made. The explicit photon block
is $f_+=-2eQF_\gamma$; an induced vector leg is
$V=\sqrt2eW(q^2)F_\gamma$, with the signed complete spectral response
of the stage-4 action. An external vector leg uses its own $1/M_i$.

The script derives the epsilon contractions with all derivative placements.
Two tensor factors each contribute a Fourier $i$; $u$ and an additional
derivative each contribute another $i$. Thus derivative, mass and
two-derivative mass vertices have phases plus, minus and plus, respectively.
The exported higher VVP action factor two and both Wick assignments are
included. Leading raw scalar columns have their quotient convention's
normalized Bose average: published $c_3$ and $d_2$ each embed with a factor
two. This is a column-normalization dictionary, not a change of physical
$c_i,d_i$ definitions.

The original VJP derivative coordinate embedding and the VVP symmetric
embedding are contracted anew. Exact symbolic comparison recovers the
established $C$ and $H$ kernels, independently of EOM constraints on the
mass traces. The original neutral $c_4$ image vanishes by the explicit
charge commutator, as do all neutral $f_-$ parents here.

# Why the EOM contact matters

Set $w=m_a^2=s+t+2h$. The LO Goldstone equation gives, after the Fourier
specialization,

$$
Y=-2wU+2Z_aI.
$$

Before making this specialization the derivative and mass sources remain
independent. Three of the thirteen action relations then have nonzero cubic
images: the two leading VJP and one leading VVP Goldstone-EOM relations.
Consequently this independent-source cubic projection does **not** descend
to the EOM-reduced action quotient. Its rank cannot be called the physical
rank of that quotient. The script exports every residual.

After retaining $Z_a$ and using the LO EOM-compatible specialization, all
thirteen relation images vanish exactly. The script checks
$M R^T=0$ with the complete parent relation matrix. This is a conditional
one-scalar source specialization, not a claim that the EOM may be imposed
on arbitrary independent external sources before differentiation.

For any supplied ALP UV profile, first evaluate the exported matrix defect

$$
\Delta_{\rm EOM}=Y_a+2m_a^2U_a-2Z_aI.
$$

Only if it vanishes within the stated LO source conventions may the
compatible maps be used. For SM mesons, the same issue distinguishes LO
eigenstate/mass inputs from an unmatched substitution of physical masses.
No topological contact coefficient is zeroed to force compatibility.

# Executed joint primitive result

The map distinguishes direct VJP and VVP kernel slots, with independent
diagonal boson generators and constant scalar vacuum. Exact SymPy row
reduction gives:

| Source prescription | Leading | Through soft first corrections | Full |
|---|---:|---:|---:|
| Independent derivative/mass sources | 7 | 32 | 73 |
| LO EOM-compatible sources with contact | 7 | 29 | 70 |

With hard-vector derivative parents promoted, the corresponding nested
ranks become $(12,43,73)$ and $(12,40,70)$. These are joint formal primitive
kernel dimensions, distinct from the action dimensions $(12,101,228)$.
They use all raw parents selected at each ceiling, not trace labels of the
previously reduced 51 combinations.

# Channel limits and validity

The companion JSON contains each normalized parent gradient for
$a\to VV$, $a\to V\gamma$ and $a\to\gamma\gamma$, with both VJP photon
assignments. It also computes on-shell neutral channel maps after imposing
$s=M_1^2,t=M_2^2$, $s=M_1^2,t=0$, or $s=t=0$, respectively.
Their coefficient-field ranks hold the Alda:2026nrd central quark spurion fixed at
$\chi=(0.013136,0.023480,0.471842)$ GeV squared, represented as exact
decimals. Vector masses and the diagonal $W(0)$ response remain coefficient-
field parameters. They vary the compatible source family
and labelled neutral external generators; they do not claim ranks for a
fixed ALP model, selected measurements or every numerical mass/sign branch.

| Neutral source family | Leading soft | Through soft first corrections | Through promoted first corrections | Full |
|---|---:|---:|---:|---:|
| $a\to VV$ | 3 | 12 | 15 | 25 |
| $a\to V\gamma$ | 6 | 19 | 23 | 27 |
| $a\to\gamma\gamma$ | 3 | 6 | 7 | 7 |

The promoted leading dimensions are four, eight and five, respectively.
For real photons the promoted counting is needed only if their scalar
invariant $h=m_a^2/2$ is hard; their zero virtualities alone do not require it.
The table does not pool target ranks: the same parent coefficients can enter
more than one channel family.

For rank arithmetic only, invertible source coordinates remove overall
mass and response scales: $m_a^2=M_1^2x$, $Z_a=M_1^2z$,
$M_2^2=M_1^2r_M$ and $W(0)=W_u\operatorname{diag}(1,r_d,r_s)$.
The script keeps the original dimensional gradients separately. The
[coordinate-verification executable](../../../scripts/verify_odd_channel_rank_coordinates.py)
derives and checks 539 exact column-scale identities. It proves that these
coordinate changes preserve every retained-parent rank on the generic
domain $M_1\ne0$ and, in photon families, $W_u\ne0$. Special parameter
loci, identical-vector restrictions, universal source response and a chosen
UV profile require their own specialization; generic ranks cannot simply
be carried over to those cases.

Photon Ward tensors vanish before projection and VVP images satisfy exact
assignment symmetry. For a single odd insertion in vacuum $a\to3P$, the
script derives the zero determinant
$\epsilon(-p_1-p_2-p_3,p_1,p_2,p_3)=0$; additional contracted derivatives
cannot add a fourth independent external scalar momentum. This result does
not apply with additional external sources or higher-order odd insertions.

# Reproduction and artifacts

```text
.venv/bin/python scripts/derive_odd_power_counted_channels.py
.venv/bin/python scripts/verify_odd_channel_rank_coordinates.py
```

- [Executed derivation](../../../data/vector_sector/odd_power_counted_channels.md) — source matrices, hard/residual split, ranks and EOM residuals.
- [Machine-readable maps](../../../data/vector_sector/odd_power_counted_channels.json) — raw parent gradients, filtered combinations, normalizations and input hashes.
- [Coordinate-scale verification](../../../data/vector_sector/odd_channel_rank_coordinate_verification.json) — 539 exact identities and the nonzero-scale domain, with a [short transcript](../../../data/vector_sector/odd_channel_rank_coordinate_verification.md).

Exact reconstruction is checked in the coefficient domain, including the
$\sqrt2$ extension and symbolic fixed-input parameters. In the latter case
an exact rational evaluation selects candidate independent rows; a symbolic
polynomial reconstruction of every column proves the matching upper bound.
That evaluation alone is not treated as rank saturation. The implementation
does not infer numerical correction sizes, coefficient values or convergence
from these ranks. Next assemble a chosen LO-consistent ALP profile and
anomaly/LSZ sector, derive the missing charged SM higher kernels, then build
fixed-scenario low-energy Jacobians before calibration or uncertainty bands.

The deterministic calculation cache under
`/private/tmp/alpaca-power-counted-channel-cache/` is keyed by the script and
full map inputs. It can be regenerated; only completed exact reconstructions
are saved. Full-map reductions are reused at identical truncation ceilings.
The established scalar-kernel regression runner also passes all 32 tests.
