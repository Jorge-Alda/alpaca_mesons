---
title: Physical neutral TFF matching obstruction from source breaking
status: mixed
created: 2026-10-13
updated: 2026-10-15
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx, Roig:2013baa, Knecht:2001qf]
provenance: analytical-derivation
---

# Physical-channel high-energy conditions

## Inputs and assumptions

Use [stage-4 action amplitudes](vector-stage4-amplitudes.md),
[physical Alda generators](alda-neutral-amplitude-inputs.md),
[stage-3 source/mass scenarios](../foundations/vector-stage3-calibration.md),
and the [stage-5 likelihood specification](vector-stage5-joint-fit-specification.md).
Tree level, fixed central physical kinematics, independent quark spurion,
single-trace odd vertices and the declared partial-isospin vector sector apply.
Self-energies remain deferred by user instruction.
The [convention dictionary](levi-civita-convention.md) fixes $(+---)$ and
$\epsilon^{0123}=+1$, hence $\epsilon_{0123}=-1$.
Below, $c_i=-c_i^R$, $d_i=-d_i^R$ and the fixed-phase scalar TFF is
the negative of its legacy value at corresponding couplings. The linked
original script/transcript/results retain source-convention coefficients;
their scalar results require this adapter. The conversion preserves
the source obstruction, residual norms and scenario counts.

## Definitions for the source mismatch

All flavor matrices are $3\times3$. The photon charge matrix is
$Q=\operatorname{diag}(2/3,-1/3,-1/3)$; its symbol is distinct from the
positive Euclidean virtuality $Q^2$. The vector index $i$ runs over
$\rho^0,\omega,\phi$ in the physical neutral basis of the selected mass action.

| Quantity | Definition, interpretation and units |
|---|---|
| $T_i$ | Dimensionless Hermitian diagonal generator multiplying the canonical physical tensor field: $V_{\mu\nu}=\sum_iT_iV_{i,\mu\nu}$ in this neutral subspace. Normalization is $\operatorname{tr}(T_iT_j)=\delta_{ij}$. The omega/phi generators are rotated consistently with the mass action, retaining the positive-strange phi phase. |
| $M_i$ | Positive tree-level physical vector pole mass in GeV for that scenario. It is distinct from the common operator-normalization scale $M_V$. |
| $R_i$ | Signed electromagnetic source overlap in GeV, $R_i=F_V\operatorname{tr}(T_iQ)+z_V\operatorname{tr}(T_i\{X,Q\})+y_V\operatorname{tr}(T_i)\operatorname{tr}(XQ)$. Electronic widths constrain its squared magnitude, not its relative sign. It is not a dimensionless photon coupling or an external-state wavefunction residue. |
| $W(s)$ | Flavor-matrix response $\sum_iR_iT_i/(M_i^2-s)$, in GeV$^{-1}$, for Minkowski photon virtuality $s$. In the conventions of [stage 4](vector-stage4-amplitudes.md), the induced tensor matrix is $V_{\mu\nu}=\sqrt2eW(s)F_{\gamma,\mu\nu}$ at leading electromagnetic order. This is an internal-line response, not the external-vector factor $1/M_i$. |
| $W_0$ | $W(0)$, the response to a real photon, in GeV$^{-1}$. |
| $S$ | First spectral source moment $\sum_iR_iT_i$, in GeV. It is the coefficient of $1/Q^2$ in $W(-Q^2)$, not a new fitted parameter independent of the $R_i$. |
| $J$ | Next moment $\sum_iR_iM_i^2T_i$, in GeV$^3$; $W(-Q^2)=S/Q^2-J/Q^4+\cdots$. |
| $q(U,V)$ | Channel-dependent symmetric bilinear trace $\operatorname{tr}(\{U,V\}T_P)$, with $\{U,V\}=UV+VU$ and $T_P$ the canonical generator of the chosen physical pseudoscalar $P$. It is neither a momentum nor a universal number. Its units are the product of the units of $U$ and $V$; $q(Q,Q)$ is dimensionless, while $q(S,Q)$ is in GeV. |
| $T_P$ | Dimensionless canonical neutral-pseudoscalar generator reconstructed from Alda arXiv v1: $T_P=F\lambda_P/(\sqrt2f_\pi)$. It retains the physical neutral admixtures and is independent of the vector mixing angle. |
| $h(U,V)$ | Companion mass trace $\operatorname{tr}(\{U,V\}\{T_P,\chi\})$, where $\chi=2B_0m_q$ has units GeV$^2$. It has two additional powers of mass compared with $q(U,V)$. |
| $X$ | Vector source's vacuum spurion $\chi_+|_{\Phi=0}=2\chi$, isospin-projected as $\operatorname{diag}(x_l,x_l,x_s)$. This $X$ is distinct from Alda's topological invariant $X_{TC}$. |
| $F_V,z_V,y_V$ | Coefficients of the common, single-trace spurion, and double-trace spurion electromagnetic source interactions, respectively. $F_V$ is in GeV; $z_V,y_V$ are in GeV$^{-1}$. |

These definitions follow the [source-action derivation](../foundations/vector-leptonic-width.md)
and [canonical amplitude insertion](vector-stage4-amplitudes.md), not an
independent effective $F_V$ assigned separately to each decay.

Only necessary vanishing conditions are imposed on the three physical
neutral-channel high-energy amplitudes. No numerical eta/eta-prime leading
coefficient is supplied: the present input record does not give the required
source-normalized axial decay-constant matrix or singlet matching. No full
finite-spurion correlator OPE is inferred from these TFF conditions.

## Executable derivation

Run [the matching script](../../../scripts/match_vector_short_distance.py).
Its [transcript](../../../data/vector_sector/short_distance_derivation.md)
and [results](../../../data/vector_sector/short_distance_results.json)
record the calculation, equations, all 64 separate scenarios and hashes.

Define the spectral moments

$$
S=\sum_i R_iT_i,\qquad J=\sum_iR_iM_i^2T_i,\qquad W_0=W(0).
$$

Then $W(-\alpha Q^2)=S/(\alpha Q^2)-J/(\alpha^2Q^4)+O(Q^{-6})$.
Insert this expansion into both action-derived VJP placements and the
bilinear VVP kernel. SymPy expands and extracts the Laurent coefficients.
It also verifies the single/double coefficients against the full numerical
kernel at $Q^2=10^5,2\times10^5$ GeV$^2$, including double-virtual ratios
$1$ and $2$. These artificial hard points are algebra checks, not data.

For compact notation define

$$
A=-c_1+c_2+c_5-2c_6,\qquad B=c_1-c_2+c_5,
$$

$$
L(U,V)=(c_1+c_2-c_5)m_P^2q(U,V)+4c_3h(U,V),
$$

$$
H_0(U,V)=(d_1-d_3)m_P^2q(U,V)+4d_2h(U,V).
$$

Adapting the script gives the canonical single-virtual coefficients

$$
a_P=-\frac{4B}{FM_V}q(W_0,Q),
$$

$$
b_P=-\frac{\sqrt2N_Cq(Q,Q)}{8\pi^2F}
-\frac{4A}{FM_V}q(S,Q)+\frac4{FM_V}L(W_0,Q)
+\frac{4\sqrt2d_3}{F}q(S,W_0),
$$

$$
k_P=\frac4{FM_V}[Aq(J,Q)+L(S,Q)]
-\frac{4\sqrt2}{F}[H_0(S,W_0)+d_3q(J,W_0)].
$$

For $s_1=-\alpha Q^2$, $s_2=-\beta Q^2$, with fixed positive
$\alpha,\beta$, the constant is

$$
b_P^{(2)}=-\frac{\sqrt2N_Cq(Q,Q)}{8\pi^2F}
-\frac4{FM_V}q(S,Q)
\left[2A+B\left(\frac\alpha\beta+\frac\beta\alpha\right)\right].
$$

The inverse coefficient is retained in the executed transcript. The script
checks that after $B=0$ the double constant is independent of the hard-ray
ratio; its symmetric cancellation tests therefore cover all fixed positive
ratios at this order.

## Source obstruction

In the retained scenarios, $q(W_0,Q)$ is nonzero for at least one physical
channel, so removing all single-virtual growth requires $B=0$.
The three double-virtual constant conditions then imply

$$
q(S,Q)=\Lambda q(Q,Q),\qquad
A\Lambda=-\frac{\sqrt2N_CM_V}{64\pi^2},
$$

with one common $\Lambda$ for all three pseudoscalars. The Alda generators
span the diagonal flavor space. Thus their three trace equations imply
$\{S,Q\}=\Lambda\{Q,Q\}$; since $S$ and $Q$ are diagonal and every
charge entry is nonzero, $S=\Lambda Q$.

The declared source action gives

$$
S=F_VQ+z_V\{X,Q\}+y_V\operatorname{tr}(XQ)I,
\qquad X=\operatorname{diag}(x_l,x_l,x_s).
$$

SymPy solves its three component equations and, for $x_l\ne x_s$, obtains
$z_V=y_V=0$, $\Lambda=F_V$. In the symmetric spurion limit the $z_V$
term is a common rescaling and the $y_V$ term vanishes. The obstruction is
therefore flavor breaking in the ultraviolet source moment without
corresponding odd interaction breaking; changing pole masses or omega–phi
basis alone cannot remove it.

## Fixed-central scenario test

Solve the nine affine zero-coefficient equations (single growing/constant
and double constant for each meson) by least squares, after normalizing
each equation by the larger of its row norm and target magnitude. This
is an algebra consistency diagnostic, **not chi-squared**, and has no
experimental error weights. The reported minimum-norm solutions are algebraic
representatives of underdetermined equations, not physical coupling fits.

| Source hypothesis | Consistent scenarios out of 16 | Relative failure of $S\propto Q$ |
|---|---:|---:|
| Universal | 16 | below $4\times10^{-16}$ |
| Single-trace positive | 0 | 0.0170–0.0392 |
| Single-trace negative | 0 | 0.3109–0.4564 |
| Single plus double trace | 0 | 0.0445–0.0564 |

All equation matrices have rank 5 out of 11. The source residual is
$\|v-\Lambda w\|/\|v\|$, with $v_P=q(S,Q)$, $w_P=q(Q,Q)$ and
the Euclidean least-squares $\Lambda$. It is not a percentage theoretical
error on a decay width. Both mass branches and all rho/K-star and phi-width
choices remain separate.

The universal scenarios can cancel the necessary zero coefficients at their
fixed physical masses; this does not establish their leading coefficients,
full correlator matching or adequacy. The universal-source leptonic test
already fails in [stage 3](../foundations/vector-stage3-calibration.md).
The [chiral-octet control](vector-chiral-octet-matching.md) separately shows
why even these necessary cancellations cannot supply the old OPE prior.

## Outcome and next derivation

The currently declared mass/source enlargement with unchanged odd vertices
has no demonstrated high-energy-consistent and leptonic-adequate joint-fit
prescription. The obstruction is established, rather than bypassed with an
arbitrary momentum cutoff or softened precision prior.

A consistent extension must first specify odd flavor-breaking/source
interactions capable of compensating $S\not\propto Q$, and enough resonance
or momentum structure to satisfy the single- and symmetric-double-virtual
leading terms. The [Roig analysis](../../extracts/Roig:2013baa/tff-correlator-compatibility.md)
provides a source-backed pseudoscalar-resonance route for correlator/falloff
compatibility; the [Knecht analysis](../../extracts/Knecht:2001qf/pion-tff-short-distance.md)
provides LMD+V as a distinct two-vector ansatz route. Neither has been
implemented or established sufficient for this finite-spurion $U(3)$ action.
Do not silently import their couplings, add independently saturated local
terms, or call the present necessary TFF equations full QCD matching.

The subsequent user-selected route is an
[odd single/double-trace extension with low-energy calibration](odd-breaking-low-energy-calibration-plan.md).
It does not require solving the ultraviolet completion before a separately
validated low-energy fit. The scope and limitations are recorded in the
[saved mismatch explanation](vector-matching-mismatch-interpretation.md).
