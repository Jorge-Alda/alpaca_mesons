---
title: Current projections and vector-pole bookkeeping for the ALP TFF
status: derived
created: 2026-10-25
updated: 2026-10-25
bibtex_keys: [Alda:2026nrd, Bauer:2020jbp, Ecker:1988te, Ecker:1989yg, Ruiz-Femenia:2003jdx, Messerli:2025rnv]
provenance: analytical-derivation
---

# Inputs, conventions, and scope

Use the independent diagonal no-mixing sources
$U_a=2c_q^A/f_a$ and $Y_a=-2\{c_q^P,\chi\}/f_a$, the ordered
$V_{i\gamma}$ and $P_{ij}$ numerators, and the scalar convention
$i\mathcal M=i e^2R F_{a\gamma^*\gamma^*}$ in the
[tree TFF record](alp-doubly-virtual-tff.md). Their source provenance is the
[LO source/LSZ assembly](../../rcht/odd-parity/odd-alp-source-lsz-assembly.md)
[@Alda:2026nrd; @Bauer:2020jbp]. No source EOM or equality of axial and
pseudoscalar traces is imposed in this calculation.

The vector generators and the absence of rho--isoscalar mixing in the
even sector are the assumptions of the
[neutral-rho record](alp-isovector-doubly-virtual-tff.md). The linear photon
source, including its single- and double-trace terms, is taken from the
[vector response record](../../rcht/odd-parity/vector-stage4-amplitudes.md)
[@Ecker:1988te; @Ecker:1989yg]. Retain stable tree poles, the original
eleven odd operators, and leading electromagnetic order
[@Ruiz-Femenia:2003jdx]. Units are GeV, the metric is $(+---)$, and
$\epsilon^{0123}=+1$. The scalar TFF has dimension $\mathrm{GeV}^{-1}$
[@Messerli:2025rnv]. This is a decomposition of the retained tree
expression, not a spectral or high-energy completion.

# Executable derivation

The paired [script](../../../scripts/derive_alp_isoscalar_tff_decomposition.py)
does the following in SymPy.

1. Interchange the up and down entries of the supplied charge matrix $Q$,
   then take its symmetric and antisymmetric parts. This derives the
   three-flavor current generators
   $$
   Q^{(0)}=\operatorname{diag}(1/6,1/6,-1/3),\qquad
   Q^{(1)}=\operatorname{diag}(1/2,-1/2,0),\qquad Q=Q^{(0)}+Q^{(1)}.
   $$
   The upper light-quark block of $Q^{(0)}$ is the existing
   $Q^{(0)}_\ell$. The strange entry must be retained for a full
   three-flavor current projection.

2. Insert each generator into the existing photon source. In the notation
   $R_i[Q^{(a)}]$ and $W_i(s;Q^{(a)})$, the extra argument only identifies
   which current is inserted in the already defined residue or response:
   $$
   \begin{aligned}
   R_i[Q^{(a)}]={}&F_V\operatorname{tr}(T_iQ^{(a)})
   +z_V\operatorname{tr}(T_i\{X,Q^{(a)}\})
   +y_V\operatorname{tr}(T_i)\operatorname{tr}(XQ^{(a)}),\\
   W_i(s;Q^{(a)})={}&\frac{R_i[Q^{(a)}]}{M_i^2-s},\qquad
   X=\operatorname{diag}(x_l,x_l,x_s).
   \end{aligned}
   $$
   The executed traces give
   $$
   R_\rho[Q^{(0)}]=0,\qquad
   R_\rho[Q^{(1)}]=\frac{F_V+2x_lz_V}{\sqrt2},\qquad
   R_\omega[Q^{(1)}]=R_\phi[Q^{(1)}]=0.
   $$
   The script also prints the full mixed omega/phi residues for $Q^{(0)}$.
   It verifies $R_i[Q^{(0)}]+R_i[Q^{(1)}]=R_i[Q]$ without a numerical fit.

3. Reuse the existing Lorentz-contraction calculation to derive each
   ordered VJP numerator with $Q$ replaced only at its explicit current
   vertex, and each VVP numerator with independent $U_a,Y_a$. Attach
   independent formal labels to the two external currents before
   assembling either photon placement. Differentiating this expression
   with respect to those formal labels derives the coefficient for
   photon-current choices $Q^{(a)},Q^{(b)}$:
   $$
   \begin{aligned}
   &2\sqrt2\sum_i\big[
   W_i(s_1;Q^{(a)})V_{i\gamma}(s_1,s_2;Q^{(b)})
   +W_i(s_2;Q^{(b)})V_{i\gamma}(s_2,s_1;Q^{(a)})\big]\\
   &\quad-2\sum_{i,j}W_i(s_1;Q^{(a)})W_j(s_2;Q^{(b)})
   P_{ij}(s_1,s_2),\qquad a,b\in\{0,1\}.
   \end{aligned}
   $$
   The notation after the semicolon specifies the explicit charge matrix
   in the existing VJP trace, not a new coupling. Summing all four
   coefficients reconstructs the complete retained resonance scalar.
   Exchanging the photon labels sends $(a,b)$ to $(b,a)$ and leaves the
   reconstructed expression unchanged. The polynomial is linear in each
   of the six independent diagonal no-mixing source entries, with no
   source-independent resonance term.

4. Independently select contributions with at least one rho response and
   contributions with at least one omega/phi response in the full retained
   tree expression. The two selections share exactly
   $$
   -2\sum_{I=\omega,\phi}\big[
   W_\rho(s_1)W_I(s_2)P_{\rho I}(s_1,s_2)
   +W_I(s_1)W_\rho(s_2)P_{I\rho}(s_1,s_2)\big].
   $$
   The script derives this overlap by the ordered sums and verifies that
   subtracting it once from their sum reconstructs the full resonance
   scalar. It is nonzero for general diagonal ALP sources. Thus adding an
   omega/phi selection to the existing rho selection without this
   subtraction double counts mixed VVP terms. Selecting a hadronic pole
   is also different from projecting both photon currents onto $Q^{(0)}$.

5. Construct the coefficient of $\varepsilon_1^\mu\varepsilon_2^\nu$
   in $R$ directly with SymPy's Levi-Civita tensor. Both contractions
   with the respective photon momenta vanish. These are Ward identities
   of the retained two-photon tensor; they do not verify an unconstructed
   off-shell multiparticle production amplitude. Finally insert the
   stable-pole responses into the assembled scalar and subtract the
   existing `assemble_scalar` result. The difference is identically zero.

# Contact contribution and interpretation

Keep the supplied $+\widehat c_\gamma/(4\pi^2f_a)$ contribution separate
from the four resonance coefficients above. The physical electromagnetic
coefficient alone does not specify an independently matched two-current
local term for all pairs of projected external currents. A microscopic
external-source completion or a declared allocation convention is needed
before dividing that total contact among such pairs. This calculation
does not choose one. The total unprojected TFF, including its existing
contact, remains fully specified within the stated retained tree
truncation.

For a three-pion unitarity cut, project the hadronic state on that cut and
retain the full electromagnetic current on the ALP-production side.
The isospin projection of the external current does not replace that
state projection. In particular, arbitrary diagonal ALP sources retain
mixed rho--isoscalar traces. No complex three-pion amplitude, phase,
subtraction function, or continuum is generated by this bookkeeping.

# Verification and remaining steps

Executed invocations:

```text
.venv/bin/python scripts/derive_alp_isoscalar_tff_decomposition.py
.venv/bin/python scripts/derive_alp_isoscalar_tff_decomposition.py --test
```

All four unit tests pass. They cover charge/residue recombination, ordered
current recombination and Bose exchange, source linearity and mixed-term
overlap, and stable-pole normalization plus both photon Ward identities.
The script prints its inputs, component residues and vertices, the four
derived current coefficients, and zero differences for the checks.

The supported objective is $2m_e<m_a<2m_\tau$ with the virtualities of the
full two-photon loop. This low-energy tree derivation supplies no uniform
description over that domain or the unbounded loop virtualities. The
[isoscalar plan](alp-isoscalar-tff-plan.md) still requires compatible
three-pion and kaon amplitudes, pole/continuum separation, higher states,
source-specific short-distance matching, and any local ALP--lepton term.
