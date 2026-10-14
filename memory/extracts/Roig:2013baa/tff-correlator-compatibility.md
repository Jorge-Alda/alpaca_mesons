---
title: Single-vector TFF and VVP high-energy compatibility
status: literature
created: 2026-10-13
updated: 2026-10-15
bibtex_keys: [Roig:2013baa]
provenance: literature
---

# Separate high-energy conditions

Displayed odd couplings are adapted to the repository metric $(+---)$ and
$\epsilon^{0123}=+1$, hence $\epsilon_{0123}=-1$. They satisfy
$c_i=-c_i^R$, $d_i=-d_i^R$ relative to the source orientation, and
$\kappa_3^{PV}=-\kappa_3^{PV,R}$ for its epsilon-linear operator;
the even-parity $d_m$ is unchanged. See the
[convention dictionary](../../rcht/odd-parity/levi-civita-convention.md).
The physical actions and vanishing tests are unchanged by this adaptation.

The archived `AcceptedVersion.tex`, section on the pion TFF, gives the
single-virtual vanishing conditions

$$
c_1-c_2+c_5=0,\qquad
c_5-c_6=-\frac{N_CM_V}{32\sqrt2\pi^2F_V}
+\frac{F_V}{\sqrt2M_V}d_3.
$$

The symmetric double-virtual vanishing condition is

$$
c_5-c_6=-\frac{N_CM_V}{64\sqrt2\pi^2F_V}.
$$

These are equations labelled `eq.TFF-constraint1` through
`eq.TFF-constraint3`, source lines 695–728 [@Roig:2013baa]. Vanishing is
distinct from determining the coefficient of the leading inverse power.

The paper distinguishes the vector-only correlator relations from a theory
including pseudoscalar resonances. In the latter, the correlator expressions
for $d_1+8d_2$ and $d_3$ contain an additional
$4\sqrt2d_m\kappa_3^{PV}/F_V$ contribution, and imposing the single-virtual
condition gives $1-32\sqrt2F_Vd_m\kappa_3^{PV}/F^2=0$ (adapted from source
`eq.constraint-Kampf-D/E`, lines 623–630, and
`eq: relation pigammagamma Kampf`, lines 737–741) [@Roig:2013baa].

The discussion explicitly identifies the wrong single-virtual high-energy
behavior obtained from vector-only VVP constraints (source lines 803–822)
[@Roig:2013baa]. Its compatible resonance relations rely on an enlarged
action; they cannot be imported into the repository's vector-only action
without implementing those interactions. The paper's vanishing tests alone
do not establish a correct nonzero symmetric $1/Q^2$ coefficient.

Source archive/navigation: [downloaded source](../../sources/Roig:2013baa/INDEX.md).
