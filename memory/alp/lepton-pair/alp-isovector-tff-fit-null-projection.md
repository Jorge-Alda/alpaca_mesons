---
title: Preliminary fit null projection onto the neutral-rho ALP TFF term
status: numerical
created: 2026-10-23
updated: 2026-10-24
bibtex_keys: [Alda:2026nrd, Ecker:1988te, Ecker:1989yg, Ruiz-Femenia:2003jdx]
provenance: numerical-computation
---

# Scope and inputs

This projects the four-dimensional fixed-nuisance null space of the selected
`central_single_trace_negative` preliminary fit onto the tree-level rho part
of the doubly virtual ALP TFF. The fit null basis and its limitations are
documented in the [preliminary-fit coefficient record](../../rcht/odd-parity/preliminary-fit-lagrangian-and-null-directions.md).
The target numerator follows the
[rho TFF derivation](alp-isovector-doubly-virtual-tff.md). Odd coefficients
use the canonical sign convention of the
[Levi-Civita conversion record](../../rcht/odd-parity/levi-civita-convention.md).

The calculation fixes the even-sector and photon-source parameters at the
preliminary fit point. The ALP axial and pseudoscalar sources remain
independent. This is a coefficient-null projection, not a joint uncertainty
profile over the fixed nuisances and not an estimate of the full hadronic TFF.

# Projection

In the ordered coefficient vector
$(c_1,\ldots,c_7,d_1,\ldots,d_4)$, the regenerated null space satisfies
$$
\delta c_3=\delta c_4=\delta d_1=\delta d_2=\delta d_3=0,\qquad
\delta c_5=\delta c_1+\delta c_2,\qquad
\delta c_6=\delta c_2,
$$
while $\delta c_1$, $\delta c_2$, $\delta c_7$ and $\delta d_4$ are free
coordinates. Insert these relations into
$$
L(s,t)=(c_1+c_2-c_5)m_a^2
+(-c_1+c_2+c_5-2c_6)s+(c_1-c_2+c_5)t.
$$
The mass term and the coefficient of $s$ cancel, leaving
$$
\delta L(s,t)=2\delta c_1 t,\qquad
\delta V_{\rho\gamma}(s,t)=\frac{4A_\rho}{M_V}\delta c_1 t.
$$
The rho-rho and rho-isoscalar VVP numerators have zero projection because
$d_1,d_2,d_3$ do not vary along the fit null space, while the neutral $d_4$
contraction vanishes. Therefore
$$
\boxed{\left.\delta F_{a\gamma^*\gamma^*}(s_1,s_2)\right|_{\text{at least one rho}}
=\frac{8\sqrt2 A_\rho}{M_V}\delta c_1
\left[s_2W_\rho(s_1)+s_1W_\rho(s_2)\right].}
$$
In the regenerated orthonormal SVD coordinates $\boldsymbol\alpha$,
$$
\delta c_1=0.623446314172\alpha_1
-0.459681627058\alpha_2+0.002737716309\alpha_3,
$$
and $\alpha_4$ is the $d_4$ direction. Thus the functional image of the
four-dimensional fit null space has rank one for generic virtualities and
generic ALP axial source. The other three null combinations are invisible to
this rho TFF projection. The image vanishes for two real photons; it is
generically nonzero for singly or doubly virtual photons. It also vanishes at
this order if $A_\rho=0$ for the chosen ALP axial source.

# Reproduction and interpretation

Run
```text
.venv/bin/python scripts/project_alp_isovector_tff_fit_nulls.py
```
The executable regenerates the null basis through
`summarize_odd_preliminary_fit.main(write_output=False)`, so it reads the fit
inputs without rewriting the production JSON. It projects the symbolic VJP,
rho-rho VVP and rho-isoscalar VVP coefficient numerators and computes the
singular values of their independent polynomial coefficients. The reported
image rank is one at absolute tolerance $10^{-9}$; the exact null identities
have maximum numerical residual $5.6\times10^{-15}$.

The central fit representative can be used as a conditional benchmark, but
its rho TFF value is not uniquely determined by the measured fit rows: the
unmeasured coordinate $\delta c_1$ changes its off-shell shape. Setting that
coordinate to zero selects the stored representative and is not a measured
constraint. The result applies to the retained tree-level operator set and
fixed fitted nuisances only. It is the neutral-rho pole projection; the
two-pion double-spectral contribution is a separate representation and must
not be added without a matching prescription, as explained in the
[spectral-construction record](alp-isovector-tff-spectral-construction.md).
The related charged-pion P-wave projection is recorded in
[the virtual-photon pion null note](alp-virtual-photon-pion-null-sensitivity.md).
