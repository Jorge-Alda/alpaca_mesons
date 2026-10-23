---
title: Odd-fit null sensitivity of the virtual-photon ALP two-pion amplitude
status: derived
created: 2026-10-23
updated: 2026-10-24
bibtex_keys: [Alda:2026nrd, Ecker:1988te, Ruiz-Femenia:2003jdx, Messerli:2025rnv]
provenance: analytical-derivation
---

# Inputs and scope

Check the four coefficient null directions of the preliminary original
seven-VJP/four-VVP fit in the tree-level charged-pion partial wave for
$a\gamma^*(y)\to\pi^+\pi^-$. The ALP sources are in the no-mixing basis,
the retained even vertex is the leading $G_V$ interaction, strong isospin is
used, and the odd coefficients use the repository's canonical
$\epsilon^{0123}=+1$ convention. The odd fit's null basis is documented in
[the fit Lagrangian record](../../rcht/odd-parity/preliminary-fit-lagrangian-and-null-directions.md);
its neutral-rho TFF projection and exact null relations are summarized in
[the rho TFF null projection](alp-isovector-tff-fit-null-projection.md).
The spectral role and required normalization are in
[the isovector TFF construction](alp-isovector-tff-spectral-construction.md).

This checks the retained odd-resonance tree amplitude. It does not include
local anomalous terms, additional odd operators, loops, or pion rescattering.

# Result

Use the crossed decay routing $k=-q_\gamma$, $k^2=y$, with
$s=(p_++p_-)^2$, $t=(k+p_-)^2$, and $u=(k+p_+)^2$. Write the amplitude as
the canonical epsilon tensor times a scalar coefficient $\mathcal A$.
Within the four-dimensional fit null space, this channel depends on the
direction represented by
$$
\delta c_1=\delta c_5=\alpha
$$
with $\delta c_2=\delta c_6=0$. Its scalar variation is
$$
\delta\mathcal A=
-\frac{8\sqrt2 eG_Vy\alpha}{3M_Vf_af_\pi^2}
\left[\frac{2(c_d^A+2c_u^A)}{s-M_\rho^2}
+(c_d^A-c_u^A)\left(\frac1{t-M_\rho^2}
+\frac1{u-M_\rho^2}\right)\right].
$$
The common $\delta c_2=\delta c_5=\delta c_6$ direction, the $c_7$
direction, and the $d_4$ direction project to zero, consistent with the
shared fit-null relations. The VVP coefficient projection is zero, so
changing the fitted photon-vector response does not remove the surviving VJP
dependence. This is the same null direction that changes the neutral-rho
TFF through its off-shell VJP numerator, although the two observables have
different kinematic weights.

For the allowed equal light-quark axial profile $c_u^A=c_d^A=c$, the crossed
terms cancel. With the projector convention
$$
f_1=\frac34\int_{-1}^{1}dz\,(1-z^2)\,\mathcal A(z),
$$
the remaining constant scalar has
$$
\delta f_1=-\frac{16\sqrt2 eG_Vcy\alpha}
{M_Vf_af_\pi^2(s-M_\rho^2)}.
$$
It is nonzero at $y\ne0$ and vanishes at $y=0$. Thus a chosen zero-null
representative is needed for a numerical virtual-photon partial wave unless
the ALP source profile or additional information fixes this coefficient
combination.

# Relation to Messerli:2025rnv

The analogous $\eta\gamma^*\to\pi\pi$ amplitude also enters the isovector
double-spectral density. Messerli:2025rnv includes its $P$ wave
$\mathcal F_{\eta\pi\pi\gamma^*}$ in the density, constructed from a
dispersive representation of $\eta\to2(\pi^+\pi^-)$, with pairwise pion
rescattering and an $a_2(1320)$ left-hand-cut contribution. It does not
obtain that amplitude from this repository's preliminary RChT fit, so the
fit's four null directions are not an uncertainty prescription used in
that construction. Its dispersive model and input uncertainties are a
separate issue. The ALP calculation needs its own source-dependent partial
wave and uncertainty treatment; replacing the ALP sources by eta-like fitted
weights would not resolve the null direction. [@Messerli:2025rnv]

# Reproducibility

The executable [`derive_alp_virtual_photon_pion_nulls.py`](../../../scripts/derive_alp_virtual_photon_pion_nulls.py)
computes the fit-null projections before and after the charged-pion flavor
sum and checks the real-photon limit, the surviving $P$ wave, and vector-basis
normalization. Run:

```text
.venv/bin/python scripts/derive_alp_virtual_photon_pion_nulls.py
.venv/bin/python scripts/derive_alp_virtual_photon_pion_nulls.py --test
```
