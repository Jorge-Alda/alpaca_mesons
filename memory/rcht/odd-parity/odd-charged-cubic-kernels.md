---
title: Normalized charged cubic kernels of the bounded odd action
status: derived
created: 2026-10-14
updated: 2026-10-14
bibtex_keys: [Ruiz-Femenia:2003jdx, Alda:2026nrd]
provenance: analytical-derivation
---

# Inputs and scope

Inputs are the [normalized neutral parent maps](odd-power-counted-channel-maps.md),
[leading SM action normalization](radiative-sm-amplitude.md), and
[bounded quotient](odd-higher-covariant-operator-count.md). The source and
field conventions follow Alda:2026nrd; the original odd action follows
Ruiz-Femenia:2003jdx. The extension and selected tree truncation are local
assumptions. Use GeV, metric $(+---)$, source $\epsilon_{0123}=+1$, incoming
momenta, canonical tensor-vector LSZ, and $\chi=2B_0m_q$.

# Calculation

Run

```text
.venv/bin/python scripts/derive_odd_charged_kernels.py
```

The script contracts the Lorentz tensors and evaluates the full flavor words
for conjugate charged legs $T_P=E_{ij}$ and $T_V=E_{ji}$. It retains
$u=-\sqrt2 T_P/F$, $Y=\sqrt2\{T_P,\chi\}/F$ and the photon commutator in
$f_-=-i\sqrt2 e[T_P,Q]F_\gamma/F$. The derivative of this last building
block has momentum

$$
p_P+q_\gamma=-q_V,
$$

because it differentiates the complete product. Assigning it only
$q_\gamma$ drops the pseudoscalar derivative. Exact tensor reduction gives
the nine higher $f_-$ parent images, in saved parent order,

$$
(0,-2s,0,0,2s,0,0,-2s,0),\qquad s=q_V^2.
$$

Their cubic kinematic image has rank one; six vanish. This is a projection
result, not an action-basis reduction. The surviving term probes the hard
on-shell vector mass and is consistent with the
[derivative promotion prescription](odd-resonance-derivative-counting.md).

The original $c_4$ gradients are

$$
\frac{8\sqrt2 e}{F M_1M_V}
(\chi_d-\chi_u,\;\chi_s-\chi_u,\;0)
$$

for $\rho^+\pi^+\gamma$, $K^{*+}K^+\gamma$, and
$K^{*0}K^0\gamma$, respectively. Here the superscripts denote conjugate
incoming legs; $M_1$ is external LSZ mass and $M_V$ the operator normalization.

# Verification and outputs

All 13 saved action relations vanish in each channel at
$m_{P,ij}^2=(\chi_i+\chi_j)/2$. This is the LO charged EOM mass, not a
physical-mass replacement. The artifact retains residuals for independent
mass inputs, the signed symbolic internal response, and all 241 raw parent
gradients. Raw parents remain redundant: the covariant quotient has 228
independent representatives.

The [JSON artifact](../../../data/vector_sector/odd_cubic_parent_kernels.json)
also exports compact flavor words, Fourier phases, Bose factors, normalization
powers and scalar polynomials for subsequent numerical assembly. The
[transcript](../../../data/vector_sector/odd_charged_cubic_kernels.md) records
the derivation. Photon Ward tensors and exact wedge reconstruction pass.
No measurement Jacobian or fit is claimed.
