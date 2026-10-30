---
title: Possible a2 exchange sensitivity in the inclusive KS K pi fit
status: mixed
created: 2026-10-30
updated: 2026-10-30
bibtex_keys: [Qin:2024ulb, BaBar:2007ceh]
provenance: numerical-computation
---

# Inputs and established structure

The full inclusive form factors in [@Qin:2024ulb, Eq. (kskpi)] contain
an $a_2(1320)$ contribution in $e^+e^-\to K_S^0K^\pm\pi^\mp$:
$$
F_{a_2}^{3}=
\frac{t-u}{M_{a_2}^2-s-iM_{a_2}\Gamma_{a_2}}
\left[
-\frac{2g_T(2c_1^T-2c_2^T-c_3^T)}{F^3}
+\frac{\sqrt2 F_Vg_T(2d_1^T-2d_2^T-d_3^T)}
{F^3(M_\rho^2-q^2-iM_\rho\Gamma_\rho(q^2))}
\right].
$$
This term is absent from the charged-kaon neutral-pion channel's direct
tensor contribution. The source's $K_2^*(1430)$ terms are separate.

For total energy below $1.78$ GeV, the largest $K\pi$ mass is about
$1.78-m_K=1.286$ GeV. Neither $K^*(1410)$ nor $K_2^*(1430)$ can therefore
occur on shell in that subchannel. The largest $KK$ mass is about
$1.78-m_\pi=1.640$ GeV, so an $a_2(1320)$ band is allowed once the total
energy exceeds $M_{a_2}+m_\pi\simeq1.458$ GeV.

# Result and interpretation

The unified inclusive fit in
[babar-kkpi-unified-vector-coefficient-fit.md](babar-kkpi-unified-vector-coefficient-fit.md)
includes both $a_2$ tensor products as nuisance parameters and fits them
together with all ten vector-sector coefficient combinations. Its full
Tables 2/3 residual is $\chi^2/\mathrm{dof}=35.716/29$. The fit does not
add the Table 7 isobar component as independent data.

A source-constrained control with two vector multiplets gives
$\chi^2=85.70/38$ without the $a_2$ term and $85.33/36$ with it. This small
change does not establish a separate $a_2$ determination. In the earlier
pole-shape control $D_T$ reached its numerical upper bound, so that
diagnostic should not be interpreted as a tensor-coupling measurement.
The unified fit's coefficient profiles remain broad and uncalibrated.
