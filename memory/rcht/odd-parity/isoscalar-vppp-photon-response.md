---
title: Tree photon-source response of the local isoscalar VPPP terms
status: derived
created: 2026-10-27
updated: 2026-10-27
bibtex_keys: [GomezDumm:2010ra, Ecker:1988te, Ecker:1989yg]
provenance: analytical-derivation
---

# Assumptions and inputs

This calculation inserts the isoscalar charged-three-pion local vertex from
the [VPPP projection](isoscalar-vppp-three-pion-projection.md) into the
canonical tree photon--vector response. It uses the quadratic antisymmetric-
tensor source term of [@Ecker:1988te; @Ecker:1989yg], the stage-four response
convention recorded in [the vector-amplitude record](vector-stage4-amplitudes.md),
and $v_\mu=-eQA_\mu$. The calculation is at tree level with exact isospin,
equal on-shell pion masses, and no vector self-energies, local source
contacts, or pion rescattering. The scalar divides out the common
$i\epsilon(p_+,p_-,p_0,\varepsilon_\gamma)$ and is written divided by $e$.
No additional external-current phase is assumed.

# Executed source insertion

Varying the quadratic tensor expression
$$
\frac{F_\gamma^2}{4}\left[(M_i^2-s)r_i^2
-2\sqrt2 eR_i r_i\right]
$$
gives $V_i/F_\gamma=\sqrt2 eR_i/(M_i^2-s)$. Inserting this response removes
the fixed external-state factor $1/M_i$ in the vector-decay projection. The
result for the local isoscalar source amplitude is
$$
\frac{\mathcal F_{\gamma,\mathrm{local}}(s)}{e}=
\frac{16\sqrt2}{M_Vf_\pi^3}
\sum_{i=\in\{\omega,\phi\}}
\frac{t_{n,i}R_i}{M_i^2-s}
\left[(g_1-g_2-g_3)(s-3m_\pi^2)
+3m_\pi^2(2g_4+g_5)\right].
$$
Here $R_i$ is the signed vector-source residue in the convention of the
stage-four record, and $t_{n,i}$ is the nonstrange diagonal entry of the
physical vector generator. With its mixing convention,
$t_{n,\omega}=\cos\theta_V/\sqrt2$ and
$t_{n,\phi}=\sin\theta_V/\sqrt2$. The ideal limit removes the $\phi$ term.
The executable [source-response script](../../../scripts/derive_vppp_photon_response.py)
also checks the pole residue against the external-vector vertex; this check
prevents retaining an erroneous $1/M_i$ factor in the photon-source result.

# Pole/contact degeneracy relevant to a fit

For one bare vector pole, the local polynomial obeys
$$
\frac{(g_1-g_2-g_3)(s-3m_\pi^2)
+3m_\pi^2(2g_4+g_5)}{M_i^2-s}
 =-(g_1-g_2-g_3)
+\frac{(g_1-g_2-g_3)(M_i^2-3m_\pi^2)
+3m_\pi^2(2g_4+g_5)}{M_i^2-s}.
$$
Thus the two combinations map to a constant and a pole residue only when
the constant source term is fixed independently. If an independent free
constant contact is added, the local combinations' profiled rank drops from
two to one: the first combination's constant shift can be absorbed, leaving
only its contribution to the pole residue combined with the second
combination. This is a structural degeneracy and must be handled in the fit
parameterization rather than hidden with a narrow prior.

# Reproducibility and limits

Run:

    .venv/bin/python scripts/derive_vppp_photon_response.py

The SymPy derivation prints the stationary source response, the physical
$\omega/\phi$ sum, the pole/contact decomposition, and the coefficient ranks.
Its focused checks cover stationarity, pole-residue normalization, the ideal
$\omega$ source limit, null directions, and transversality. This determines
the bare tree response only. A timelike dispersive fit still needs a declared
self-energy and local-source prescription plus the rescattered complex
$\gamma^*\to3\pi$ amplitude.
