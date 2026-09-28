---
title: Alda quark rotations and the photon anomaly normalization
status: mixed
created: 2026-09-28
updated: 2026-09-28
bibtex_keys: [Alda:2026bar, Bai:2025fvl]
provenance: literature
---

# Photon anomaly convention

## Evidence

Alda defines the ALP-dependent quark-field redefinition (Eq. (II.2), p. 6) by

$$
q_{L,R}\to g_{L,R}q_{L,R},\qquad
g_{L,R}=\exp\!\left[\frac{i a}{f_a}(\delta\mp\kappa)\right].
$$

Thus the upper sign is for $q_L$ and the lower sign for $q_R$:

$$
q_L\to e^{+ia(\delta-\kappa)/f_a}q_L,
\qquad
q_R\to e^{+ia(\delta+\kappa)/f_a}q_R.
$$

For the purely axial part, this is equivalently

$$
q\to \exp\!\left(+i\frac{a}{f_a}\kappa\gamma_5\right)q.
$$

Here the standard eigenvalues $\gamma_5q_L=-q_L$ and $\gamma_5q_R=+q_R$
reproduce the two chiral phases above.

With Alda's quark-level normalization

$$
\mathcal L\supset \frac{\partial_\mu a}{f_a}
\bar q\gamma^\mu(c_q^V+c_q^A\gamma_5)q
+\frac{\alpha_s}{4\pi}c_G\frac{a}{f_a}G\widetilde G,
$$

the coefficients transform as (Eq. (II.3))

$$
c_q^A\to c_q^A+\kappa,\qquad c_G\to c_G-\operatorname{tr}\kappa.
$$

These statements are directly evidenced by `Alda:2026bar`, whose available local copy is a manuscript and is not currently indexed by INSPIRE-HEP. [@Alda:2026bar]

## Derived photon coefficient under a stated normalization assumption

Assume the explicit photon term is normalized as

$$
\mathcal L_{\gamma\gamma}=\frac{\alpha_{\rm em}}{4\pi}
c_\gamma\frac{a}{f_a}F_{\mu\nu}\widetilde F^{\mu\nu}.
$$

Using the same $\widetilde F^{\mu\nu}=\tfrac12\epsilon^{\mu\nu\rho\sigma}F_{\rho\sigma}$ convention and $N_c$ colored copies, assume the axial Jacobian under Alda's rotation is normalized as

$$
c_\gamma\to c_\gamma-2N_c\,\operatorname{tr}(Q^2\kappa).
$$

This is the photon-only specialization of the anomaly bookkeeping used in the WZW construction. Bai writes the explicit gauge-field term as $(a/f_a)c_{\mathcal A_1\mathcal A_2}F_{\mathcal A_1}\widetilde F_{\mathcal A_2}$ and embeds the derivative couplings as left/right background-field shifts, Eqs. (2.1) and (2.4); that literature evidence motivates the bookkeeping but does not by itself fix the displayed $\alpha_{\rm em}/(4\pi)$ normalization. [@Bai:2025fvl]

It then follows algebraically that the basis-invariant photon combination is

$$
c_\gamma^{\rm inv}=c_\gamma+2N_c\operatorname{tr}(Q^2c_q^A).
$$

For $N_c=3$ and $Q=\operatorname{diag}(2/3,-1/3,-1/3)$,

$$
c_\gamma^{\rm inv}
=c_\gamma+\frac{8}{3}c_{q,u}^A
+\frac{2}{3}c_{q,d}^A
+\frac{2}{3}c_{q,s}^A.
$$

## Scope and normalization caveat

Alda's paper explicitly treats amplitudes outside the odd-intrinsic-parity WZW sector; it does not itself define a $c_\gamma$ term. The photon shift above therefore requires the displayed anomaly normalization assumption, while the rotation signs and the $c_q^A$ shift are literature evidence from Alda. Bai's WZW treatment is used only for the photon/anomaly bookkeeping; no vector-meson or resonance terms are included here. [@Alda:2026bar; @Bai:2025fvl]
