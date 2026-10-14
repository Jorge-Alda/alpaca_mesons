---
title: Alda quark rotations and the photon anomaly normalization
status: mixed
created: 2026-09-28
updated: 2026-10-15
bibtex_keys: [Alda:2026nrd, Bai:2025fvl, Bauer:2020jbp]
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

These statements are directly evidenced by the
[archived arXiv-v1 source](../../sources/Alda:2026nrd/INDEX.md)
of `Alda:2026nrd`; its bibliography record is already catalogued. [@Alda:2026nrd]

## Source-backed photon operator and separate coefficient rotation

Use the explicit photon-operator normalization of
[@Bauer:2020jbp, Eqs. (39)/(54)], as recorded in the
[source digest](../Bauer:2020jbp/photon-operator-normalization.md):

$$
\mathcal L_{\gamma\gamma}=\frac{\alpha_{\rm em}}{4\pi}
c_\gamma\frac{a}{f_a}F_{\mu\nu}\widetilde F^{\mu\nu}.
$$

Using $\widetilde F^{\mu\nu}=\tfrac12\epsilon^{\mu\nu\rho\sigma}_\gamma F_{\rho\sigma}$
and $N_c$ colored copies, the separately declared axial-Jacobian input under
Alda's rotation is

$$
c_\gamma\to c_\gamma-2N_c\,\operatorname{tr}(Q^2\kappa).
$$

This is the photon-only specialization of the anomaly bookkeeping used in the WZW construction. Bai writes the explicit gauge-field term as $(a/f_a)c_{\mathcal A_1\mathcal A_2}F_{\mathcal A_1}\widetilde F_{\mathcal A_2}$ and embeds the derivative couplings as left/right background-field shifts, Eqs. (2.1) and (2.4). The displayed $\alpha_{\rm em}/(4\pi)$ coefficient convention is fixed by Bauer rather than inferred from Bai's model. Eq. (39) defines the operator coefficient; it alone does not derive the Jacobian or the repository's meson mixing invariant. [@Bai:2025fvl; @Bauer:2020jbp]

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

The initially unspecified photon-dual orientation is completed in the
[odd-action matching record](../../alp/photon/alp-photon-odd-action-phase-matching.md).
Its inherited positive WZW coefficient uses the single canonical tensor
$\epsilon^{0123}=+1$, hence $\epsilon_{0123}=-1$ for $g=(+---)$.
The former calculation used an oppositely oriented RChT source tensor;
the [canonical conversion audit](../../rcht/odd-parity/levi-civita-convention.md)
documents the action-preserving odd-coupling and signed-amplitude conversion.
This orientation is a convention choice rather than an Alda result; no
equality with an unspecified gluon-dual orientation is inferred. The canonical
photon coefficient $c_\gamma$ is unchanged by the source-convention conversion.

Alda's paper explicitly treats amplitudes outside the odd-intrinsic-parity WZW sector; it does not itself define a $c_\gamma$ term. Bauer supplies the photon-operator normalization, while the rotation signs and the $c_q^A$ shift are literature evidence from Alda. The photon Jacobian and leading-nonet WZW meson contribution remain separate from that operator definition. Bai's WZW treatment is used only for the photon/anomaly bookkeeping; no vector-meson or resonance terms are included here. [@Alda:2026nrd; @Bai:2025fvl; @Bauer:2020jbp]
