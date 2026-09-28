---
title: One-ALP photon coefficient from WZW and meson mixing
status: derived
created: 2026-09-28
updated: 2026-09-28
bibtex_keys: [Alda:2026bar, Bai:2025fvl]
provenance: analytical-derivation
---

# One-ALP Photon Coefficient from WZW and Meson Mixing

## Scope and provenance

This record derives the result with
[`scripts/derive_alp_photon_invariant.py`](../scripts/derive_alp_photon_invariant.py).
The literature inputs are Alda's field redefinition, LSZ coefficients, inverse
operator, and generator completeness [@Alda:2026bar, Eqs. (II.2), (II.3),
(III.25), (III.28), and (B.1)--(B.2)], as recorded in
[`appendix-a-expansion.md`](extracts/Alda:2026bar/appendix-a-expansion.md) and
[`photon-anomaly-normalization.md`](extracts/Alda:2026bar/photon-anomaly-normalization.md).
Bai's background-field construction establishes that the complete gauged WZW
functional and its counterterm must be combined for gauge invariance and
rotation-parameter cancellation; this is literature evidence recorded in
[`background-fields-wzw.md`](extracts/Bai:2025fvl/background-fields-wzw.md)
[@Bai:2025fvl]. The photon formulas below are derived claims, not claims made in
Alda, whose stated scope excludes the WZW sector.

## Conventions, units, and assumptions

All $c$ coefficients are dimensionless, $a$ and $f_a$ have mass dimension one,
and

$$
\mathcal L_{\gamma\gamma}
=\frac{\alpha_{\rm em}}{4\pi}c_\gamma\frac{a}{f_a}
F_{\mu\nu}\widetilde F^{\mu\nu},
\qquad
\widetilde F^{\mu\nu}=\frac12\epsilon^{\mu\nu\rho\sigma}F_{\rho\sigma}.
$$

The dimensionful coefficient multiplying $aF\widetilde F$ is therefore
$\alpha_{\rm em}C_{a\gamma\gamma}/(4\pi f_a)$; omitting the displayed $1/f_a$
would change the coefficient's mass dimension and normalization.

Alda's redefinition is

$$
q_L\to e^{ia(\delta-\kappa)/f_a}q_L,
\qquad
q_R\to e^{ia(\delta+\kappa)/f_a}q_R.
$$

For $\delta=0$ and the standard eigenvalues $\gamma_5q_L=-q_L$ and
$\gamma_5q_R=+q_R$, this is

$$
q\to e^{+ia\kappa\gamma_5/f_a}q.
$$

The associated transformations used here are

$$
c_q^A\to c_q^A+\kappa,
\qquad
c_\gamma\to c_\gamma-2N_c\operatorname{tr}(Q^2\kappa).
$$

The photon shift is an explicit anomaly-normalization assumption consistent
with the displayed photon operator. We retain one ALP at $O(1/f_a)$, the
leading pseudoscalar WZW photon vertex, and no vector mesons or resonance
terms. The inverse $\mathbb M_a^{-1}$ is assumed to exist, so the derivation is
not applicable at an ALP--meson degeneracy. No flavour-diagonality assumption
is imposed on $c_q^A$, $C$, or $\mathbb M_a^{-1}(C)$.

## Step 1: explicit-contact invariant

Applying the two coefficient shifts gives

$$
\begin{aligned}
c_\gamma' +2N_c\operatorname{tr}(Q^2c_q^{A\prime})
&=c_\gamma-2N_c\operatorname{tr}(Q^2\kappa)
  +2N_c\operatorname{tr}[Q^2(c_q^A+\kappa)]\\
&=c_\gamma+2N_c\operatorname{tr}(Q^2c_q^A).
\end{aligned}
$$

Hence

$$
\boxed{c_\gamma^{\rm inv}
=c_\gamma+2N_c\operatorname{tr}(Q^2c_q^A)}.
$$

## Step 2: LSZ components and the WZW contribution

Alda's physical-ALP meson components are

$$
w_i^{(a)}=\frac{f_\pi}{f_a}\left[
\frac{\operatorname{tr}(C\lambda_i^\dagger)}{m_i^2-m_a^2}
+\operatorname{tr}(c_q^A\lambda_i^\dagger)\right].
$$

The spectral definition of $\mathbb M_a^{-1}$ makes the bracket

$$
\operatorname{tr}\!\left[
\left(c_q^A+\mathbb M_a^{-1}(C)\right)\lambda_i^\dagger
\right].
$$

With $\operatorname{tr}(\lambda_i\lambda_j^\dagger)=2\delta_{ij}$,
completeness states

$$
A=\frac12\sum_i\operatorname{tr}(A\lambda_i^\dagger)\lambda_i.
$$

The no-vector pseudoscalar WZW vertex is taken with the sign required by the
Alda anomaly convention and the coefficient shifts above:

$$
\mathcal L_{\rm WZW}\supset
\frac{\alpha_{\rm em}}{4\pi}\frac{N_c}{f_\pi}
\operatorname{tr}(Q^2\Phi)F\widetilde F.
$$

Replacing each meson field by its ALP LSZ component gives

$$
\begin{aligned}
\Delta C_{a\gamma\gamma}^{\rm WZW/mix}
&=N_c\frac{f_a}{f_\pi}
\sum_i w_i^{(a)}\operatorname{tr}(Q^2\lambda_i)\\
&=2N_c\operatorname{tr}\!\left[
Q^2\left(c_q^A+\mathbb M_a^{-1}(C)\right)\right].
\end{aligned}
$$

Adding the explicit contact term yields

$$
\boxed{
C_{a\gamma\gamma}(m_a)
=c_\gamma^{\rm inv}
+2N_c\operatorname{tr}[Q^2\mathbb M_a^{-1}(C)]
=c_\gamma+2N_c\operatorname{tr}\!\left[
Q^2(c_q^A+\mathbb M_a^{-1}(C))\right]}.
$$

Here $C$ and therefore $\mathbb M_a^{-1}(C)$ are invariant under Alda's field
redefinition. The split between the explicit contact term and the
WZW/meson-pole term is basis dependent; only their displayed sum is invariant.

## Step 3: Alda no-mixing basis

Choose

$$
\kappa=-c_q^A-\mathbb M_a^{-1}(C).
$$

Then

$$
\boxed{\widehat c_q^A=-\mathbb M_a^{-1}(C)}
$$

and the transformed explicit photon coefficient is

$$
\begin{aligned}
\widehat c_\gamma
&=c_\gamma-2N_c\operatorname{tr}(Q^2\kappa)\\
&=c_\gamma+2N_c\operatorname{tr}\!\left[
Q^2(c_q^A+\mathbb M_a^{-1}(C))\right]\\
&=\boxed{C_{a\gamma\gamma}(m_a)}.
\end{aligned}
$$

Thus the meson-pole contribution vanishes in the no-mixing variables and its
effect is carried entirely by the transformed explicit contact coefficient.

## Three-flavour electromagnetic specialization

For $N_c=3$ and
$Q=\operatorname{diag}(2/3,-1/3,-1/3)$,

$$
c_\gamma^{\rm inv}=c_\gamma+\frac83(c_q^A)_{uu}
+\frac23(c_q^A)_{dd}+\frac23(c_q^A)_{ss},
$$

and

$$
\begin{aligned}
C_{a\gamma\gamma}(m_a)=c_\gamma
&+\frac83\left[(c_q^A)_{uu}+(\mathbb M_a^{-1}(C))_{uu}\right]\\
&+\frac23\left[(c_q^A)_{dd}+(\mathbb M_a^{-1}(C))_{dd}\right]\\
&+\frac23\left[(c_q^A)_{ss}+(\mathbb M_a^{-1}(C))_{ss}\right].
\end{aligned}
$$

The flavour matrices remain fully general: off-diagonal entries are retained
throughout the completeness derivation, but only diagonal entries can enter
$\operatorname{tr}(Q^2A)$ because $Q^2$ is diagonal.

## Verification

Run

```text
.venv/bin/python scripts/derive_alp_photon_invariant.py
```

The script constructs a complete normalized $U(3)$ generator basis, applies
completeness to fully general $3\times3$ matrices, derives the WZW sum, and
asserts the contact invariant, the full coefficient, both no-mixing identities,
and the absence of off-diagonal entries after the diagonal-$Q$ trace.
