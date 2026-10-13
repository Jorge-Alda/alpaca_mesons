---
title: Flavor factors for space-like pseudoscalar transition form factors
status: derived
created: 2026-10-13
updated: 2026-10-13
bibtex_keys: [Ruiz-Femenia:2003jdx, Aloni:2018vki]
provenance: analytical-derivation
---

# Flavor factors for $P\to\gamma^*\gamma$

## Assumptions and inputs

Use the nonet unit fields and $η$--$η'$ rotation in
[the SM flavor algebra derivation](radiative-sm-flavor-algebra.md):
$T_{\pi^0}=\mathrm{diag}(1,-1,0)/\sqrt2$,
$T_8=\mathrm{diag}(1,1,-2)/\sqrt{12}$, and
$T_1=\mathrm{diag}(1,1,1)/\sqrt3$. Set
$T_\eta=\cos\theta T_8+\sin\theta T_1$ and
$T_{\eta'}=-\sin\theta T_8+\cos\theta T_1$. The electromagnetic matrix is
$Q=\mathrm{diag}(2/3,-1/3,-1/3)$. The source identifies the one-off-shell-photon
$t_2$ flavor structure as $\mathrm{tr}(Q\{Q,T_P\})$
[@Ruiz-Femenia:2003jdx, source lines 440--444]. The illustrative numerical
mixing uses $\sin\theta=-1/3$, $\cos\theta=2\sqrt2/3$
[@Aloni:2018vki], as in the parent derivation.

## Derivation

Since all matrices are diagonal, $\mathrm{tr}(Q\{Q,T_P\})=2\sum_iQ_i^2(T_P)_{ii}$.
Direct substitution gives

$$
F_{\pi^0}=\frac{\sqrt2}{3},\qquad
F_\eta=\frac{\cos\theta+4\sin\theta}{3\sqrt3},\qquad
F_{\eta'}=\frac{4\cos\theta-\sin\theta}{3\sqrt3}.
$$

For the stated illustrative angle these are $\sqrt2/3$,
$-\sqrt3/18$, and $5\sqrt6/18$, respectively. The overall signs follow the
declared meson-field convention; measured widths alone do not determine them.
The executable derivation is
[`scripts/derive_radiative_tff_flavor.py`](../../../scripts/derive_radiative_tff_flavor.py).

## Scope and next input needed

These are flavor factors only. They do not define a physical $Q^2F_{P\gamma^*\gamma}$
predictor or an absolute normalization: the cited RChT paper gives the on-shell
$\pi^0\to\gamma\gamma$ correction and says $t_2$ contributes with one photon
off shell, but does not print an off-shell form-factor formula
[@Ruiz-Femenia:2003jdx; see the source digest
](../../extracts/Ruiz-Femenia:2003jdx/pseudoscalar-two-photon-form-factor.md).
Constructing observations therefore still requires the full tree-level
off-shell amplitude and its WZW/resonance normalization, including the relevant
VVP exchange, to be derived consistently with the fitted coupling convention.
