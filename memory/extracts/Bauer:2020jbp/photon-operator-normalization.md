---
title: Source-backed low-energy ALP photon operator normalization
status: literature
created: 2026-10-15
updated: 2026-10-15
bibtex_keys: [Bauer:2020jbp]
provenance: literature
---

# Evidence and coefficient convention

In [@Bauer:2020jbp, Eq. (39)], the weak-scale ALP Lagrangian includes

$$
\mathcal L\supset c_{\gamma\gamma}\frac{\alpha}{4\pi}
\frac{a}{f}F_{\mu\nu}\widetilde F^{\mu\nu}.
$$

Eq. (54) uses the same positive coefficient normalization below the weak
scale. Immediately after Eq. (1), the source defines
$\widetilde F^{\mu\nu}=\epsilon^{\mu\nu\rho\sigma}F_{\rho\sigma}/2$
with $\epsilon^{0123}=+1$. These are directly evidenced by the archived
[source index](../../sources/Bauer:2020jbp/INDEX.md), `paper_v2.tex:541--558`
and `707--714`, with the dual at line 95.

The repository adopts this coefficient convention, writing $c_\gamma$ for
the local photon Wilson coefficient and $f_a$ for the suppression scale.
Thus $c_\gamma/f_a$ denotes the coefficient ratio corresponding to
$c_{\gamma\gamma}/f$ at the selected scale and in the selected quark basis.
Here the project's $f_a$ labels Bauer's EFT suppression scale $f$; it is
not automatically Bauer's separately defined QCD-axion decay constant.
This is a notation identification; no numerical UV coefficient or decay
constant is inferred from it. The source explicitly distinguishes coefficients
above and below weak-scale matching thresholds.

# Scope

The factor $+\alpha_{\rm em}/(4\pi)$ and tensor orientation are literature
definitions, not an additional photon-normalization assumption. Eq. (39)
alone does not derive the repository's quark-rotation Jacobian, leading-nonet
WZW meson contribution or the invariant $C_{a\gamma\gamma}$. Those remain
separate results and declared inputs in the
[photon invariant record](../../alp/photon/alp-photon-invariant-from-mixing.md).
Independent local/anomaly matching for the enlarged RChT action is a further
physical task; it does not make the photon Wilson-coefficient normalization
uncertain.
