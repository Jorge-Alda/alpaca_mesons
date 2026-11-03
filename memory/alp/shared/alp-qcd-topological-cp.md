---
title: QCD topological ALP couplings under chiral redefinitions and CP
status: derived
created: 2026-11-03
updated: 2026-11-03
bibtex_keys: [Alda:2026nrd, Bonnefoy:2022rik]
provenance: analytical-derivation
---

# QCD topological ALP couplings under chiral redefinitions and CP

## Scope and conventions

Use the partonic QCD truncation and $x=a/f_a$ from the
[general chiral-redefinition record](alp-general-chiral-redefinition-invariants.md).
The ALP is CP odd. Normalize the dimensionless topological coefficient
$\Theta(x)$ so its term linear in $x$ is Alda:2026nrd's $c_G$ convention
[@Alda:2026nrd]. CP statements below assume a CP-adapted field basis and
quark masses whose phases give a CP-conserving physical vacuum angle.

## Effect of a chiral field redefinition

For $q_{R,L}'=g_{R,L}(x)q_{R,L}$, the anomalous Jacobian shifts the
topological coefficient by the singlet axial angle. In Alda:2026nrd's
normalization,
$$
\Theta'(x)=\Theta(x)-\frac{1}{2i}
\log\!\left(\frac{\det g_R(x)}{\det g_L(x)}\right),
$$
with the logarithm chosen continuously from $x=0$. This depends only on the
singlet axial part; traceless flavor rotations have no QCD topological shift.
For $g_R=e^{ix(\delta+\kappa)}$ and
$g_L=e^{ix(\delta-\kappa)}$, the determinant ratio gives
$$
\Theta'(x)=\Theta(x)-x\operatorname{tr}\kappa.
$$
Thus only the linear coefficient shifts for this linear-in-$x$ rotation:
$c_G'=c_G-\operatorname{tr}\kappa$. This is the shift that leaves
$c_{GA}=c_G+\operatorname{tr}c_q^A$ invariant when
$c_q^{A\prime}=c_q^A+\kappa$ [@Alda:2026nrd].

For a general nonlinear ALP-dependent redefinition, write the singlet axial
angle as a Taylor series in $x$. Its order-$n$ trace shifts $\Theta_n$ at
the same Taylor order. A nonlinear rotation can therefore generate higher
powers $x^nG\widetilde G$, but their coefficients are fixed by that field
redefinition; they are not new Jacobian structures independent of the
rotation. The anomaly produces the topological density $G\widetilde G$; it
does not produce the CP-even density $G^2$.

## CP selection of gluonic operators

Under CP, $x\to -x$, $G^A_{\mu\nu}\widetilde G^{A\mu\nu}\to
-G^A_{\mu\nu}\widetilde G^{A\mu\nu}$, while
$G^A_{\mu\nu}G^{A\mu\nu}$ is even. Therefore
$$
\operatorname{CP}\!\left[x^nG\widetilde G\right]=(-1)^{n+1},\qquad
\operatorname{CP}\!\left[x^nG^2\right]=(-1)^n.
$$
If the ALP interactions conserve CP, the coefficient functions must have
the corresponding parity:
$$
\Theta(-x)=-\Theta(x),\qquad \Phi(-x)=\Phi(x),
$$
where $\Theta(x)G\widetilde G$ is the topological term and $\Phi(x)G^2$ is
the scalar gluon coupling. Thus odd powers $xG\widetilde G$,
$x^3G\widetilde G$, $x^5G\widetilde G$, and so on are CP allowed; even
powers multiplying $G\widetilde G$ violate CP in this basis. Conversely,
$x^2G^2$, $x^4G^2$, and so on are CP allowed, while odd powers multiplying
$G^2$ violate CP.

The constant term in $\Theta$ combines with the phase of $\det M(0)$ to form
the physical QCD vacuum angle. CP requires that physical angle to be at a
CP-conserving value; $\Theta(0)$ by itself depends on the quark phase
convention.

CP alone permits these higher powers. Exact continuous ALP shift symmetry
restricts the topological coefficient to the linear shift coupling, up to
the usual topological total derivative; higher powers require explicit
shift breaking or spurions that transform to compensate it
[@Bonnefoy:2022rik]. The scalar gluon terms $x^{2m}G^2$ are independent
shift-breaking Wilson operators, not consequences of the axial anomaly.
They fit the current truncation because they contain no ALP derivative.

If a field redefinition does not preserve the simple CP assignment of the
fields, the transformed coefficient functions need not look even or odd by
the naive rule; CP is then realized with the transformed field convention.
The parity conditions above apply in a CP-adapted basis.

## Reproducible sign and parity bookkeeping

Run:

```sh
.venv/bin/python scripts/derive_alp_qcd_topological_cp.py --order 5
```

The script expands the anomaly shift using the stated Alda normalization and
checks the CP signs for arbitrary ALP power. It does not recalculate the
normalization of the QCD anomaly.
