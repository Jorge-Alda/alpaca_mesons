---
title: Canonical Levi-Civita tensor and opposite-convention source conversion
status: mixed
created: 2026-10-15
updated: 2026-10-15
bibtex_keys: [Bauer:2020jbp, Ruiz-Femenia:2003jdx, Alda:2026nrd]
provenance: analytical-derivation
---

# Inputs and assumptions

The user selects $\epsilon^{0123}=+1$ as in [@Bauer:2020jbp] immediately
after Eq. (1). The archived [source](../../sources/Bauer:2020jbp/INDEX.md)
has the dual definition at `paper_v2.tex:95`. With the repository metric
$(+---)$ there is one metric-related tensor, hence $\epsilon_{0123}=-1$.
This fixes a tensor convention, not additional anomaly matching.

The original Ruiz action uses the opposite lower-index orientation
[@Ruiz-Femenia:2003jdx, source `vvp.tex:352--354`]. Its
[operator digest](../../extracts/Ruiz-Femenia:2003jdx/vector-odd-parity.md)
and [two-photon digest](../../extracts/Ruiz-Femenia:2003jdx/pseudoscalar-two-photon-form-factor.md)
are now adapted; original archived source text remains literal literature.
The positive photon-dual normalization is the source-backed definition in
[@Bauer:2020jbp, Eqs. (39)/(54)], recorded in the
[photon-operator digest](../../extracts/Bauer:2020jbp/photon-operator-normalization.md).
The coefficient rotation and leading-nonet WZW meson contribution are separate
inputs/results in the [photon record](../../alp/photon/alp-photon-invariant-from-mixing.md).
Alda's independent axial/mass/topological source definitions are unchanged;
its WZW-excluding meson action does not supply the photon operator or WZW meson contribution.

# Executed conversion

Run the deterministic SymPy derivation and meaningful tensor tests:

```text
.venv/bin/python scripts/derive_levi_civita_convention.py
.venv/bin/python -m unittest scripts.implementation.test_levi_civita_convention
```

The executable constructs the metric, lowers the all-upper tensor, contracts
basis vectors, substitutes the coefficient dictionary in a generic action
and signed anomaly-plus-resonance amplitude, and evaluates both photon Wick
assignments. Its [transcript](../../../data/vector_sector/levi_civita_convention.md)
and [exact artifact](../../../data/vector_sector/levi_civita_convention.json)
derive the action-preserving dictionary, with $R$ denoting reference values:

$$
\epsilon_{\mu\nu\rho\sigma}=-\epsilon^R_{\mu\nu\rho\sigma},\qquad
c_i=-c_i^R,\qquad d_i=-d_i^R,\qquad g_j=-g_j^R.
$$

The same rule applies to every freely normalized single-epsilon odd parent.
Even-sector coefficients, masses, source weights, Fourier convention and
external-state phases stay fixed. The explicit WZW derivative coefficient
reverses; both parts of each inhomogeneous OPE right-hand side reverse.
Homogeneous relations and rank/null dimensions are unchanged.

For a signed scalar multiplying a fixed canonical epsilon contraction,

$$
F(c,d)=-F^R(-c,-d).
$$

Thus the TFF anomaly constant becomes negative while the linear resonance
slopes retain their signs. This preserves the full tensor amplitude and its
anomaly--resonance interference, rather than independently flipping its
anomaly. A positive-pion scalar convention would require a separate meson
rephasing and the corresponding ALP LSZ rephasing; it is not used here.

With $E=\epsilon(q_1,\varepsilon_1,q_2,\varepsilon_2)$, the executed photon
ratios are $\langle F\widetilde F\rangle/E=-4$ and
$\langle\epsilon\partial A\partial A\rangle/E=-2$. Therefore

$$
\mathcal L_{a\gamma\gamma}
=\frac{e^2}{16\pi^2f_a}C_{a\gamma\gamma}aF\widetilde F
=\frac{e^2}{8\pi^2f_a}C_{a\gamma\gamma}a\epsilon\partial A\partial A,
\qquad G_{\rm contact}=-\frac{C_{a\gamma\gamma}}{4\pi^2f_a}.
$$

The positive coefficient $C_{a\gamma\gamma}$ and its source invariant are
unchanged. The old photon-phase transcript had a second, oppositely
oriented source tensor; its $\eta_\gamma=-1$ is now historical conversion
information, not a live convention. The manuscript states only canonical
physics; all reference-conversion clarifications reside in memory.

# Saved outputs and reproducibility boundary

Existing derivation solvers and `data/vector_sector/` caches produced before
this conversion remain **reference-backend outputs** unless an output
explicitly records this canonical convention. They were not all regenerated.
Their fitted odd coefficients and covariance samples are $g^R,c^R,d^R$;
convert their signs before combining them with canonical operators. A saved
signed scalar must also reverse when its epsilon contraction reverses.
For a linear gradient expressed per canonical odd coefficient, the two
reversals cancel, so its scalar matrix is unchanged. Contact/anomaly constants
must still reverse in that scalar basis. Widths, ranks, objective values and
coefficient-covariance matrices alone cannot certify signed interference.
Stored signed target/profile values therefore require the adapter above;
they are not relabeled as newly regenerated canonical predictors.

The shared `levi_civita_contract` helper now contracts lowered arguments with
the canonical upper tensor. `legacy_epsilon_determinant` explicitly retains
the old contravariant determinant for the historical kinematics probe.
The probe and all three convention tests pass. The conversion artifact
checks action and interference-sensitive total squared-amplitude invariance.

# Planning

This completes the requested manuscript/memory convention adaptation. The
next predictor step is to use the canonical coefficient/tensor adapter at
every legacy-backend boundary, or regenerate an explicitly canonical output
when needed. No previous calibration/profile cache has been silently
regenerated. The broader completion gates remain in [next-steps.md](next-steps.md).
