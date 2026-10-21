---
title: Conditional four-row fit of the leading vector-pseudoscalar coupling
status: numerical
created: 2026-10-21
updated: 2026-10-21
bibtex_keys: [Ecker:1988te, Ecker:1989yg, ParticleDataGroup:2024cfk, Alda:2026nrd]
provenance: numerical-computation
---

# Conditional four-row $G_V$ fit from $K^*\to K\pi$ and $\phi\to K\bar K$

## Scope and inputs

The selected provisional odd-parity fit does not determine $G_V$. Its central
negative single-trace solution does fit $F_V=0.168952777554148$ GeV and
$z_V=-0.170269957284197$ GeV$^{-1}$ together with vector-sector inputs; these
are fitted parameters, not values fixed by chiral symmetry. Other theoretical
prescriptions can give different source couplings: for example, the Aloni
mass-mixing convention is matched to $F_V=f_\pi/\sqrt2$ in the separate
prescription record. This is not the value used by the provisional fit.

The leading even-parity RChT interaction is

$$
\mathcal L_{VPP}=\frac{iG_V}{\sqrt2}\langle V_{\mu\nu}u^\mu u^\nu\rangle,
$$

which gives the generic-flavor $VPP$ vertex in
[the interaction expansion](rcht-vpp-vpa-photon-pa-expansion.md)
[@Ecker:1988te]. The fit holds every parameter of
`central_single_trace_negative` in
[`odd_preliminary_fit.json`](../../../data/vector_sector/odd_preliminary_fit.json)
fixed, including the fitted tensor masses, neutral mixing inputs, and shared
total-width nuisances. It uses the charged and neutral $K^*$ partial widths
$51.4\pm0.8$ and $47.3\pm0.5$ MeV, plus the $\phi\to K^+K^-$ and
$\phi\to K_LK_S$ branching fractions $(49.1\pm0.5)\%$ and
$(33.9\pm0.4)\%$ [@ParticleDataGroup:2024cfk]. The PDG constrained-fit
correlation coefficient between the two $\phi$ fractions is $-0.72$ and is
included. The partial-width entries come from PDG constrained fits that also
use total widths. The $\phi$ branching fractions use the fixed provisional
$\phi$ total-width nuisance to compare predicted widths with measured
fractions. Therefore this supplemental four-row fit overlaps the archived
provisional inputs; its $\chi^2$ must not be added to the old fit as
independent data.

## Width formula and fit

For each charge mode, tensor-to-Proca conversion gives

$$
\frac{\Gamma(V\to P_1P_2)}{G_V^2}
=\left|\operatorname{Tr}\!\left(T_V[T_{P_1},T_{P_2}]\right)\right|^2
\frac{M_{\rm phys}^2p^3}{12\pi F^4M_{\rm tensor}^2},
\qquad
p=\frac{\sqrt{[M_{\rm phys}^2-(m_1+m_2)^2]
\,[M_{\rm phys}^2-(m_1-m_2)^2]}}{2M_{\rm phys}}.
$$

The charged $K^*$ row sums the $K^+\pi^0$ and $K^0\pi^+$ widths; the neutral
row sums $K^0\pi^0$ and $K^+\pi^-$. Both use the published neutral-pion
generator. The two $\phi$ rows use charged and neutral physical kaon masses;
negligible CP violation identifies the $K_LK_S$ width with the corresponding
$K^0\bar K^0$ width, with no identical-particle factor. Physical masses enter
phase space, while the fixed provisional tensor masses enter the vertex
normalization. The $\phi$ predictions are divided by the fitted total-width
nuisance, held fixed. The script minimizes the one-parameter Gaussian
likelihood exactly in $G_V^2\geq0$. The positive branch is adopted by
convention because widths alone determine only $|G_V|$.

Run:

```text
.venv/bin/python scripts/fit_vector_gv.py data/vector_sector/gv_measurements.json data/vector_sector/gv_conditional_fit.json
```

The four-row common-coupling fit gives $|G_V|=0.0592514$ GeV with a
conditional $\Delta\chi^2=1$ interval $[0.0591415,0.0593611]$ GeV and
$\chi^2=95.9921$ for three nominal degrees of freedom. The marginal pulls at
the common best fit are $-8.02$ (charged $K^*$), $-3.99$ (neutral $K^*$),
$+3.37$ ($\phi\to K^+K^-$), and $-0.98$ ($\phi\to K_LK_S$). Fitting each
row separately gives $|G_V|=63.338$, $60.543$, $58.261$, and $59.597$ MeV,
respectively. Under this tree-level model and the adopted summary covariance,
one universal coupling does not describe the selected rows. The nominal
goodness-of-fit does not account for fixed-input uncertainty, omitted model
uncertainty, or unavailable cross-correlations; it is not a precision
exclusion of the RChT operator.

The approximately $100\%$ $\rho\to\pi\pi$ branching fraction is not a further
independent width datum: the provisional fit already uses the rho total-width
constraint, and PDG supplies no uncertainty on that saturation fraction. It is
therefore excluded rather than counted again. The charged and neutral $K^*$
partial-width rows also overlap the PDG total-width information in the
provisional calibration. The fit is a supplemental conditional diagnostic,
not an independent update of its likelihood.

## Isospin-violating $\omega,\phi\to\pi^+\pi^-$ modes

PDG lists $\mathcal B(\phi\to\pi^+\pi^-)=(7.3\pm1.3)\times10^{-5}$ and
$\mathcal B(\omega\to\pi^+\pi^-)=(1.53\pm0.12)\%$
[@ParticleDataGroup:2024cfk]; their source details are in the
[PDG decay-input extract](../../extracts/ParticleDataGroup:2024cfk/vector-to-two-pseudoscalar-inputs.md).
They are deliberately not added as $G_V$ fit observations. In the exact-isospin
LO action used above, the charged-pion commutator is proportional to
$\operatorname{diag}(1,-1,0)$. Both ideal isoscalar generators,
$T_{\omega_I}=\operatorname{diag}(1,1,0)/\sqrt2$ and
$T_{\phi_I}=\operatorname{diag}(0,0,1)$, have zero trace against it. Any
rotation confined to the isoscalar $\omega$--$\phi$ subspace preserves that
zero, so the predicted $V\to\pi^+\pi^-$ amplitudes are exactly zero for every
$G_V$. The nonzero PDG rates require isospin-violating contributions, such as
$\rho$ admixtures or electromagnetic effects, absent from the frozen
provisional model. They cannot constrain its universal $G_V$ without adding
and fitting those effects.

## Verification

[`scripts/fit_vector_gv.py`](../../../scripts/fit_vector_gv.py) records the
input hashes, frozen parameter values, summed flavor traces, prediction, and
fit result in
[`gv_conditional_fit.json`](../../../data/vector_sector/gv_conditional_fit.json).
The accompanying unit checks verify the equal-mass $\rho\to\pi\pi$
normalization, flavor factors, threshold closure, correlated-fit stationarity,
and the fixed provisional input selection. The result is conditional on the
leading tree-level RChT vertex and the quoted-summary Gaussian likelihood.
