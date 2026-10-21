---
title: PDG 2024 vector-to-two-pseudoscalar inputs for the conditional G_V fit
status: literature
created: 2026-10-21
updated: 2026-10-21
bibtex_keys: [ParticleDataGroup:2024cfk]
provenance: literature
---

# PDG vector-to-two-pseudoscalar inputs

The PDG 2024 charged and neutral $K^*$ listings give the partial widths

$$
\Gamma(K^{*+}\to K\pi)=51.4\pm0.8\;\mathrm{MeV},\qquad
\Gamma(K^{*0}\to K\pi)=47.3\pm0.5\;\mathrm{MeV}.
$$

The charged value is the page-5 partial-width entry, ``OUR FIT''; the neutral
value is the page-6 partial-width entry, ``OUR FIT'', whose quoted error already
includes scale factor $S=2.0$. These arise from PDG constrained fits that use
total and partial widths, and are not independent of those total-width inputs
[@ParticleDataGroup:2024cfk]. The corresponding listings give $K\pi$ branching
fractions $99.902\pm0.009\%$ and $99.754\pm0.021\%$, respectively. They are
not included together with the partial-width entries in the repository fit,
to avoid counting the same constrained decay-width information twice.

The PDG $\phi(1020)$ listing gives

$$
\mathcal B(\phi\to K^+K^-)=(49.1\pm0.5)\%,\qquad
\mathcal B(\phi\to K_LK_S)=(33.9\pm0.4)\%,
$$

with scale factors $1.3$ and $1.2$ stated in the listing. Its constrained-fit
correlation table gives correlation coefficient $-0.72$ for these two
branching fractions. The percentages and quoted errors are retained as listed;
no extra multiplication by the displayed scale factors is made. They are fit
fractions subject to the PDG branching-fraction sum constraint. In a
conditional $G_V$ extraction, convert predicted partial widths to branching
fractions using the already fixed provisional-fit $φ$ total-width nuisance;
do not refit that width or count its constraint again.

The same listings also quote

$$
\mathcal B(\phi\to\pi^+\pi^-)=(7.3\pm1.3)\times10^{-5},\qquad
\mathcal B(\omega\to\pi^+\pi^-)=(1.53\pm0.12)\%.
$$

The $\phi$ value is listed on page 3 and the $\omega$ value on page 2. The
$\omega$ listing reports scale factor $S=1.2$; the quoted error is used as
printed without another rescaling. These are nonzero PDG decay summaries, but
they are not inputs to the current leading exact-isospin $G_V$ fit. Their
model status is explained in
[`vector-gv-fit.md`](../../../rcht/foundations/vector-gv-fit.md).

The listing assigns approximately $100\%$ of the $\rho$ width to $\pi\pi$ but
does not quote an independent uncertainty for that branching fraction. The
rho total-width constraint already used in the provisional fit therefore is
not an additional independent $\rho\to\pi\pi$ datum for the supplemental fit.

Source pages are the PDG 2024 charged $K^*$ listing pages 5 and 8, neutral
$K^*$ listing pages 6 and 8, and $\phi$ listing pages 3 and 4. Local text and
PDF files are archived under
`memory/sources/ParticleDataGroup:2024cfk/`; the corresponding source URLs and
checksums are indexed in
[`vector-sector-calibration-inputs.md`](vector-sector-calibration-inputs.md)
[@ParticleDataGroup:2024cfk].
