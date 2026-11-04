---
title: Shared even-sector and conditional G_V constraint audit
status: numerical
created: 2026-11-04
updated: 2026-11-04
bibtex_keys: [Ecker:1988te, Ecker:1989yg, ParticleDataGroup:2024cfk, Alda:2026nrd, BaBar:2007ceh, Qin:2024ulb, Ruiz-Femenia:2003jdx]
provenance: numerical-computation
---

# Shared even-sector and conditional $G_V$ constraint audit

## Scope and inputs

This record audits the current stored global vector-sector likelihood and its conditional $G_V$ decay constraints. The baseline is not a converged fit: it has `success=false`, stopped at 1500 function evaluations, and reports $\chi^2=307.9848522$ for 43 degrees of freedom. All comparisons below are diagnostics against that stored parameter point; none establishes a converged or global minimum. The block-level baseline is also catalogued in the [joint-fit baseline record](babar-kkpi-global-vector-sector-fit-baseline.md).

The executable calculation is [`scripts/audit_shared_even_gv.py`](../../../scripts/audit_shared_even_gv.py). It reads the stored joint result, its input data, the conditional $G_V$ measurements, and the source-backed [PDG decay-input extract](../../extracts/ParticleDataGroup:2024cfk/vector-to-two-pseudoscalar-inputs.md). Its machine-readable result is [`data/vector_sector/shared_even_gv_audit.json`](../../../data/vector_sector/shared_even_gv_audit.json). Run from the repository root:

```text
.venv/bin/python scripts/audit_shared_even_gv.py --output data/vector_sector/shared_even_gv_audit.json
```

The six input hashes recorded by the joint-fit baseline match their current files. The output records hashes for the audit inputs and scripts. The charged/neutral $K^*$, $\phi$, and $\omega$ listing PDFs have been restored to the ignored local source archive and all three match the SHA-256 values in the data metadata. The $K^*$ widths and saturation fractions, plus the $\phi$ branching fractions and their $-0.72$ correlation, were checked against the downloaded listing text. The source files are indexed under [`memory/sources/ParticleDataGroup:2024cfk/`](../../sources/ParticleDataGroup:2024cfk/INDEX.md), and their values are summarized in the [PDG extract](../../extracts/ParticleDataGroup:2024cfk/vector-to-two-pseudoscalar-inputs.md) [@ParticleDataGroup:2024cfk]. The relevant bibliography entries are already present in `bibliography.bib`.

The conditional rows are the charged and neutral $K^*\to K\pi$ partial widths, $51.4\pm0.8$ MeV and $47.3\pm0.5$ MeV, and the fractions $\mathcal B(\phi\to K^+K^-)=(49.1\pm0.5)\%$ and $\mathcal B(\phi\to K_LK_S)=(33.9\pm0.4)\%$ [@ParticleDataGroup:2024cfk]. The two $\phi$ fractions have correlation coefficient $-0.72$; the quoted $K^*$ rows are treated as uncorrelated because their cross-covariance is unavailable. Other missing cross-correlations are set to zero as an explicit summary-covariance approximation.

## Mass, width, and overlap conventions

The PDG physical charged and neutral $K^*$ masses are $0.89167$ GeV and $0.89555$ GeV. Physical vector and daughter masses enter the external two-body phase space. The baseline vertex normalization uses the fitted common tensor mass $m_{K^*}$; the internal $K^*$ propagator uses the charge-mean shared width. The charged and neutral modes are summed separately with the published neutral-pion generator. The $\phi$ branching predictions use physical kaon masses and the fitted shared $\phi$ total-width nuisance. These conventions follow the existing [conditional $G_V$ calculation](../foundations/vector-gv-fit.md), the global fit code, and the input mass data [@Ecker:1988te; @Ecker:1989yg; @Alda:2026nrd; @ParticleDataGroup:2024cfk].

The global fit has seven shared-width rows. Each charged and neutral $K^*$ width appears once in the `shared_widths` block and is predicted from $G_V$ by imposing $K^*\to K\pi$ saturation. The two $\phi\to K\bar K$ fractions appear once in `phi_KK`. The separate four-row conditional fit reuses these same measurements; its $\chi^2$ is therefore not added to the joint likelihood. The conditional fit itself gives $|G_V|=0.05925138$ GeV and $\chi^2=95.99208$ for three nominal degrees of freedom with the provisional masses, mixing angle, and $\phi$ width fixed. This is a supplemental conditional diagnostic, not independent evidence to combine with the joint result.

The PDG extract gives $K\pi$ fractions of $99.902\pm0.009\%$ for charged $K^*$ and $99.754\pm0.021\%$ for neutral $K^*$. Their central deficits from unit saturation are $0.098\%$ and $0.246\%$. These fractions and the width summaries arise from overlapping constrained fits. Because their cross-covariance is unavailable, a non-saturated likelihood alternative is not evaluated and the fractions are not added as independent observations.

## Block and subset results

The stored block contributions are reproduced by recomputing all 75 joint residuals; the largest residual-level difference is $7.5\times10^{-14}$.

| Joint-fit block | Rows | $\chi^2$ |
| --- | ---: | ---: |
| Odd low-energy summaries | 12 | 11.15417 |
| Masses and leptonic widths | 7 | 10.02971 |
| Shared total widths | 7 | 102.45464 |
| $\phi\to K\bar K$ fractions | 2 | 16.73735 |
| Inclusive BaBar $KK\pi$ bins | 45 | 156.63438 |
| Normalization constraints | 2 | 10.97459 |

The selected 12-row even/$G_V$ diagnostic contains the seven mass/leptonic rows, the charged and neutral $K^*$ width rows, the $\phi$ total-width row, and the two $\phi$ branching rows. It refits the six even-sector coordinates, $G_V$, and the $\phi$ total-width pull for diagnosis only. Re-evaluation gives $\chi^2=125.02918430$, compared with $125.02918430$ in the stored result; the largest individual residual difference is $1.61\times10^{-5}$ from repeating this small conditional minimization. This subset is not another contribution to the joint total.

For a two-row model $y_i=\lambda a_i$ with covariance $C$, the common scale is profiled as

$$
\widehat\lambda=\max\!\left(0,\frac{a^T C^{-1}y}{a^T C^{-1}a}\right),\qquad
\chi^2_{\min}=(\widehat\lambda a-y)^T C^{-1}(\widehat\lambda a-y).
$$

Using the documented row covariances gives a charged/neutral $K^*$ pair minimum $\chi^2=22.21066$ for one nominal degree of freedom. The predicted neutral-to-charged coefficient ratio at the stored shared-width point is $1.00718$, while the observed width ratio is $0.92023$. The correlated $\phi$ pair minimum is $\chi^2=4.92390$ for one nominal degree of freedom. These are separate one-scale diagnostics, not independent likelihood terms.

## Prescription sensitivity and retained assumptions

As a limited external-width sensitivity, replace the common fitted tensor mass in the two external $K^*\to K\pi$ width normalizations by the respective PDG physical $K^*$ masses, while retaining physical phase space and the charge-mean internal propagator prescription. The profiled $K^*$ pair minimum changes from $22.21066$ to $18.21070$ ($\Delta\chi^2=-3.99997$). Propagating the altered shared widths through the stored parameter point, with every fit coordinate fixed, changes the full residual sum from $307.98485$ to $320.37901$ ($\Delta\chi^2=+12.39416$): the shared-width block changes by $+10.96925$, the BaBar block by $+1.59710$, and the odd block by $-0.17219$. This is a fixed-candidate sensitivity only. It is neither a refit nor a complete isospin-breaking vector model, and its full-sum difference must not be interpreted as a change between fitted minima.

The retained baseline assumptions are the negative single-trace source, near-ideal inverse-mass branch, $K^*\to K\pi$ saturation, one count for each duplicated constrained $K^*$ width, diagonal PDG summary errors apart from the stated $\phi$ correlation, physical phase-space masses, fitted tensor vertex masses, and the charge-mean internal $K^*$ width. No uncertainty from frozen provisional parameters or unavailable cross-covariances is propagated. The nominal goodness-of-fit is not calibrated for those omissions.

## Verification and provenance

The script checks all six input hashes stored in the global fit JSON before recomputing its residuals. Two runs to separate output files produced byte-identical JSON. The SHA-256 of the recorded audit result, including the now-verified source archive hashes, is `7ea5a4a581f08297189848bd57883eb229644bad5af4bbd4cf0bfc2b039e103a`; the SHA-256 of `scripts/audit_shared_even_gv.py` is recorded in the result JSON. The analytic pair-fit checks pass through:

```text
.venv/bin/python scripts/audit_shared_even_gv.py --self-test
```

The fit formulas use the leading even-parity interaction $\mathcal L_{VPP}=iG_V\langle V_{\mu\nu}u^\mu u^\nu\rangle/\sqrt2$ [@Ecker:1988te]. The stored global likelihood and BaBar data scope are documented in the [baseline record](babar-kkpi-global-vector-sector-fit-baseline.md) [@BaBar:2007ceh; @Qin:2024ulb; @Ruiz-Femenia:2003jdx].
