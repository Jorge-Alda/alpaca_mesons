---
title: Isospin-breaking VPP vertex as a switchable option of the joint vector likelihood
status: mixed
created: 2026-11-09
updated: 2026-11-09
bibtex_keys: [Alda:2026nrd, ParticleDataGroup:2024cfk, BaBar:2007ceh]
provenance: analytical-derivation
---

# Isospin-breaking $VPP$ vertex as a likelihood option

This record describes how the three isospin-breaking $VPP$ operators with coefficients $r_1,r_2,r_3$ enter the joint vector likelihood as a switchable option, so that the full likelihood can be evaluated at points of the $(r_1,r_2,r_3)$ scan. **No fit or refit has been run with the option.** The baseline likelihood and its earlier variations are unchanged.

Script: [`scripts/joint_vector_isospin_breaking_vpp_option.py`](../../../scripts/joint_vector_isospin_breaking_vpp_option.py) (class `IsospinBreakingVppLikelihood`, a subclass of `VariedLikelihood`); tests in [`scripts/test_joint_vector_isospin_breaking_vpp_option.py`](../../../scripts/test_joint_vector_isospin_breaking_vpp_option.py).

```bash
.venv/bin/python scripts/joint_vector_isospin_breaking_vpp_option.py
.venv/bin/python -m unittest scripts.test_joint_vector_isospin_breaking_vpp_option
```

Minimal use at a scan point (evaluation only):

```python
model = IsospinBreakingVppLikelihood(2, 16, 16, bin_block_size=4, isospin_breaking_vpp=True)
residual = model.residual(model.extend(baseline_values, (-0.350, 1.855, 0.729)))
```

## Conventions and units

GeV for masses, widths, energies, $G_V$ and $F=0.093$ GeV; GeV$^2$ for invariants and for $X=\chi_+|_{\Phi=0}=4B_0\,\mathrm{diag}(m_u,m_d,m_s)$; GeV$^{-2}$ for $r_i=\sqrt2\lambda^V_{\chi,i}/G_V$; residuals dimensionless. Conventions of Alda:2026nrd and of the baseline likelihood [@Alda:2026nrd]. The notation $r_i$ and $\lambda^V_{\chi,i}$ is the one proposed (awaiting approval) in [joint-vector-isospin-breaking-vpp-operators.md](joint-vector-isospin-breaking-vpp-operators.md).

## Inputs

- Effective flavor trace of each mode,
$$T(r)=\mathrm{tr}(V[P_1,P_2])+r_1\,\mathrm{tr}(V\{[P_1,P_2],X\})+r_2\,\mathrm{tr}\big(V(P_1XP_2-P_2XP_1)\big)+r_3\,\mathrm{tr}(V)\,\mathrm{tr}([P_1,P_2]X),$$
  with the flavor matrices of `fit_vector_gv.flavor_factors`: [joint-vector-isospin-breaking-vpp-operators.md](joint-vector-isospin-breaking-vpp-operators.md) (Step 2) and `scripts/derive_isospin_breaking_vpp_operators.py`.
- Per-mode factors $T(r)/T(0)$ and the scan point $r=(-0.350,\,1.855,\,0.729)$ GeV$^{-2}$: `scripts/scan_isospin_breaking_vpp_simultaneous.py` (no memory record yet).
- Substitution rule $G_V\to\sqrt{\kappa_q}G_V$ on the BaBar terms with $k_t,k_u$, and the term-by-term charge assignment: [joint-vector-kstar-kpi-saturation-freedom.md](joint-vector-kstar-kpi-saturation-freedom.md) (Step 3).
- $X=(0.026272,\,0.046960,\,0.943684)$ GeV$^2$ and the published $\pi^0$ diagonals: `data/vector_sector/alda_2026nrd_arxiv_v1.json` ([alda-neutral-amplitude-inputs.md](alda-neutral-amplitude-inputs.md)) [@Alda:2026nrd].
- Baseline likelihood rows and BaBar amplitude: [babar-kkpi-global-vector-sector-fit-baseline.md](babar-kkpi-global-vector-sector-fit-baseline.md), [joint-vector-likelihood-variation-options.md](joint-vector-likelihood-variation-options.md); PDG rows [@ParticleDataGroup:2024cfk]; BaBar $KK\pi$ data [@BaBar:2007ceh].

## Assumptions and choices (labelled)

1. **Vertex only** (inherited from the operator record): the $r_i$ change only $VPP$ couplings. $K^*$ masses and running-width shapes are unchanged; only the on-shell widths are rescaled.
2. **One factor per decay mode.** Since $r_2$ gives different factors to the two $K\pi$ modes of one $K^*$ charge state, the per-charge $\sqrt{\kappa_q}$ of the saturation record is replaced by the factor of the $K^*\to K\pi$ mode at each BaBar vertex.
3. **$\phi$ flavor matrix.** The $\phi$ factors use the fitted $\phi$ generator `sector.generators[2]`, the same matrix whose LO trace the baseline $\phi\to K\bar K$ rows use (choice). The scan script uses the fixed provisional mixing angle $0.0824$; at the accepted endpoint the fitted angle is $0.0952$.
4. **Excited $K^*$ terms.** Option `vpp_kstar_excitation_scaling`: `'ground_only'` (default) multiplies only the $K^*(892)$ pole by the mode factor; `'whole_propagator'` multiplies the whole propagator including the $\beta$-weighted $K^*(1410)$, $K^*(1680)$ terms. This is the open choice of the saturation record; the default is an assumption.
5. **Charge-mean propagator width** of the baseline is kept (isospin-limit assumption of the baseline), now built from the modified $K^*$ total widths.
6. **Sign.** The factors are signed: the BaBar terms use $T(r)/T(0)$, the widths $(T(r)/T(0))^2$.

## Method

1. Option `isospin_breaking_vpp=True` appends the coordinates `isospin_breaking_vpp_r1_GeV_minus2`, `..._r2_...`, `..._r3_...` after the 32 baseline coordinates. Start $(0,0,0)$; bounds $[-5,5]$ GeV$^{-2}$ (the scan range); outside the bounds `unpack` raises and the residual becomes the baseline penalty vector, as for the baseline pole-residue domain.
2. $T(r)/T(0)$ is evaluated with numpy for the four $K^*$ modes, $\phi\to K^+K^-$, $\phi\to K^0\bar K^0$, $\rho^0\to\pi^+\pi^-$ and $\rho^+\to\pi^+\pi^0$.
3. Rows where the factors enter (implemented):
   - $K^{*\pm}$ and $K^{*0}$ total-width rows: each mode of $\Gamma=G_V^2\sum$`width_coefficient` is multiplied by its squared factor ($K^{*\pm}$: $K^+\pi^0$, $K^0\pi^+$; $K^{*0}$: $K^0\pi^0$, $K^+\pi^-$).
   - $K^{*0}\to K^0\gamma$ and $K^{*\pm}\to K^\pm\gamma$ branching rows: through their total-width denominators. The radiative amplitudes contain no $VPP$ vertex.
   - BaBar $K^*$ propagator widths: through the charge-mean width.
   - BaBar $K^*$-channel terms (`direct` and `exchange` times $k_t$ or $k_u$): the $K^*$ propagator is multiplied by the factor of its vertex mode, which equals $G_V\to G_V\,T(r)/T(0)$ on exactly these terms (saturation record, Step 3). $K^+K^-\pi^0$, $t$ and $u$: $K^{*\pm}\to K^\pm\pi^0$. $K_SK^\pm\pi^\mp$, $t$ (neutral $K^*$): $K^{*0}\to K^+\pi^-$. $K_SK^\pm\pi^\mp$, $u$ (charged $K^*$): $K^{*\pm}\to K^0\pi^\pm$. $C$-conjugate modes have equal factors (tested).
   - $\phi\to K^+K^-$ and $\phi\to K_LK_S$ branching rows: predictions multiplied by the squared $\phi$ factors; the $\phi$ total width is a free pull in the baseline and is not changed.
4. Jacobian: inherited columns unchanged; the three added columns by forward differences with step $\sqrt{\epsilon}\max(1,|r_i|)$, as in the $K^*$ charge option.
5. Degrees of freedom: $75$ residuals ($73$ measurements, $2$ normalization constraints) minus $35$ coordinates $=40$ (baseline $43$).

## Not implemented (left at the baseline)

- **BaBar $s$-channel terms** with `rs`, `lls` ($\rho,\omega,\phi\to K\bar K$ vertices, including the charged $\rho\to K\bar K$ vertex in $K_SK\pi$). The baseline writes these vertices in the light/strange basis of the rotated neutral propagator; $r_1,r_2$ act differently on its light and strange parts and $r_3$ adds a $\mathrm{tr}(V)$ term. The operator record gives factors only for the physical $\rho^0,\omega\to K\bar K$ modes at ideal mixing and none for the charged $\rho$, so the substitution is not unambiguous. The operator record notes that these factors are large ($\rho,\omega\to K\bar K$ relative to $K^*$ up to $3.85$ for $r_2$ alone), so this omission can matter for the BaBar rows.
- Running-width shapes, $K^*$ masses, $\rho$, $\omega$, $\eta'$ widths (free pulls), mass and leptonic-width rows: no $VPP$ vertex or not affected by assumption 1.
- No likelihood row measures $\rho^\pm/\rho^0$ couplings; the $\rho$ factors are computed only for the isospin-limit check.

## Results at the scan point (numerical; evaluation of individual rows, not a fit)

At $r=(-0.350,1.855,0.729)$ GeV$^{-2}$, accepted endpoint `basin_probes/even_shared_widths`, orders $(1,3,4)$:

| Quantity | Value |
|---|---|
| $T(r)/T(0)$: $K^+\pi^0$, $K^0\pi^+$, $K^0\pi^0$, $K^+\pi^-$ | $0.72811,\ 0.74763,\ 0.72086,\ 0.70201$ |
| $\Gamma_{K^{*\pm}}(r)/\Gamma_{K^{*\pm}}(0)$, $\Gamma_{K^{*0}}(r)/\Gamma_{K^{*0}}(0)$ | $0.54872,\ 0.50135$ (ratio $1.094483$, the scan target) |
| $T(r)/T(0)$ for $\phi\to K^+K^-$, $K^0\bar K^0$ (fitted angle) | $1.03249,\ 1.05634$ |
| change of the $\phi$ rate ratio $K^+K^-/K^0\bar K^0$ | $-4.46\%$ (scan, fixed angle: $-4.43\%$) |

Both $K^*$ widths drop by about half at this point: with $G_V$ held fixed, the $K^*$ width rows move far from the PDG values unless $G_V$ readjusts, which would in turn change the $\phi$ and BaBar $s$-channel terms. This is a statement about the residual rows at fixed baseline coordinates, not about the profiled likelihood.

## Verification (run 2026-11-09)

`scripts.test_joint_vector_isospin_breaking_vpp_option`: 18 tests pass.
- Option off: residual and Jacobian bitwise equal to `GlobalLikelihood(2,16,16,bin_block_size=4)` at the accepted endpoint; option on at $r=0$ (both excitation scalings): residual and first 32 Jacobian columns bitwise equal (orders $(1,3,4)$).
- Factors equal the scan script's `trace_coefficients` to $10^{-12}$ ($K^*$, $\rho$, and $\phi$ at the fixed angle); $\phi$ factors at the fitted angle equal the SymPy traces of `derive_isospin_breaking_vpp_operators.effective_traces`.
- At the scan point: $K^*$ width ratios equal the scan's phase-space-weighted squared traces; $\Gamma_{K^{*\pm}}/\Gamma_{K^{*0}}$ rises by the scan's $R$; $K^*\to K\gamma$ predictions scale by the inverse width ratio; $\phi$ rows move in the scan's direction; the BaBar rows equal `babar_amplitudes` of `derive_kstar_kpi_saturation_freedom.py` with per-channel $\kappa$ set to the squared mode factors, for both excitation scalings ($10^{-11}$); only the expected rows change; $r_3$ alone changes only the two $\phi$ rows.
- $\rho^+/\rho^0$ factor ratio equal to 1 to $10^{-14}$ with the isospin-limit $\pi^0$ at three points; with the published $\pi^0$ it equals the scan's value ($|\,\cdot-1|<2\times10^{-3}$).
- Jacobian at the scan point finite, shape $(75,35)$; the $r_3$ column vanishes outside the $\phi$ rows.
- Invalid options, wrong coordinate counts and $|r_i|>5$ raise (or give the penalty residual).

Also rerun and passing: `test_joint_vector_likelihood_options`, `test_joint_vector_kstar_charge_options`, `test_global_kkpi_running_propagator`, `test_derive_kstar_kpi_saturation_freedom`, `test_derive_isospin_breaking_vpp_operators`, `test_scan_isospin_breaking_vpp_simultaneous` (51 tests). Not checked: any fit or refit with the option, high quadrature orders at nonzero $r$, the BaBar $s$-channel factors.
