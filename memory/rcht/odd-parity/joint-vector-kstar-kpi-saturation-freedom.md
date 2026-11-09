---
title: Freeing the K pi saturation of the K-star couplings in the joint vector likelihood
status: derived
created: 2026-11-09
updated: 2026-11-09
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Freeing the $K\pi$ saturation of the $K^*(892)$ couplings

This record derives how to replace the single $G_V$ that fixes both the $K^*(892)$ widths and the $K^*K\pi$ vertex in the BaBar $KK\pi$ amplitude by one free coupling per $K^*$ charge state. No refit is done and the baseline likelihood is unchanged.

Script: [`scripts/derive_kstar_kpi_saturation_freedom.py`](../../../scripts/derive_kstar_kpi_saturation_freedom.py), tests in [`scripts/test_derive_kstar_kpi_saturation_freedom.py`](../../../scripts/test_derive_kstar_kpi_saturation_freedom.py).

```bash
.venv/bin/python scripts/derive_kstar_kpi_saturation_freedom.py
.venv/bin/python -m unittest scripts.test_derive_kstar_kpi_saturation_freedom
```

## Conventions and units

Natural units, GeV for masses, widths, $G_V$ and $F=0.093$ GeV; GeV$^2$ for the invariants $q^2,s,t,u$. Conventions of Alda:2026nrd and of the baseline likelihood [@Alda:2026nrd]. In $K^+K^-\pi^0$ both $t$ and $u$ are $K^\mp\pi^0$ invariants (charged $K^*$); in $K_SK^\pm\pi^\mp$ the integrated $t=(p_{K^\pm}+p_{\pi^\mp})^2$ is the neutral $K^*$ and $u$ the charged $K^*$ (assignment read off the phase-space construction, as in [`joint_vector_kstar_charge_options.py`](../../../scripts/joint_vector_kstar_charge_options.py) line 195).

## Proposed notation (awaiting user approval)

For $q\in\{\pm,0\}$, the **$K^{*}_q$ $K\pi$ saturation factor**
$$\kappa_{K^{*\pm}}=\left(\frac{g_{K^{*\pm}}}{g^{\rm sat}_{K^{*\pm}}}\right)^2,\qquad \kappa_{K^{*0}}=\left(\frac{g_{K^{*0}}}{g^{\rm sat}_{K^{*0}}}\right)^2,$$
where $g_{K^{*}_q}$ is the $K^*_q\to K\pi$ coupling used in the likelihood and $g^{\rm sat}_{K^*_q}$ the value implied by the common $G_V$ of the baseline. $\kappa=1$ is the baseline. The symbols and names are proposed, not approved.

## Inputs

- Two-body width of the leading VPP interaction and its flavor traces: [conditional $G_V$ record](../foundations/vector-gv-fit.md); code `scripts/fit_vector_gv.py` line 66 (`width_coefficient`) and `flavor_factors`.
- Baseline saturated widths: `scripts/fit_babar_kkpi_global_vector_sector.py` lines 245–256 (`totals[name]=gv**2*sum(width_coefficient(...))`), with physical parent masses $0.89167$, $0.89555$ GeV and the fitted tensor mass $m_{K^*}$ in the vertex; baseline description in [babar-kkpi-global-vector-sector-fit-baseline.md](babar-kkpi-global-vector-sector-fit-baseline.md).
- Where the widths enter the likelihood: total-width rows (line 276), $K^*\to K\gamma$ branching-fraction denominators (line 270), BaBar $K^*$ propagators through the charge-mean width (line 358) and the running-width rescaling (line 81).
- BaBar amplitude terms: lines 385–393 (`direct`, `exchange`, `kt`, `ku`).
- Charge-specific propagator option: [`joint_vector_kstar_charge_options.py`](../../../scripts/joint_vector_kstar_charge_options.py).
- Evaluation point: accepted candidate `basin_probes/even_shared_widths` from `data/vector_sector/joint_vector_endpoint_comparison.json` ([endpoint comparison](joint-vector-fit-endpoint-comparison.md)); reduced quadrature orders $(1,3,4)$ without energy-bin blocks.
- PDG 2024 $K^*$ total widths $0.0514$ and $0.0473$ GeV from the likelihood rows ([width audit](shared-even-sector-gv-constraint-audit.md)).

## Step 1: width and vertex scale together (SymPy, derived)

The baseline on-shell coupling is
$$g^{\rm sat}_{VPP}=\frac{G_V\,M_{\rm phys}^2\,\sqrt{T^2}}{\sqrt2\,F^2\,M_T},\qquad \Gamma=\frac{g^2\,p^3}{6\pi M_{\rm phys}^2},$$
with $T^2$ the squared flavor trace, $M_T$ the tensor mass and $p$ the daughter momentum. SymPy confirms $\Gamma^{\rm sat}/G_V^2=T^2M_{\rm phys}^2p^3/(12\pi F^4M_T^2)$, identical to `width_coefficient`. With $g_q=\sqrt{\kappa_q}\,g^{\rm sat}_q$,
$$\Gamma_{K^*_q}=\kappa_q\,\Gamma^{\rm sat}_{K^*_q}(G_V,m_{K^*},\text{masses})$$
exactly, for each $K\pi$ mode and hence for their sum. The equation $\Gamma(r\,g^{\rm sat})=\kappa\Gamma^{\rm sat}$ has the two roots $r=\pm\sqrt\kappa$: the width fixes only $|g_q|$. The parametrization by $\kappa_q>0$ keeps the baseline sign ($r=+\sqrt{\kappa_q}$); the sign-flipped branch is a separate discrete alternative and is **not** included (assumption, stated).

Exact squared traces: in the isospin limit $\pi^0=\mathrm{diag}(1,-1,0)/\sqrt2$ they are $1/2$ ($K^+\pi^0$, $K^0\pi^0$) and $1$ ($K^0\pi^+$, $K^+\pi^-$), summing to $3/2$ for both charges. With the published Alda:2026nrd $\pi^0$ diagonals used in the code, $K^+\pi^0$ gives $0.52888$ and $K^0\pi^0$ gives $0.47159$; the sums are $1.52888$ (charged) and $1.47159$ (neutral). These exact values match `kaon_flavors` to $10^{-15}$.

## Step 2: saturated widths at the accepted endpoint (numerical check)

At `basin_probes/even_shared_widths` ($G_V=0.05906781503$ GeV, $m_{K^*}=0.8908385838$ GeV) the SymPy formula gives $\Gamma^{\rm sat}_{K^{*\pm}}=0.0447846333258$ GeV and $\Gamma^{\rm sat}_{K^{*0}}=0.0451061788348$ GeV, equal to the baseline `totals` to $1.1\times10^{-16}$ relative. The factors that would reproduce the PDG widths are $\kappa_{K^{*\pm}}=1.147715$ and $\kappa_{K^{*0}}=1.048637$.

## Step 3: term-by-term charge assignment in the BaBar amplitude (derived, numerically checked)

`direct` $=2\sqrt2\,G_V/(3F^3M_V)$ and `exchange` $=2G_V/F^3$ each contain one $VPP$ vertex. In the terms with an $s$-channel $\rho,\omega,\phi$ propagator (`rs`, `lls`) this vertex is $V\to K\bar K$; in the terms with `kt` or `ku` it is $K^*\to K\pi$. The substitution is therefore:

| Channel | Term (code line) | $K^*$ charge | Factor |
|---|---|---|---|
| $K^+K^-\pi^0$ | $-$direct$\,[\mathrm{crk}(t)\,k_t+\mathrm{crk}(u)\,k_u]$ (388) | $K^{*\pm}$ | $\sqrt{\kappa_{K^{*\pm}}}$ |
| $K^+K^-\pi^0$ | $+$exchange$\,(\ldots)[\mathrm{dr}(t)\,k_t+\mathrm{dr}(u)\,k_u]$ (389) | $K^{*\pm}$ | $\sqrt{\kappa_{K^{*\pm}}}$ |
| $K_SK\pi$ | $-2\,$direct$\,\mathrm{cr}(t)\,k_t$; $+$exchange$\,(\ldots)\mathrm{dr}(t)\,k_t$ (392–393) | $K^{*0}$ | $\sqrt{\kappa_{K^{*0}}}$ |
| $K_SK\pi$ | $+$direct$\,\mathrm{crk}(u)\,k_u$; $-$exchange$\,(\ldots)\mathrm{dr}(u)\,k_u$ (392–393) | $K^{*\pm}$ | $\sqrt{\kappa_{K^{*\pm}}}$ |
| both | direct/exchange terms with `rs`, `lls` | none ($\rho,\omega,\phi\to K\bar K$) | 1 |
| both | contact, WZW constant, $a_2$ terms | none (no $G_V$) | 1 |

Since every $k_t$, $k_u$ appears only multiplied by `direct` or `exchange`, multiplying each $K^*$ propagator by $\sqrt{\kappa_q}$ of its charge state is the same as $G_V\to\sqrt{\kappa_q}\,G_V$ on exactly these terms. No single term contains both charge states. The only place where both enter together is the baseline charge-mean width $(\Gamma_{K^{*0}}+\Gamma_{K^{*\pm}})/2$ used in every $K^*$ propagator (line 358): with $\kappa_{K^{*\pm}}\neq\kappa_{K^{*0}}$ each propagator width depends on both factors. Two consistent treatments exist: keep the mean width (inherited isospin-limit assumption), or use the subclass option `separate_kstar_charge_states`, where each propagator carries the width of its own charge state.

Numerical checks (all run, reduced orders):
- $\kappa=(1,1)$ with the mean width: BaBar residuals bitwise equal to the baseline; forcing the code path through the ground/excited split at $\kappa=1+1$ ulp changes amplitudes by $1.9\times10^{-16}$ relative.
- $G_V=0$: amplitudes with the $K^*$ propagators set to zero are bitwise equal to the full ones, so no $G_V$-independent term contains a $K^*$.
- Widths held fixed, every $K^*$ propagator scaled by $1.37$ equals $G_V\to1.37\,G_V$ on the $K^*$ terms only (difference $3.3\times10^{-16}$).
- Widths held fixed, the amplitude is affine in $(\sqrt{\kappa_{K^{*\pm}}},\sqrt{\kappa_{K^{*0}}})$: cross term $1.7\times10^{-16}$, second difference $5.1\times10^{-16}$. $\kappa_{K^{*0}}$ has exactly zero effect on $K^+K^-\pi^0$; for $\sqrt\kappa$ changing $1\to1.5$, $K_SK\pi$ changes by $10.3\%$ (neutral) and $6.8\%$ (charged) of its largest amplitude.
- With charge-specific widths, $\kappa=1$ and $\Delta m=0$, the BaBar residuals agree exactly with `KstarChargeLikelihood(separate_kstar_charge_states=True, free_kstar_mass_difference=True)`.

### Open choice: excited $K^*$ terms

`kt`, `ku` contain the $K^*(892)$ pole plus $\beta$-weighted $K^*(1410)$, $K^*(1680)$ terms, and $\beta$ is shared by all vector species within a channel (Qin:2024ulb convention in `fit_babar_kkpi_g_coefficients.propagator`). Two options are not equivalent, because a $K^*$-only rescaling of $\beta$ is not available:
- (a) $\sqrt{\kappa_q}$ multiplies only the $K^*(892)$ pole (the factor is then strictly a $K^*(892)K\pi$ coupling);
- (b) $\sqrt{\kappa_q}$ multiplies the whole $K^*$ propagator (keeps the excited-to-ground ratios $\beta$ fixed).

At the endpoint with $\kappa_q$ set to the PDG-matching values and charge-specific widths, the BaBar $\chi^2$ (statistical residuals only, reduced orders, other parameters not refitted) is $271.306$ for (a) and $271.341$ for (b). The choice is left to the user; it is not resolved here.

## Step 4: consequence

$\Gamma_{K^*_q}$ and the $K^*_q$ vertex depend on the same coupling, so $\kappa_q$ must change both: width rows, $K^*\to K\gamma$ denominators and propagator widths scale by $\kappa_q$, and the $K^*$ vertex by $\sqrt{\kappa_q}$. Freeing only the total width (as the subclass option `kstar_width_treatment='free_total_widths'` does) describes a different model in which the $K^*$ total width is not produced by the $K\pi$ coupling that appears in the BaBar amplitude; since the measured $K\pi$ fractions are $99.90\%$ and $99.75\%$ ([width audit](shared-even-sector-gv-constraint-audit.md)), that model is not consistent with $K^*\to K\pi$ dominance. Numerically, at the PDG-matching factors with charge-specific widths, width-only and width-plus-vertex amplitudes differ by up to $2.0\%$ relative, and the reduced-order BaBar $\chi^2$ values are: baseline $269.723$, width only $268.817$, width and vertex (a) $271.306$, (b) $271.341$, width and vertex with the mean width $273.126$. These are evaluations at a fixed point, not fit results.

## Status of each ingredient of the baseline relation

- $\Gamma_q\propto g_q^2$ and vertex $\propto g_q$ at tree level: **derived** (exact).
- $\Gamma_{\rm tot}(K^*)=\Gamma(K\pi)$, i.e. $K\gamma$ and other modes neglected: **assumption** (deficits $0.10\%$ and $0.25\%$).
- One $G_V$ for $K^*$ and $\rho,\omega,\phi$ (SU(3) relation of the single VPP operator): **assumption**; this is exactly what $\kappa_q$ relaxes for the $K^*$.
- Physical parent and daughter masses in phase space, fitted tensor mass $M_T=m_{K^*}$ in the $1/M_T^2$ vertex factor: **assumption** (inherited baseline convention).
- The BaBar `direct` normalization uses the fixed $M_V=0.770$ GeV, while the width uses $M_T=m_{K^*}$: **inherited assumption, not audited here**. Hence $\kappa=1$ means "same $G_V$", not literally the same on-shell $K^*K\pi$ coupling in both places.
- Running $K^*$ width = on-shell width times a fixed reference shape that includes the $K\eta$ channel (`fit_babar_kkpi_g_coefficients.py` line 74): **assumption**; with $\kappa_q$ the $K\eta$ part is scaled too.
- Charge-mean propagator width: **isospin-limit assumption** of the baseline, removable with `separate_kstar_charge_states`.
- Published $\pi^0$ diagonals (isospin breaking in the traces, sums $1.52888$ vs $1.47159$): **literature input** [@Alda:2026nrd].
- Positive sign branch $g_q/g_q^{\rm sat}=+\sqrt{\kappa_q}$: **assumption**.

## Added coordinates and bounds

Two coordinates, $\kappa_{K^{*\pm}}$ and $\kappa_{K^{*0}}$; one common $\kappa$ if the isospin limit is imposed (then the mean-width treatment is exact and no charge mixing arises). Natural bound: $\kappa_q>0$ (the baseline rejects nonpositive widths); no natural upper bound. Start: $\kappa_q=1$ (baseline). With $\kappa_q$ free, the two $K^*$ width rows no longer constrain $G_V$ (each factor can absorb its row up to the $K\gamma$ and BaBar correlations), so $G_V$ is then fixed by the $\phi\to K\bar K$ rows and the $s$-channel BaBar terms; the nominal degrees of freedom decrease by the number of added factors. An equivalent parametrization uses the couplings $\sqrt{\kappa_q}\,G_V$ directly; the choice of optimizer coordinate is not made here.

## Verification summary

Script assertions (Steps 1–3) and five unit tests pass. Not checked: high quadrature orders, any refit, the $M_V$ versus $m_{K^*}$ normalization of `direct`.
