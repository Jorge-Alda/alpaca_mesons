---
title: Double-trace singlet-source structure y_V as a switchable option of the joint vector likelihood
status: mixed
created: 2026-11-17
updated: 2026-11-17
bibtex_keys: [Alda:2026nrd, Ecker:1988te, BaBar:2007ceh]
provenance: analytical-derivation
---

# Double-trace singlet-source structure $y_V$ as a likelihood option

This record describes how the coefficient $y_V$ of the double-trace singlet-source structure enters the joint vector likelihood as a switchable coordinate (task 4.2 of `extend-joint-vector-model-diagnostics`). The physics input is the [4.1 derivation](joint-vector-double-trace-source-mass-and-kkpi.md). The baseline likelihood and its earlier variations are unchanged. Fit plan and dry run: [joint-vector-double-trace-singlet-source-fit-plan.md](joint-vector-double-trace-singlet-source-fit-plan.md).

Script: [`scripts/joint_vector_double_trace_singlet_source_option.py`](../../../scripts/joint_vector_double_trace_singlet_source_option.py) (class `DoubleTraceSingletSourceLikelihood`, a subclass of `VariedLikelihood`); tests in [`scripts/test_joint_vector_double_trace_singlet_source_option.py`](../../../scripts/test_joint_vector_double_trace_singlet_source_option.py).

```bash
.venv/bin/python scripts/joint_vector_double_trace_singlet_source_option.py
.venv/bin/python -m unittest scripts.test_joint_vector_double_trace_singlet_source_option
```

Minimal use:

```python
model = DoubleTraceSingletSourceLikelihood(2, 16, 16, bin_block_size=4, double_trace_singlet_source=True)
residual = model.residual(model.extend(baseline_values, 0.01))   # y_V in GeV^-1
```

## Conventions and units

Conventions of Alda:2026nrd [@Alda:2026nrd]; antisymmetric-tensor vector fields [@Ecker:1988te]. GeV for masses, widths, energies and $F_V$; GeV$^2$ for invariants and $x_q=(\chi_+)_{qq}$ at the vacuum; GeV$^{-1}$ for $z_V$ and $y_V$; residuals dimensionless. The photon source matrix is
$$S_\gamma=F_VQ+z_V\{\chi_+,Q\}+y_V\langle\chi_+Q\rangle\mathbf 1,$$
and $F_l=F_V+2z_Vx_l$, $F_s=F_V+2z_Vx_s$ ([4.1 record](joint-vector-double-trace-source-mass-and-kkpi.md)). The code name of the coordinate is `double_trace_singlet_source_y_V_GeV_inverse`; the user approved $y_V$ and the name "double-trace singlet-source structure".

## Inputs

- Residue shifts $\delta R_{\omega_I}=\sqrt2y_V(x_l-x_s)/3$, $\delta R_{\phi_I}=y_V(x_l-x_s)/3$, $\delta R_\rho=0$, and the $KK\pi$ replacements: [4.1 record](joint-vector-double-trace-source-mass-and-kkpi.md), `scripts/derive_double_trace_vector_mass_and_kkpi.py` (functions `photon_residues`, `kkpi_source_factors`).
- `source_residues(generators, chi_plus, FV, z, y)` of `scripts/implementation/vector_sector.py` and the `double_trace` hypothesis of `action_sector` in `scripts/derive_odd_low_energy_jacobians.py` (uses `point[6]` as $y_V$).
- Baseline likelihood and BaBar amplitude: [babar-kkpi-global-vector-sector-fit-baseline.md](babar-kkpi-global-vector-sector-fit-baseline.md), [joint-vector-likelihood-variation-options.md](joint-vector-likelihood-variation-options.md); BaBar $KK\pi$ data [@BaBar:2007ceh].

## Method

1. Option `double_trace_singlet_source=True` appends $y_V$ after the 32 baseline coordinates. Start $0$; bounds $[-0.5,0.5]$ GeV$^{-1}$ (reasons in the [fit plan](joint-vector-double-trace-singlet-source-fit-plan.md)); outside the bounds `unpack` raises and the residual is the baseline penalty vector.
2. **Residues.** `unpack` runs the inherited unpack and then rebuilds the vector sector with `action_sector((point, y_V), 'double_trace')`. The mass matrix, masses and mixing angle do not depend on $y_V$ (4.1 result (a)); only the residues change. They enter the $\rho,\omega,\phi$ leptonic widths and, through the neutral response $W(s)=\sum_iR_iT_i/(M_i^2-s)$, the photon-insertion parts of the 12 odd radiative rows.
3. **BaBar current.** In the baseline the factor `light` $=F_l$ multiplies both the $\rho^0$ insertion (`source_rho`) and the $\omega_I$ insertion (inside `source_light`, `source_strange`). With the option on, `_babar_blocks_full` is a copy of the baseline method in which
   - `source_rho` $=F_l\,D_\rho(q)$ (isovector factor, also used by the $a_2$ tensor term);
   - `source_light` $=F_l^{(y)}\Pi_{ll}-\sqrt2F_s^{(y)}\Pi_{ls}$ and `source_strange` $=F_s^{(y)}\Pi_{ss}-F_l^{(y)}\Pi_{ls}/\sqrt2$, with the isoscalar factors
   $$F_l^{(y)}=F_l+2y_V(x_l-x_s),\qquad F_s^{(y)}=F_s-y_V(x_l-x_s)$$
   (function `isoscalar_source_factors`), and $\Pi$ the ideal-basis neutral propagators at the $e^+e^-$ invariant;
   - the propagator cache key also contains $y_V$.
   Everything else in the amplitude is the baseline code.
4. Jacobian: inherited 32 columns unchanged; the $y_V$ column by forward differences with step $\sqrt\epsilon\max(1,|y_V|)$, or the step supplied by the quadrature audit.
5. Degrees of freedom: $75-33=42$ (baseline 43).

## Assumptions (labelled)

- Inherited from the 4.1 record: tree level, order-$e^2$ electromagnetic mass shifts omitted for every source term, exact strong isospin $x_u=x_d$, excited vector poles of the BaBar current unchanged.
- The baseline source-branch check ($F_l>0$, $F_s<0$) is kept on $F_V,z_V$; the $y_V$-shifted isoscalar factors are not constrained in sign (choice).
- Only the baseline source hypothesis, mass branch, regulator and fixed poles are supported with the option on; other `VariedLikelihood` options that act on `unpack` or on the BaBar amplitude are refused. Row-selection and error options (excluded blocks, energy window, $\phi$ correlation, BaBar systematic) remain available.

## Rows affected

At $y_V\ne0$ and fixed other coordinates: the 12 odd rows, the $\omega$ and $\phi$ leptonic widths and the 45 BaBar rows change. The four mass rows, the seven shared-width rows, the two $\phi\to K\bar K$ rows and the normalization constraints are bitwise unchanged. The $\rho$ leptonic width changes only by rounding ($<10^{-12}$ in the pull) because $\langle T_\rho\rangle=0$.

At the best accepted endpoint, $y_V=0.022$ GeV$^{-1}$ moves the physical residues from $(R_\rho,R_\omega,R_\phi)=(0.1106,\,0.0318,\,0.0545)$ GeV to $(0.1106,\,0.0231,\,0.0470)$ GeV and the isoscalar factors from $(F_l,F_s)=(0.1564,-0.1537)$ GeV to $(0.1165,-0.1338)$ GeV (printed by the script).

## Verification (run 2026-11-17)

`scripts.test_joint_vector_double_trace_singlet_source_option`: 17 tests pass (orders $(1,3,4)$, best accepted endpoint).
- Option off: residual and Jacobian bitwise equal to `GlobalLikelihood`. Option on at $y_V=0$: residual and first 32 Jacobian columns bitwise equal, with four-bin blocks and with full grids; the BaBar derivative block is bitwise equal.
- `isoscalar_source_factors` equals the SymPy `shifted_light`, `shifted_strange` of the 4.1 script to $10^{-14}$; the 4.1 `kkpi_source_factors` differences vanish.
- Ideal residues from `source_residues` equal the SymPy `photon_residues` at $y_V=0.013$ GeV$^{-1}$; the sector's physical residues are the rotated ideal residues.
- Each cached BaBar insertion equals propagator times residue, assembled independently in the test: `source_rho`$/\sqrt2=D_\rho R_\rho$, `source_light`$/(3\sqrt2)=\Pi_{ll}R_{\omega_I}+\Pi_{ls}R_{\phi_I}$, $-$`source_strange`$/3=\Pi_{ls}R_{\omega_I}+\Pi_{ss}R_{\phi_I}$, for both channels ($10^{-11}$). This ties the BaBar insertions to the same residues as the leptonic widths.
- `rs`, `lls`, `kt`, `ku` and `source_rho` are bitwise unchanged by $y_V$; only the expected rows change; the cache follows $y_V$.
- Names, bounds, start, length check, penalty outside the bounds, refused option combinations, finite Jacobian and supplied finite-difference steps.

Also rerun and passing: `scripts.test_derive_double_trace_vector_mass_and_kkpi` (9 tests), `scripts.test_joint_vector_likelihood_options`.
