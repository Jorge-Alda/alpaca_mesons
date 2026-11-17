---
title: rho0-omega mixing and the phi to K+K- versus K0 K0bar split at first order in m_u - m_d
status: mixed
created: 2026-11-17
updated: 2026-11-17
bibtex_keys: [Alda:2026nrd, Ecker:1988te, ParticleDataGroup:2024cfk]
provenance: analytical-derivation
---

# $\rho^0$–$\omega$ mixing and the $\phi\to K^+K^-$ versus $K^0\bar K^0$ split

This record settles task 5.1 of the change `extend-joint-vector-model-diagnostics` (item 4 of the [candidate causes](joint-vector-fit-inadequacy-candidate-causes.md)). The algebra is derived with SymPy. The numbers are evaluated at the fixed accepted endpoint `basin_probes/even_shared_widths`; no fit, refit or scan was run.

**Chosen parameterization.** Extension 4 needs no new coordinate: at the order of the joint fit, $\rho^0$–$\omega$ mixing and the $\phi$ charged/neutral split are fixed by the existing coordinates $m_\rho$, $m_{K^*}$, $\kappa_V$, $z_V$ and the Alda:2026nrd quark-mass spurion, and the task-1.2 diagnostic shows a $\phi\to K\bar K$ residual mainly in the sum of the two channels, which isospin breaking cannot change at first order.

Script: [`scripts/derive_rho_omega_mixing_and_phi_channel_split.py`](../../../scripts/derive_rho_omega_mixing_and_phi_channel_split.py); tests: [`scripts/test_derive_rho_omega_mixing_and_phi_channel_split.py`](../../../scripts/test_derive_rho_omega_mixing_and_phi_channel_split.py).

```bash
.venv/bin/python scripts/derive_rho_omega_mixing_and_phi_channel_split.py
.venv/bin/python -m unittest scripts.test_derive_rho_omega_mixing_and_phi_channel_split
```

## Operators and order

No operator is added. The isospin breaking comes from the operators already in the joint fit, evaluated with $x_u\ne x_d$ in the vacuum spurion $\chi_+=\mathrm{diag}(x_u,x_d,x_s)$, $x_q=4B_0m_q$:

$$
\mathcal L_{\rm mass}=\frac{M_V^2}{4}\langle V_{\mu\nu}V^{\mu\nu}\rangle-e_m^V\langle V_{\mu\nu}V^{\mu\nu}\chi_+\rangle+\frac{\Delta_V}{12}\langle V_{\mu\nu}\rangle\langle V^{\mu\nu}\rangle+\frac{\kappa_V}{4}\langle V_{\mu\nu}\rangle\langle V^{\mu\nu}\chi_+\rangle,
$$

from the [nonideal mass record](../foundations/vector-nonideal-mass-action.md) and the [mass extension](../foundations/vector-mass-extension.md), and the photon source term

$$
\frac{F_V}{2\sqrt2}\langle V_{\mu\nu}f_+^{\mu\nu}\rangle+\frac{z_V}{2\sqrt2}\langle V_{\mu\nu}\{\chi_+,f_+^{\mu\nu}\}\rangle+\frac{y_V}{2\sqrt2}\langle V_{\mu\nu}\rangle\langle\chi_+f_+^{\mu\nu}\rangle
$$

from the [photon–vector diagonalization](../foundations/rcht-photon-vector-diagonalization.md) ($y_V=0$ in the baseline).

**Order.** Tree level; first order in $x_u-x_d$ (that is, in $m_u-m_d$); the mass terms at order $e^0$ and linear in $\chi_+$; the photon couplings at order $e$. This is the order of the joint fit, which so far sets $x_u=x_d=x_l$. The order-$e^2$ photon-exchange mixing is evaluated only as a size estimate (Step 6) and is not part of the result.

## Conventions and units

Natural units, GeV. $e_m^V,\kappa_V$ are dimensionless; $\Delta_V$, $x_q$ and mass-matrix entries are in GeV$^2$; $F_V$ in GeV; $z_V,y_V$ in GeV$^{-1}$. Generators $T_\rho=\mathrm{diag}(1,-1,0)/\sqrt2$, $T_{\omega_I}=\mathrm{diag}(1,1,0)/\sqrt2$, $T_{\phi_I}=\mathrm{diag}(0,0,1)$, physical $T_\omega=\cos\theta_VT_{\omega_I}-\sin\theta_VT_{\phi_I}$, $T_\phi=\sin\theta_VT_{\omega_I}+\cos\theta_VT_{\phi_I}$ ([nonideal mass record](../foundations/vector-nonideal-mass-action.md)). For real neutral components $V=\sum_iv_iT_i$ the mass matrix is $(\mathcal M_V^2)_{ij}=2\,\partial^2\mathcal L/\partial v_i\partial v_j$; for a complex charged field it is $2\,\partial^2\mathcal L/\partial V\partial V^\dagger$. Photon residues $R_i=\langle T_iS_\gamma\rangle$. The $\phi\to K\bar K$ flavor traces and widths are those of `flavor_factors` and `width_coefficient` in [`scripts/fit_vector_gv.py`](../../../scripts/fit_vector_gv.py), as used by the joint fit.

## Assumptions and approximations

1. **First order in $m_u-m_d$** (approximation). Products of two isospin-breaking entries are dropped. The unit test compares with exact complex eigenvectors: the difference is below $5\%$ of the first-order value.
2. **Constant on-shell widths in the mixing denominators** (approximation). Each neutral state enters the propagator matrix $(\mathcal M^2-s)^{-1}$ with $m_V^2\to m_V^2-im_V\Gamma_V$, using the endpoint total widths. Running widths are not used. A real, zero-width diagonalization is **not** valid for $\rho^0$–$\omega$ (Step 2).
3. **Spurion values** (input). $x_u,x_d,x_s=(0.026272,\,0.046960,\,0.943684)$ GeV$^2$, the Alda:2026nrd arXiv-v1 central values; their mean equals the fit's $x_l=0.036616$ GeV$^2$ (asserted).
4. **Electromagnetic terms omitted** (approximation, as in the fit and in the [double-trace record](joint-vector-double-trace-source-mass-and-kkpi.md)): order-$e^2$ photon exchange between vectors, local $\langle VQVQ\rangle$-type mass terms, the direct $\gamma^*\to K^+K^-$ term and the Coulomb interaction of $K^+K^-$ in $\phi$ decays. The last one is a loop effect beyond tree level. Its size is not computed here; for slow kaons it is expected to be of the same order as the residual ratio change (unsourced expectation, to be checked against a source before use).
5. **No direct isospin-breaking $VPP$ coupling.** The order-$p^4$ operators with coefficients $r_1,r_2,r_3$ belong to the separate [isospin-breaking $VPP$ option](joint-vector-isospin-breaking-vpp-likelihood-option.md); they are only used for comparison in Step 5.
6. **Leptonic widths** (assumption). The PDG $\Gamma(V\to e^+e^-)$ is identified with the squared modulus of the full pole residue, including the admixture of the other isospin state.

## Inputs

- Mass terms and isospin-limit matrix entries: [nonideal mass record](../foundations/vector-nonideal-mass-action.md), [mass extension](../foundations/vector-mass-extension.md) [@Ecker:1988te].
- Source matrix $S_\gamma$, residues, and the order-$e^2$ kinetic matrix $K_V=\mathbf 1-\boldsymbol\epsilon\boldsymbol\epsilon^T/Z$ with $\epsilon_i=\sqrt2eR_i/M_i$: [photon–vector diagonalization](../foundations/rcht-photon-vector-diagonalization.md).
- Spurion: `data/vector_sector/alda_2026nrd_arxiv_v1.json` ([Alda:2026nrd inputs](alda-neutral-amplitude-inputs.md)) [@Alda:2026nrd].
- Endpoint coordinates, total-width pulls and residues: `data/vector_sector/joint_vector_fit_adequacy.json` and the neutral-sector construction used by the fit ([goodness of fit](joint-vector-fit-goodness-of-fit.md), [likelihood baseline](babar-kkpi-global-vector-sector-fit-baseline.md)). At the endpoint $m_\rho=0.77539$ GeV, $m_{K^*}=0.89084$ GeV, $\Delta_V=0.017829$ GeV$^2$, $\kappa_V=0.046146$, $F_V=0.16896$ GeV, $z_V=-0.17098$ GeV$^{-1}$, $G_V=0.059068$ GeV, $\theta_V=0.09524$ rad.
- $\phi\to K\bar K$ measured fractions, covariance and endpoint predictions: the [task-1.2 channel diagnostic](joint-vector-phi-kaon-residual-diagnostic.md) (its script is imported) [@ParticleDataGroup:2024cfk].
- PDG $\mathcal B(\omega\to\pi^+\pi^-)=(1.53\pm0.12)\%$: [PDG extract](../../extracts/ParticleDataGroup:2024cfk/vector-to-two-pseudoscalar-inputs.md); $K^*$ masses and leptonic widths: `data/vector_sector/pdg2024.json`; radiative branching fractions and total widths: `data/odd_parity/odd_low_energy_measurements.json` [@ParticleDataGroup:2024cfk].
- $r_2-r_3$ split of the $\phi$ coupling and the task-3.2 best-fit $r$: [operators record](joint-vector-isospin-breaking-vpp-operators.md), [fit outcome](joint-vector-isospin-breaking-vpp-fit-outcome.md).

## Method and results, step by step

### Step 1: mass matrix with $x_u\ne x_d$ (derived)

SymPy differentiates $\mathcal L_{\rm mass}$. In the $(\rho^0,\omega_I,\phi_I)$ basis the isospin-limit entries reproduce the mass-extension matrix exactly, and the new entries are

$$
(\mathcal M_V^2)_{\rho\omega_I}=(x_u-x_d)\Big(\frac{\kappa_V}{2}-2e_m^V\Big),\qquad
(\mathcal M_V^2)_{\rho\phi_I}=\frac{\kappa_V(x_u-x_d)}{2\sqrt2}.
$$

$\Delta_V$ does not enter, because $\langle T_\rho\rangle=0$. The diagonal entries keep their isospin-average form, so $m_{\rho^0}=m_{\rho^\pm}$ at this order. The charged masses are $m_{\rho^\pm}^2=M_V^2-2e_m^V(x_u+x_d)$, $m_{K^{*+}}^2=M_V^2-2e_m^V(x_u+x_s)$, $m_{K^{*0}}^2=M_V^2-2e_m^V(x_d+x_s)$. Hence

$$
(\mathcal M_V^2)_{\rho\omega_I}-\frac{\kappa_V}{2}(x_u-x_d)=m_{K^{*+}}^2-m_{K^{*0}}^2 .
$$

If the fit coordinate $m_{K^*}$ is read as the $K^{*+}$ mass, $e_m^V=(m_\rho^2-m_{K^{*+}}^2)/[2(x_s-x_d)]$, which is $1.15\%$ larger than the isospin-average calibration.

In the physical basis, $(\mathcal M_V^2)_{\rho\omega}=\cos\theta_V(\mathcal M_V^2)_{\rho\omega_I}-\sin\theta_V(\mathcal M_V^2)_{\rho\phi_I}$ and $(\mathcal M_V^2)_{\rho\phi}=\sin\theta_V(\mathcal M_V^2)_{\rho\omega_I}+\cos\theta_V(\mathcal M_V^2)_{\rho\phi_I}$.

At the endpoint, $e_m^V=-0.10603$:

| entry | GeV$^2$ |
|---|---|
| $(\mathcal M_V^2)_{\rho\omega_I}$ | $-4.865\times10^{-3}$ |
| $(\mathcal M_V^2)_{\rho\phi_I}$ | $-3.38\times10^{-4}$ |
| $(\mathcal M_V^2)_{\rho\omega}$ | $-4.810\times10^{-3}$ |
| $(\mathcal M_V^2)_{\rho\phi}$ | $-7.99\times10^{-4}$ |
| $m_{K^{*+}}^2-m_{K^{*0}}^2$ | $-4.387\times10^{-3}$ ($K^{*+}$ lighter by $2.46$ MeV; PDG $3.88$ MeV) |

### Step 2: first-order pole residue (derived)

For an inverse propagator $\begin{pmatrix}a-s&\delta\\\delta&b-s\end{pmatrix}$, SymPy finds that the residue at the pole $s=b$ factorizes as $w\,w^T$ with $w=(\delta/(b-a),\,1)$ at first order in $\delta$. The physical $\omega$ therefore couples as $T_\omega$ plus a $\rho^0$ component

$$
\frac{(\mathcal M_V^2)_{\rho\omega}}{(m_\omega^2-im_\omega\Gamma_\omega)-(m_\rho^2-im_\rho\Gamma_\rho)}=-0.0047+0.0443\,i\quad(\text{modulus }0.0445),
$$

and the physical $\rho^0$ carries the opposite $\omega$ component. The same formula with $\phi$ gives the $\rho^0$ component of $\phi$: $-1.72\times10^{-3}+0.43\times10^{-3}\,i$; without widths it is $-1.82\times10^{-3}$.

For $\rho^0$–$\omega$ the zero-width value would be $-0.42$, with a rotation angle of $-0.35$ rad. This is far outside first order, because $m_\omega^2-m_\rho^2=0.0114$ GeV$^2$ is ten times smaller than $m_\rho\Gamma_\rho\approx0.114$ GeV$^2$. A real diagonalization of the $3\times3$ matrix, like the one used for $\omega$–$\phi$, must therefore not be used for $\rho^0$–$\omega$.

### Step 3: isospin breaking in the photon residues (derived)

With $x_u\ne x_d$ the source term gives, relative to the isospin-average residues,

$$
\delta R_\rho=\frac{z_V(x_u-x_d)}{3\sqrt2},\qquad
\delta R_{\omega_I}=\frac{(z_V+y_V)(x_u-x_d)}{\sqrt2},\qquad
\delta R_{\phi_I}=\frac{y_V(x_u-x_d)}{2}.
$$

$F_V$ gives no isospin breaking. The $\omega$ shift is large compared with $R_\omega$, because the isoscalar charge factor suppresses $R_\omega$ while the $z_V$ term weights the $u$ quark with its full charge. The symbolic residues agree with `source_residues` of `scripts/implementation/vector_sector.py` (unit test).

### Step 4: $\phi\to K\bar K$ traces (derived)

With $T_\phi+c\,T_\rho$, where $c$ is the $\rho^0$ component of $\phi$, the traces used by the fit are

$$
\langle T[K^+,K^-]\rangle=a+\frac{c}{\sqrt2},\qquad
\langle T[K^0,\bar K^0]\rangle=a-\frac{c}{\sqrt2},\qquad
a=\frac{\sin\theta_V}{\sqrt2}-\cos\theta_V .
$$

The ratio of the two squared traces changes by $2\sqrt2\,c/a$ at first order, and their sum does not change at first order. **First-order isospin breaking can move only the charged-to-neutral ratio, not the sum.**

### Step 5: effect on the task-1.2 residuals (numerical, fixed endpoint)

| quantity | value |
|---|---|
| ratio change $\Gamma(K^+K^-)/\Gamma(K^0\bar K^0)$ from the $\rho^0$ component of $\phi$ | $+0.52\%$ (linear $+0.523\%$) |
| same, adding the order-$e^2$ photon-exchange estimate of Step 6 | $-0.43\%$ |
| ratio change preferred by the residual alone | $-5.9\%$ |
| predicted fractions $(K^+K^-,K^0\bar K^0)$ | $(0.50850,0.33551)\to(0.50982,0.33464)$ |
| Cholesky-whitened residuals | $(3.4991,2.3738)\to(3.7650,2.3341)$ |
| symmetric-whitened residuals (task-1.2 rotation) | $(4.1057,1.0108)\to(4.3413,0.8812)$ |
| block $\chi^2$ | $17.879\to19.623$ |

For comparison, the $VPP$ option gives a $\phi$ coupling split $(r_2-r_3)(x_u-x_d)$ ([operators record](joint-vector-isospin-breaking-vpp-operators.md)). At the task-3.2 best fit $r_2-r_3=1.24$ GeV$^{-2}$ this is a ratio change of about $-5.1\%$.

**Diagnostic decomposition.** Let $d$ be prediction minus measurement and $W$ the inverse covariance. Freeing one linear direction $v$ leaves $d^TWd-(v^TWd)^2/v^TWv$:

| direction freed | $\chi^2$ removed | $\chi^2$ left |
|---|---|---|
| common normalization of both channels (isospin symmetric: $G_V$, $\Gamma_\phi$, $SU(3)$ breaking) | $12.95$ | $4.92$ |
| charged/neutral split only (isospin breaking) | $9.64$ | $8.23$ |
| both | $17.88$ | $0$ |

The two directions are correlated through the $-0.72$ measurement correlation, so the removed amounts do not add up.

### Step 6: checks and size estimates (numerical, fixed endpoint)

- **$\omega\to\pi^+\pi^-$** (not a fit row). This decay proceeds only through the $\rho^0$ component at this order. Predicted $\mathcal B=2.8\%$ with the fit's $G_V$ for $\rho\to\pi\pi$ (tree width $117$ MeV), or $3.6\%$ with the measured $\rho^0$ width; PDG gives $(1.53\pm0.12)\%$. The order-$m_q$ mixing entry is therefore about $1.4$ to $1.5$ times too large in modulus. The omitted pieces are electromagnetic terms, a direct $\omega\pi\pi$ coupling and the energy dependence of the widths.
- **$K^*$ mass split.** The strong part predicts $m_{K^{*+}}-m_{K^{*0}}=-2.46$ MeV; PDG gives $-3.88$ MeV. The electromagnetic part is not included.
- **Order-$e^2$ photon exchange** (estimate, not included). From $\det(sK_V-D)=0$, one photon exchange adds $s\,\epsilon_i\epsilon_j=2e^2R_iR_js/(M_iM_j)$ to the off-diagonal entry. This gives $+6.5\times10^{-4}$ GeV$^2$ for $\rho\omega$ (13% of the strong entry, opposite sign) and $+1.45\times10^{-3}$ GeV$^2$ for $\rho\phi$, which is larger than the strong entry and has the opposite sign. With the $\rho\omega$ term added, $\mathcal B(\omega\to\pi^+\pi^-)$ becomes $2.1\%$. The direct $\gamma^*\to K^+K^-$ term and the Coulomb interaction of the same order are not computed. For the $\phi$ split only the size, $|\Delta|\lesssim0.5\%$, is therefore meaningful; its sign is not determined at this order.
- **Leptonic rows**, relative change at fixed coordinates:
  - $\Gamma(\omega\to e^+e^-)$: $+15.1\%$, of which $+16.3\%$ comes from $\delta R_\omega$ ($z_V$) and $-0.9\%$ from the $\rho^0$ component. The row error is $3.3\%$.
  - $\Gamma(\rho\to e^+e^-)$: $+1.8\%$. The row error is $0.85\%$.

  A refit would partly absorb these changes in $F_V$ and $z_V$; that is not tested here.
- **Radiative rows.** The bound is $2r|\mathrm{Re}\,c|+r^2|c|^2$, where $r$ is the amplitude ratio taken from the measured widths:

  | row | $r$ | bound | in row errors |
  |---|---|---|---|
  | $\rho^0\to\pi^0\gamma$ | $3.18$ | $5.0\%$ | $0.29$ |
  | $\omega\to\eta\gamma$ | $2.98$ | $4.5\%$ | $0.22$ |
  | $\omega\to\pi^0\gamma$ | | $0.3\%$ | $0.10$ |
  | $\rho^0\to\eta\gamma$ | | $0.3\%$ | $0.05$ |

  The $\rho^0\to\eta\gamma$ tension (pull $-3.27$) is not affected.
- **BaBar $KK\pi$.** At $\sqrt s=1.37$ GeV, the off-diagonal $\rho\omega$ propagator is $0.4\%$ of the diagonal one. The change is negligible compared with the bin errors (estimate).

## Choice of the minimal parameterization

1. **Pattern in the data.** The task-1.2 residual is mainly in the **sum** of the two channels. A common normalization removes $12.95$ of $17.88$. What is left in the ratio direction is $\chi^2=4.92$: nominal $p=0.026$ for one degree of freedom, and $p=0.085$ for the two-row block, which is above the block threshold $10^{-3}$.
2. **What isospin breaking can do.** By Step 4, first-order isospin breaking acts only on the ratio. The ratio part fixed by the fit's own operators is $+0.5\%$ (or $-0.4\%$ with the photon-exchange estimate), against the preferred $-4$ to $-6\%$. This is an order of magnitude too small. With the sign of the order-$m_q$ result it raises the block $\chi^2$ by $1.7$ at the fixed endpoint.
3. **Decision.**
   - The $\phi$ split gets **no new coordinate**. The diagnostic does not require a ratio coordinate for the block criterion.
   - If a ratio coordinate is wanted later, the existing combination $r_2-r_3$ of the $VPP$ option already spans that direction. Task 3.2 used it: $\phi$ block $2.5$.
   - $\rho^0$–$\omega$ mixing gets **no new coordinate**. Its entries are fixed by $e_m^V$, from $m_\rho$ and $m_{K^*}$, together with $\kappa_V$ and $x_u-x_d$.
   - The sum residual belongs to candidate 2, the $SU(3)$-breaking $VPP$ coupling, not to extension 4.
4. **Consequence for task 5.2.** A zero-coordinate option ($x_u\ne x_d$ in the existing mass and source terms, with the complex first-order mixing of Step 2) can be implemented without changing the degrees of freedom.
   - Its largest effect at fixed coordinates is on $\Gamma(\omega\to e^+e^-)$, through the $z_V$ residue shift, and not on $\phi\to K\bar K$.
   - The $\omega\to\pi^+\pi^-$ check shows that the order-$m_q$ mixing is uncertain by about 50% in modulus.
   - Whether this option is worth a refit is a decision for the lead agent and the user. No refit was run.

## Verification

Run on 2026-11-17:
- The script assertions pass:
  - the isospin limit reproduces the mass-extension matrix and the charged masses;
  - the off-diagonal entries and the $K^*$ identity hold;
  - the pole-residue component is $\delta/(b-a)$;
  - the residue shifts hold;
  - the first-order $\phi$ sum change is zero and the ratio formula holds;
  - the symmetric-limit endpoint matrix reproduces the fit's $\omega,\phi$ masses and $\theta_V$ to $10^{-12}$;
  - the spurion mean equals $x_l$.
- The 13 unit tests pass:
  - the symbolic results;
  - agreement with `source_residues`;
  - first order against exact complex eigenvectors (within $5\%$);
  - the zero-width $\rho\omega$ value exceeding first order;
  - the $K^*$ identity in numbers;
  - the rotated residuals reproduce the covariance quadratic form;
  - the endpoint block $\chi^2=17.8785$ of the [task-1.2 record](joint-vector-phi-kaon-residual-diagnostic.md);
  - the projection remainders, checked by a brute-force scan and by the two-direction solve;
  - the bracket on $\mathcal B(\omega\to\pi^+\pi^-)$.

Not checked: any refit; running widths in the mixing denominators; the complete order-$e^2$ isospin breaking (Coulomb, direct photon, local electromagnetic terms); the BaBar likelihood with the mixing included.
