---
title: Isospin-breaking vector-pseudoscalar-pseudoscalar operators and the charged/neutral K-star split
status: mixed
created: 2026-11-09
updated: 2026-11-09
bibtex_keys: [Alda:2026nrd, Ecker:1988te, ParticleDataGroup:2024cfk]
provenance: analytical-derivation
---

# Isospin-breaking $VPP$ operators and the charged/neutral $K^*$ split

This record asks whether an isospin-breaking even-parity $VPP$ operator can produce the charged/neutral $K^*(892)$ coupling asymmetry that the joint vector likelihood prefers, and what the same operator predicts for $\phi\to K\bar K$, $\rho\to\pi\pi$ and $\rho,\omega\to K\bar K$. No refit is done. The algebra is derived; the coefficient values and size comparisons are numerical estimates at fixed inputs.

Script: [`scripts/derive_isospin_breaking_vpp_operators.py`](../../../scripts/derive_isospin_breaking_vpp_operators.py), tests in [`scripts/test_derive_isospin_breaking_vpp_operators.py`](../../../scripts/test_derive_isospin_breaking_vpp_operators.py).

```bash
.venv/bin/python scripts/derive_isospin_breaking_vpp_operators.py
.venv/bin/python -m unittest scripts.test_derive_isospin_breaking_vpp_operators
```

## Conventions and units

Natural units; GeV for masses, widths, $G_V$ and $F=0.093$ GeV; GeV$^{-2}$ for the ratios $r_i$ below. Chiral building blocks, $U=u^2=\exp(i\sqrt2\Phi/F)$, $u_\mu=-\sqrt2\,\partial_\mu\Phi/F+O(\Phi^3)$ and $\chi_\pm$ as in [rcht-basics.md](../foundations/rcht-basics.md) [@Ecker:1988te]. The LO even-parity interaction is $\frac{iG_V}{\sqrt2}\langle V_{\mu\nu}u^\mu u^\nu\rangle$ ([vector-gv-fit.md](../foundations/vector-gv-fit.md)). Flavor matrices of the decaying vector and of the two pseudoscalars are exactly those of `flavor_factors` in [`scripts/fit_vector_gv.py`](../../../scripts/fit_vector_gv.py) (line 34), including the published $\pi^0$ diagonals of Alda:2026nrd (line 43) and the nonideal $\omega$–$\phi$ angle of the provisional fit (`fixed_parameters`, line 102) [@Alda:2026nrd]. Charge conjugation: $V_{\mu\nu}\to-V_{\mu\nu}^T$, $u_\mu\to u_\mu^T$, $\chi_\pm\to\chi_\pm^T$ [@Ecker:1988te].

## Proposed notation (awaiting user approval)

Coefficients of the three surviving operators: $\lambda^V_{\chi,1}$, $\lambda^V_{\chi,2}$, $\lambda^V_{\chi,3}$ (mass dimension $-1$), and their ratios to the LO coupling
$$r_i=\frac{\sqrt2\,\lambda^V_{\chi,i}}{G_V}\quad(\mathrm{GeV}^{-2}).$$
In the code: `coefficient_anticommutator`, `coefficient_sandwich`, `coefficient_vector_singlet` and `r1`, `r2`, `r3`. The per-charge factors $\kappa_{K^{*\pm}},\kappa_{K^{*0}}$ are those of the [saturation-freedom record](joint-vector-kstar-kpi-saturation-freedom.md), also awaiting approval.

## Inputs

- Required asymmetry: $\kappa_{K^{*\pm}}=1.147715$, $\kappa_{K^{*0}}=1.048637$ ([saturation-freedom record](joint-vector-kstar-kpi-saturation-freedom.md), Step 2). The script recomputes $\kappa_{K^{*\pm}}/\kappa_{K^{*0}}$ from the PDG $K\pi$ widths $51.4$ and $47.3$ MeV and parent masses $0.89167$, $0.89555$ GeV in `data/vector_sector/gv_measurements.json` [@ParticleDataGroup:2024cfk] and asserts agreement to $2\times10^{-6}$.
- Quark-mass spurion: $B_0m_q=(0.006568,\,0.01174,\,0.235921)$ GeV$^2$ for $u,d,s$ from `data/vector_sector/alda_2026nrd_arxiv_v1.json` ([alda-neutral-amplitude-inputs.md](alda-neutral-amplitude-inputs.md)) [@Alda:2026nrd]. Hence $\chi_+|_{\Phi=0}=4B_0\,\mathrm{diag}(m_u,m_d,m_s)=(0.026272,\,0.046960,\,0.943684)$ GeV$^2$, consistent with $4B_0\hat m=0.036616$ GeV$^2$ of [vector-stage3-calibration.md](../foundations/vector-stage3-calibration.md), and $(m_d-m_u)/m_s=0.02192$.
- $\phi\to K\bar K$: measured branching fractions $0.491\pm0.005$, $0.339\pm0.004$, correlation $-0.72$ (same JSON) [@ParticleDataGroup:2024cfk]; endpoint predictions $(0.5084953,\,0.3355123)$ from [joint-vector-phi-kaon-residual-diagnostic.md](joint-vector-phi-kaon-residual-diagnostic.md).
- Candidate motivation: candidates 2 and 4 of [joint-vector-fit-inadequacy-candidate-causes.md](joint-vector-fit-inadequacy-candidate-causes.md).

## Assumptions (none of these is established)

1. **Chiral counting.** $u_\mu=O(p)$, $\chi_\pm=O(p^2)$ ($m_q=O(p^2)$), $e^2=O(p^2)$; operators are classified by the order of the chiral tensor multiplying $V_{\mu\nu}$, as for the LO $O(p^2)$ term. Only terms linear in one quark-mass insertion, or bilinear in the charge matrix, are kept.
2. **Spurion values.** Constant quark masses (no source derivatives); $\chi_+$ at zero fields equals $4B_0m_q$ with the Alda:2026nrd central values; electromagnetic spurions $Q_L=Q_R=Q=\mathrm{diag}(2/3,-1/3,-1/3)$, with photon loops not computed.
3. **Vertex only.** Isospin breaking is placed in the $VPP$ vertex. The propagators and the $K^{*\pm}$, $K^{*0}$ masses are untouched; physical parent and daughter masses enter phase space as in the baseline; the tensor mass cancels in every ratio used here.
4. **$\pi^0$–$\eta$ mixing.** The published $\pi^0$ diagonals $(1.01707,-0.98257,-0.0114)/\sqrt2$ already contain the $\pi^0$–$\eta$–$\eta'$ admixture at LO [@Alda:2026nrd]. They are used in the LO trace and, unchanged, inside the new operators (no separate correction of the mixing by the new terms).
5. **Large $N_c$.** Single traces are leading in $1/N_c$; double traces are kept in the enumeration but only the one with $\langle V_{\mu\nu}\rangle$ is carried to the channels, because the others do not split $K^*$ or $\phi$ charge states (see Step 1). For the electromagnetic part only single traces are enumerated.
6. **Positive LO sign branch** and tree level, as in the baseline.

## Step 1: allowed structures (derived)

At $O(p^2)$ the only antisymmetric chiral tensors are $u^\mu u^\nu$ and $f_+^{\mu\nu}$; $\chi_\pm$ carry no Lorentz index and the antisymmetric part of $\nabla^\mu u^\nu$ is $f_-^{\mu\nu}/2$, which needs external fields. **No isospin-breaking $VPP$ operator exists at the order of the $G_V$ term.** The first ones are $O(p^4)$, i.e. one quark-mass insertion ($m_q$) or two charge matrices ($e^2$) on top of the two derivatives.

Expanding to two pseudoscalars, every such operator gives the LO momentum structure $V_{\mu\nu}\partial^\mu\Phi\,\partial^\nu\Phi$ times a flavor function of the vector matrix $V$, the two pseudoscalar matrices $P_1$ (index $\mu$), $P_2$ (index $\nu$), and the diagonal spurion $S$. The script enumerates all products of traces of $\{V,P_1,P_2,S\}$, antisymmetrizes in $P_1\leftrightarrow P_2$, projects on even $C$, and computes the rank on random matrices: **rank 5** for both $3\times3$ and $5\times5$ matrices (so no Cayley–Hamilton relation acts in this sector), of which **2 single traces**. An independent set is

| Structure (with its $C$ image) | Operator | Fate |
|---|---|---|
| $\langle V\{[P_1,P_2],S\}\rangle$ | $i\lambda^V_{\chi,1}\langle V_{\mu\nu}\{u^\mu u^\nu,\chi_+\}\rangle$ | kept |
| $\langle V(P_1SP_2-P_2SP_1)\rangle$ | $i\lambda^V_{\chi,2}\langle V_{\mu\nu}u^\mu\chi_+u^\nu\rangle$ | kept |
| $\langle V\rangle\langle[P_1,P_2]S\rangle$ | $i\lambda^V_{\chi,3}\langle V_{\mu\nu}\rangle\langle u^\mu u^\nu\chi_+\rangle$ | kept (only $\omega,\phi$ via their singlet part) |
| $\langle S\rangle\langle V[P_1,P_2]\rangle$ | $\langle\chi_+\rangle\langle V_{\mu\nu}u^\mu u^\nu\rangle$ | discarded: isospin symmetric, absorbed in $G_V$ |
| $\langle V[P_1,S]\rangle\langle P_2\rangle$ | $\langle V_{\mu\nu}[u^\mu,\chi_+]\rangle\langle u^\nu\rangle$ | discarded: needs a singlet pseudoscalar ($\eta'$); zero for $K\pi$, $K\bar K$, $\pi^+\pi^-$ |

$\langle V_{\mu\nu}u^\mu\rangle\langle u^\nu\chi_+\rangle$ (which would split $K^{*\pm}\to K^\pm\pi^0$ from $K^{*0}\to K^0\pi^0$ and $\rho^\pm$ from $\rho^0$) is $C$-odd and excluded. The $\chi_-$ operator $\langle V_{\mu\nu}[u^\mu,\nabla^\nu\chi_-]\rangle$ is not independent at two pseudoscalars: SymPy verifies $[P_1,\{P_2,S\}]-[P_2,\{P_1,S\}]=\{[P_1,P_2],S\}+2(P_1SP_2-P_2SP_1)$, so it equals the first plus twice the second structure. $\langle V_{\mu\nu}\{f_+^{\mu\nu},\chi_+\}\rangle$ and $\langle V_{\mu\nu}[f_-^{\mu\nu},\chi_-]\rangle$ need a photon and do not contribute to $VPP$ (they change $V$–$\gamma$ couplings, not considered here).

All kept operators are isospin breaking through the $\lambda_3$ part of $\chi_+$ and, unavoidably, $SU(3)$ breaking through its $\lambda_8$ and singlet parts: the same coefficient multiplies $m_u$, $m_d$ and $m_s$.

## Step 2: channel couplings (derived, SymPy)

With $u^\mu\to-\sqrt2\partial^\mu\Phi/F$ and $\chi_+\to4B_0m_q$, the vertex keeps the LO form $g_{VPP}=G_VM_{\rm phys}^2\,T/(\sqrt2F^2M_T)$ with the flavor trace replaced by
$$T=\mathrm{tr}(V[P_1,P_2])+r_1\,\mathrm{tr}\big(V\{[P_1,P_2],4B_0m_q\}\big)+r_2\,\mathrm{tr}\big(V(P_1\,4B_0m_q\,P_2-P_2\,4B_0m_q\,P_1)\big)+r_3\,\mathrm{tr}(V)\,\mathrm{tr}\big([P_1,P_2]\,4B_0m_q\big).$$
Relative shifts $\delta T/T$ per channel (isospin-limit $\pi^0$, ideal mixing; exact, linear in $r_i$):

| Channel | $r_1$ | $r_2$ | $r_3$ |
|---|---|---|---|
| $K^{*+}\to K^+\pi^0$ | $4B_0(m_u+m_s)$ | $4B_0m_u$ | 0 |
| $K^{*+}\to K^0\pi^+$ | $4B_0(m_u+m_s)$ | $4B_0m_d$ | 0 |
| $K^{*0}\to K^0\pi^0$ | $4B_0(m_d+m_s)$ | $4B_0m_d$ | 0 |
| $K^{*0}\to K^+\pi^-$ | $4B_0(m_d+m_s)$ | $4B_0m_u$ | 0 |
| $\phi\to K^+K^-$ | $8B_0m_s$ | $4B_0m_u$ | $4B_0(m_s-m_u)$ |
| $\phi\to K^0\bar K^0$ | $8B_0m_s$ | $4B_0m_d$ | $4B_0(m_s-m_d)$ |
| $\rho^0\to\pi^+\pi^-$, $\rho^+\to\pi^+\pi^0$ | $4B_0(m_u+m_d)$ | $2B_0(m_u+m_d)$ | 0 |
| $\rho^0\to K^+K^-$ / $K^0\bar K^0$ | $8B_0m_u$ / $8B_0m_d$ | $4B_0m_s$ | 0 |
| $\omega\to K^+K^-$ / $K^0\bar K^0$ | $8B_0m_u$ / $8B_0m_d$ | $4B_0m_s$ | $8B_0(m_u-m_s)$ / $8B_0(m_d-m_s)$ |

The $r_2$ structure weights each mode by the mass of the quark pair created between the two pseudoscalars ($u\bar u$ for $K^+\pi^0$ and $K^+\pi^-$, $d\bar d$ for $K^0\pi^+$ and $K^0\pi^0$, $s\bar s$ for $\rho,\omega\to K\bar K$). LO $\omega\to\pi^+\pi^-$ vanishes; the new terms give it a direct coupling proportional to $4B_0(m_u-m_d)(r_1-r_2/2)$ in units of the $\rho^0$ LO trace.

**Zero-coefficient reduction (checked).** At $r_i=0$, with the published $\pi^0$ diagonals and the fit angle, $T^2$ equals every entry of `flavor_factors` to $1.1\times10^{-16}$; in the isospin limit $T^2=1/2$ ($K^+\pi^0$) and $1$ ($K^+\pi^-$) exactly.

**Linear splits in the isospin limit (SymPy, equal-phase-space weights $1/3$, $2/3$):**
$$\delta_{K^{*\pm}}-\delta_{K^{*0}}=\Big(r_1-\frac{r_2}{3}\Big)4B_0(m_u-m_d),\qquad \delta_{\phi\to K^+K^-}-\delta_{\phi\to K^0\bar K^0}=(r_2-r_3)\,4B_0(m_u-m_d),\qquad \delta_{\rho^+}-\delta_{\rho^0}=0.$$
The vanishing $\rho^\pm/\rho^0$ split is the expected G-parity selection rule: a single $\lambda_3$ insertion is G-odd and $\rho\to\pi\pi$ is G-even.

## Step 3: coefficient that gives the required $K^*$ asymmetry (numerical, exact in $r$)

Exact ratio $R(r)=[\Gamma_{K^{*\pm}}(r)/\Gamma_{K^{*\pm}}(0)]/[\Gamma_{K^{*0}}(r)/\Gamma_{K^{*0}}(0)]$ with `width_coefficient` phase space (line 66), published $\pi^0$ and Alda spurions. Required value $R=\kappa_{K^{*\pm}}/\kappa_{K^{*0}}=1.094483$, i.e. a coupling ratio $1.046175$ ($4.6\%$). $r_3$ does not enter $R$ (checked symbolically).

| Scenario | Root nearest zero | $r\cdot4B_0(m_d-m_u)$ | $r\cdot4B_0m_s$ |
|---|---|---|---|
| $r_1$ only | $-0.6951$ GeV$^{-2}$ (other root $-1.0199$) | $-0.0144$ | $-0.656$ |
| $r_2$ only | $+3.7362$ GeV$^{-2}$ (other root $-27.59$) | $+0.0773$ | $+3.53$ |

**Size estimate (labelled as such).** In units of the relative coupling shift, the required isospin asymmetry $\simeq0.045$ is $2.1\times(m_d-m_u)/m_s$ and $19\times\alpha/\pi$. Because the operator coefficient multiplies $m_s$ with the same strength as $m_d-m_u$, the corresponding $SU(3)$-breaking shift is $|r\,4B_0m_s|=0.66$ ($r_1$) or $3.5$ ($r_2$). Assuming (estimate, not sourced) that natural $SU(3)$ breaking of vector couplings is $|r\,4B_0m_s|\lesssim0.3$, the isospin asymmetry these operators can give is $\lesssim0.3\times0.022\times(1\text{ or }1/3)\approx0.7\%$ or $0.2\%$ in the coupling, against $4.6\%$ required. For $r_1$ the exact solution only reaches $4.6\%$ because the $K^*$ coupling itself is reduced to $0.33$ of its LO value, which divides the isospin shift by $0.31$; the linear approximation does not apply there.

## Step 4: predictions at those coefficients (numerical)

| Quantity | $r_1$ only | $r_2$ only | Present rows |
|---|---|---|---|
| $\Gamma(\phi\to K^+K^-)/\Gamma(\phi\to K^0\bar K^0)$, relative change | $+0.9\%$ | $-16.3\%$ | preferred $-4.4\pm2.0\%$; pull after: $+2.75$ / $-6.06$ |
| $g(K^{*\pm}\to K^0\pi^+)/g(\rho^0\to\pi\pi)$ | $0.343$ | $1.034$ | $\rho\to\pi\pi$ not a $G_V$ row of the joint fit |
| $g(\phi\to K^0\bar K^0)/g(K^{*\pm}\to K^0\pi^+)$ | $-1.19$ (sign flip) | $0.823$ | tested by $\phi$ and $K^*$ rows together |
| $g(\rho^0\to K\bar K)/g(K^{*\pm}\to K^0\pi^+)$ | $2.96$ | $3.85$ | enters BaBar $s$-channel terms; not evaluated |
| $g(\rho^+)/g(\rho^0)$ | $1.000000$ | $0.99941$ (from the published $\pi^0$ only) | no $\rho^\pm/\rho^0$ coupling row |
| $\omega\to\pi^+\pi^-$ direct, in units of $\rho^0$ LO | $+0.0143$ | $+0.0385$ | diagnostic only; interference with $\rho$–$\omega$ mixing not computed |
| $K^*$ mode factors $T(r)/T(0)$: $K^+\pi^0$, $K^0\pi^+$, $K^0\pi^0$, $K^+\pi^-$ | $0.326,0.326,0.311,0.311$ | $1.136,1.175,1.136,1.098$ | — |

With $r_2$ the factor is not the same for the two modes of one charge state, so the single per-charge $\sqrt{\kappa_q}$ of the saturation record is not exact: the BaBar $K^+K^-\pi^0$ terms use $K^{*\pm}\to K^\pm\pi^0$ ($u\bar u$ created), while the charged width is dominated by $K^0\pi^+$ ($d\bar d$ created).

**Two-coefficient solutions with the $\phi$ ratio fixed at its preferred change $-4.43\%$:**
- $r_2=+3.7362$, $r_3=+1.2482$ GeV$^{-2}$: $\phi$ split repaired, but $g(\phi\to K^0\bar K^0)/g(\rho^0)$ becomes $2.02$ and $g(\phi)/g(K^*)=1.95$.
- $r_1=-0.8415$, $r_2=-0.7871$ GeV$^{-2}$: $g(K^*)/g(\rho^0)=0.16$, $g(\phi)/g(K^*)=-4.6$.

**Sign versus the candidate-causes motivation.** With $r_2$ the $\phi$ ratio moves in the preferred direction (lower $K^+K^-$ relative to $K^0\bar K^0$), with about four times the preferred size. The motivation of candidate 2 concerns the *sum* $\phi\to K\bar K$ against the $K^*$ widths: $r_2$ changes that sum only through the $\omega$–$\phi$ mixing angle ($g(\phi)/g(K^*)=0.82$), and $r_1$ changes it by an amount far beyond the residual. In every scenario the $4.6\%$ $K^*$ split comes with an $SU(3)$-breaking change of the $\phi$, $\rho$ or $K^*$ couplings that is excluded by the $\phi$ and $K^*$ rows (or a $\phi$ charge split excluded at $6\sigma$); the $\rho,\omega\to K\bar K$ factors would also enter the BaBar amplitudes, but this is not evaluated.

## Step 5: electromagnetic structures (numerical, estimate)

Single traces of $\{V,P_1,P_2,Q,Q\}$, antisymmetrized and $C$-even: rank 4. A representative set, with the linear charged-minus-neutral $K^*$ split, the $\phi$ split and the $\rho^+-\rho^0$ split per unit coefficient (relative to the LO trace, $Q$ in units of $e$), and the coefficient needed for the $4.6\%$ split:

| Structure | $K^*$ split | $\phi$ split | $\rho^+-\rho^0$ | Needed coefficient |
|---|---|---|---|---|
| $\langle VP_1P_2QQ\rangle$ | $+0.167$ | 0 | 0 | $+0.271=117\,\alpha/\pi$ |
| $\langle VP_1QQP_2\rangle$ | $-0.109$ | $+0.333$ | 0 | $-0.414=-178\,\alpha/\pi$ |
| $\langle VP_1QP_2Q\rangle$ | $+0.113$ | $-0.333$ | $+0.25$ | $+0.398=171\,\alpha/\pi$ |
| $\langle VQP_1P_2Q\rangle$ | $-0.333$ | 0 | $-0.5$ | $-0.135=-58\,\alpha/\pi$ |

Taking $\alpha/\pi$ as the natural size of such a coefficient (estimate, assumption), the electromagnetic operators fall short by a factor 60 to 180. Unlike the quark-mass operators, they can split $\rho^\pm$ from $\rho^0$ (the product $Q\otimes Q$ contains an $I=2$ part).

## Conclusion

The leading isospin-breaking $VPP$ operators are $O(p^4)$: two single traces, $i\lambda^V_{\chi,1}\langle V_{\mu\nu}\{u^\mu u^\nu,\chi_+\}\rangle$ and $i\lambda^V_{\chi,2}\langle V_{\mu\nu}u^\mu\chi_+u^\nu\rangle$, plus the $1/N_c$-suppressed $i\lambda^V_{\chi,3}\langle V_{\mu\nu}\rangle\langle u^\mu u^\nu\chi_+\rangle$ that acts only on $\omega,\phi$. Within the stated assumptions, none of them, alone or in pairs, gives the $4.6\%$ charged/neutral $K^*$ coupling asymmetry without an $SU(3)$-breaking change of order one in the $K^*$, $\phi$ or $\rho$ couplings or a $\phi\to K^+K^-/K^0\bar K^0$ shift excluded by the present $\phi$ rows. Electromagnetic $O(e^2p^2)$ operators would need coefficients 60–180 times $\alpha/\pi$. This disfavours candidate 4 as the explanation of the $K^*$ width tension when it is realized through local $VPP$ operators; it does not test photon loops, isospin breaking in the $K^*$ propagators, or the $K^*$ width measurements themselves.

## Verification

Run on 2026-11-09: script assertions (zero-coefficient reduction against `flavor_factors`, $\chi_-$ identity, linear-split formulas, recomputed $\kappa$ ratio, LO $\phi$ ratio equal to the endpoint ratio $1.515579$) and 8 unit tests pass. Not checked: any refit; the BaBar likelihood with the $\rho,\omega\to K\bar K$ factors; $\omega,\phi\to\pi\pi$ widths including $\rho$–$\omega$ mixing; double traces with two charge matrices; uncertainties of $B_0m_q$.
