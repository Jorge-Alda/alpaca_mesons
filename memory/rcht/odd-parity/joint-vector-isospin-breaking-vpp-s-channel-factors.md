---
title: Isospin-breaking VPP factors of the s-channel vector to K Kbar vertices in the BaBar amplitude
status: mixed
created: 2026-11-09
updated: 2026-11-09
bibtex_keys: [Alda:2026nrd, Ecker:1988te, BaBar:2007ceh]
provenance: analytical-derivation
---

# Isospin-breaking $VPP$ factors of the $s$-channel $V\to K\bar K$ vertices

This record derives the factors $T(r)/T(0)$ that the three isospin-breaking $VPP$ operators (coefficients $r_1,r_2,r_3$) give to the $s$-channel $\rho,\omega,\phi\to K\bar K$ vertices of the BaBar $KK\pi$ amplitude of the joint vector likelihood. The algebra is derived (SymPy); the numbers at the fitted angle are numerical evaluations. No fit or refit is run.

Script: [`scripts/derive_isospin_breaking_vpp_s_channel_factors.py`](../../../scripts/derive_isospin_breaking_vpp_s_channel_factors.py); checks also in [`scripts/test_joint_vector_isospin_breaking_vpp_option.py`](../../../scripts/test_joint_vector_isospin_breaking_vpp_option.py) (classes `SChannelDerivation`, `SChannelTerms`).

```bash
.venv/bin/python scripts/derive_isospin_breaking_vpp_s_channel_factors.py
```

## Conventions and units

GeV$^2$ for $X=\chi_+|_{\Phi=0}=4B_0\,\mathrm{diag}(m_u,m_d,m_s)\equiv\mathrm{diag}(X_u,X_d,X_s)$; GeV$^{-2}$ for $r_i$; the factors are dimensionless. Conventions of Alda:2026nrd and of the operator record [@Alda:2026nrd; @Ecker:1988te]. Flavor matrices: $K^+=E_{13}$, $K^0=E_{23}$, $K^-=E_{31}$, $\bar K^0=E_{32}$ ($E_{ij}$ the matrix unit, 1-based); a decaying $\rho^+$ is represented by $E_{21}$ and $\rho^-$ by $E_{12}$, as $\rho^+\to\pi^+\pi^0$ in `fit_vector_gv.flavor_factors`; $\rho^0=\mathrm{diag}(1,-1,0)/\sqrt2$; light $L=\mathrm{diag}(1,1,0)/\sqrt2$ and strange $S=\mathrm{diag}(0,0,1)$, as in the baseline BaBar amplitude; $\omega=\cos\theta\,L-\sin\theta\,S$, $\phi=\sin\theta\,L+\cos\theta\,S$. The notation $r_i$ is the one proposed (awaiting approval) in the operator record.

## Inputs

- Vertex of each mode, linear in the vector matrix $V$:
$$T(r)=\mathrm{tr}(V[P_1,P_2])+r_1\,\mathrm{tr}(V\{[P_1,P_2],X\})+r_2\,\mathrm{tr}\big(V(P_1XP_2-P_2XP_1)\big)+r_3\,\mathrm{tr}(V)\,\mathrm{tr}([P_1,P_2]X),$$
  each $VPP$ coupling multiplied by $T(r)/T(0)$: [joint-vector-isospin-breaking-vpp-operators.md](joint-vector-isospin-breaking-vpp-operators.md) (Step 2).
- BaBar amplitude and its $s$-channel terms `rs`, `lls`: `_babar_blocks_full` of `scripts/fit_babar_kkpi_global_vector_sector.py`, [babar-kkpi-global-vector-sector-fit-baseline.md](babar-kkpi-global-vector-sector-fit-baseline.md) [@BaBar:2007ceh].
- Fitted $\omega$–$\phi$ generators `sector.generators[1]`, `[2]` at the accepted endpoint `basin_probes/even_shared_widths` ([joint-vector-fit-endpoint-comparison.md](joint-vector-fit-endpoint-comparison.md)); the same $\phi$ matrix used for the $\phi$ rows in [joint-vector-isospin-breaking-vpp-likelihood-option.md](joint-vector-isospin-breaking-vpp-likelihood-option.md).
- $X=(0.026272,\,0.046960,\,0.943684)$ GeV$^2$ from `data/vector_sector/alda_2026nrd_arxiv_v1.json` ([alda-neutral-amplitude-inputs.md](alda-neutral-amplitude-inputs.md)) [@Alda:2026nrd].

## Assumptions (labelled)

1. **Vertex only**, as in the operator record: propagators, masses and widths of $\rho,\omega,\phi$ are unchanged.
2. **Ground states only (assumption).** The $r_i$ are coefficients of operators built with the ground-state nonet field. The factors multiply the $\rho(770)$ pole of `rs` and the $\omega(782)$, $\phi(1020)$ poles of `lls`; the $\beta$-weighted excited $\rho$, $\omega$ terms of the baseline are not changed, because their couplings are separate parameters not fixed by these operators. This is the same choice as the default `'ground_only'` treatment of the $K^*$ terms; no `'whole_propagator'` alternative is provided for the $s$-channel.
3. **Production on the light end** (inherited from the baseline): the $\pi^0$ production vertices of $K^+K^-\pi^0$ couple to $L$ only; the $\omega$–$\phi$ propagator is closed on the $K^+K^-$ vertex. The $r_i$ do not change the production vertices (odd-parity, no $VPP$ vertex).
4. **Charged $\rho$ width.** In $K_SK\pi$ the baseline uses the $\rho^0$ width for the charged $\rho$ (isospin limit); unchanged.

## Where the $s$-channel vertices sit (derived from the baseline code)

- $K^+K^-\pi^0$: $s=m^2(K^+K^-)$. `rs` is the $\rho^0(s)$ propagator with the $\rho^0\to K^+K^-$ vertex. `lls = ll - sqrt2 ls` is the light-produced neutral propagator closed on the $K^+K^-$ vertex: $\mathrm{ll}=\sum_nR_{Ln}^2D_n$, $\mathrm{ls}=\sum_nR_{Ln}R_{Sn}D_n$ with $R_{Ln}=\mathrm{tr}(V_nL)$, $R_{Sn}=\mathrm{tr}(V_nS)$ and $D_n$ the $\omega,\phi$ pole propagators. The coefficient $-\sqrt2$ is $T_S(0)/T_L(0)$ for $K^+K^-$ (checked, Step 3).
- $K_SK^\pm\pi^\mp$: $s=m^2(K_SK^\pm)$ (masses $0.497611$, $0.493677$ GeV in `phase_points`) is a charged pair; `rs` is the charged $\rho$ with $\rho^+\to K^+\bar K^0$ (or $\rho^-\to K^-K^0$). The $K_S$ is the $\bar K^0$ ($K^0$) component, since $K^+K^0\pi^-$ has strangeness 2. No $\omega,\phi$ appear.
- Every appearance of `rs`, `lls` in the `direct` and `exchange` terms carries exactly one such vertex; `rq` and the source terms at $q$ carry none.

## Method and result

**Step 1–2 (SymPy).** Factors $f=T(r)/T(0)$; all equal 1 at $r=0$ (asserted):
$$f_{\rho^0\to K^+K^-}=1+2X_ur_1+X_sr_2,\qquad f_{\rho^0\to K^0\bar K^0}=1+2X_dr_1+X_sr_2,$$
$$f_{\rho^\pm\to K\bar K}=1+(X_u+X_d)\,r_1+X_sr_2,$$
$$f_L=1+2X_ur_1+X_sr_2+2(X_u-X_s)\,r_3,\qquad f_S=1+2X_sr_1+X_ur_2+(X_s-X_u)\,r_3\quad(K^+K^-),$$
and for the physical states
$$f_{\omega\to K^+K^-}=\frac{\cos\theta\,f_L/\sqrt2+\sin\theta\,f_S}{\cos\theta/\sqrt2+\sin\theta},\qquad f_{\phi\to K^+K^-}=\frac{\sin\theta\,f_L/\sqrt2-\cos\theta\,f_S}{\sin\theta/\sqrt2-\cos\theta}.$$
$\rho$ factors contain no $r_3$ ($\mathrm{tr}\,\rho=0$).

**Step 3.** $T_S(0)/T_L(0)=-\sqrt2$ for $K^+K^-$, reproducing the baseline `lls`.

**Step 4 (SymPy, all $\theta$, $r$).** Physical-basis and light/strange assemblies agree:
$$\sqrt2\sum_{n=\omega,\phi}R_{Ln}\,T_n(r)\,D_n=f_L\,\mathrm{ll}-\sqrt2\,f_S\,\mathrm{ls},\qquad \sqrt2\,T_n(0)=R_{Ln}-\sqrt2R_{Sn},$$
and at $r=0$ both equal the baseline $\mathrm{ll}-\sqrt2\,\mathrm{ls}$. Hence the substitution in the likelihood is
$$\texttt{rs}\to\texttt{rs}+(f_\rho-1)\,D_\rho^{\rm ground},\qquad \texttt{lls}\to\texttt{lls}+\sum_nR_{Ln}(R_{Ln}-\sqrt2R_{Sn})(f_n-1)D_n,$$
with $f_\rho=f_{\rho^0\to K^+K^-}$ in $K^+K^-\pi^0$ and $f_{\rho^\pm\to K\bar K}$ in $K_SK\pi$.

**Step 5 (charge structure, SymPy).** $f_{\rho^+\to K^+\bar K^0}=f_{\rho^-\to K^-K^0}$ ($C$); the charged-$\rho$ factor is the mean of the two $\rho^0$ factors; for $X_u=X_d$ all three $\rho$ factors coincide.

**Step 6 (numerical).** The fitted angle is $\theta=0.0952440$ rad, read from `sector.generators[2]` of `GlobalLikelihood.unpack` at `basin_probes/even_shared_widths` (both generators lie in the $L,S$ plane to $10^{-14}$); the provisional fixed angle of the scan is $0.0824135$. At $r=(-0.350,1.855,0.729)$ GeV$^{-2}$:

| Mode | $\rho^0\to K^+K^-$ | $\rho^0\to K^0\bar K^0$ | $\rho^\pm\to K\bar K$ | $f_L$ | $f_S$ | $\omega\to K^+K^-$ | $\phi\to K^+K^-$ |
|---|---|---|---|---|---|---|---|
| $T(r)/T(0)$ | $2.7321$ | $2.7177$ | $2.7249$ | $1.3946$ | $1.0569$ | $1.3544$ | $1.0325$ |

The $\phi\to K^+K^-$ value equals the factor of the $\phi$ rows of the likelihood option. At this point the $\rho\to K\bar K$ couplings are multiplied by about $2.7$ (the $r_2X_s$ term), while the $K^*\to K\pi$ factors are about $0.7$: the $SU(3)$-breaking side effect of the operator record, now entering the BaBar rows.

## Verification (run 2026-11-09)

- Script assertions pass (closed forms, $r=0$ limit, $-\sqrt2$ ratio, basis equivalence, charge structure, generator plane).
- Unit tests in `scripts.test_joint_vector_isospin_breaking_vpp_option`: the option's numpy factors equal the SymPy factors to $10^{-12}$ at $r=(0.8,-1.3,2.1)$ GeV$^{-2}$ with the fitted angle; the cached `rs`, `lls` of the likelihood equal an independent light/strange assembly to $10^{-12}$.

Not checked: any fit or refit; high quadrature orders at nonzero $r$; the size of the BaBar $\chi^2$ change.
