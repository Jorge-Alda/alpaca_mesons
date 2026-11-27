---
title: ALP decays to octet baryon pairs in the only-up model with pointlike couplings
status: mixed
created: 2026-11-19
updated: 2026-11-19
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh, ParticleDataGroup:2024cfk]
provenance: analytical-derivation
---

# $a\to B_1\bar B_2$ in the only-up model

**Result.** For the vertex $\mathcal L_{aB}=\frac{\partial_\mu a_{\rm phys}}{f_a}\sum_{B_1,B_2}\bar B_1\gamma^\mu(c^V_{B_1B_2}+c^A_{B_1B_2}\gamma_5)B_2$,
$$
\Gamma(a\to B_1\bar B_2)=\frac{\lambda^{1/2}(m_a^2,m_{B_1}^2,m_{B_2}^2)}{8\pi f_a^2m_a^3}
\Big[|c^A_{B_1B_2}|^2(m_{B_1}+m_{B_2})^2\big(m_a^2-(m_{B_1}-m_{B_2})^2\big)
+|c^V_{B_1B_2}|^2(m_{B_1}-m_{B_2})^2\big(m_a^2-(m_{B_1}+m_{B_2})^2\big)\Big],
$$
which reduces to Eq. (ALP_baryon_pair_width) of Alda:2026nrd for $B_1=B_2$. In the only-up model the couplings to the canonically normalized ALP are the Section V couplings at the physical pole mass divided by $\sqrt{N(m_a^2)}$; $c^V_{B_1B_2}=0$ for all neutral pairs. For $\kappa_u\le1$ the open channels are $p\bar p$, $n\bar n$, $\Lambda\bar\Lambda$, $\Lambda\bar\Sigma^0$ and $\Sigma^0\bar\Lambda$.

## Conventions and units

Conventions of Alda:2026nrd [@Alda:2026nrd]: $\bm\Phi=\sum_kP_k\bm\lambda_k$ with $\langle\bm\lambda_k\bm\lambda_l\rangle=2\delta_{kl}$, $\bm U=\bm u^2$, $\bm u=e^{i\bm\Phi/(2f_\pi)}$, $\bm u_\mu=i(\bm u^\dagger\partial_\mu\bm u-\bm u\partial_\mu\bm u^\dagger)$, meson kinetic term $\tfrac14\langle\partial\bm\Phi\partial\bm\Phi\rangle$, ALP kinetic term $\tfrac12(\partial a)^2$. Octet matrix $\bm{\mathcal B}=\sum_B\bm\lambda_{B}B$ of Section V. Masses, $f_a$, $f_\pi$ and widths in GeV; $c^V$, $c^A$ dimensionless. $\lambda(x,y,z)=x^2+y^2+z^2-2xy-2yz-2zx$. The final state $B_1\bar B_2$ receives only the $(B_1,B_2)$ term of the sum; $c_{B_2B_1}=c_{B_1B_2}^*$.

## Assumptions

- Tree level, LO $U(3)$ meson Lagrangian and LO octet baryon Lagrangian with $D_B,F_B,S_B,\widetilde S_B$; radial field frozen ($\sigma=0$); no expansion in $f_\pi/f_a$.
- **Approximation (user choice): baryon form factors set to one.** The vertex is the pointlike LO one. For $m_a\simeq2$ GeV this lies outside the domain of the chiral expansion; Alda:2026nrd uses a form-factor continuation there (Sec. baryon_GeV_continuation), which is not included here.
- $\Sigma^0$–$\Lambda$ mixing from isospin breaking is not included; physical isospin-split baryon masses are used in the phase space.
- The ALP-like pole is simple and $\mathbb M_{m_a^2}$ is invertible (true in the region where baryon channels are open).

## Inputs

- Only-up model, $f_a=\sqrt2m_u/\kappa_u$, quark masses, $f_\pi$, $B_0$, $m_1^2$: [only-up inputs](../meson-mixing/finite-mixing/spectrum-and-poles/quadratic-mixing/only-up-finite-mixing-inputs.md) and `scripts/implementation/only_up_inputs.py`.
- Source coefficients $\widehat P_{01}=-m_u\,\mathrm{diag}(1,0,0)$, $\widehat S_{02}=-\tfrac{m_u}{2}\mathrm{diag}(1,0,0)$, $c_q^{AP}=\mathrm{diag}(1/2,0,0)$, $c_{GA}=0$; $\mathbb M_s$, pole equation, $N(m^2)$ and the basis without ALP–meson mixing: [radial-mode record](../radial-mode/couplings-and-decays/finite-mixing/murayama-sigma-two-pseudoscalar-decays.md), [quartic-vertex record, section Check](../three-pseudoscalar/finite-mixing/murayama-unmixed-basis-quartic-vertex.md).
- Physical ALP mass: `spectrum` of `scripts/implementation/only_up_decays.py` (finite-mixing pole).
- Section V couplings and building blocks: [Section V extract](../../extracts/Alda:2026nrd/baryon-couplings.md), [finite-mixing baryon couplings](finite-mixing/finite-mixing-baryon-couplings.md).
- $D_B=0.815$, $F_B=0.434$, $S_B=-0.417$ (Alda:2026nrd Eq. baryon_LECs); $\widetilde S_B=0.025$, central value of Eq. stilde_value (naive uncertainty $0.031$, conservative $0.062$) [@Alda:2026nrd].
- Baryon masses (GeV), PDG 2024 [@ParticleDataGroup:2024cfk]: $p$ 0.93827209, $n$ 0.93956542, $\Lambda$ 1.115683, $\Sigma^+$ 1.18937, $\Sigma^0$ 1.192642, $\Sigma^-$ 1.197449, $\Xi^0$ 1.31486, $\Xi^-$ 1.32171. These values were entered from the PDG summary tables and are not yet in a memory extract.

## Method

**Step 1: invariant source.** With $\chi_0=2B_0m_q$, $\mathcal C=-\tfrac12\{\chi_0,c_q^{AP}\}-\tfrac{m_1^2}{3}c_{GA}\mathbf1=-B_0m_u\,\mathrm{diag}(1,0,0)=B_0\widehat P_{01}$ (SymPy). The $\mathbb M_a$ operator of Section V coincides with $\mathbb M_s$ of the radial-mode record at $s=m_a^2$, so $\mathbb M_a^{-1}(\mathcal C)=B_0\mathbb M^{-1}_{m_a^2}(\widehat P_{01})$.

**Step 2: building blocks in the basis without ALP–meson mixing.** In the polar fields the ALP enters only through the quark mass source, so the original baryon Lagrangian has no direct ALP coupling. After $\bm U=\bm g_R\bm U'\bm g_R$, $\bm g_R=\exp\!\big(\tfrac{iaB_0}{f_a}\mathbb M^{-1}_{m_a^2}(\widehat P_{01})\big)$, the terms linear in $a$ without mesons have $\bm u=\bm g_R$, giving (SymPy for diagonal shift; numerical check with a non-diagonal Hermitian shift)
$$
\bm u_\mu\big|_a=-\frac{2\partial_\mu a}{f_a}\mathbb M_a^{-1}(\mathcal C),\qquad
\bm\Gamma_\mu\big|_a=0,\qquad
\partial_\mu X_{TC}\big|_a=-\frac{2\partial_\mu a}{f_a}\langle\mathbb M_a^{-1}(\mathcal C)\rangle,
$$
the last from $X_{TC}=-\langle\bm\Phi\rangle/f_\pi$ with $\bm\Phi=\bm\Phi'+\tfrac{2f_\pi}{f_a}\mathbb M_a^{-1}(\mathcal C)\,a$. These are Eq. (baryon_ALP_building_blocks) of Alda:2026nrd with $c_q^{VS}=0$, $c_{GA}=0$, evaluated at the physical $m_a$. The further linear change $a'=a''-\sum_iw_i^{(a)}P_i'/N(m_a^2)$ of the draft adds only terms with a meson field and does not change this vertex.

**Step 3: normalization (derived, included).** Inserting $\bm\Phi=\bm\Phi'+\tfrac{2f_\pi}{f_a}\mathbb M_a^{-1}(\mathcal C)a$ into $\tfrac14\langle\partial\bm\Phi\partial\bm\Phi\rangle+\tfrac12(\partial a)^2$ gives the coefficient of $\tfrac12(\partial a)^2$ (SymPy)
$$
K'_{aa}=1+\frac{2f_\pi^2}{f_a^2}\big\langle[\mathbb M_a^{-1}(\mathcal C)]^2\big\rangle=N(m_a^2).
$$
The ALP column of the quadratic form vanishes at $p^2=m_a^2$, so the propagator residue is $1/N(m_a^2)$ and $a_{\rm phys}=\sqrt{N(m_a^2)}\,a''$. All couplings are divided by $\sqrt{N(m_a^2)}$. This is not a small correction to be dropped silently: it is $O(f_\pi^2/f_a^2)\times\langle[\mathbb M_a^{-1}(\mathcal C)]^2\rangle$; numerically $N-1\approx2.5\times10^{-3}$–$3.7\times10^{-3}$ where the channels are open, while on the full grid its largest value is $N=1.996$ at $\kappa_u=0.218$ ($m_a=0.486$ GeV, near the $\eta$ crossing).

**Step 4: couplings.** Substitution into the Section V Lagrangian gives
$$
c^V_{B_1B_2}=-\frac{\langle c_q^{VS}[\bm\lambda_{B_1}^\dagger,\bm\lambda_{B_2}]\rangle}{\sqrt{N(m_a^2)}}=0,\qquad
c^A_{B_1B_2}=-\frac{(D_B+F_B)\langle\bm\lambda_{B_1}^\dagger X\bm\lambda_{B_2}\rangle+(D_B-F_B)\langle\bm\lambda_{B_1}^\dagger\bm\lambda_{B_2}X\rangle+(S_B+\widetilde S_B)\langle X\rangle\langle\bm\lambda_{B_1}^\dagger\bm\lambda_{B_2}\rangle}{\sqrt{N(m_a^2)}},
$$
with $X=\mathbb M_a^{-1}(\mathcal C)$. This is the pole-basis result of the [finite-mixing record](finite-mixing/finite-mixing-baryon-couplings.md); the two bases agree. In the only-up model $X=\mathrm{diag}(X_u,X_d,X_s)$ with ($s=m_a^2$, SymPy)
$$
X_q=\frac{\mathcal C_{qq}-\tfrac{m_1^2}{3}\langle X\rangle}{2B_0m_q-s},\qquad
\langle X\rangle=\frac{-B_0m_u/(2B_0m_u-s)}{1+\tfrac{m_1^2}{3}\sum_{q'}(2B_0m_{q'}-s)^{-1}}.
$$
Hence, e.g., $c^A_{pp}\sqrt N=-(D_B+F_B)X_u-(D_B-F_B)X_s-(S_B+\widetilde S_B)\langle X\rangle$ and $c^A_{\Sigma^0\Lambda}\sqrt N=c^A_{\Lambda\Sigma^0}\sqrt N=-D_B(X_u-X_d)/\sqrt3$, independent of $\widetilde S_B$. Neutral flavour-conserving pairs (equal charge and strangeness): the eight diagonal pairs plus $\Lambda\bar\Sigma^0$ and $\Sigma^0\bar\Lambda$. The vector coupling vanishes for all of them because $c_q^{VS}=0$; for $\Sigma^0\Lambda$ it would vanish also for any diagonal $c_q^{VS}$, since the two generators commute.

**Step 5: amplitude and width.** With $\partial_\mu a\to-iq_\mu$, $\mathcal M=\tfrac{1}{f_a}\bar u_1\not q(c^V+c^A\gamma_5)v_2$. The spin sum $\mathrm{tr}[(\not p_1+m_1)\Gamma(\not p_2-m_2)\gamma^0\Gamma^\dagger\gamma^0]/f_a^2$, evaluated with explicit Dirac matrices in the ALP rest frame (SymPy), is
$$
\sum|\mathcal M|^2=\frac{2}{f_a^2}\Big[|c^A|^2(m_1+m_2)^2\big(m_a^2-(m_1-m_2)^2\big)+|c^V|^2(m_1-m_2)^2\big(m_a^2-(m_1+m_2)^2\big)\Big],
$$
with no $c^Vc^{A*}$ interference. On shell, $\bar u_1\not q\gamma_5v_2=(m_1+m_2)\bar u_1\gamma_5v_2$ and $\bar u_1\not qv_2=(m_1-m_2)\bar u_1v_2$. $\Gamma=|\vec p|\sum|\mathcal M|^2/(8\pi m_a^2)$ gives the Result. Near threshold the axial term behaves as $|\vec p|$ (S wave) and the vector term as $|\vec p|^3$ (P wave). No symmetry factor ($B_1$ and $\bar B_2$ are distinct).

## Numerical results

Physical-pole ALP mass from `spectrum`; $\kappa_u$ grid `logspace(-3,0,500)` of `murayama.py`.

| channel | threshold (GeV) | opens at $\kappa_u$ | open grid points |
|---|---|---|---|
| $p\bar p$ | 1.876544 | 0.800502 | 17 |
| $n\bar n$ | 1.879131 | 0.801610 | 16 |
| $\Lambda\bar\Lambda$ | 2.231366 | 0.952438 | 4 |
| $\Lambda\bar\Sigma^0$, $\Sigma^0\bar\Lambda$ | 2.308325 | 0.985377 | 2 each |
| $\Sigma^{\pm,0}\bar\Sigma^{\pm,0}$, $\Xi\bar\Xi$ | $\ge2.378740$ | closed ($m_a(\kappa_u{=}1)=2.342491$ GeV) | 0 |

Widths with $\widetilde S_B=0.025$ (range for $\widetilde S_B=0.025\pm0.031$ in brackets):

| $\kappa_u$ | $m_a$ (GeV) | $f_a$ (GeV) | $N(m_a^2)$ | channel | $c^A$ | $\Gamma$ (GeV) |
|---|---|---|---|---|---|---|
| 0.85 | 1.992105 | $3.6603\times10^{-3}$ | 1.003684 | $p\bar p$ | $-1.40\times10^{-3}$ | $1.368\times10^{-2}$ [$1.253$, $1.489$]$\times10^{-2}$ |
| | | | | $n\bar n$ | $+5.6\times10^{-4}$ | $2.180\times10^{-3}$ [$1.736$, $2.675$]$\times10^{-3}$ |
| 0.90 | 2.108874 | $3.4570\times10^{-3}$ | 1.003216 | $p\bar p$ | $-1.24\times10^{-3}$ | $1.740\times10^{-2}$ |
| | | | | $n\bar n$ | $+5.1\times10^{-4}$ | $2.889\times10^{-3}$ |
| 1.00 | 2.342491 | $3.1113\times10^{-3}$ | 1.002526 | $p\bar p$ | $-1.00\times10^{-3}$ | $2.025\times10^{-2}$ [$1.863$, $2.194$]$\times10^{-2}$ |
| | | | | $n\bar n$ | $+4.2\times10^{-4}$ | $3.550\times10^{-3}$ |
| | | | | $\Lambda\bar\Lambda$ | $+1.1\times10^{-4}$ | $1.634\times10^{-4}$ [$0.62$, $3.14$]$\times10^{-4}$ |
| | | | | $\Lambda\bar\Sigma^0$, $\Sigma^0\bar\Lambda$ | $-5.3\times10^{-4}$ | $2.484\times10^{-3}$ each |

The couplings $c^A$ are small because $X\sim B_0m_u/m_a^2\sim10^{-3}$ at $m_a\simeq2$ GeV, but $f_a\simeq3$ MeV, so $c^A/f_a=O(0.3)$ GeV$^{-1}$. $\Lambda\bar\Lambda$ is the most sensitive channel to $\widetilde S_B$ because there the $D_B$ term, $-D_B(X_u+X_d+4X_s)/3$, and the $S_B$ term partly cancel.

## Verification

[Unit tests](../../../scripts/test_derive_murayama_baryon_decays.py), 16 tests, all pass (2026-11-19):

- SymPy: $\mathcal C=B_0\widehat P_{01}$; $\bm u_\mu$, $\bm\Gamma_\mu$, $\partial X_{TC}$ from $\bm u=\bm g_R$; $K'_{aa}=N$; closed form of $X_q$; Dirac trace equals the closed $\sum|\mathcal M|^2$; diagonal width equals Alda:2026nrd Eq. (ALP_baryon_pair_width).
- Section V trace formulas against the Lagrangian evaluated with explicit $3\times3$ matrices, random non-commuting Hermitian $\mathbb M_a^{-1}(\mathcal C)$ and $c_q^{VS}$ (seed 2026), all 64 octet pairs, and hermiticity $c_{B_2B_1}=c_{B_1B_2}^*$: agreement to $10^{-13}$.
- Chiral shift with a non-diagonal Hermitian generator (finite differences of `expm`): $\bm u_\mu=-2A/f_a$, $\bm\Gamma_\mu=0$ to $10^{-8}$.
- $\mathbb M_0^{-1}(\mathcal C)$ equals Alda:2026nrd Eq. (Ma_light_limit) exactly; $\mathbb M_s^{-1}$ inverts $\mathbb M_s$ on a non-commuting source.
- **Basis independence:** couplings in the basis without mixing (Section V$/\sqrt N$) equal those obtained in the original polar fields from the canonically normalized eigenvector of the $(\pi_3,\pi_8,\eta_0,a)$ quadratic form, with $\bm u_\mu=-\partial\bm\Phi/f_\pi$, $X_{TC}=-\langle\bm\Phi\rangle/f_\pi$, at $\kappa_u=0.003,0.05,0.3,0.9,1$ (relative $10^{-8}$); eigenvalue equals `spectrum` $m_a^2$ to $10^{-9}$.
- Small mixing: $N-1\propto\kappa_u^2$; at $\kappa_u=10^{-3}$ the couplings equal the Section V expression with the light limit of $\mathbb M_a^{-1}(\mathcal C)$ to $10^{-3}$, and $c^A_{pp}$, $c^A_{nn}$ agree with Eq. (effcoupl_fc) within the expected 10% (large-$m_0$ and $m_{u,d}/m_s$ approximations of that equation).
- Width against explicit Dirac spinors and a numerical spin sum, for $p\bar p$, $\Lambda\bar\Sigma^0$, $\Sigma^0\bar\Lambda$ with complex test couplings including $c^V\ne0$, and for the model couplings at $\kappa_u=1$: relative $10^{-12}$.
- Threshold: zero width at and below threshold; $\Gamma/|\vec p|$ (axial) and $\Gamma/|\vec p|^3$ (vector) tend to finite limits; the $p\bar p$ width switches on at the computed threshold $\kappa_u$.

- Implementation: the numerical part ($D_B$, $F_B$, $S_B$, $\widetilde S_B$, baryon masses, octet generators, $\mathbb M_s^{-1}$, Section V couplings, width formula) lives in `scripts/implementation/only_up_decays.py` as `baryon_pair_width(inputs, baryons, topological_singlet=0.025)` and `baryon_pair_widths(inputs, topological_singlet=0.025, channels=None)`, with `baryons = (B_1, B_2)` for $a\to B_1\bar B_2$; `scripts/derive_murayama_baryon_decays.py` keeps the SymPy derivations and the cross-checks and imports these functions. The widths at $\kappa_u=0.85$, $0.9$, $0.95$, $1$ (central and lower $\widetilde S_B$) are bit-identical before and after the move. [Implementation tests](../../../scripts/implementation/test_only_up_decays.py) (`BaryonPairWidthTests`, 3 tests, pass 2026-11-19): table values at $\kappa_u=0.9$, $1$ to relative $10^{-12}$, zero widths below threshold ($\kappa_u=0.3$), dependence on $\kappa_u$ only through `inputs['f_a']`, $\widetilde S_B$ argument.

Open points: the baryon masses are not yet backed by a memory extract of ParticleDataGroup:2024cfk; the pointlike LO vertex at $m_a\simeq2$ GeV is a user-fixed approximation and differs from the form-factor continuation used in Alda:2026nrd.

## Commands

```sh
.venv/bin/python scripts/derive_murayama_baryon_decays.py
.venv/bin/python scripts/derive_murayama_baryon_decays.py --test
.venv/bin/python -m unittest scripts.test_derive_murayama_baryon_decays
.venv/bin/python -m unittest scripts.implementation.test_only_up_decays
```
