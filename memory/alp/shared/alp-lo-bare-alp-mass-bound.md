---
title: Lower bound on the bare ALP mass from vacuum stability at LO
status: derived
created: 2026-11-12
updated: 2026-11-12
bibtex_keys: [Alda:2026nrd, ParticleDataGroup:2024cfk]
provenance: analytical-derivation
---

# Lower bound on the bare ALP mass from vacuum stability at LO

## Assumptions and inputs

- LO blocks $K$, $M^2$ of the [finite-mixing record](finite-alp-meson-mixing.md) [@Alda:2026nrd], with $\epsilon=f_\pi/f_a$ (notation pending approval, see the [weak-mixing record](alp-lo-weak-mixing-expansion.md)).
- The dressed mass function $\widehat M(x)$, $x=a/f_a$, and the invariants $c_q^{AP}=c_q^A-c_q^P$, $c_q^{VS}=c_q^V-c_q^S$, $c_{GA}$ of the [dressed-mass record](alp-dressed-mass-first-order-invariants.md) and the [dressed topological record](alp-dressed-topological-invariant.md). Here $\widehat M_n\equiv d^n\widehat M/dx^n|_{x=0}$.
- $m_{a,0}^2$ multiplies a term that contains no quark field, so chiral redefinitions of the quarks do not change it.
- The quark masses are positive and $m_0^2\ge0$, so $M^2_{PP}>0$. Flavour-diagonal sources $\chi_0=\mathrm{diag}(\chi_u,\chi_d,\chi_s)$; a d-s entry is covered by the same formula (symbolic), and charge-violating sources only at random exact rational points.
- Constant-generator redefinitions were tested explicitly (noncommuting Hermitian vector and axial generators, complex d-s entries). For a general $g_{R,L}(x)$ the invariance of the final statement follows from the invariance of the minimum of the potential over $U$ under $U\to g_RUg_L^\dagger$, but a dedicated test was not made.

## Result

**Frame-independent form.** With $N_f=3$,
$$
m_{a,0}^2\;\ge\;\frac{f_\pi^2}{f_a^2}\left[B_0\,\mathrm{Re}\,\mathrm{tr}\,\widehat M_2+2\,\mathrm{tr}(c_q^{AP}\chi_0c_q^{AP})-\frac{2\,(m_0^2/N_f)\,\big(c_{GA}-\mathrm{tr}\,c_q^{AP}\big)^2}{1+(m_0^2/N_f)\sum_q\chi_q^{-1}}\right].
$$
Equivalently, $B_0\,\mathrm{Re}\,\mathrm{tr}\,\widehat M_2+2\,\mathrm{tr}(c_q^{AP}\chi_0c_q^{AP})=B_0\,d^2\,\mathrm{tr}\big(\widehat M^\dagger\widehat M\big)^{1/2}/dx^2|_0$, the second derivative of the sum of the singular values of $M(x)$ (equal to those of $\widehat M(x)$), and $2\,\mathrm{tr}(c_q^{AP}\chi_0c_q^{AP})=2B_0\sum_{q,q'}|\widehat P_{1,qq'}|^2/(m_q+m_{q'})$ with $c_{GA}-\mathrm{tr}\,c_q^{AP}=\widehat\Theta_1+\tfrac12\mathrm{tr}(m_q^{-1}\widehat P_1)$. A d-s entry adds $2(\chi_d+\chi_s)|c^{AP}_{ds}|^2$ to the first two terms. The combination $c_q^{VS}$ does not enter. The bound is necessary and sufficient for $M^2\ge0$, and equals the condition that the second derivative of the vacuum energy of the LO chiral Lagrangian, minimized over the mesons, plus $m_{a,0}^2$, is nonnegative.

**What is not fixed by first-order invariants.** The real number $\mathrm{Re}\,\mathrm{tr}\,\widehat M_2$ is the only second-order datum that enters. $\widehat\Theta_2$ and the imaginary part of $\mathrm{tr}\,\widehat M_2$ do not.

**Linear-spurion convention** (the convention of the finite-mixing record: the quark spurion is linear in $a$, and $C_{R,L}$ constant). Then $\mathrm{Re}\,\mathrm{tr}\,\widehat M_2=4\,\mathrm{tr}(c_q^Pm_qc_q^P)-4\,\mathrm{tr}(c_q^{AP}m_qc_q^{AP})$, independent of $c_q^V$, $c_q^S$ and of the $x$-derivatives of the connections, and the bound reduces to
$$
m_{a,0}^2\;\ge\;\frac{f_\pi^2}{f_a^2}\left[2\,\mathrm{tr}(c_q^P\chi_0c_q^P)-\frac{2\,(m_0^2/N_f)\,\big(\mathrm{tr}\,c_q^P+c_G\big)^2}{1+(m_0^2/N_f)\sum_q\chi_q^{-1}}\right].
$$
This expression is **not** frame independent by itself. Under a redefinition with axial generator $G_A$ the sources shift as $c_q^{P}\to c_q^P+G_A$, $c_q^A\to c_q^A+G_A$, $c_G\to c_G-\mathrm{tr}\,G_A$. Then $c_q^{AP}$, $c_q^{VS}$, $c_{GA}$ and $\mathrm{tr}\,c_q^P+c_G$ are unchanged, but $2\,\mathrm{tr}(c_q^P\chi_0c_q^P)$ changes by $2\,\mathrm{tr}(G_A\chi_0G_A)+2\,\mathrm{tr}(G_A\{c_q^P,\chi_0\})$. The change is compensated by the $a^2$ term that the transformed spurion acquires, which is $\mathrm{Re}\,\mathrm{tr}\,\chi'_{(2)}=2\,\mathrm{tr}(c_q^P\chi_0c_q^P)-2\,\mathrm{tr}(c_q^{P\prime}\chi_0c_q^{P\prime})$. In the transformed frame this term must be kept in $M^2_{aa}$. A vector generator $G_V$ alone changes nothing. Physically, the positive term appears because a linear spurion gives singular values $m_q\sqrt{1+4c_q^2x^2}$, whose moduli grow at order $a^2$.

## Derivation (linear-spurion convention; the frame-independent form follows by the replacement stated above)

1. The topological part of $M^2$ is $(m_0^2/2N_f)\,u\,u^{\mathsf T}$ with $u=(t_i,\,2\epsilon c_G)$, $t_i=\mathrm{tr}\,\lambda_i$. It is positive semidefinite.
2. Eliminate the mesons: $M^2\ge0$ iff $M^2_{aa}-M^2_{aP}(M^2_{PP})^{-1}M^2_{Pa}\ge0$, with $M^2_{aa}=m_{a,0}^2+2\epsilon^2m_0^2c_G^2/N_f$.
3. In the quark basis the meson block is $\mathrm{diag}(\chi_q)/2+(m_0^2/2N_f)\mathbb 1\mathbb 1^{\mathsf T}$ and the mixing is $\epsilon\,(-c_q\chi_q+m_0^2c_G/N_f)$. Invert the meson block with the Sherman-Morrison formula and collect the terms.

## Limits and special points (linear-spurion convention)

- **Sufficient bound:** dropping the (positive semidefinite) topological square gives $m_{a,0}^2\ge2\epsilon^2\,\mathrm{tr}(c_q^P\chi_0c_q^P)$, which is never smaller than the exact bound.
- **$m_0^2\to\infty$:** the bracket tends to $2\sum_qc_q^2\chi_q-2(\mathrm{tr}\,c^P+c_G)^2/\sum_q\chi_q^{-1}$ and decreases monotonically with $m_0^2$.
- **$c_G=0$:** the bracket is positive for every nonzero diagonal $c_q^P$, and for $m_0^2\to\infty$ it vanishes only if $c_q\propto1/\chi_q$.
- **Negative bound:** the bracket is smallest, $-2m_0^2c_G^2/N_f$, at $c_q=m_0^2c_G/(N_f\chi_q)$ (vanishing mass mixing). For $c_q^P=0$ it is $-2(m_0^2/N_f)c_G^2/(1+(m_0^2/N_f)\sum_q\chi_q^{-1})$, so $m_{a,0}^2$ may then be zero or negative. In that case, with $m_{a,0}^2=0$, the lightest pole divided by the $\epsilon^2$ bracket tends to $1$ as $\epsilon\to0$.
- **Flavour-violating d-s entry:** it adds $2(\chi_d+\chi_s)|c^P_{ds}|^2$ to the bracket.

## Meaning

For $f_a\gg f_\pi$ the bound is a correction of relative size $f_\pi^2/f_a^2$ to the ALP mass squared, and it does not restrict ALPs lighter than the hadronic scale. It becomes comparable to the pion mass only for $f_a$ a few times $f_\pi$ (order-of-magnitude values below, linear-spurion convention). Since the sources of Alda:2026nrd are defined up to the redefinitions above, the statement "the bracket is positive" has physical content only in the frame-independent form: its value depends on $\mathrm{Re}\,\mathrm{tr}\,\widehat M_2$, which a different choice of spurion convention (that is, a different $a^2$ term) changes. A coupling with $c_q^{AP}=0$ and $\widehat M_2=0$ leaves only the non-positive topological term.

## Order-of-magnitude illustration (assumptions, not a result)

Inputs: $f_\pi=0.0924$ GeV; $2B_0m_{ud}=m_{\pi^0}^2$ with $m_{\pi^0}=0.1349768$ GeV; $m_u/m_d=0.474$ and $m_s/m_{ud}=27.33$ (rounded values, ParticleDataGroup:2024cfk); $m_0^2=m_\eta^2+m_{\eta'}^2-2m_K^2$ at LO, $=0.726$ GeV$^2$ with $m_K=0.4956$ GeV; $N_f=3$. This gives $\chi_q=(0.0117,0.0247,0.498)$ GeV$^2$.

| $c_q^P$, $c_G$ | bracket (GeV$^2$) | $\sqrt{\text{bound}}$ at $f_a=1,10,100$ GeV | $f_a/f_\pi$ at which bound $=m_\pi^2$ |
|---|---|---|---|
| $\mathbb 1$, $0$ | $0.932$ | $89,\ 8.9,\ 0.89$ MeV | $7.2$ |
| $\mathbb 1$, $1$ | $0.826$ | $84,\ 8.4,\ 0.84$ MeV | $6.7$ |
| $\mathrm{diag}(1,-\tfrac12,-\tfrac12)$, $1$ | $0.270$ | $48,\ 4.8,\ 0.48$ MeV | $3.8$ |
| $0$, $1$ | $-0.015$ | none (negative) | — |

The bracket is dominated by $\chi_s(c_s^P)^2$. These values illustrate the scaling for order-one couplings; they have no meaning outside those inputs. Values quoted by the diagnostic script for its own test inputs ($f_a=0.2$ GeV, so $\epsilon=23/50$) are not physical bounds.

## Script

[`scripts/derive_alp_lo_stability_bound.py`](../../../scripts/derive_alp_lo_stability_bound.py) verifies the independence of $M^2$ from $c_q^A$, the Schur complement in the generator basis against the Sherman-Morrison form, the limits above, and the sign of the smallest eigenvalue of $M^2$ just above and just below the bound at random points. All checks pass.

[`scripts/derive_alp_lo_stability_bound_invariance.py`](../../../scripts/derive_alp_lo_stability_bound_invariance.py) verifies the frame independence: $K'=T^{\mathsf T}KT$ and $M^{2\prime}=T^{\mathsf T}M^2T$ exactly under constant-generator redefinitions (the $a^2$ spurion term kept), invariance of $m_{a,0}^2$ and of the Schur complement, the source shifts above, the failure of the linear-spurion bracket in the other frames, identical $\widehat M_{0,1,2}$, $\widehat\Theta_1$ and bound in three frames (linear spurion, generic generators, frame without kinetic mixing; 5 and 9 mesons), a sign change of the smallest eigenvalue of $M^2$ at the invariant threshold in each frame, and agreement of its second derivative with that of the nonlinear potential minimized over $U(3)$ (relative $10^{-6}$).

```bash
.venv/bin/python scripts/derive_alp_lo_stability_bound.py
.venv/bin/python scripts/derive_alp_lo_stability_bound_invariance.py
```
