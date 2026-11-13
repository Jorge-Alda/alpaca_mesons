---
title: Scalar and pseudoscalar quark-mass sources of an ALP carried by a Peccei-Quinn phase
status: derived
created: 2026-11-13
updated: 2026-11-13
bibtex_keys: [DiLuzio:2026poh, Alda:2026nrd]
provenance: analytical-derivation
---

# Quark-mass sources from a Peccei-Quinn phase

## Conventions and units

Natural units, $q=(u,d,s)$, $x=a/f_a$ real and dimensionless. The partonic mass term follows the [general chiral-redefinition record](alp-general-chiral-redefinition-invariants.md):
$$
\mathcal L\supset-\bar q_RXq_L-\bar q_LX^\dagger q_R=-\bar q\,(S-iP\gamma_5)\,q,\qquad X=S+iP,\quad X(0)=m_q,
$$
with $S$, $P$ Hermitian and of mass dimension one. Taylor coefficients are derivatives at $x=0$: $F(x)=\sum_k\frac{x^k}{k!}F_k$.

## Inputs

- Mass term of DiLuzio:2026poh [@DiLuzio:2026poh, Eqs. (3), (7)], $-m_q(\sqrt2\varphi)^{n_q}\bar q_Lq_R+\text{h.c.}$, with $\varphi=\frac{1+\sigma}{\sqrt2}e^{ia}$ (Eq. (9)); see [model and field normalization](../../extracts/DiLuzio:2026poh/model-lagrangian-and-field-normalization.md), [angular-mode parametrization](../../extracts/DiLuzio:2026poh/angular-mode-parametrization.md) and [angular-mode couplings](../../extracts/DiLuzio:2026poh/angular-mode-couplings.md).
- $m_q=\mathrm{diag}(m_u,m_d,m_s)$ real and positive; PQ charges $N=\mathrm{diag}(n_u,n_d,n_s)$ real (integers in the model, kept symbolic).
- Source transformations, transport $\partial_xW_\eta=iC_\eta W_\eta$, $W_\eta(0)=\mathbf 1$, and the dressed source $\widehat M=W_R^\dagger MW_L$: [chiral source transformations](alp-chiral-source-transformations.md).
- First-order parametrization $S_1=i[c_q^S,m_q]$, $P_1=\{c_q^P,m_q\}$ of Alda:2026nrd [@Alda:2026nrd]: [Appendix-A extract](../../extracts/Alda:2026nrd/appendix-a-expansion.md), [dressed-mass first-order record](alp-dressed-mass-first-order-invariants.md).

## Assumptions

- **Assumption:** the radial mode is set to zero, $\sigma=0$.
- **Assumption:** $x$ is identified with the paper's dimensionless phase $a=a_{\rm can}/f_\varphi$, i.e. $f_a=f_\varphi$. With the light-axion definition $f_a=f_\varphi/\sum_qn_q$ [@DiLuzio:2026poh, Eq. (53)], every $N$ below must be read as $N/\sum_qn_q$.
- In this basis the phase is carried entirely by the mass term: $C_R=C_L=0$ and $\Theta$ does not depend on $x$.

## Method

1. At $\sigma=0$, $(\sqrt2\varphi)^{n_q}=e^{in_qx}$, so the paper's $\bar q_Lq_R$ coefficient is $m_qe^{iNx}$. Comparing with $-\bar q_LX^\dagger q_R$ gives $X^\dagger=m_qe^{iNx}$, hence
$$
X(x)=m_q\,e^{-iNx},\qquad S(x)=m_q\cos(Nx),\qquad P(x)=-m_q\sin(Nx).
$$
2. Sign check at Dirac level, with $\bar q_{R}q_{L}=\bar q\frac{1-\gamma_5}2q$ and $\bar q_Lq_R=\bar q\frac{1+\gamma_5}2q$: the paper's term is $-m_q[\cos(Nx)\bar qq+i\sin(Nx)\bar q\gamma_5q]$, and $-\bar q(S-iP\gamma_5)q$ has coefficients $-S$ for $\bar qq$ and $+iP$ for $\bar q\gamma_5q$. Both coefficients agree.
3. Differentiate $S$ and $P$ at $x=0$.
4. Match the first order to $S_1=i[c_q^S,m_q]$, $P_1=\{c_q^P,m_q\}$; all matrices are diagonal.
5. Invariance check: apply $g_R=e^{icx}$, $g_L=\mathbf 1$ with Hermitian constant $c$. Then $M'=g_RM$, $C_R'=i\,g_R\partial_xg_R^\dagger=c$, $W_R'=e^{icx}$, $W_L'=\mathbf 1$. The script builds the Taylor coefficients of $M'$ and $C_R'$ by the Leibniz rule, $W_R'$ by the transport recursion $W_{k+1}=i\sum_j\binom kjC_jW_{k-j}$, and compares $\widehat M_k=\sum_j\binom kj(W_R'^\dagger)_jM'_{k-j}$ with $M_k$.

## Result

For all $j\ge0$:
$$
S_{2j}=(-1)^j\,m_qN^{2j},\qquad S_{2j+1}=0,\qquad P_{2j}=0,\qquad P_{2j+1}=-(-1)^j\,m_qN^{2j+1}.
$$
Explicitly, $S_0=m_q$, $P_1=-m_qN$, $S_2=-m_qN^2$, $P_3=m_qN^3$. At first order $c_q^S=0$ (any diagonal $c_q^S$ gives the same $S_1=0$) and
$$
c_q^P=-\frac N2 .
$$
With $C_A=0$ this gives $c_q^{AP}=N/2$ and $\widehat P_1=-\{c_q^{AP},m_q\}=-Nm_q$, consistent with the [dressed-mass first-order record](alp-dressed-mass-first-order-invariants.md).

## Verification

[`derive_alp_pq_phase_mass_sources.py`](../../../scripts/derive_alp_pq_phase_mass_sources.py) asserts, with exact SymPy arithmetic:

- the closed forms for $k=0,\dots,6$ (default; `--max-order` changes $K$);
- the Dirac-level sign mapping between the paper's form and $-\bar q(S-iP\gamma_5)q$, for both $\bar qq$ and $\bar q\gamma_5q$;
- the first-order matching $S_1=0$, $P_1=-Nm_q$, $c_q^P=-N/2$;
- $C_R'_0=c$, $C_R'_{k\ge1}=0$ and $\widehat M_k=M_k$ for $k\le6$, for symbolic diagonal $c$ and for a non-diagonal rational Hermitian $c$;
- the closed forms, sign mapping and first-order matching with exact rational masses and $N=\mathrm{diag}(1,0,0)$, $\mathrm{diag}(2,-1,1)$.

The unit tests in `scripts/test_derive_alp_pq_phase_mass_sources.py` repeat these checks and confirm that a wrong sign of $P_1$ is rejected. All checks pass.

```sh
.venv/bin/python scripts/derive_alp_pq_phase_mass_sources.py
.venv/bin/python -m unittest scripts.test_derive_alp_pq_phase_mass_sources
```
