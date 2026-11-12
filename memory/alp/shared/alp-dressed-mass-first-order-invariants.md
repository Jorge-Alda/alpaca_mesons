---
title: First-order expansion of the dressed mass function and comparison with the Alda:2026nrd invariants
status: derived
created: 2026-11-11
updated: 2026-11-12
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# First-order expansion of the dressed mass function and comparison with the Alda:2026nrd invariants

## Scope and inputs

This record expands $\widehat M(x)=W_R^\dagger(x)M(x)W_L(x)$ to first order in $x=a/f_a$ and compares the result with the invariants $c_q^{VS}$ and $c_q^{AP}$ of Alda:2026nrd [@Alda:2026nrd]. Definitions of $W_{R,L}$, $M=S+iP$, $C_{R,L}=C_V\pm C_A$ and the restricted partonic Lagrangian are in the [general chiral-redefinition record](alp-general-chiral-redefinition-invariants.md). The source parametrization is that of the [Appendix-A extract](../../extracts/Alda:2026nrd/appendix-a-expansion.md):

$$
S(x)=m_q+ix[c_q^S,m_q]+\cdots,\qquad
P(x)=x\{c_q^P,m_q\}+\cdots,\qquad
C_V(0)=c_q^V,\quad C_A(0)=c_q^A,
$$

with $m_q$, $c_q^{V,A,S,P}$ Hermitian. The identification of $C_{V,A}(0)$ with $c_q^{V,A}$ and of $M$ with Alda:2026nrd's $\chi/(2B_0)$ was checked against the vielbein of Eq. (A.13): for $\Phi=0$, $u_\mu=2(\partial_\mu x)C_A$, and the term linear in $\phi$ is $-2i(\partial_\mu x)[\phi,C_V]$. Only the order-$x$ coefficient is studied; the anomalous $\Theta$ shift is not part of $\widehat M$ (see the [topological-sector record](alp-qcd-topological-cp.md)); the singlet invariant $c_{GA}$ is recovered in the [dressed topological record](alp-dressed-topological-invariant.md).

## Derivation

Integrating $\partial_xW_\eta=iC_\eta W_\eta$ with $W_\eta(0)=\mathbf1$ gives $W_\eta=\mathbf1+ixC_{\eta,0}+O(x^2)$. Hence

$$
\widehat M(x)=M_0+x\left[M_1-i\left(C_{R,0}M_0-M_0C_{L,0}\right)\right]+O(x^2),
\qquad M_0=m_q .
$$

Truncating $W_\eta$ at order $x$ is exact for this coefficient. Using $C_{R,0}M_0-M_0C_{L,0}=[C_{V,0},M_0]+\{C_{A,0},M_0\}$ and $M_1=i[c_q^S,m_q]+i\{c_q^P,m_q\}$,

$$
\boxed{\widehat M_1=-i\,[c_q^V-c_q^S,\,m_q]-i\,\{c_q^A-c_q^P,\,m_q\}
=-i[c_q^{VS},m_q]-i\{c_q^{AP},m_q\}.}
$$

The commutator term is Hermitian and the anticommutator term is anti-Hermitian, so

$$
\widehat S_1=-i[c_q^{VS},m_q],\qquad \widehat P_1=-\{c_q^{AP},m_q\},
\qquad \widehat S_0=m_q,\ \widehat P_0=0 .
$$

In the CCWZ spurion, $\chi(x)=2B_0\widehat M(x)$, so the $O(1/f_a)$ source has $\chi_{S,1}=2B_0\widehat S_1$ and $\chi_{P,1}=2B_0\widehat P_1$.

## Comparison with Alda:2026nrd

1. **Invariance.** Under $c_q^{V,S}\to c_q^{V,S}+\delta$ and $c_q^{A,P}\to c_q^{A,P}+\kappa$, i.e. $g_{R,L}=e^{ix(\delta\pm\kappa)}$, the combinations $c_q^{VS}$ and $c_q^{AP}$ do not change, so $\widehat M_1$ does not change. This reproduces the transformation laws and invariants of Alda:2026nrd [@Alda:2026nrd], and the first-order laws $S_1'$, $P_1'$ of the general record.
2. **Alda:2026nrd's prescription.** The choice $\delta=-c_q^S$, $\kappa=-c_q^A-\mathbb M_a^{-1}(C)$ gives $\widehat c_q^S=0$, $\widehat c_q^V=c_q^{VS}$, $\widehat c_q^P-\widehat c_q^A=-c_q^{AP}$. Dressing this representative gives the same $\widehat M_1$, independent of $\mathbb M_a^{-1}(C)$. The dressed representative has instead $\widehat C_V=\widehat C_A=0$ and $\widehat S_1,\widehat P_1$ as above; the two representatives are related by the constant redefinition $\delta=-c_q^V$, $\kappa=-c_q^A$.
3. **Information content.** The anticommutator with a positive-definite $m_q$ is invertible (for $m_q=\mathrm{diag}(m_u,m_d,m_s)$ its determinant on $3\times3$ matrices is $8m_um_dm_s(m_u+m_d)^2(m_u+m_s)^2(m_d+m_s)^2$), so $\widehat P_1$ determines $c_q^{AP}$ completely. The commutator with $m_q$ is not invertible: its kernel is the set of matrices commuting with $m_q$ (dimension 3 for nondegenerate $m_q$, 5 for $m_u=m_d\neq m_s$). The part of $c_q^{VS}$ commuting with $m_q$ does not enter $\widehat M_1$; for example, a diagonal $c_q^{VS}$ with diagonal $m_q$ gives $\widehat S_1=0$. That part is still an invariant, and it appears in the vector-current coupling, e.g. $[c_q^{VS},\mathcal B]$ in Alda:2026nrd's baryon couplings. So $\widehat M_1$ carries all of $c_q^{AP}$ but only the off-commutant part of $c_q^{VS}$. The same limitation holds for $c_q^S$ itself, since $S_1=i[c_q^S,m_q]$ only defines $c_q^S$ up to this kernel.
4. **Higher orders.** This reduction is first order only. At order $x^n$, $\widehat M_n$ contains products and commutators of lower-order coefficients, as stated in the general record, and is not a pairwise difference. The explicit second-order coefficients $\widehat S_2$, $\widehat P_2$ are in the [second-order record](alp-dressed-mass-second-order-invariants.md).

## Reproducible calculation

The script [`derive_alp_dressed_mass_first_order.py`](../../../scripts/derive_alp_dressed_mass_first_order.py) uses noncommuting symbols, with $x$ real, and checks: the general first-order formula; the Alda:2026nrd form; the Hermitian and anti-Hermitian split; invariance under the first-order redefinition; an exact rational $3\times3$ check; the invertibility and kernel statements above; and that Alda:2026nrd's representative gives the same $\widehat M_1$. All checks pass. A unit-test file, `scripts/test_derive_alp_dressed_mass_first_order.py`, accompanies it.

```sh
.venv/bin/python scripts/derive_alp_dressed_mass_first_order.py
```
