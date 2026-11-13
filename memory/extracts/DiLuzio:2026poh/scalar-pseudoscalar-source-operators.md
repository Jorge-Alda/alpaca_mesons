---
title: DiLuzio:2026poh scalar and pseudoscalar QCD operators sourced by varphi
status: literature
created: 2026-11-13
updated: 2026-11-13
bibtex_keys: [DiLuzio:2026poh]
provenance: literature
---

# Scalar and pseudoscalar QCD operators coupled to $\varphi$

The paper's appendix "Effective pion potential" treats the PQ scalar as an external source in the chiral limit. The paper uses the operator names $\mathcal O_{\mathbf 1,S}$, $\mathcal O_{\mathbf 1,P}$, $\mathcal O^a_{\mathbf{Adj},S}$, $\mathcal O^a_{\mathbf{Adj},P}$; it does not use the symbols $S_n$ or $P_n$.

## Source and decomposition

Upgrade $\kappa_q\varphi\to\Phi$ with $\Phi\to L\Phi R^\dagger$. For quadratic order in $\Phi$ and the chiral limit [@DiLuzio:2026poh, Eq. (70)]:
$$
\mathscr L_{\rm eff}=\Pi_A(p)\,\mathrm{tr}[\Phi\Phi^\dagger]+\Pi_B(p)\,\mathrm{tr}[\Phi U^\dagger]\,\mathrm{tr}[\Phi^\dagger U]+\Pi_0\,\mathrm{tr}[\Phi U^\dagger+\text{h.c.}].
$$
The operators $\mathcal O_{LR}^{ij}\equiv\bar q_L^iq_R^j$ are decomposed [@DiLuzio:2026poh, Eq. (71)]:
$$
\mathcal O_{LR}^{ij}=\frac{\delta^{ij}}{\sqrt{2N_F}}\mathcal O_{\mathbf1}+T^a_{ij}\mathcal O^a_{\mathbf{Adj}},\quad
\mathcal O_{\mathbf1}=\mathcal O_{\mathbf1,S}-i\mathcal O_{\mathbf1,P},\quad
\mathcal O^a_{\mathbf{Adj}}=\mathcal O^a_{\mathbf{Adj},S}-i\mathcal O^a_{\mathbf{Adj},P},
$$
with $T^a$ the $\mathrm{SU}(N_F)$ generators, $\mathrm{tr}[T^aT^b]=\delta^{ab}/2$ and $\Phi=\Phi_1/\sqrt{2N_F}+\Phi^a_{\mathbf{Adj}}T^a$. The singlet operators are $\mathcal O_{\mathbf1,S}\equiv\sqrt{2/N_F}\,\bar qq$ and $\mathcal O_{\mathbf1,P}\equiv\sqrt{2/N_F}\,i\bar q\gamma_5q$ (and "similarly" for the adjoint). The paper takes $N_F=3$ but writes general-$N_F$ formulas.

## Two-point functions

$$
\langle\mathcal O_{\mathbf1,(S,P)}(p)\mathcal O_{\mathbf1,(S,P)}(-p)\rangle=i\Pi_{\mathbf1,(S,P)}(p),\qquad
\langle\mathcal O^a_{\mathbf{Adj},(S,P)}(p)\mathcal O^b_{\mathbf{Adj},(S,P)}(-p)\rangle=i\Pi_{\mathbf{Adj},(S,P)}(p)\delta^{ab}\quad[\text{Eq. (72)}].
$$
Single-pole (large-$N$) saturation [@DiLuzio:2026poh, Eq. (73)]:
$$
\Pi_{\mathbf1,P}=\frac{f_{\eta'}^2B_0^2}{p^2-m_{\eta'}^2}+\cdots,\quad
\Pi_{\mathbf{Adj},P}=\frac{f_\pi^2B_0^2}{p^2-m_\pi^2}+\cdots,\quad
\Pi_{\mathbf1,S}=\frac{f_{\sigma_0}^2B_0^2}{p^2-m_{\sigma_0}^2}+\cdots,\quad
\Pi_{\mathbf{Adj},S}=\frac{f_{\sigma_8}^2B_0^2}{p^2-m_{\sigma_8}^2}+\cdots .
$$
Matching on $U=1$ [@DiLuzio:2026poh, Eq. (74)]:
$$
\Pi_A=-\frac{\Pi_{\mathbf{Adj},S}+\Pi_{\mathbf{Adj},P}}2,\qquad
\Pi_B=\frac{\Pi_{\mathbf{Adj},S}+\Pi_{\mathbf{Adj},P}-\Pi_{\mathbf1,S}-\Pi_{\mathbf1,P}}{2N_F},\qquad\Pi_0=\frac{f_\pi^2B_0}2 .
$$

## Pion splitting from the one-loop effective potential

The effective potential from integrating $\varphi$ gives, at leading order in $\kappa_{u,d}$ [@DiLuzio:2026poh, Eq. (77)]:
$$
\Delta(m_{\pi^0}^2-m_{\pi^\pm}^2)\approx(\kappa_un_u-\kappa_dn_d)^2\frac2{f_\pi^2}\int\frac{d^4Q}{(2\pi)^4}\frac{\Pi_B(Q^2)}{-Q^2-m_\varphi^2},
$$
and a closed form for $\int\frac{d^4Q}{(2\pi)^4}\Pi_B$ in terms of $m_{\eta'},m_{a_0},m_\sigma,f_{\eta'},f_\pi$ is given in Eq. (78).

## Assumptions and ambiguities

- Single-pole dominance; the scalar decay constants $f_{\sigma_{0,8}}$ are unknown and are fixed by requiring the integral to be finite for heavy $m_\varphi$.
- Treats $\kappa_q$ as a perturbation; the paper itself says the result is uncertain for $\kappa_q\sim1$.
- Possible normalization inconsistency (extractor): the stated $\mathcal O_{\mathbf1,S}-i\mathcal O_{\mathbf1,P}=\sqrt{2/N_F}(\bar qq+\bar q\gamma_5q)$ equals $2\sqrt{2/N_F}\sum_i\bar q_L^iq_R^i$ with $\bar q_Lq_R=\bar q\frac{1+\gamma_5}2q$, which gives $\mathrm{tr}\,\mathcal O_{LR}$ twice the sum of $\bar q_L^iq_R^i$ in Eq. (71). The factor 2 is not resolved in the paper and is not needed for the quoted results; verify before using the operator normalization quantitatively.
- Eq. (73) normalizes poles by $f^2B_0^2$ without stating the matching to $\mathcal O_{\mathbf1,P}$ normalization.
