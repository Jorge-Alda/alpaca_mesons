---
title: One-loop RChT generating functional with S and P resonances generates O(p^2) counterterms
status: literature
created: 2026-11-05
updated: 2026-11-05
bibtex_keys: [Rosell:2005ai]
provenance: literature
---

# O(p^2) divergences in the one-loop RChT generating functional

Source: `memory/sources/Rosell:2005ai/rct.tex`. [@Rosell:2005ai]

## Setup (assumptions of the source)

- $U(3)_L\otimes U(3)_R$ with one nonet each of scalar $S$ and pseudoscalar $P$ resonances; $SU(3)$ external currents (lines 421-424). Vector and axial-vector resonances are not included.
- Lagrangian (Eq. `eq:lagr1`, lines 425-460): $\mathcal{L}_\chi^{(2)}=\frac{F^2}{4}\langle u_\mu u^\mu+\chi_+\rangle$, kinetic terms of $S,P$ (with covariant derivatives), and $\mathcal{L}_2(S),\mathcal{L}_2(P),\mathcal{L}_2(S,P)$: couplings $c_d,c_m,d_m$ (one resonance) and $\lambda_i^{SS},\lambda_i^{PP},\lambda_i^{SP}$ (two resonances). No $\mathcal{L}_\chi^{(4)}$ is included (lines 476-480). Operators with three or more resonances are dropped (lines 1046-1062).
- Counting: $F,c_d,c_m,d_m\sim\mathcal{O}(\sqrt{N_C})$, $\lambda_i^{RR}\sim\mathcal{O}(1)$ (lines 531-535).
- The results are given without imposing short-distance constraints (lines 555-557).
- Method: background-field expansion, heat kernel, dimensional regularization; divergence $S_1^{\rm div}$ in Eq. `eq:s1div` (lines 789-791); result Eq. `eq:main` (line 898): $\mathcal{L}_1=\sum_{i=1}^{18}\alpha_i\mathcal{O}_i+\sum_{i=1}^{68}\beta_i^R\mathcal{O}_i^R+\sum_{i=1}^{383}\beta_i^{RR}\mathcal{O}_i^{RR}$, with $\mu\,d\alpha_i^r/d\mu=-\gamma_i/(16\pi^2)$ (Eq. `eq:rgeq`).

## Result: $\mathcal{O}(p^2)$ operators

Table `tab:e1` (appendix, lines 1954-1965) lists the Goldstone-only operators and their coefficients $\gamma_i$ (both $S$ and $P$ included). The $\mathcal{O}(p^2)$ operators are $\mathcal{O}_1=\langle u\cdot u\rangle$, $\mathcal{O}_2=\langle\chi_+\rangle$ and $\mathcal{O}_3=-\langle u_\mu\rangle^2$, with $N$ the number of flavours:

$$
\gamma_1=-2N\lambda_1^{PP}M_P^2+\tfrac12 N M_P^2(\lambda_1^{SP})^2-2N\lambda_1^{SS}M_S^2+\tfrac12 N M_S^2(\lambda_1^{SP})^2+N\frac{c_d^2}{F^2}M_S^2,
$$
$$
\gamma_2=-2N\lambda_3^{PP}M_P^2-2N\lambda_3^{SS}M_S^2,
$$
$$
\gamma_3=2\lambda_2^{PP}M_P^2-\tfrac12M_P^2(\lambda_1^{SP})^2+2\lambda_2^{SS}M_S^2-\tfrac12M_S^2(\lambda_1^{SP})^2-\frac{c_d^2}{F^2}M_S^2 .
$$

Since $\mathcal{L}_\chi^{(2)}=\frac{F^2}{4}\langle u\cdot u+\chi_+\rangle$, a nonzero $\gamma_1$ is a one-loop divergence of the coefficient of $\langle u\cdot u\rangle$ (i.e. of $F^2$), and a nonzero $\gamma_2$ one of the coefficient of $\langle\chi_+\rangle$ (i.e. of $F^2B_0$). This reading of the table is mine; the source lists the coefficients without this wording. Every term is $M_R^2\times\mathcal{O}(1)$ couplings, or $c_d^2M_S^2/F^2$, so relative to $F^2\sim N_C$ the correction is of relative order $1/N_C$. Interpretation of the diagrams (not stated by the source): $\lambda^{RR}$ terms come from resonance loops with one resonance propagator closing on a two-resonance vertex; $(\lambda_1^{SP})^2$ terms from the $SP$ bubble; $c_d^2M_S^2/F^2$ from the $S\pi$ bubble.

## Other statements

- Lines 943-966: Table `tab:e1` collects the full basis of $\mathcal{O}(p^2)$ and $\mathcal{O}(p^4)$ $U(3)$ $\chi$PT operators generated in the integration; it agrees with Herrera-Siklody:1996tqr once resonance couplings are switched off. $\mathcal{L}_{\rm kin}(S,P)$ without resonance couplings contributes to $\gamma_4,\gamma_7,\gamma_{13},\gamma_{14},\gamma_{15}$ only.
- Running of $c_d^r,c_m^r,d_m^r$ (Eq. `eq:evolution`, lines 995-1015) is a next-to-leading order effect in $1/N_C$, from $M_R^2/F^2$ factors (lines 1016-1020). Using Eqs. `eq:largecd`, `eq:largecd2` (lines 539-550; $\lambda_3^{SS}=\lambda_3^{PP}=0$, labelled "controversial" by the source) the operators with $\chi_\pm$, including $\mathcal{O}_2$, do not run (lines 1032-1043). For $\gamma_1$ the source makes no such statement.
- Worked example (lines 1085-1105): scalar-current two-point function in the chiral limit with only $S$; counterterms $\mathcal{O}_{12}+2\mathcal{O}_{16}=\langle\chi_+^2\rangle$, $\mathcal{O}_4^S=\langle S\chi_+\rangle$, $\mathcal{O}_1^{SS}=\langle SS\rangle$.

## Verdict

- O(p^2) affected: **yes**, divergences renormalizing $F^2$ ($\gamma_1$) and $F^2B_0$ ($\gamma_2$), plus the singlet $\mathcal{O}_3$ ($\gamma_3$), all $1/N_C$ suppressed. $\mathcal{O}(p^4)$ operators are generated as well.
- Conditions and caveats: (i) no short-distance constraints imposed; with $\lambda_3^{RR}=0$ the $B_0$ term $\gamma_2$ vanishes, but $\gamma_1$ is not shown to vanish; (ii) $U(3)$ with $S,P$ only, no $V,A$; the $V\pi$ contribution is in Rosell:2004mn (see `../Rosell:2004mn/pion-self-energy-o-p2-divergence.md`); (iii) only divergent parts are computed; finite local $\mathcal{O}(p^2)$ remainders are scheme dependent and not given; (iv) the operator basis was reduced with equations of motion (appendix `app:tres`) and the source says possible $U(3)$ algebraic relations were not used, so the basis is not necessarily minimal; (v) $N$ is the number of flavours.
- Thresholds: this source gives no kinematic discussion.
