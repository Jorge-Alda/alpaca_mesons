---
title: Included and excluded loop classes in the U(3) NLO theory with explicit resonances
status: mixed
created: 2026-11-05
updated: 2026-11-05
bibtex_keys: [Kaiser:1998ds, Herrera-Siklody:1996tqr, Rosell:2004mn, Rosell:2005ai, Rosell:2006dt, Pich:2010sm, Pich:2008jm, ParticleDataGroup:2024cfk]
provenance: analytical-derivation
---

# Loop classes

## Assumptions and inputs

- The fields are those of [field-content.md](field-content.md). Orders follow the rule $\delta^{\,n-k}$ of [power-counting.md](power-counting.md), which assigns that order to a term $\propto N_c^k p^{2n}$ normalized with $F^2$.
- The declared calculation is $\mathcal O(\delta^2)$ in the Goldstone/ALP sector, at first order in $1/f_a$.
- **User decision:** only Goldstone loops are included. Loops with an internal resonance line are excluded by a declared approximation.
- **Counting assumption:** a loop with only Goldstone internal lines has loop momenta of the order of the external momenta and meson masses. It multiplies the tree graph built from the same vertices by $N_c^{-1}p^2$, that is, by $p^2/(16\pi^2F^2)$.
- **Counting assumption:** a loop with an internal resonance line has loop momenta of order $M_R$. It multiplies by $N_c^{-1}$ without a factor $p^2$, because the vertex momentum factors are supplied by $M_R$.
- Masses for thresholds are PDG values in GeV [@ParticleDataGroup:2024cfk].

The executable classification is in [`scripts/derive_nlo_u3_resonance_counting.py`](../../../scripts/derive_nlo_u3_resonance_counting.py):

```bash
.venv/bin/python scripts/derive_nlo_u3_resonance_counting.py --loops
```

## Classification

The classification covers contributions to Goldstone/ALP two-point functions and to ALP decay vertices. Orders are the low-energy limit unless stated otherwise.

| Class | Graphs | $k$ | $p^{2n}$ | Order | Status |
|---|---|---|---|---|---|
| G-LO | one loop with only Goldstone lines, including $\eta'$, all vertices from $\mathcal L^{(0)}$ | 0 | $p^4$ | $\delta^2$ | included |
| G-ALP | as G-LO, with one $\mathcal L^{(0)}$ vertex containing an ALP source ($c_q^A,c_q^V,c_q^S,c_q^P,c_G$) | 0 | $p^4$ | $\delta^2$ | included |
| G-NLO | Goldstone loop with one $L_5$ or $L_8$ vertex | 0 | $p^6$ | $\delta^3$ | excluded: beyond the declared order |
| G-NLO | Goldstone loop with one $\Lambda_1$ or $\Lambda_2$ vertex | $-1$ | $p^4$ | $\delta^3$ | excluded: beyond the declared order |
| R-int, $\mathcal O(p^2)$ part | one loop with at least one internal resonance line | 0 | $p^2$ | $\delta$ | excluded: declared approximation |
| R-int, $\mathcal O(p^4)$ part | same | 0 | $p^4$ | $\delta^2$ | excluded: declared approximation |
| R-vertex | Goldstone loop on a resonance vertex of a tree exchange graph, including tadpoles | 0 | $p^6$ | $\delta^3$ | excluded: beyond the declared order |
| R-dress | Goldstone loop in a resonance propagator, $\Sigma_R$ | 0 | $p^8$ | $\delta^4$ | excluded: beyond the declared order |
| R-dress near the pole | $\Sigma_R/M_R^2\sim\Gamma_R/M_R\sim N_c^{-1}$ at $s\simeq M_R^2$ | – | – | – | outside the strict counting; line-shape prescription |
| A-int | loop with an internal ALP line | – | – | – | excluded: first order in $1/f_a$ |
| G-2loop | two Goldstone loops with $\mathcal L^{(0)}$ vertices | $-1$ | $p^6$ | $\delta^4$ | excluded: beyond the declared order |

The included classes are exactly G-LO and G-ALP. They have the order of $\mathcal L^{(2)}$, so their UV divergences are absorbed by the $\delta^2$ constants $L_4$, $L_6$, $L_7$, $L_{18}$, $L_{25}$, and by $L_5$ and $L_8$ with their $U(3)$ coefficients [@Kaiser:1998ds; @Herrera-Siklody:1996tqr]. The counterterm basis is given in the counterterm record of this folder.

## Step-by-step derivation of each order

1. **G-LO, G-ALP.** The tree graph from $\mathcal L^{(0)}$ is $N_cp^2$. The Goldstone loop gives $N_c^{-1}p^2$, so the product is $N_c^0p^4$, which is $\delta^2$. An ALP source in one vertex changes neither power, only the order in $1/f_a$. The $\eta'$ runs in these loops, because $m_0^2=\mathcal O(\delta)$ makes it a light field of the counting ([power-counting.md](power-counting.md)).
2. **G-NLO.** A vertex from $\mathcal L^{(1)}$ is already $\delta$, either $N_cp^4$ or $N_c^0p^2$. The loop adds $\delta^2$, giving $\delta^3$.
3. **R-int.** The loop momenta are of order $M_R$, so the loop factor is $N_c^{-1}$ alone. Chiral symmetry still requires a Goldstone two-point function to start at $p^2$ or $m_q$. The lowest contribution is therefore $N_c^{-1}$ times the $F^2$-normalized term $N_cp^2$, which gives $N_c^0p^2$, that is, $\delta$. Contributions at $p^4$ then give $\delta^2$. The sources confirm the $\mathcal O(p^2)$ piece:
   - In the pion self-energy, the $V\pi$ and $S\pi$ loops give an $\mathcal O(p^2)$ divergence proportional to $6G_V^2M_V^2/(16\pi^2F^4)$, together with a $c_d^2M_S^2$ term. It is absorbed in the pion wave-function renormalization ([extract](../../extracts/Rosell:2004mn/pion-self-energy-o-p2-divergence.md)) [@Rosell:2004mn].
   - In the one-loop generating functional with $U(3)$ scalar and pseudoscalar resonances, the divergences multiply $\langle u_\mu u^\mu\rangle$ and $\langle\chi_+\rangle$, that is, they renormalize $F^2$ and $F^2B_0$ ([extract](../../extracts/Rosell:2005ai/o-p2-counterterms-one-loop-generating-functional.md)) [@Rosell:2005ai]. Reading the wave-function renormalization of Rosell:2004mn as a shift of $F^2$ uses this result.
4. **R-vertex.** Tree exchange is $\delta$ ([power-counting.md](power-counting.md), Step 5). A Goldstone loop on one of its vertices adds $N_c^{-1}p^2$, giving $\delta^3$.
5. **R-dress.** Each vertex of $\Sigma_R$ is $\sqrt{N_c}\,p^2/F^2$, so $\Sigma_R/M_R^2\sim N_c^{-1}p^4$. In the low-energy limit the result is $\delta\times\delta^3=\delta^4$. Near the pole, $s\simeq M_R^2$, the factor is $N_c^{-1}$ with no momentum suppression, so it is not small.
6. **A-int.** The loop factor is $1/(16\pi^2f_a^2)$. It is second order in $1/f_a$ and has no $\delta$ order.
7. **G-2loop.** Two Goldstone loops give $N_c^{-2}p^4$ times $N_cp^2$, which is $N_c^{-1}p^6$, that is, $\delta^4$.

## The excluded resonance loops

R-int starts at $\delta$, below the declared order. Its exclusion is therefore an approximation, not a statement of the counting. It is controlled as follows.

- **Below the $R\pi$ thresholds (assumption).** R-int is analytic in $s$ and $m_q$, so its effect has the form of local terms at $\mathcal O(p^2)$ and $\mathcal O(p^4)$. The sources absorb these terms in two places:
  - in the Goldstone constants: the pion wave function and $F^2$, $F^2B_0$ in Rosell:2004mn and Rosell:2005ai;
  - in the resonance parameters: $F_V$, $F_A$, $M_V$, $M_A$ through the NLO Weinberg sum rules, $c_m$, $d_m$, $G_V$, and $L_{10}^r$, $C_{87}^r$. In Pich:2008jm, $F$ itself is an input without a loop shift ([extract](../../extracts/Pich:2008jm/nlo-weinberg-sum-rules-and-resonance-couplings.md)) [@Pich:2008jm; @Rosell:2006dt].

  The exclusion is therefore consistent only if $F$, $B_0$, $\tau$ (that is, $m_0^2$), $\Lambda_1$, $\Lambda_2$, $L_i^r(\mu)$ **and the resonance masses and couplings** are fixed from data, either as physical values or as fits. Estimates from large $N_c$, resonance saturation or leading-order short-distance relations do not contain the resonance-loop part. This condition binds the matching in alp-gii.3 and the inputs in alp-gii.5.
- **Short-distance condition.** The loop corrections are local and their polynomial divergences cancel only if the spectral functions vanish at high energy. This is imposed through the short-distance constraints, for example the Weinberg sum rules at NLO ([extract](../../extracts/Pich:2008jm/two-particle-cuts-vv-minus-aa-correlator.md)) [@Pich:2008jm]. Without them, absorption into constants fixed from data is not guaranteed by the sources. Here $L_i^r(\mu)$ is the constant that remains after explicit resonance exchange is removed (see the design decision on notation).
- **Thresholds.** Resonance–resonance cuts lie at $(M_R+M_{R'})^2\gtrsim2\,$GeV$^2$, outside the validity range ([extract](../../extracts/Pich:2010sm/resonance-loop-cuts-and-divergence-structure.md)) [@Pich:2010sm]. Only resonance–Goldstone cuts can lie below $1\,$GeV. In the pseudoscalar channel the $V\pi$ and $S\pi$ intermediate states contribute ([extract](../../extracts/Rosell:2006dt/resonance-pion-cuts-scalar-pseudoscalar-correlator.md)) [@Rosell:2006dt]. The lowest vector–pseudoscalar thresholds are:

| State | Threshold (GeV) |
|---|---|
| $\rho\pi$ | 0.9148 |
| $\omega\pi$ | 0.9222 |
| $\rho\eta$ | 1.3231 |
| $\omega\eta$ | 1.3305 |
| $K^*K$ | 1.3853 |

- **Consequence for the validity range.** Between about $0.91\,$GeV and $1\,$GeV, the non-analytic part of R-int is open in channels where a vector–pseudoscalar state is allowed. It is not included there, and it enters the truncation uncertainty. Selection rules for specific channels are not applied in this record.
- **Size.** For the relative size of the excluded loops, the one-loop RChT analyses of the vector form factor and of $\Pi_{VV}-\Pi_{AA}$ provide the reference [@Rosell:2004mn; @Pich:2010sm; @Pich:2008jm]. In Pich:2008jm, at $\mu=770\,$MeV, $L_{10}^r$ goes from $-5.3\times10^{-3}$ (LO) to $(-4.4\pm0.9)\times10^{-3}$ (NLO in $1/N_c$), and $C_{87}^r$ from $4.3\times10^{-5}$ to $(3.1\pm1.1)\times10^{-5}$. The dominant error comes from the second resonance multiplets $V'$, $A'$ ([extract](../../extracts/Pich:2008jm/nlo-size-L10-C87-numerics.md)). The central values shift by about 17% ($L_{10}^r$) and 28% ($C_{87}^r$), consistent with a $1/N_c$ correction. The numerical estimate used here is set in the validity-range record of this folder.
- **Flavor content.** In Pich:2008jm, $\eta_1$ enters neither at tree level nor at one loop in $\Pi_{VV}-\Pi_{AA}$, so $L_{10}^r$ and $C_{87}^r$ are the same in $U(3)$ and $SU(3)$ at NLO in $1/N_c$ ([extract](../../extracts/Pich:2008jm/u3-su3-matching-singlet.md)). This does not hold for the scalar and pseudoscalar correlators that feed the quadratic sector ([extract](../../extracts/Rosell:2006dt/resonance-pion-cuts-scalar-pseudoscalar-correlator.md)).

## Resonance widths

R-dress is beyond the declared order at low energy. Near a resonance pole it is included only through the line-shape prescription of [vector-width-lineshape-prescriptions.md](../../rcht/foundations/vector-width-lineshape-prescriptions.md). This is a resummation outside the strict counting. Rosell:2004mn and Pich:2010sm show that a consistent dressed propagator requires the corresponding loop-induced vertex and counterterm pieces ([extract](../../extracts/Rosell:2004mn/vector-self-energy.md), [extract](../../extracts/Pich:2010sm/vector-form-factor-nlo.md)). Downstream calculations therefore report both the strict-order result with a narrow propagator and the result with a line shape.

## Verification

The script asserts the following, and all runs exit with code 0 (`--loops`, `--orders`, the default run and `--test`):

- the included set is exactly {G-LO, G-ALP};
- every class labelled "beyond the declared order" is above $\delta^2$;
- R-int has orders $\{\delta,\delta^2\}$, with the lowest below the declared order;
- the $\rho\pi$ threshold is below $1\,$GeV.

The tables in this record are copied from its output.
