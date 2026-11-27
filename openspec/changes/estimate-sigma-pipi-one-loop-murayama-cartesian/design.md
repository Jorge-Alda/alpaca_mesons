# Design

## Context

See proposal.md. The polar estimate (`estimate-sigma-pipi-one-loop-murayama`) fixes the assumptions and the size rule to be reused: D1 size rule, D2 products, D3 loop scales, D4 exclusion, D10 regime, D11 CP and flavour-diagonality, D12 routing, and $M_1^{\rm est}=\sum_\Gamma S_\Gamma\,(m_\Gamma^2)^{2-r_\Gamma}T_\Gamma/(16\pi^2)$. These are rules and assumptions, not results. The vertices, graph list, traces and powers are derived anew in the Cartesian fields. The existing notes are used for one thing: the tree amplitude $M_0$ and the equivalence of the Cartesian and polar trees (`radial-decay-polar-cartesian-general-sources.md`).

Cartesian frame: $\varphi=f_a+s+ip$, kinetic term $\tfrac12(\partial s)^2+\tfrac12(\partial p)^2$ (no derivative vertex of $s$ or $p$), potential from $V_\varphi(s,p)$ and from the source $\widehat M(s,p)=\sum_{r,s}(\widehat S^C_{r,s}+i\widehat P^C_{r,s})(s/f_a)^r(p/f_a)^s$ coupled to $U$ in the chiral Lagrangian. The invariants $\widehat S^C_{r,s}$, $\widehat P^C_{r,s}$ are Hermitian $3\times3$ matrices of mass dimension one.

## Goals / Non-Goals

**Goals:**
- Closed symbolic $M_1^{\rm est}$ in $\widehat S^C_{r,s}$, $\widehat P^C_{r,s}$, with the same rule and assumptions as the polar estimate and an independent derivation.
- A diagnostic comparison with the polar estimate, with the order-of-magnitude caveat.

**Non-Goals:**
- Computing any loop integral, scalar function, coefficient or amplitude; numerical scans; absorptive parts; renormalization; counterterm analysis.
- Reproducing the polar result or the polar up-only prefactor, or deciding which estimate is better.
- Changing the polar records, the draft or the main specs.

## Decisions

**C1. Rules and assumptions reused, results not.** D1–D5, D7, D9–D12 are applied as stated. The vertices, graph list, traces and $\epsilon$ powers are derived in the Cartesian fields and are not copied. Alternative rejected: translating the polar vertices through the field map, which would inherit the polar structure.

**C2. Assumptions in the Cartesian invariants.** Each polar assumption is rewritten for $\widehat S^C$, $\widehat P^C$ (no tadpoles, no quadratic $s$–$p$ or $s$–meson mixing, $\widehat\Theta_1=0$, CP, flavour-diagonality, $f_\pi$ independent of the radial fields), with the condition on the coefficients derived from the Cartesian quadratic Lagrangian, and linked to the polar assumption it corresponds to. Legs are $s$ (radial), $p$ (ALP field) and mesons; external $\sigma$ is the pole of $s$.

**C3. Vertices from the Cartesian Lagrangian.** For $n=3,4,5$ expand $\mathcal L(s,p,\Phi)$ noncommutatively in formal legs and extract the coefficient of $\phi_1\cdots\phi_n$ (construction 1). Write a closed rule for the potential vertices in terms of $\widehat S^C_{r,s}$, $\widehat P^C_{r,s}$, the $V_\varphi(s,p)$ polynomial and the meson kinetic vertices (construction 2). The two agree as exact identities of cyclic words. No polar vertex or field map enters.

**C4. Checks of the vertices.** The checks are: (a) relabelling symmetry of each vertex; (b) dimensional scaling $c^{4-n}$; (c) the up-only limit; (d) **tree equivalence**, the only comparison with the polar results: the on-shell cubic vertex $V(s;\pi^+,\pi^-)$ equals the polar tree $M_0$ of the existing notes, using the relation between the Cartesian and polar coefficients, stationarity of $V_\varphi$ and the pole equation, as in `radial-decay-polar-cartesian-general-sources.md`. No other comparison with polar vertices (off-shell vertices, quartic or quintic vertices) is used.

**C5. Graph list, traces, powers derived anew.** The one-loop 1PI graphs and external-leg corrections are enumerated from the Cartesian vertices by script: topologies (tadpole $(5)$, bubble $(3,4)$, triangle $(3,3,3)$; box and 1PI $(2,5)$ checked, not assumed), symmetry factors, $r_\Gamma$, the D3 split of the meson line, D4 terms, D12 routing and the large-$\epsilon$ powers at fixed sources and fixed $\lambda$. Only afterwards is the graph list compared with the polar one.

**C6. Diagnostic comparison.** After the derivation, a table compares the classes, leading powers and coefficients (generic, up-only) of the two estimates. It is a diagnostic: results are orders of magnitude, graphs with the same $\epsilon$ scaling can cancel in the exact loop in one parametrization and not be visible in the other, so no agreement is expected and neither result is authoritative. The record states, per class, whether same-scaling graphs have leading coefficients that could cancel.

**C7. Check, not scan.** Axial-redefinition invariance on the traces, completeness identities and trace equalities for generic noncommuting symbols, all exact; no loop is evaluated.

## Risks / Trade-offs

- [The size rule gives only an order of magnitude per class and can miss cancellations between same-scaling graphs, differently in each parametrization] → Stated in the records; neither estimate is authoritative; the comparison is a diagnostic.
- [Many source-coefficient combinations at five legs] → Generated by script, verified by relabelling symmetry and the leg-pattern count.
- [$\lambda$ enters the $p$-vertices directly through $V_\varphi(s,p)$] → Keep $\lambda$ symbolic and eliminate it only through $m_\sigma^2$ as in D10, stating the substitution.
- [Thresholds, absorptive parts, naturalness without counterterm analysis] → Restated limits.
