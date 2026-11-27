---
title: Cartesian frame - limits and assumptions of the one-loop size estimate of sigma to pi+ pi-, including the vacuum condition
status: mixed
created: 2026-11-27
updated: 2026-11-27
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: assumption
---

# Cartesian frame: limits of the size estimate of $\sigma\to\pi^+\pi^-$ at one loop

Task 4.5 of the OpenSpec change `estimate-sigma-pipi-one-loop-murayama-cartesian`. This record lists what the Cartesian estimate $M_1^{\rm est}$ and the ratio $R=|M_1^{\rm est}|/|M_0|$ of the [Cartesian estimate record](cartesian-sigma-pipi-loop-estimate.md) do not contain, and the assumptions on which they rest. The limits of the polar estimate ([polar limits record](../one-loop-estimate/murayama-sigma-pipi-loop-limits.md), used here for its list of limits only) are restated for the Cartesian fields. **No loop integral, scalar loop function, loop coefficient or loop amplitude is computed**, and this record contains no value of any of them. $R$ is a size ratio, not a computed correction.

## Caveat

The Cartesian and polar estimates are order-of-magnitude estimates per graph class. Graphs with the same $\epsilon$ scaling can cancel in the exact loop, and the cancellations can be different in each parametrization, because the size rule acts on different vertices and lines in each. Agreement between the Cartesian and the polar estimate, including the polar up-only prefactor $-21/(64\pi^2)$, is not expected, and neither estimate is authoritative. The comparison of the two is a diagnostic (task 4.4).

## Assumptions entering the graph list

- **CP and flavour-diagonality (D11, assumption).** CP with $s$ even and $p$, $\Phi$ odd, and flavour-diagonal sources and quark masses. They remove the one-particle-reducible graphs in which an external leg mixes at one loop with another state ($s$ with $p$ or a neutral meson, $\pi^\pm$ with $K^\pm$); without them these graphs would have to be added. The trace algebra itself keeps generic noncommuting sources. Rule R7 of the [Cartesian selection record](cartesian-sigma-pipi-loop-graph-selection.md).
- **Tadpoles (assumption).** The one-loop tadpoles of $\pi^0$, $\eta$, $\eta'$ and $p$ vanish under CP; the one-loop $s$ tadpole is set to zero (rule R8, assumption 9 of the [Cartesian frame and tree record](cartesian-sigma-pipi-loop-frame-and-tree.md)). This is a vacuum condition; its difference from the polar one is treated below.
- **Residue prescription (assumption).** The external legs enter through $\delta Z=-\Pi'(m^2)$, with the derivative acting on the explicit $p^2$ dependence of the numerator and, for $r_\Gamma\ge2$, on one propagator, each replaced by its size (rule R10, [Cartesian external-leg record](cartesian-sigma-pipi-loop-traces-external-legs.md)). The overall sign follows $M=V$; the size rule does not fix relative signs.
- **Other assumptions** of the frame record: no tree tadpoles, no quadratic $s$–$p$ or $s$–meson mixing, $\mathcal C=B_0\widehat P^C_{0,1}$ (Cartesian form of $\widehat\Theta_1=0$), no derivative coupling of $s$, $p$ to the mesons at leading $U(3)$ chiral order, $f_\pi$ independent of $s$ and $p$, regular simple poles. The traces are invariant under a constant axial redefinition, but redefinitions that depend on $s$ are not covered ([Cartesian frame-invariance record](cartesian-sigma-pipi-loop-frame-invariance.md)).

## Vacuum condition

The Cartesian estimate sets the one-loop $s$ tadpole to zero, $\langle{\rm Re}\,\varphi\rangle=f_a$. The polar estimate sets the one-loop $\sigma$ tadpole to zero, $\langle|\varphi|\rangle=f_a$. These are two different definitions of $f_a$ at one loop. With $\varphi=f_a+s+ip$ (script output),
$$
|\varphi|=f_a+s+\frac{p^2}{2f_a}+\dots,\qquad \langle|\varphi|\rangle-\langle{\rm Re}\,\varphi\rangle=\frac{\langle p^2\rangle}{2f_a}\quad\text{(one loop, with }\langle s\rangle=0\text{)} .
$$
The tree is $M_0=-\frac{B_0}{2f_a}\langle\widehat S^C_{1,0}\{\lambda_{\pi^+},\lambda_{\pi^-}\}\rangle$ in both parametrizations, so $M_0\propto1/f_a$ changes when one definition is replaced by the other. **Size (assumption: size rule D1 with the D3 scale $m_a$ of the bare $p$ line).** $\langle p^2\rangle$ is a scalar moment with one denominator, so $\langle p^2\rangle\to m_a^2/(16\pi^2)$; only its size is meaningful. The one-loop relative change of $M_0$ is the term of first order in the loop counting,
$$
\frac{\delta M_0}{M_0}=-\frac{\delta f_a}{f_a}=-\frac{m_a^2}{32\pi^2f_a^2}\ \to\ -\frac{A}{32\pi^2f_\pi^2}\,\epsilon^4\qquad(\text{large }\epsilon,\ m_a^2\simeq\epsilon^2A,\ f_a=f_\pi/\epsilon),
$$
and $-B_0m_u\epsilon^4/(32\pi^2f_\pi^2)$ for the up-only sources. The higher orders in $\delta f_a$ belong to two and more loops and are not used; the script shows that keeping them would not give a one-loop quantity (the factor $f_a/(f_a+\delta f_a)$ saturates at large $\epsilon$).

- **Generic sources.** The leading power $\epsilon^4$ is the leading power of $R$ ([Cartesian large-$\epsilon$ graph record](cartesian-sigma-pipi-loop-large-epsilon-graphs.md)). The choice of vacuum condition therefore reaches the leading power of $R$: the two definitions of $f_a$ change $M_1^{\rm est}/M_0$ at the same order as the loop estimate itself.
- **Up-only sources.** $R\sim\epsilon^2$ ([Cartesian up-only record](cartesian-sigma-pipi-loop-up-only-estimate.md)), so the change from the vacuum condition is larger than $R$ by $\frac{\epsilon^2m_P^2}{4B_0m_u}$ at large $\epsilon$.

A comparison of the Cartesian and polar values of $M_1^{\rm est}/M_0$ must therefore state which definition of $f_a$ is used.

## Excluded chiral background (D4)

The all-generator-sum terms of $C3$, $C15$ and $C\pi3$ are not estimated (rule R6, [Cartesian graph-list record](cartesian-sigma-pipi-loop-graph-list.md)). They have the scale $m_P$, contain no ALP-sector letter, and their ratio to $M_0$ depends on $\epsilon$ only through $p_\sigma^2=m_\sigma^2$. At large $\epsilon$ their powers relative to $M_0$ are $\epsilon^0$, $\epsilon^2$ and $\epsilon^0$, below the leading power $\epsilon^4$ of the generic $R$ ([large-$\epsilon$ graph record](cartesian-sigma-pipi-loop-large-epsilon-graphs.md)). For the up-only sources, whose $R$ is of order $\epsilon^2$, the term of $C15$ is of the same power. Their size relative to the tree is set by $\max(m_P^2,m_\sigma^2)/(16\pi^2f_\pi^2)$; this is a scale, not a computed value.

## Size rule, loop scales and routing

- **Size rule (D1, D2).** Each scalar moment with $r$ denominators is replaced by $(m_{\rm term}^2)^{2-r}/(16\pi^2)$, and each product with the loop momentum by $m_{\rm term}^2$. No dependence on the external momenta beyond the polynomial numerator is kept.
- **On-shell cancellations.** Setting $k\cdot k\to m_{\rm term}^2$ cannot resolve cancellations that occur only on shell. In the Cartesian frame the $s\,p\,p$ vertex is a constant, with no momentum dependence: $V(s;p,p)=-2\lambda f_a+2B_0\epsilon^3\langle\widehat S^C_{1,2}\rangle/f_\pi$, made of $\lambda f_a=\epsilon\big(m_\sigma^2-\epsilon^2\kappa_\sigma\big)/(2f_\pi)$ and a source trace ([Cartesian vertex record](cartesian-sigma-pipi-loop-vertex-rules.md)). The only momentum-dependent vertex is the four-meson kinetic vertex, and it does not enter the leading power: the leading numerators contain neither $\widehat P^C_{0,1}$ nor $\widehat S^C_{0,0}$, which the meson components of the ALP-like lines would bring (checked in the large-$\epsilon$ tests). So the polar mechanism, a momentum-dependent radial–ALP vertex against a source vertex near the mass shell, has no counterpart at the leading power of the Cartesian estimate. The leading graphs can still cancel in the full loop: $C1$, $C2$, $C10$, $C11$, $C12$, $C\sigma4$, $C\sigma5$, $C\sigma6$ reach $\epsilon^4$ with numerators whose relative signs depend on the sources ("depends on sources" in the large-$\epsilon$ record). In particular $C2$, $C12$ and $C\sigma6$ all contain $\widehat S^C_{1,2}$, $C1$, $C10$ and $C\sigma4$ contain $\widehat S^C_{3,0}$, and $C11$, $C\sigma5$ contain $\widehat S^C_{2,1}$. Their loop functions differ (tadpoles at $m_\sigma$ or $m_a$, bubbles at $p_\sigma^2=m_\sigma^2$, residue derivatives), so the size rule cannot decide whether they cancel. This record does not assert that they do.
- **Routing (D12, assumption).** For each independent product of traces the routing with the largest $|C_t^\rho|$ is used, sign kept, with the tie rule of R9. This bounds the routing ambiguity of the size rule per graph, $|\sum_tC_t^\rho t|\le\sum_t|C_t^{\rho_t^*}||t|$ for every routing $\rho$, but not the true loop, and not the full $M_1^{\rm est}$, because graphs of opposite sign can cancel. Only $C12$, $C13$, $C15$ depend on the routing; 34 products of $M_1^{\rm est}$ have mass-dependent branches.
- **Loop scale (D3, assumption).** Each term uses the largest scale among its lines: $m_\sigma$ (radial), $m_a$ (ALP-like line, and pieces (ii)+(iii)+(iv) of a meson line), $m_P$ (bare generator sum). One symbolic $m_P$ for the nine meson-like poles is an approximation. Assigning $m_a$ to the rest of a meson line neglects the ALP–meson mass mixing inside the loops; it gives $\delta Z_{C\pi2}+\delta Z_{C\pi3}=0$ and the cancellation of $\widehat S^C_{0,0}$ in the sum of $C12$, $C13$, $C15$ ([estimate record](cartesian-sigma-pipi-loop-estimate.md)). No lower limit on the scale is imposed: the all-generator term of $C\sigma9$ has the scale $m_P<m_\sigma$, and for the up-only sources it alone sets the leading term $-B_0^2m_u^2\epsilon^2/(8\pi^2f_\pi^2m_P^2)$, which depends on this choice. For triangles the largest mass gives the smallest factor $(m^2)^{-1}$, so the choice is not conservative there.
- **Residues.** $\delta Z=-\Pi'$ with the size rule applied after the derivative; the prescription is an assumption, and the relative sign of the external-leg term and the 1PI graphs is not fixed by the size rule.

## Near degeneracy

$\mathbb M^{-1}_{m_a^2}(\mathcal C)$ has poles at the unmixed meson masses. Near an ALP–meson degeneracy the pole vector and $N(m_a^2)$ change rapidly; the formulas stay exact, but one loop scale per line is not adequate there ([finite-mixing limits record](../../../../meson-mixing/finite-mixing/pole-limits/alp-lo-finite-mixing-limits.md)). The boundary $r=\kappa_\sigma/A=1$ of the large-$\epsilon$ regions ($m_a\simeq m_\sigma$ at order $\epsilon^2$), where the up-only sources sit, is a degenerate mass relation that the size rule does not resolve.

## Validity of leading-order chiral perturbation theory and of $R\ll1$

- At large $\epsilon$ with fixed sources, $m_a^2\simeq\epsilon^2A$ and $m_\sigma^2\simeq\epsilon^2\kappa_\sigma$ grow, and the leading-order chiral Lagrangian is used at $m_a$, $m_\sigma$ beyond its usual range ($m_\sigma\lesssim4\pi f_\pi$). This is a limit of the estimate.
- For generic sources $R\sim\epsilon^4X/(16\pi^2f_\pi^2)$ with $X$ a combination of $A$, $\kappa_\sigma$ and $B_0$ times source traces, which equals $\epsilon^2X/(16\pi^2f_a^2)$, i.e. $m_a^2/(16\pi^2f_a^2)$ or $m_\sigma^2/(16\pi^2f_a^2)$ up to ratios of source traces (identity checked by the script). $R\ll1$ requires $m_a,m_\sigma\ll4\pi f_a$, which bounds $\epsilon$ from above; the large-$\epsilon$ results apply only below that bound. Since $f_a<f_\pi$ for $\epsilon>1$, $R\ll1$ fails before the condition of leading-order chiral perturbation theory.
- For the up-only sources $R\sim\epsilon^2B_0^2m_u^2/(8\pi^2f_\pi^2m_P^2)$ gives a weaker bound, but the vacuum-condition difference above is then larger than $R$.

## Absorptive parts and thresholds

At $p_\sigma^2=m_\sigma^2$ the $\pi\pi$ channel is open and, depending on the masses, also $K\bar K$, $\eta\eta$ and channels with ALP-like states. The size rule gives a magnitude only and carries no absorptive part. It can miss an enhancement near thresholds or from rescattering of the final pions, in particular the $S$-wave $\pi\pi$ rescattering at $s=m_\sigma^2$.

## Cancellations between graphs

$M_1^{\rm est}$ is a signed sum of graphs with exact relative signs from the vertex rules, but the size rule does not fix the relative magnitudes of different loop functions. Graphs with the same size factor can cancel or add in the full calculation, and the large-$\epsilon$ powers are given for orientation only; the estimate itself is the full exact sum, valid for any $\epsilon$ ([cancellation record](../../../../loop-corrections/finite-mixing/power-counting-and-validity/crossing-and-validity/alp-loop-estimate-cancellations.md)).

## Naturalness and counterterms

The terms not computed, including the counterterms needed to absorb the divergences of the loops, are assumed to be of the natural size given by the size rule ([naturalness record](../../../../loop-corrections/finite-mixing/power-counting-and-validity/loop-counting/alp-loop-naturalness-assumption.md), [radial naturalness record](../../loop-estimates/power-counting/murayama-radial-naturalness-assumption.md)). No counterterm basis is written and no renormalization-scale analysis is done. A counterterm of unnaturally large or small size would change $R$.

## Verification

The record states limits and assumptions; it contains no value of a loop (tested: no loop-function name and no component notation). The vacuum-condition statements are the output of the script and are tested: the expansion of $|\varphi|$, $\delta f_a=\langle p^2\rangle/(2f_a)$, the one-loop relative change $-m_a^2/(32\pi^2f_a^2)$, its large-$\epsilon$ term $-A\epsilon^4/(32\pi^2f_\pi^2)$ and the up-only term, the comparison with the leading powers of $R$ (generic $\epsilon^4$, up-only $\epsilon^2$ with the factor $\epsilon^2m_P^2/(4B_0m_u)$), the saturation of the all-orders factor (negative control), the identity behind the scale of $R$, and the powers of the D4 terms for the up-only sources ($\epsilon^0$, $\epsilon^2$, $\epsilon^0$ for $C3$, $C15$, $C\pi3$). The other statements are checked by the scripts of the linked records.

## Reproduction

```sh
PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_cartesian_loop_limits.py
PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_cartesian_loop_limits
```

[Script](../../../../../../scripts/derive_sigma_pipi_cartesian_loop_limits.py), [tests](../../../../../../scripts/test_derive_sigma_pipi_cartesian_loop_limits.py).
