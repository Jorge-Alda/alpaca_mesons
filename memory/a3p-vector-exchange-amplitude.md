---
title: Canonical vector-pole contribution to ALP three-pseudoscalar amplitudes
status: mixed
created: 2026-09-24
updated: 2026-09-24
bibtex_keys: [Alda:2026bar, Ecker:1988te, Ecker:1989yg, Ruiz-Femenia:2003jdx, ParticleDataGroup:2024cfk]
provenance: analytical-derivation
---

# Canonical Vector Exchange in $a\to3P$

## Scope, conventions, and inputs

The executable derivation is [`scripts/derive_a3p_vector_exchange_amplitude.py`](../scripts/derive_a3p_vector_exchange_amplitude.py), and the numerical implementation is [`scripts/implementation/vector_meson_exchange.py`](../scripts/implementation/vector_meson_exchange.py). Natural units are used; all masses, widths, $f_\pi$, $G_V$, and $f_a$ are in GeV.

The hatted Hermitian $c_q^A$ and $c_q^P$ are the no-ALP-meson-mixing inputs defined in [the Alda Appendix A extract](extracts/Alda:2026bar/appendix-a-expansion.md) [@Alda:2026bar]. The leading $G_V$ pole depends on $c_q^A$. The implementation validates $c_q^P$ to retain the public interface of the existing three-meson amplitude, but $c_q^P$ does not enter this contribution.

The pseudoscalar generators, masses, aliases, and coherent $K_S/K_L$ convention are imported directly from [`three_meson_width.py`](../scripts/implementation/three_meson_width.py), whose conventions are recorded in [the generic-generator derivation](alp-three-meson-generic-generators.md). In particular, the Alda-normalized matrix has off-diagonal generators $\sqrt2E_{ij}$ and its fitted neutral generators, including $\Lambda_{\pi^0}\simeq\operatorname{diag}(1,-1,0)$. Here

$$
u=\exp\!\left(\frac{i\Phi_{\rm Alda}}{2f_\pi}\right),\qquad
U=u^2=\exp\!\left(\frac{i\Phi_{\rm Alda}}{f_\pi}\right).
$$

This is equivalent to the RChT convention

$$
U=\exp\!\left(\frac{i\sqrt2\Phi_R}{F}\right),\qquad
\Phi_{\rm Alda}=\sqrt2\Phi_R,\qquad F=f_\pi,
$$

used for the leading $G_V$ operator in [the RChT record](rcht-basics.md) [@Ecker:1988te].

## Explicit model assumptions

1. **Ideal vector mixing:** the physical-vector generators are
   $T_{\rho^0}=\operatorname{diag}(1,-1,0)/\sqrt2$,
   $T_\omega=\operatorname{diag}(1,1,0)/\sqrt2$,
   $T_\phi=\operatorname{diag}(0,0,1)$, and unit-normalized $E_{ij}$ for charged $\rho$ and $K^*$ states. Ideal $\omega$--$\phi$ mixing is a model approximation, not a symmetry requirement [@Ecker:1988te]. The conjugate vector uses $T_{\bar V}=T_V^\dagger$.
2. **Single-resonance short-distance value:** by default
   $$G_V=f_\pi/\sqrt2,$$
   following the single-resonance short-distance constraints summarized in the existing Ruiz-Femenia extract [@Ruiz-Femenia:2003jdx]. `g_v` can override it.
3. **Constant widths:** the propagator uses $M_V^2-s-iM_V\Gamma_V$. Default $M_V,\Gamma_V$ values are the PDG listing central values transcribed in [the light-vector input record](extracts/ParticleDataGroup:2024cfk/light-vector-pole-parameters.md) [@ParticleDataGroup:2024cfk]. A constant-width line shape is an explicit implementation assumption, not a PDG result. `pole_parameters` can override each named state.
4. **Pole-only scope:** the calculation retains only the canonical Proca propagator between the $WPa$ and $WPP$ vertices inherited from the $G_V$ operator. Eliminating the original antisymmetric tensor also induces a local cross term between its $Pa$ and $PP$ currents, schematically $J_{Pa}^{\mu\nu}J^{PP}_{\mu\nu}$. This is a direct, non-pole $aP^3$ interaction of $O(p^4)$ in Goldstone chiral counting. It is mandatory if the Proca formulation is meant to reproduce the complete antisymmetric-tensor result, rather than an additional optional resonance coupling [@Ecker:1989yg]. The present module omits it because it implements only the resonant pole contribution. The flavor-resolved coefficient of this matching contact has not yet been derived in the repository conventions.

   This matching contact is distinct from the leading $O(p^2)$ Goldstone-only ChPT amplitude implemented in `three_meson_width.py`; adding both does not by itself double count. Double counting would arise if one also used $O(p^4)$ ChPT low-energy constants that already contain the same vector-resonance contribution without first matching or subtracting it. Conversely, one should use either the original antisymmetric-tensor exchange or the equivalent package “Proca pole plus induced contact,” not add all three pieces independently.

   Further effects not included here are genuinely model-dependent vector-sector corrections: higher-derivative and explicit-$SU(3)$-breaking resonance operators, nonideal neutral-vector mixing, heavier vector multiplets, multi-resonance operators, loops and final-state rescattering, and improved energy-dependent or dispersive line shapes. These are logically separate from the mandatory representation-matching contact. At strong leading order, $F_V$ does not supply another $VPP$ or $VPa$ vertex; its photon--vector mixing effects are electromagnetic and begin beyond the $e^0$ two-$G_V$ contribution considered here.

## Alda-normalization conversion and sign

The RChT interactions in the $\Phi_R$ convention are derived in [the flavor-expansion record](rcht-vpp-vpa-photon-pa-expansion.md). Substituting each $\lambda_R=\Lambda_{\rm Alda}/\sqrt2$ gives, for distinct labeled components and one fixed pair $i<j$,

$$
\mathcal L_{WPP}
=\frac{iG_V}{\sqrt2f_\pi^2M_V}
W_{\mu\nu}\partial^\mu P_i\partial^\nu P_j
\operatorname{Tr}\!\left(T_{\bar V}[\Lambda_i,\Lambda_j]\right),
$$

$$
\mathcal L_{WPa}
=-\frac{\sqrt2iG_V}{f_\pi f_aM_V}
W_{\mu\nu}\partial^\mu P_k\partial^\nu a
\operatorname{Tr}\!\left(T_V[\Lambda_k,c_q^A]\right).
$$

The SymPy script performs this substitution rather than assuming the converted coefficients. Their product is

$$
\frac{G_V^2}{f_\pi^3f_aM_V^2}.
$$

With the script's $e^{-ip\cdot x}$ derivative convention, the scalar factors from the two Feynman vertices are $-ia$ and $+ib$ for positive coefficient magnitudes $a,b$. Their product is $+ab$. Together with the propagator convention this yields the positive amplitude sign below. This audit fixes only the executable scalar $i$ bookkeeping; field-wide phase conventions remain conventional when transformed consistently.

## Kinematic derivation

For spectator $k$, complementary ordered pair $i<j$, $q=p_i+p_j$, and $s=q^2$, define

$$
J_a=(q\cdot p_k)p_a-(q\cdot p_a)p_k,
\qquad
J_{PP}=(q\cdot p_i)p_j-(q\cdot p_j)p_i.
$$

The SymPy Gram-matrix calculation gives

$$
q\cdot J_a=q\cdot J_{PP}=0
$$

and

$$
J_a\cdot J_{PP}=\frac{s}{4}\left[
(s+m_j^2-m_i^2)(s_{ik}-m_i^2-m_k^2)
-(s+m_i^2-m_j^2)(s_{jk}-m_j^2-m_k^2)
\right].
$$

Thus one labeled contribution is

$$
\mathcal M_{V,k;ij}=
+\frac{G_V^2}{f_\pi^3f_aM_V^2}
\frac{X_{V;ij}Y_{V;k}\,J_a\cdot J_{PP}}
{M_V^2-s-iM_V\Gamma_V},
$$

where

$$
X_{V;ij}=\operatorname{Tr}\!\left(T_{\bar V}[\Lambda_i,\Lambda_j]\right),
\qquad
Y_{V;k}=\operatorname{Tr}\!\left(T_V[\Lambda_k,c_q^A]\right).
$$

Reversing $i,j$ changes the signs of both $X$ and $J_{PP}$, leaving the contribution invariant. The script asserts this, transversality, linear $c_q^A$ and inverse-$f_a$ scaling, and zero total mass dimension of the four-point amplitude.

## Identical particles and Mandelstam dependence

This pole amplitude does **not** inherit the one-variable reduction of the leading chiral contact amplitude recorded in [the three-meson amplitude derivation](alp-three-meson-amplitudes-widths.md). The executable script derives the following statements.

For three identical pseudoscalars, every spectator assignment leaves an identical pair at the $VPP$ vertex. Therefore

$$
[\Lambda_P,\Lambda_P]=0
$$

and the complete canonical vector-pole contribution vanishes. It is consequently independent of all Mandelstam variables only in the trivial sense that it is zero.

For

$$
a\to P(p_1)P(p_2)Q(p_3),
$$

the $Q$-spectator graph also vanishes because its $VPP$ pair is $P(p_1)P(p_2)$. The two labeled $P$-spectator graphs survive together or vanish together for a given vector flavor factor. Their pair invariants and pole denominators are respectively

$$
s_{23},\quad \Delta_V-s_{23},
\qquad\text{and}\qquad
s_{13},\quad \Delta_V-s_{13},
$$

where $\Delta_V=M_V^2-iM_V\Gamma_V$. Writing their kinematic numerators as $N_{23}$ and $N_{13}$, one vector contribution is proportional to

$$
\frac{N_{23}}{\Delta_V-s_{23}}
+\frac{N_{13}}{\Delta_V-s_{13}}.
$$

This is Bose symmetric under $p_1\leftrightarrow p_2$, or $s_{13}\leftrightarrow s_{23}$, but Bose symmetry does not imply dependence on $s_{12}$ alone. Define

$$
\Sigma=m_a^2+2m_P^2+m_Q^2,
\qquad
S=s_{13}+s_{23}=\Sigma-s_{12},
\qquad
R=s_{13}s_{23}.
$$

SymPy reduces the numerator sum to

$$
N_{13}+N_{23}=F(s_{12};m_a^2,m_P^2,m_Q^2)
-\frac{\Sigma+s_{12}}{4}R,
$$

for a mass-dependent polynomial $F$. Thus even replacing both propagators by the same constant does not produce a function of $s_{12}$ alone. In the physical decay region $\Sigma+s_{12}>0$, so the coefficient of $R$ cannot vanish. Equal final masses and $\Gamma_V=0$ do not alter this conclusion. With finite pole denominators, the common denominator is

$$
(\Delta_V-s_{13})(\Delta_V-s_{23})
=\Delta_V^2-\Delta_V S+R,
$$

which adds further $R$ dependence. The script finds a nonzero, parameter-independent quadratic coefficient in the numerator of $\partial\mathcal M_V/\partial R$, proving that a nonzero finite single-vector pole contribution cannot become a function of $s_{12}$ alone by choosing masses or $\Delta_V$.

The precise exceptions are: an exact zero of the coupling or either flavor trace; exact cancellation of complete residues among degenerate exchanges with tuned flavor weights; isolated kinematic zeros; or restriction to a one-dimensional locus such as $s_{13}=s_{23}$. The algebraic condition $\Sigma+s_{12}=0$ removes the $R$ term from the numerator sum but lies outside physical decay kinematics and still does not remove finite-propagator $R$ dependence. Cancellations among nondegenerate poles are not a generic symmetry reduction and require all relevant pole residues to cancel.

This statement concerns the **amplitude**. A total width may always be represented as an outer one-dimensional integral after numerically performing the inner Dalitz integration over $s_{23}$ at fixed $s_{12}$. That does not make the pole amplitude one-variable, and unlike the chiral contact $PPQ$ amplitude, the inner integration is not merely multiplication by the Dalitz-interval width.

## Selection and implementation details

For every labeled spectator index, the implementation retains the complementary pair in increasing index order. If $c_q^A\propto\mathbf1$, it returns an empty contribution list before exchange-specific work. For diagonal $c_q^A$, diagonal $\rho^0,\omega,\phi$ directions are skipped before either flavor trace according to [the neutral-vector selection rule](flavor-conserving-neutral-vector-selection-rule.md). For each remaining candidate it tests $X$ first and rejects exact or scale-relative roundoff zeros before computing $Y$ or kinematics. Repeated final fields are not compressed into multiplicities, because their labeled pair invariants differ. Importing the existing $K_S/K_L$ generators preserves coherent neutral-kaon sums automatically.

For an isospin-symmetric but nonidentity diagonal spurion, $c_q^A=\operatorname{diag}(c_\ell,c_\ell,c_s)$, the complete $\pi^+\pi^-\eta$ contribution vanishes: the charged-pion spectator commutators vanish, while the diagonal $\eta$ spectator cannot couple through a diagonal spurion. For a generic diagonal spurion with $c_u\ne c_d$, the fitted Alda $\eta$ generator's unequal $u,d$ entries produce small charged-$\rho$ isospin-breaking traces. Setting those entries equal would be an additional isospin-limit approximation and is not silently imposed.

The public model is called as `model(s12, s23)` with exactly

$$
s_{12}=(p_1+p_2)^2,\qquad s_{23}=(p_2+p_3)^2,
$$

and derives $s_{13}=m_a^2+\sum_i m_i^2-s_{12}-s_{23}$. `model.contributions` exposes each exchange name, spectator index, pair indices, traces, and pole inputs; `model.exchange_names` is a concise name tuple.

## Invocation and verification

```text
.venv/bin/python scripts/derive_a3p_vector_exchange_amplitude.py
.venv/bin/python -m unittest scripts.implementation.test_vector_meson_exchange
```

The implementation tests cover identity and diagonal-vector zeros, the full $\pi^+\pi^-\eta$ zero for an isospin-symmetric nonidentity diagonal spurion, allowed charged-$\rho$ exchange in $\pi^+\pi^-\pi^0$, $c_q^P$ validation and independence, $c_q^A/f_a$ scaling, ordered permutation covariance including repeated fields, pair reversal, finite default-width behavior at a pole, and contribution metadata.
