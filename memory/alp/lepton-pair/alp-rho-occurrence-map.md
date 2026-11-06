---
title: Rho occurrence and cut map for the low-energy ALP TFF
status: mixed
created: 2026-11-06
updated: 2026-11-06
bibtex_keys: [Alda:2026nrd, Messerli:2025rnv, Ecker:1988te, Ecker:1989yg, Ruiz-Femenia:2003jdx, GomezDumm:2013sib, Watson:1952ji, Omnes:1958hv]
provenance: literature
---

# Scope and conventions

This record maps the rho-dependent terms already present in the indexed ALP
TFF derivations. It distinguishes the pion-current form factor, the
production partial wave used on a two-pion unitarity cut, the photon-leg
response, and separate TFF pole terms. The map records their analytic
variables and phase conventions; it does not select a new running width or
derive a matching prescription. Use the unnormalized scalar TFF, the
no-mixing ALP sources, natural units, metric $(+---)$, and
$\epsilon^{0123}=+1$ [@Alda:2026nrd; @Messerli:2025rnv].

# Occurrence map

| Object | Rho dependence and invariant | Cut or phase content | Status and source |
|---|---|---|---|
| Pion vector form factor $F_\pi^V(s)$ | Enters the vacuum-to-$\pi^+\pi^-$ current matrix element with $s=(p_++p_-)^2$. The elastic representation is $F_\pi^V(s)=P(s)\Omega_1^1(s)$; it contains rho dynamics through the measured/modelled $\pi\pi$ P-wave phase, not an added stable rho pole. | In the elastic isovector region, $\arg F_\pi^V(s+i0)=\delta_1^1(s)$. The Omnès right-hand cut begins at $4M_\pi^2$; phase continuation, subtractions, and inelastic effects are inputs. | Dispersive current representation. A stable tree rho pole is an alternative approximation and needs matching before it can be combined with this current form factor. See [pion current and production wave](alp-pion-vector-and-production-form-factors.md). |
| Dipion $s$-channel in $a\gamma^*\to\pi^+\pi^-$ | The retained tree production scalar contains $N_s(s,y)/(s-M_\rho^2+i0)$, with $s=(p_++p_-)^2$ and photon virtuality $y$. | This stable tree pole is not the physical $\pi\pi$ discontinuity. In the full production partial wave, rho enhancement and elastic rescattering belong to the dipion right-hand-cut treatment. The tree expression does not supply that rescattering phase. | Symbolic source-resolved tree kernel. Its width or replacement by an Omnès/rescattered representation is not selected. See [pion current and production wave](alp-pion-vector-and-production-form-factors.md) and [two-pion spectral construction](alp-isovector-tff-spectral-construction.md). |
| Crossed $t/u$ terms in the production wave | The retained scalar contains $N_t(t,y)/(t-M_\rho^2+i0)+N_u(u,y)/(u-M_\rho^2+i0)$, where $t=(k+p_-)^2$ and $u=(k+p_+)^2$. Angular projection yields logarithms with the stated $i0$. | These are crossed-channel singularities of the production kernel, not the dipion $s$-channel unitarity cut. Their analytic continuation and phases must be kept distinct from the elastic $F_\pi^V(s)$ phase. | Stable tree poles with an explicit boundary prescription. The records do not justify replacing every crossed denominator by the same running-width expression used on the dipion cut. See [pion current and production wave](alp-pion-vector-and-production-form-factors.md). |
| Photon response $W_\gamma(y)$ | The production amplitude also contains the full neutral-vector response $W_\gamma(y)=\sum_i R_iT_i/(M_i^2-y-i0)$, including a rho response at the photon virtuality $y=k^2$. | This is a photon-leg vector response. It is distinct from the two-pion cut in $s$ and from the crossed $t/u$ poles. Its tree expression does not itself provide a $\pi\pi$ spectral density or the production-wave rescattering phase. | Full neutral-vector matrix at the retained tree level; no rho-specific finite-width prescription is selected for this slot. See [pion current and production wave](alp-pion-vector-and-production-form-factors.md). |
| Two-pion TFF discontinuity | The charged two-pion density is proportional to $[F_\pi^V(x)]^*\mathcal F_{a\pi\pi\gamma^*}(x,y)$, with $\mathcal F_{a\pi\pi\gamma^*}=f_{1,\mathrm{current}}=f_1/e$ in the declared real-epsilon convention. | The single-cut coefficient is $x\sigma_\pi^3/(96\pi)$; the per-ordering double-spectral density uses $x\sigma_\pi^3/(192\pi)$ and the exchanged photon ordering. The production wave is an input to this unitarity product, not an independent additive TFF term. Its phase must be compatible with the current phase. | Source-specific low-energy isovector construction. The double-spectral term may require additional single-spectral/subtraction pieces; those are not fixed by this map. See [two-pion spectral construction](alp-isovector-tff-spectral-construction.md). |
| Stable neutral-rho TFF pole | The retained tree TFF contains VJP terms with one $W_\rho(s_i)$ and VVP terms with $W_\rho(s_1)W_j(s_2)$, including $\rho\rho$ and mixed $\rho I$ assignments. | This is a pole approximation to rho-mediated TFF strength. The selected low-energy isovector representation uses the two-pion spectrum instead of an additional stable rho pole for the same contribution. | The tree rho-pole record is a diagnostic component, not an additive term on top of the two-pion spectral component without matching. See [neutral-rho pole contribution](alp-isovector-doubly-virtual-tff.md). |
| Mixed dressed $F^{(\rho I)}$, $I=\omega,\phi$ | The mixed VVP term has a rho response on one photon leg and a narrow isoscalar pole on the other, summed over both photon orderings. The rho response is $R_\rho/[M_{0\rho}^2-s-\Delta_\rho^{\mathrm{self}}(s)]$. | The conditional self-energy convention is $\operatorname{Im}\Delta_\rho^{\mathrm{self}}(s+i0)=+\sqrt{s}\,\Gamma_{\rho\to\pi\pi}(s)$, with the real part defined by subtractions. A unit-projector width and alternatives are derived, but no numerical parameters or unique off-shell continuation are selected. This response is not the full two-pion double-spectral TFF. | The pure $\omega/\phi$ pole selection does not contain this mixed VVP term. However, the mixed term and the two-pion double-spectral reconstruction can share rho-mediated strength; their source-specific matching is open. See [mixed rho-isoscalar VVP](alp-tff-mixed-rho-isoscalar-vvp.md), [conditional vector self-energy](alp-isoscalar-vector-self-energy.md), and [rho-to-two-pion self-energy](alp-rho-pipi-self-energy.md). |

# Phase and overlap conclusions

The two-pion density depends on the relative phase in
$[F_\pi^V]^*\mathcal F_{a\pi\pi\gamma^*}$. In a strictly elastic model
with matching Watson phases, these phases cancel in the product; the
retained tree production kernel alone does not establish this physical
phase. The crossed $t/u$ logarithms and the photon response at $y$ carry
their own boundary prescriptions and must not be assigned the dipion phase
by analogy [@Watson:1952ji; @Omnes:1958hv].

The indexed derivations establish three separate statements: the stable rho
pole is not an extra term on top of the same two-pion spectral contribution;
the rho-containing and pure $\omega/\phi$ selections of the retained tree
TFF intersect in the mixed $\rho I$ VVP terms; and the mixed dressed term
may share rho strength with the two-pion double-spectral reconstruction.
The first two are explicit component-bookkeeping results. The third is a
possible continuum overlap whose quantitative decomposition and subtraction
are not established. No common running-width prescription across $s,t,u,y$
follows from these records [@Messerli:2025rnv]. The conditional width options
and their subtraction limits are recorded in
[the rho-to-two-pion self-energy derivation](alp-rho-pipi-self-energy.md).

## Source records

- [Pion current and production wave](alp-pion-vector-and-production-form-factors.md)
- [Two-pion spectral construction](alp-isovector-tff-spectral-construction.md)
- [Neutral-rho pole contribution](alp-isovector-doubly-virtual-tff.md)
- [Mixed rho-isoscalar VVP term](alp-tff-mixed-rho-isoscalar-vvp.md)
- [Isoscalar source decomposition](alp-isoscalar-tff-source-decomposition.md)
- [Conditional vector self-energy](alp-isoscalar-vector-self-energy.md)

These records are all indexed in this directory's `INDEX.md`. This map
introduces no new numerical width, subtraction constant, or overlap model.
