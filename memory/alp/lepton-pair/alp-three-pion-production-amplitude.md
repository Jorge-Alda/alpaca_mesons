---
title: Retained source-resolved ALP production current to three pions
status: derived
created: 2026-10-26
updated: 2026-10-26
bibtex_keys: [Alda:2026nrd, Ecker:1989yg, Ruiz-Femenia:2003jdx, Niecknig:2012sj, Hoferichter:2025lcz]
provenance: analytical-derivation
---

# Scope and conventions

This record assembles the gauged leading-order Born current and the strict
spin-one $ω/φ$ pole residues for
$a(p)+\gamma^*(q)\to\pi^+(p_+)\pi^-(p_-)\pi^0(p_0)$. It uses the independent
no-mixing sources and odd couplings of the
[source-resolved TFF](alp-doubly-virtual-tff.md), with natural units, metric
$(+---)$, $\epsilon^{0123}=+1$, $S=1+iT$, and incoming-photon routing
$p+q=k=p_++p_-+p_0$ [@Alda:2026nrd; @Ecker:1989yg;
@Ruiz-Femenia:2003jdx]. The photon coupling $e$ is included.

The result is the retained Born term plus pole residues. It is not a complete
off-shell five-point amplitude: a three-pion $I=0$, $J^{PC}=1^{--}$ projection
of the full-charge Born current, the normalized complex decay functions,
rescattering, and regular production terms have not been supplied. The SM
$\gamma^*\to3\pi$ current amplitude and the ALP production amplitude are
distinct matrix elements [@Hoferichter:2025lcz].

# Executable derivation

Run

```text
.venv/bin/python scripts/derive_alp_three_pion_production_amplitude.py
```

The script reuses the executed source traces in the
[LO Born-current derivation](alp-three-pion-current-born-amplitude.md) and
the ordered VJP/VVP contractions in the
[spin-one pole-factorization derivation](alp-three-pion-vector-pole-factorization.md).
For $D_\pm=q^2-2q\cdot p_\pm$, the Born term is

$$
J^\mu_{\rm Born}=e\,\mathcal V_{a\pi^+\pi^-\pi^0}
\left[\frac{(2p_+-q)^\mu}{D_++i0}
-\frac{(2p_--q)^\mu}{D_-+i0}\right],
$$

where $\mathcal V_{a\pi^+\pi^-\pi^0}$ retains the two independent no-mixing
ALP source matrices as given in the linked Born-current record. The gauged
contact coefficient and the difference of the two shifted source vertices
vanish after the diagonal-source traces are evaluated. SymPy contracts the
current with $q_\mu$ and obtains zero.

Let $s=k^2$, $y=q^2$, and use the established ordered production numerator
$N_i(s,y)$ from the linked pole-factorization record, retaining the full
charge matrix and all intermediate $j=\rho,\omega,\phi$ photon responses.
Define
$$
\mathcal E_i^{\alpha\mu}=\epsilon^{\alpha\mu\rho\sigma}k_\rho q_\sigma,
\qquad
J_{3\pi}^\beta=\epsilon^{\beta\lambda\kappa\delta}
p_{+\lambda}p_{-\kappa}p_{0\delta},
\qquad
\Pi^{(i)}_{\alpha\beta}=-g_{\alpha\beta}
+\frac{k_\alpha k_\beta}{M_i^2}.
$$

The strict pole terms, with $i=\omega,\phi$, are
$$
J^\mu_{i,\rm pole}=
-\frac{iF_i^{\rm on}(t,u)N_i(M_i^2,y)}
{M_i(M_i^2-s-i0)}
\left[\mathcal E_i^Tg\Pi^{(i)}gJ_{3\pi}\right]^\mu.
$$

The factor $\zeta_i=-1$ is the conventional incoming-photon crossing/LSZ
factor already matched in
[the pole-phase derivation](alp-three-pion-pole-phase.md). The physical
complex phase and normalization of $F_i^{\rm on}$ remain inputs. Since
$k\cdot J_{3\pi}=0$, the longitudinal part of $\Pi^{(i)}$ drops. SymPy
verifies the photon and vector Ward identities for $\mathcal E_i$, the
three-pion-current transversality, and the Born and pole photon Ward
contractions. Thus the retained sum
$$
J^\mu_{\rm retained}=J^\mu_{\rm Born}
+J^\mu_{\omega,\rm pole}+J^\mu_{\phi,\rm pole}
$$
is transverse in the photon index.

# Isospin projection and limitations

The Born current is constructed with the full electromagnetic generator.
The result $[Q^{(0)}_\ell,T_\pm]=0$ only removes the minimal light-quark
isoscalar-photon coupling; it does not project the three-pion state. The
full-charge Born contribution to the $I=0$, $J^{PC}=1^{--}$ three-pion state
therefore remains unassigned here. Do not discard it or call it zero before
that state projection is derived.

Only the strict vector-pole residues are included. The charged-pion Born
singularities are explicit, but their unitarity dressing and compatible
subtractions are not derived here. The result also does not determine
$F_i^{\rm on}$, its correlated variations, a physical phase from rate-only
data, finite widths, or regular off-shell production terms. Consequently it
is not yet a numerical three-pion spectral contribution to the doubly
virtual TFF. The next physical input is the rescattered, complex, normalized
$V_i\to3\pi$ decay function together with the source-side state projection
and a consistent pole/continuum prescription [@Niecknig:2012sj].
