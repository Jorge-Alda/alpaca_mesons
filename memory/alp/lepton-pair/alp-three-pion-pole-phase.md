---
title: Crossing and LSZ phase of the strict ALP three-pion vector pole
status: derived
created: 2026-10-25
updated: 2026-10-25
bibtex_keys: [Alda:2026nrd, Ecker:1989yg, Ruiz-Femenia:2003jdx, Niecknig:2012sj]
provenance: analytical-derivation
---

# Inputs and declared conventions

Use the signed no-mixing production amplitude in the
[$a\to V_i\gamma$ record](../../rcht/odd-parity/alp-vgamma-amplitude.md),
the ordered numerator and tensor normalization in the
[strict-pole record](alp-three-pion-vector-pole-factorization.md), and the
decay convention in the [partial-wave record](alp-three-pion-vector-partial-wave.md).
The original odd operators and vector tensor normalization are those already
used there [@Alda:2026nrd; @Ecker:1989yg; @Ruiz-Femenia:2003jdx]. The
decay convention is the one recorded from Niecknig:2012sj
[@Niecknig:2012sj]. No new literature formula or physical fit input is used.

Natural units and GeV are used, with metric $(+---)$,
$\epsilon^{0123}=+1$, hence $\epsilon_{0123}=-1$, and $S=1+iT$.
Choose the same ordinary neutral-photon LSZ phase in production and its
crossed amplitude. Define the incoming photon momentum by the existing
scattering routing $k=p_a+q$, $s=k^2$, $y=q^2$. Its outgoing momentum
in the crossed decay amplitude is $-q$. No additional bosonic crossing
phase is introduced. These explicit choices specify the convention matched
below; a simultaneous rephasing of external states must also be applied to
their source matrix elements.

# Executed numerator and crossing match

The paired [script](../../../scripts/derive_alp_three_pion_pole_phase.py)
reuses the original Lorentz contractions, inserts the independent no-mixing
sources and full neutral photon response, and obtains for both
$i=\omega,\phi$
$$
\frac{N_i(M_i^2,0)}{M_i}-G_{V_i\gamma}=0.
$$
This comparison is with the existing signed physical amplitude, not with a
decay width. It retains independent axial and pseudoscalar sources and the
external-vector normalization $M_i$, distinct from the odd-operator scale
$M_V$. No fit representative or equality of source traces is imposed.

The stored decay amplitude has
$$
\mathcal M(a\to V_i\gamma)
=G_{V_i\gamma}\epsilon_{\mu\nu\rho\sigma}
n^{*\mu}\varepsilon^\nu k^\rho q_\gamma^\sigma.
$$
SymPy evaluates the same epsilon contraction with $q_\gamma=-q$ and
with $q$ and derives their ratio $-1$. Thus, in the incoming-photon
scattering routing, the coefficient multiplying
$\epsilon(n^*,\varepsilon,k,q)$ is
$$
-\frac{N_i(M_i^2,y)}{M_i}.
$$
The real-photon normalization fixes the sign; the virtual-current ordered
numerator is continued with the same convention. This is a continuation of
the production residue, not a determination of regular off-shell decay
terms.

# Propagator composition and conventional phase

Retain the declared signed decay amplitude
$$
\mathcal M(V_i\to3\pi)=iF_i^{\rm on}(t,u)\,n_\beta J_{3\pi}^\beta.
$$
The script multiplies $i\mathcal M_{\rm prod}$, the vector propagator
$i\Pi^{(i)}/(s-M_i^2+i0)$, and $i\mathcal M_{\rm dec}$, and divides
by $i$ to recover $\mathcal M$. It also constructs the epsilon tensor
and spin sum explicitly. Both production Ward contractions vanish and
the longitudinal projector term drops. In the notation of the strict-pole
record, the result is
$$
\mathcal M^\mu_{i,\rm pole}
=-\frac{iF_i^{\rm on}(t,u)N_i(M_i^2,y)}
{M_i(M_i^2-s-i0)}
\left[\mathcal E_i^Tg\Pi^{(i)}gJ_{3\pi}\right]^\mu.
$$
Comparison with that record's explicit $i\zeta_i$ coefficient gives
$$
\zeta_i=-1
$$
for the specified incoming-photon routing and signed decay convention.
Changing to an outgoing photon tensor reverses this displayed crossing
factor; it does not change the full amplitude when all definitions are
converted together. The script returns zero for the difference between
the composed amplitude and the recorded strict-pole expression with this
matched value.

# Physical decay normalization still needed

Fixing this conventional phase does not fix the normalization phase of a
decay solution obtained only from widths and Dalitz intensities. The script
checks invariance of the rate under multiplication of $F_i^{\rm on}$ by
an arbitrary unit-modulus phase. Such a multiplication changes the
amplitude relative to the signed production source, although the decay
intensity is unchanged. A common signed decay Lagrangian, an interference
measurement, or a specified matched current-amplitude solution is needed
to determine that remaining relative phase. It must be carried in the
physical decay function rather than silently removed after setting the
conventional $\zeta_i$.

No complex decay grid, correlated variations, direct contact terms or
regular off-shell terms are supplied. The user's scope
$2m_e<m_a<2m_\tau$ and the unbounded lepton-loop virtualities therefore
still require the additional amplitudes and high-energy continuation in the
[isoscalar plan](alp-isoscalar-tff-plan.md).

# Reproducibility and checks

Executed:

```text
.venv/bin/python scripts/derive_alp_three_pion_pole_phase.py
.venv/bin/python scripts/derive_alp_three_pion_pole_phase.py --test
```

Four unit tests pass: signed real-photon numerator matching for omega and
phi, explicit boson crossing and both Ward contractions, spin-sum and
propagator composition, and invariance of the decay intensity under its
global normalization phase. The script prints intermediate and final
expressions through SymPy's LaTeX renderer. These checks establish the
declared conventional phase only; they do not measure the physical decay
normalization phase.
