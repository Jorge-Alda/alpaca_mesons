---
title: ELI5 explanation of unitarity cuts and vector dispersive treatments in the ALP TFF
status: mixed
created: 2026-10-27
updated: 2026-10-27
bibtex_keys: [Alda:2026nrd, Messerli:2025rnv]
provenance: assumption
---

# The short picture

A unitarity cut is what appears when an intermediate set of real particles
can go on shell. Think of a virtual photon with enough energy to make two
pions: those pions can exist as real intermediate particles before the
process continues. Above their threshold, the amplitude develops a
continuous branch cut and an imaginary part, rather than only an isolated
pole.

Unitarity says how to calculate this contribution. Multiply the amplitude
for one side of the process to make the intermediate particles by the
amplitude for the other side to absorb them, and sum over their allowed
momenta and quantum numbers. In the TFF, one side is the electromagnetic
current making hadrons from the vacuum; the other is the ALP transition
amplitude connecting those hadrons to the other photon. A dispersion
relation then uses this spectral strength to reconstruct the TFF away from
the cut.

# Why the first isovector and isoscalar thresholds differ

In the exact-isospin electromagnetic-current decomposition used here, the
isovector vector current has the quantum numbers to create two pions in a
relative $P$ wave. The threshold is therefore $s=4M_\pi^2$. The isoscalar
vector current cannot create this two-pion state in the isospin-symmetric
limit; its first allowed light-hadron state is three pions, with threshold
$s=9M_\pi^2$. The current's quantum numbers are like an entry requirement
for intermediate states, and the lightest allowed state sets the first
threshold. Higher states can still contribute above their thresholds.
These statements are for the stated isospin limit; isospin breaking and
additional channels can modify the simple picture.

# How poles and cuts fit together

In general, a TFF contains both poles and cuts. A cut is the continuum of
allowed on-shell states. A vector resonance such as the $\rho$ lies within
that continuum and enhances the two-pion amplitude around its resonance
mass. In an exact treatment, the resonance and the two-pion cut are parts
of the same analytic structure.

A stable tree-level vector pole with denominator
$1/(M_\rho^2-s)$ is an approximation because it omits the $\pi\pi$
self-energy that dresses the propagator, along with other rescattering and
nonresonant contributions. Its virtue is that it isolates the retained
vector-exchange terms and their source and flavor structure. It is not the
full isovector spectral function. Adding it directly to a dispersive
two-pion contribution can count some $\rho$ strength twice.

The $\rho$ width comes mainly from its on-shell $\pi\pi$ channel. Dressing
the propagator schematically replaces the stable denominator by
$$
\frac{1}{M_{\rho,0}^2-s-\Delta_\rho^{\mathrm{self}}(s)},
\qquad
\operatorname{Im}\Delta_\rho^{\mathrm{self}}(s+i0)
=\sqrt{s}\,\Gamma_\rho(s).
$$
The imaginary part spreads the stable pole into a resonance-shaped
contribution on the two-pion cut. A constant-width Breit--Wigner can mimic
the line shape near the peak, while an energy-dependent dispersive
treatment also ties the real and imaginary parts together and can include
the measured two-pion phase.

# When each picture is useful

A pole model is a compact qualitative description when the dominant
resonance scale and residue are what matter, especially away from the peak
or when an observable averages broadly over virtualities. It is least
reliable near the $\rho$ peak, where the width and pion rescattering shape
the line, and wherever threshold behavior or precise phases matter. There
is no universal boundary: usefulness depends on the observable and target
accuracy. A quantitative comparison needs the same source amplitudes,
matching conventions, and a prescription that avoids counting resonance
strength once as a pole and again inside a continuum.

The current first-deliverable choice is therefore a dispersive two-pion
treatment for the broad $\rho$, with narrow-pole $\omega$ and $\phi$ terms.
The latter may still need finite-width treatment near their timelike peaks
or when phase precision matters. The repository has not yet quantified the
error of this mixed strategy; that and a full dispersive treatment of all
vectors are future extensions. See the
[mixed-vector strategy record](alp-tff-mixed-vector-strategy.md).
