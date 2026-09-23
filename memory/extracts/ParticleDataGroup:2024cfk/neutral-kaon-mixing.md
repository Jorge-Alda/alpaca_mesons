---
title: Neutral-kaon mixing and indirect CP violation
status: literature
created: 2026-09-23
updated: 2026-09-23
bibtex_keys: [ParticleDataGroup:2024cfk]
provenance: literature
---

# Neutral-Kaon Mixing and Indirect CP Violation

The Particle Data Group writes the neutral-meson mass eigenstates in terms of flavour states with complex coefficients $p$ and $q$, normalized by $|p|^2+|q|^2=1$ when CPT violation in mixing is neglected. For kaons, the states are conventionally labeled by lifetime as $K_S$ and $K_L$ rather than by light and heavy mass [Section 13.1.2, Eq. (13.7), @ParticleDataGroup:2024cfk].

The implementation adopts the commonly used kaon convention

$$
|K_S\rangle=p|K^0\rangle+q|\bar K^0\rangle,
\qquad
|K_L\rangle=p|K^0\rangle-q|\bar K^0\rangle,
$$

with

$$
p=\frac{1+\epsilon_K}{\sqrt{2(1+|\epsilon_K|^2)}},
\qquad
q=\frac{1-\epsilon_K}{\sqrt{2(1+|\epsilon_K|^2)}}.
$$

The numerical default uses

$$
|\epsilon_K|=(2.228\pm0.011)\times10^{-3},
\qquad
\phi_{\epsilon_K}=43.52^\circ,
$$

from the PDG kaon discussion [Section 13.4, Eqs. (13.42) and (13.46a), @ParticleDataGroup:2024cfk]. The phase assignment uses the CPT-invariant two-pion phase in Eq. (13.42) as the standard phenomenological phase of $\epsilon_K$. This identification is an explicit approximation: the PDG notes after Eq. (13.45) that the convention-dependent mixing parameter $\widetilde\epsilon=(1-q/p)/(1+q/p)$ is in general different from the decay parameter $\epsilon$. The approximation neglects that distinction, direct CP violation, and CPT violation in mixing; callers may override the complex `epsilon_k` value.

The mass eigenstates are not orthogonal when $|p|\ne|q|$:

$$
\langle K_S|K_L\rangle=|p|^2-|q|^2
=\frac{2\operatorname{Re}\epsilon_K}{1+|\epsilon_K|^2}.
$$

This last equation is a derived consequence of the stated convention, not a separately quoted PDG result.
