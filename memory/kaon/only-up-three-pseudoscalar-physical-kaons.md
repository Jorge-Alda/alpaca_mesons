---
title: Physical neutral kaons in the only-up three-pseudoscalar amplitude
status: numerical
created: 2026-11-19
updated: 2026-11-19
bibtex_keys: [ParticleDataGroup:2024cfk, DiLuzio:2026poh]
provenance: numerical-computation
---

# Physical neutral kaons in the only-up amplitude

## Inputs and assumptions

The only-up source and radial exchange are described in [the three-pseudoscalar record](../alp/three-pseudoscalar/murayama-only-up-three-pseudoscalar-amplitude.md). The outgoing-state convention and its complex conjugation are in [the kaon projection](neutral-kaon-amplitude-projection.md). The input $|\epsilon_K|=2.228\times10^{-3}$ and phase $43.52^\circ$ follow [the PDG extract](../extracts/ParticleDataGroup:2024cfk/neutral-kaon-mixing.md) [@ParticleDataGroup:2024cfk]. Identifying the phenomenological decay parameter with the convention-dependent mixing parameter is an approximation; direct CP violation and CPT violation in mixing are neglected. The two neutral kaons use the common leading-order mass; their measured mass splitting and unstable propagation are omitted.

## Method and output

The real neutral-kaon generators used by the amplitude calculator are $T_1=E_{ds}+E_{sd}$ and $T_2=-i(E_{ds}-E_{sd})$. For the only-up diagonal source they are degenerate and unmixed. The outgoing $K^0$ and $\bar K^0$ generators are $(T_1-iT_2)/\sqrt2$ and $(T_1+iT_2)/\sqrt2$. Apply the conjugated coefficients $p^*$ and $\pm q^*$ from the projection note to each physical kaon leg and sum the complex flavor amplitudes before squaring. This also applies to the optional radial exchange. The returned amplitude is dimensionless; its two invariant arguments are in GeV$^2$.

## Verification

Run `.venv/bin/python -m unittest scripts.implementation.test_only_up_kaon_projection`. The tests check the complex coefficients on two outgoing legs and execute the finite-mixing amplitude calculator with physical kaon labels. The integrated width uses the same coherent amplitude, common leading-order neutral-kaon mass, and a factorial symmetry factor for identical physical labels. It checks each flavor component for accessible stable poles before integration; this conservative check can reject a pole even if its residue cancels in the full physical amplitude. Its output is in GeV. Tests also cover closed channels, both radial settings, and the identical-particle factor; the open-channel test replaces numerical integration with a fixed integral to isolate these conventions.

## Charge-conjugation zero

For the diagonal only-up source, the strong and radial interactions conserve charge conjugation. For a neutral spectator $\pi^0$, $\eta$, or $\eta'$ and common neutral-kaon masses, the two flavor amplitudes with $K^0\bar K^0$ and $\bar K^0K^0$ at the same momenta are equal, while same-strangeness amplitudes vanish. The projection in [the earlier diagnostic](pi0-kl-ks-vector-width-slowdown-diagnostic.md) therefore cancels exactly for $K_SK_L$, for any complex $\epsilon_K$. This applies to both the contact and radial-exchange terms and at finite neutral ALP mixing. It assumes the only-up source; it does not hold for a generic charge-conjugation violating source or unequal kaon propagation.

The implementation returns an exact zero before pole checks and quadrature for these channels. This removes roundoff rather than setting a numerical amplitude threshold. Other physical-kaon channels retain coherent projection and integration. The focused test checks all three spectator labels, both radial settings, and that quadrature is not called; it also retains a nonzero identical-kaon width check with mocked quadrature.
