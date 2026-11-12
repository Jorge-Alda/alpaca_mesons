# Design

## Context

The LO spectrum is solved: `scripts/solve_alp_meson_spectrum.py` returns $m_n^2$, $v^{(n)}$ ($v^{(n)\mathsf T}Kv^{(n)}=1$) and residues of $D^{-1}$ for the neutral basis $(P_3,P_8,P_0,P_6,P_7,a)$, with the exact $K^0$ pole carrying no ALP component. Covariance under constant field changes is in [the congruence record](../../../memory/alp/shared/alp-inverse-propagator-congruence.md), and finite redefinitions of the CCWZ blocks in [the redefinition record](../../../memory/alp/shared/alp-ccwz-finite-redefinitions.md). Existing channel amplitudes use the pole-aligned route for one decaying ALP and unmixed outgoing mesons. See proposal.md for scope.

## Goals / Non-Goals

**Goals:**
- One construction valid for all LO channels, with outgoing mesons also projected on their poles.
- Cross-check: pole-projection route versus pole-aligned transformed vertices.
- A process inventory usable as the process list of alp-sjp.5. The change is purely analytical: no numerical inputs are used and nothing depends on alp-sjp.4.

**Non-Goals:**
- Resonance exchange, loops, NLO terms, widths, numerical phenomenology, phase-space integration (see proposal.md).

## Decisions

1. **Primary route: contraction in the original fields.** Expand $\mathcal L_{\rm LO}(P_i,a)$ jointly to the vertex order, contract each leg with $v^{(n)}$, divide nothing further because $v^{(n)}$ already carries the LSZ normalization ($v^{\mathsf T}Kv=1$). This treats parent and daughters on the same footing. Alternative (shift to align only the parent, daughters unmixed) rejected: it misses the ALP admixture of outgoing mesons and breaks near degenerate poles.
2. **Cross-check route: transformed vertices.** Apply the exact chiral map of the pole-aligned redefinition ([single-pole record](../../../memory/alp/shared/alp-single-pole-alignment.md), [generator record](../../../memory/alp/shared/alp-pole-alignment-chiral-generator.md)) and take vertices in the new coordinates, including all nonlinear terms. Equality with route 1 is both the completeness check and the statement required by the acceptance criteria. Alternative of using only one route rejected: the issue asks explicitly for complete transformed vertices, and the two routes fail differently.
3. **Vertex generation by a joint polynomial expansion.** Expand $U=\exp(i\Phi/f_\pi)$ and the $a$-dependent sources to the needed total order (up to four legs, for $a\to3P$ with the parent ALP-component leg) with SymPy on the $6\times6$ neutral basis plus charged fields where the channel needs them. Reuse `derive_alp_ccwz_arbitrary_order.py`. Charged final states use the exact decoupling of the charged mesons from $a$ ([decoupling record](../../../memory/alp/shared/alp-lo-charged-meson-decoupling.md)), so their legs have no ALP component.
4. **ALP dependence of sources.** **Full ALP-dependent sources, no linearity assumption.** The ALP enters only through the functions $M(x)$, $C_{R,L}(x)$ (one $\partial_\mu a$) and $\Theta(x)$ of [the invariants record](../../../memory/alp/shared/alp-general-chiral-redefinition-invariants.md), with $x=a/f_a$, kept as arbitrary functions with Taylor coefficients $M^{(k)}(0)$, $C^{(k)}_{R,L}(0)$, $\Theta^{(k)}(0)$ carried symbolically. A vertex with $k$ ALP legs uses the $k$-th Taylor coefficient, so $a^kP^m$ contact terms follow from the exact functions; the dressed $\widehat M$ and transport matrices of [the source-transformation record](../../../memory/alp/shared/alp-chiral-source-transformations.md) give the transformed sources for route 2. The linear-spurion forms of [the finite-mixing record](../../../memory/alp/shared/finite-alp-meson-mixing.md) and Alda:2026nrd are recovered as the special case $M''=0$ etc. and are used only as a check. The quadratic blocks $K$, $M^2$ depend on these functions only through their values and first and second Taylor coefficients at $x=0$; the first task verifies that the archived solver blocks agree with this general form and, where they assumed linearity, supplies the general blocks through the solver's $(K,M^2)$ interface.
5. **Exchange diagrams with the full $D^{-1}(s)$.** Required for correct off-pole behaviour of the mixed system; checked against the pole sum on shell.
6. **Photon channels.** LO is the WZW sector plus the anomaly coefficients in the Bauer normalization of Alda:2026nrd ([photon invariant record](../../../memory/alp/photon/alp-photon-invariant-from-mixing.md), [WZW $\gamma PP$ record](../../../memory/alp/photon-two-pseudoscalar/wzw-alp-photon-pp-amplitude.md)). The finite-mixing $\gamma\gamma$ amplitude sums over all meson components of the parent pole; the vector-conversion and resonance parts are deferred.
7. **Baryon channels.** Use the LO CCWZ vertex of [the baryon record](../../../memory/alp/baryon/finite-mixing-baryon-couplings.md). The decay $a\to B\bar B$ is kinematically open only above $2m_B$ and a LO treatment for such masses is outside the chiral domain; the inventory states this and gives the amplitude structure, with the open/closed test in terms of the physical parent mass.
8. **Selection rules from the vectors.** Derive them by tracing flavour structure and symmetries of the sources, not by assuming first-order mixing; verify numerically on random CP-conserving and CP-violating sources with the solver.
9. **Generic pseudoscalar amplitudes.** Derive $a\to P_iP_j$ and $a\to P_iP_jP_k$ with the generators $\lambda_i$ and indices left open, following the generic-generator approach of [the three-meson generator record](../../../memory/alp/three-pseudoscalar/alp-three-meson-generic-generators.md); specific channels follow by substitution. Alternative of one derivation per channel rejected: it multiplies work, hides common structure and invites inconsistencies. $a\to4P$ is excluded from this change.
10. **Records.** One concept per file: construction (shared), inventory (shared), selection rules (shared), and one record per channel family in its existing folder. Notation follows Alda:2026nrd; no new shorthand without approval.

## Risks / Trade-offs

- [Pole vector $v^{(n)}$ has $O(1)$ admixture, so the LO model is not controlled] → carry the spectrum diagnostics into each amplitude record and flag it.
- [Higher Taylor coefficients of $M(x)$, $C_\eta(x)$, $\Theta(x)$ are free parameters, so results carry many symbols] → keep them symbolic, report the linear-spurion special case as a check, and state which coefficients each channel depends on.
- [Archived solver blocks may assume linear sources] → task 1.1 compares them with the blocks from the general functions and extends them if they differ; any change to the spectrum is reported against `lo-finite-mixing-spectrum`.
- [Off-shell dependence of vertices with derivative couplings] → define amputated vertices as the Lagrangian expansion in the stated field coordinates and verify through the route-1/route-2 equality and field-redefinition invariance of the on-shell amplitude.
- [Generic-index amplitudes involve flavour traces of several generators] → derive them symbolically in terms of $\operatorname{Tr}$ of products of $\lambda_i$ and the pole-vector contractions, and evaluate specific assignments only as checks.
- [Overlap with existing channel records] → reproduce them as checks; only update an existing record if a discrepancy is found, and report it.
