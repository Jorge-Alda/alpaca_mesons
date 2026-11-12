---
title: Construction of leading-order amplitudes with finite ALP-meson mixing
status: derived
created: 2026-11-12
updated: 2026-11-12
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Construction of leading-order amplitudes with finite ALP-meson mixing

## Inputs and assumptions

Use the conventions of [@Alda:2026nrd] and the $U(3)$ CCWZ building blocks $u_\mu$ and $\chi_+$ in the [joint vertex derivation](alp-lo-joint-vertex-expansion.md). The undressed functions $M(x)$, $C_{R,L}(x)$ and $\Theta(x)$, with $x=a/f_a$, are independent inputs as described in the [general-source record](alp-general-chiral-redefinition-invariants.md). They are expanded locally about $x=0$ using their Taylor derivatives, with the factorials included; no linear dependence on $a$ is assumed. The currents contain one derivative of $a$. The bare ALP potential is assumed quadratic.

The expansion assumes an aligned stationary vacuum $U=1$, $a=0$, Hermitian $\chi(0)$ and $\Theta(0)=0$; tadpoles are absent or separately cancelled. Charge-conserving sources allow the neutral and charged sectors to be treated as in the [charged decoupling record](alp-lo-charged-meson-decoupling.md). The general mass-source second derivative corrects the ALP mass entry, as derived in the joint vertex record; the corrected matrices are inputs to the [spectrum solver](alp-lo-pole-solver.md). Pole vectors are normalized with $K$ and no expansion in $f_\pi/f_a$ is made. The explicit overlap in the [operator-inverse record](alp-lo-finite-mixing-operator-inverse.md) has the qualifications stated in the [pole-projected derivation](alp-lo-pole-projected-tree-amplitudes.md).

## Two constructions

1. The [joint expansion](alp-lo-joint-vertex-expansion.md) gives symmetrized, amputated vertices with incoming momenta. The [pole-projected construction](alp-lo-pole-projected-tree-amplitudes.md) contracts every external leg with its physical pole vector, including each neutral daughter's ALP component. Four-leg amplitudes include the contact vertex and all three cubic-vertex exchange partitions, using the full inverse propagator at the exchanged momentum. The pole vectors already contain the residue normalization.
2. The [pole-aligned construction](alp-lo-pole-aligned-tree-vertices.md) applies the exact symmetric axial chiral map, locally expanded to the field order required by the vertex. It transforms daughter directions and the quadratic matrices as well as the vertices. For a stationary vacuum its cubic coordinate terms suffice through four external legs. They produce a nonlinear quartic contact whose kinetic and mass pieces need not vanish separately, although their sum vanishes for on-shell external poles. The shift aligns one parent pole and does not diagonalize propagation globally. The map and its relation to the CCWZ blocks are described in the [finite-redefinition record](alp-ccwz-finite-redefinitions.md).

The two routes give the same on-shell tree amplitude in the checks documented in the pole-aligned record. The source-consistent mixed-pole comparison uses noncommuting $U(2)$ sources at two algebraic on-shell points, with nonzero exchange; those points are not asserted to be in real physical phase space. A separate $U(3)$ charged-vertex check retains both real and imaginary parts of the Hermitian $d$–$s$ vector source, but uses an illustrative shift rather than a same-source $U(3)$ spectrum. These checks are not an exhaustive verification for all sources or channels.

## Validity

The calculation is LO and tree level, with real symmetric constant $K>0$ and $M^2$, simple stable external poles and an aligned stationary vacuum. The pole-aligned map additionally requires a nonzero ALP component of the selected parent, an invertible meson block at its pole, and a local logarithm near $U=1$; the direct projection can instead use a zero-ALP pole. Degenerate poles require a subspace treatment. A genuine internal pole requires a prescription beyond this stable tree-level calculation. The exact diagonalization is exact within the LO model, not a higher-order chiral result. Large admixtures, near degeneracy and other conditions listed in the [truncation record](alp-lo-truncation-limits.md) are heuristic warning criteria, not derived bounds. Loops, NLO terms, resonance mixing, widths and complex poles are excluded.

## Reproducibility

From the repository root, the minimal examples for the [vertex generator](../../../scripts/derive_alp_lo_joint_vertices.py), [pole-projected evaluator](../../../scripts/derive_alp_pole_projected_amplitudes.py), and [pole-aligned evaluator](../../../scripts/derive_alp_pole_aligned_vertices.py) are:

```sh
.venv/bin/python scripts/derive_alp_lo_joint_vertices.py
.venv/bin/python scripts/derive_alp_pole_projected_amplitudes.py
.venv/bin/python scripts/derive_alp_pole_aligned_vertices.py
```

All three invocations completed successfully on 2026-11-12. Their outputs are illustrative symbolic examples, not phenomenological predictions. The individual derivation records link the associated test suites and describe precisely which comparisons they perform.
