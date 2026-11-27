---
title: Computation of the cubic scalar moments and crossing for three mesons
status: mixed
created: 2026-11-23
updated: 2026-11-23
bibtex_keys: [Alda:2026nrd, Manohar:1983md]
provenance: analytical-derivation
---

# Computation of the cubic scalar moments and crossing

The [on-shell loop estimate](alp-three-meson-one-loop-estimate.md) defines the diagram contributions. The [on-shell assembly](alp-three-meson-onshell-assembly.md) gives the vacuum, fixed-mass, exchange, and external-leg corrections. Supply the dressed sources, full tree pole basis and masses, the three daughter labels, two independent pair invariants, and the exact tree amplitude from the [tree pole projection](../../../decay-amplitudes/vertices-and-kinematics/tree-vertices/finite-mixing/alp-finite-mixing-three-meson-tree.md). Use pole projections analytic in $\epsilon$ on the chosen branch, with the ALP mass held fixed. The loop reference scales are inputs to the [naturalness prescription](../power-counting-and-validity/loop-counting/alp-loop-naturalness-assumption.md).

1. Expand each routed projected vertex product, pole mass, and propagator argument through cubic order. Expand the denominators before momentum reduction. In particular, the mass shift of a meson line in a term linear in $\epsilon$ gives
   $$
   \frac{\epsilon}{k^2-m_n^2}
   =
   \frac{\epsilon}{k^2-m_{n,0}^2}
   +\frac{\epsilon^3\langle\lambda_n\mathcal C\rangle^2}
    {(m_{n,0}^2-m_a^2)(k^2-m_{n,0}^2)^2}
   +O(\epsilon^5).
   $$
   The [cycle expansion script](../../../../../scripts/derive_alp_three_meson_loop_cubic_cycles.py) constructs each order needed for one- through four-point functions. It retains the powers of raised denominators and the [momentum reduction](alp-loop-momentum-reduction.md).
2. Use the [one-loop callback](../../../../../scripts/derive_alp_loop_1pi_callback.py) to sum routed cycles with their Hessian symmetry factors. Continue the external Gram products to Euclidean signature for the cycle calculation and convert the 1PI coefficients back to the Minkowski vertex conventions. Keep the scalar moments as functions of the shifted denominator, so differentiating the two-point function for the residue retains the derivative of the scalar moment.
3. The [connected driver](../../../../../scripts/derive_alp_three_meson_onshell_driver.py) passes those 1PI coefficients to the on-shell assembly and expands the resulting amplitude through order $\epsilon^3$. It groups the cubic scalar-moment contributions before applying the chosen magnitude weights. In the notation of the manuscript, the crossing satisfies
   $$
   \frac{\epsilon^3}{16\pi^2}
   \sum_{r,\alpha}(m_\alpha^2)^{2-r}\big|C_{r\alpha}^{(3)}\big|
   = R_{\max}\big|\mathcal M^{(0)}(\epsilon)\big|.
   $$
   Here $C_{r\alpha}^{(3)}$ is calculated by the cubic expansion and subsequent scalar-moment grouping, $m_\alpha$ is the reference scale of that moment, and $\mathcal M^{(0)}$ retains the exact tree pole projections and exchange denominators. The [crossing record](../power-counting-and-validity/crossing-and-validity/alp-loop-crossing-exact-tree-condition.md) states the domain and threshold convention.

The executable checks cover geometric propagator expansion, moment differentiation, on-shell assembly and a nonzero end-to-end symbolic tadpole projection through the crossing equation. The end-to-end example uses labeled meson directions sharing a flavour matrix and checks the algebra. A full channel with physical quark masses and sources has yet to be evaluated.

## Invariant-input solver

The [crossing CLI](../../../../../scripts/calculate_alp_loop_crossing.py) takes a JSON input file. Its help lists the field names and units. Supply the dressed Taylor derivatives $\widehat S_n$, $\widehat P_n$ and $\widehat\Theta_n$ for $n=0,\ldots,6$, without factorials. Each mass-source derivative is a real $3\times3$ matrix in GeV. Each topological coefficient is real and dimensionless. Exact rational values can be entered as strings.

The supplied reference scale applies to non-tadpole scalar moments. Standalone tadpoles use their internal masses. Daughter labels are zero-based indices in the sorted unmixed meson spectrum. Supply a positive $R_{\max}$ and a bracket in $f_\pi/f_a$ that spans a crossing on a continuous physical pole branch. The output reports the root, the cubic loop-size coefficient and the equation residual. The [pole adapter](alp-loop-crossing-pole-adapter.md) defines the exact tree evaluation. The full symbolic cycle expansion can be expensive.

```sh
PYTHONPATH=scripts .venv/bin/python -m unittest derive_alp_three_meson_loop_cubic_cycles derive_alp_loop_1pi_callback derive_alp_three_meson_onshell_assembly derive_alp_three_meson_onshell_driver
```
