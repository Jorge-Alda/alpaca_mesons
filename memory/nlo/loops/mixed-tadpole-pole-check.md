---
title: Mixed ALP–meson tadpoles and local pole check
status: derived
created: 2026-11-10
updated: 2026-11-10
bibtex_keys: [Alda:2026nrd, Kaiser:1998ds]
provenance: analytical-derivation
---

# Mixed tadpoles and local poles

## Inputs and assumptions

The calculation uses the $a\Phi^3$ vertices in [vertices and topologies](vertices-and-topologies.md), the tadpole contraction and subtraction in [loop functions and subtraction](loop-functions-and-subtraction.md), and the local terms in [quadratic Lagrangian](../quadratic-lagrangian.md). The loop classes are those in [loop classes](../counting/loop-classes.md). The subtraction coefficients are those in [counterterm basis](../counting/counterterm-basis.md), with the source-matched $B$ subtraction in [subtraction scheme and scales](../counting/subtraction-scheme-and-scales.md). The flavour source matrices do not commute in general. The calculation is at first order in $1/f_a$, with no CP-violating sources or constant $\theta$ term, and $Q=sK-M^2+\Sigma$. Masses and decay constants are in GeV and $s$ in GeV$^2$.

## Method and result

For each external generator, the script contracts the two remaining meson fields in each of the three possible placements in the $a\Phi^3$ vertices. It prints the resulting traces for $c_q^A$ and $c_q^P$, keeping the source matrices arbitrary. The $c_q^V$, $c_q^S$ and $c_G$ tadpoles vanish at this order, as explained in the vertex record. The finite tadpole matrix is the matrix function of the leading mass matrix described in the subtraction record. The pole is obtained by replacing the tadpole matrix by twice the leading mass matrix times $\lambda$.

The following is a diagnostic, not a completed renormalized mixed function. With the source-matched subtraction and only the stated local terms, the mixed pole for external $\pi_3$ and $c_q^A=\lambda_3$ is

$$
-\lambda(\mu)\,s\,\frac{2B_0(5m_u+5m_d+2m_s)}{3f_af_\pi}.
$$

It does not vanish for generic quark masses. The off-shell meson sector likewise required a divergent field redefinition, as documented in [meson self-energy](meson-self-energy.md). Its mixed-block counterpart has not yet been derived; cancellation cannot be asserted from the local terms alone.

## Verification

The unit test `test_mixed_contraction_against_vertex_script` compares the $c_q^A$ and $c_q^P$ contractions with the independently expanded vertex script for a selected noncommuting component. `test_topological_normalization` checks the $B$ normalization. Run

```bash
PYTHONPATH=scripts .venv/bin/python -m unittest derive_nlo_loop_self_energies.LoopSelfEnergyTests.test_mixed_contraction_against_vertex_script derive_nlo_loop_self_energies.LoopSelfEnergyTests.test_topological_normalization
.venv/bin/python scripts/derive_nlo_loop_self_energies.py --mixed
```

The first command passes. The second deliberately exits with a failed assertion because the local terms alone do not cancel every mixed pole. The mixed calculation is incomplete until the field redefinition and any resulting source terms are derived and tested.
