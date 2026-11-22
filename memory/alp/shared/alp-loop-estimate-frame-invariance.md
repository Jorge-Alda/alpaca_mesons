---
title: Constant axial redefinition of the finite-mixing loop-size estimate
status: derived
created: 2026-11-22
updated: 2026-11-22
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Frame invariance of the loop-size estimate

Consider $g_R=g_L^\dagger=\exp(ix\kappa)$, with $x=a/f_a$ and a constant Hermitian axial generator $\kappa$, as in the [invariant meson-pole calculation](alp-invariant-loop-crossing-meson-poles.md). The generator preserves the source restrictions of the dressed-frame theory. The comparison is at the same on-shell kinematics for a nondegenerate pole and a nonzero tree amplitude.

1. **Dress the mass source.** The [source transformation laws](alp-chiral-source-transformations.md) give $M'=g_R M g_L^\dagger$ and $W_\eta'=g_\eta W_\eta$ for the right and left transports. Hence
   $$
   (W_R')^\dagger M'W_L'=W_R^\dagger M W_L=\widehat M(x).
   $$
   Taylor differentiation and the Hermitian split give the same $\widehat S_n$ and $\widehat P_n$ in both frames at every order.
2. **Dress the topological source.** The anomaly changes $\Theta(x)$ by $-x\operatorname{tr}\kappa$. The trace of the axial connection changes by $\operatorname{tr}\kappa$, since conjugation preserves the trace and the right and left derivative terms have opposite signs. In the [dressed topological coefficient](alp-dressed-topological-invariant.md) these two terms cancel:
   $$
   \Theta'(x)+\int_0^x\operatorname{tr}C_A'(x')\,dx'
   =
   \Theta(x)+\int_0^x\operatorname{tr}C_A(x')\,dx'.
   $$
   All $\widehat\Theta_n$ are therefore unchanged.
3. **Construct the pole and vertex inputs.** With invariant $\widehat S_0$, $\widehat P_1$ and $\widehat\Theta_1$, the [dressed pole equations](alp-lo-exact-pole-equations.md) use the same $\mathcal C=B_0\widehat P_1-(m_0^2/N_f)\widehat\Theta_1\mathbf 1$, the same meson mass operator and the same squared pole masses. In particular, writing $m_a$ for the tree-level ALP-like pole mass, $\mathbb{M}_{m_a^2}^{-1}(\mathcal C)$ and $N(m_a^2)$ are unchanged where the inverse exists. Use the same pole labels and external-state conventions in both frames.
4. **Compare loop size to tree.** In [Results 1–4 of the dressed vertex counting](alp-dressed-vertex-power-counting-invariants.md), forming vertices, multiplying pole overlaps, and taking flavour traces preserve the equality of the dressed inputs. Each class therefore has the same invariant contraction, factor from the ALP propagators and loop scale in both frames. This includes the mixing insertions on meson-like lines, which use projections of $\mathcal C$. The one-point source is also invariant, with the same vacuum prescription understood in both frames if its effect is included. The [naturalness estimate](alp-loop-naturalness-assumption.md) and the [exact tree amplitude](alp-finite-mixing-three-meson-tree.md) use those same inputs, so their ratio is unchanged. A crossing defined by that ratio reaching the same criterion consequently has the same location in both frames whenever such a crossing exists.

## Symbolic verification

The [executable check](../../../scripts/derive_alp_loop_frame_invariance.py) expands the mass transformation and its inverse through sixth order in $x$ using noncommuting source coefficients. The difference vanishes at each order. The script also checks explicitly that the first topological coefficient is recovered. The transport and trace identities above give the corresponding all-order statements.

For a diagonal projection of the commuting class-2 quartic contraction, the script constructs the dressed first pseudoscalar and second scalar coefficients independently in each frame. It then forms the contraction and its quotient after including the tadpole scale and propagator factor. It leaves the exact tree amplitude as a function of its invariant inputs. The resulting contraction and quotient differences vanish exactly. The general proof for the other contractions follows from the common dressed inputs in steps 3 and 4.

For the same diagonal projection, the untransported second scalar coefficient changes relative to $\widehat S_2$ by $-4\kappa\widehat P_1-4\kappa^2m_q$, with the quantities interpreted as their commuting entries. Substituting that raw coefficient in the transformed-frame contraction gives a generically nonzero quotient difference. The script also reuses the [raw tadpole test](../../../scripts/derive_alp_invariant_loop_crossing.py) to verify the distinction between the frame-dependent linear-spurion term and the invariant combination.

These checks leave all physical parameters symbolic. The representative quotient verifies the algebra of the estimate, with the process-dependent projection and loop functions left unevaluated.

```sh
.venv/bin/python scripts/derive_alp_loop_frame_invariance.py
PYTHONPATH=scripts .venv/bin/python -m unittest derive_alp_loop_frame_invariance
```
