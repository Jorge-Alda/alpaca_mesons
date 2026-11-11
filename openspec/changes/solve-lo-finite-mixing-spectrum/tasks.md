# Tasks

Scripts are written by the `coder` agent; each memory record has the front matter required by AGENTS.md and is added to `memory/alp/shared/INDEX.md` (and parent indexes). Track with `bd` (alp-sjp.2, LO part only). Do not close alp-sjp.2.

## 1. LO solver

- [x] 1.1 Write a script that builds $K$ and $M^2$ for $(P_i,a)$, $N_f=3$, from the closed forms of `finite-alp-meson-mixing.md`, and verifies them against the expansion of the LO Lagrangian in `derive_alp_lo_finite_mixing.py`.
- [x] 1.2 Implement the general-$(K,M^2)$ solver (Cholesky reduction): poles, $K$-orthonormal eigenvectors, compositions, residues of $D^{-1}$; verify completeness $\sum_nv^{(n)}v^{(n)\mathsf T}=K^{-1}$ and $D^{-1}(s)$ against a direct inverse, with exact rational data and with random symmetric $(K,M^2)$.
- [x] 1.3 Verify agreement with the single-pole eigenvector $w^{(a)}_i(m_a^2)$ and the LSZ factor of `finite-alp-meson-mixing.md` and with `derive_alp_finite_mixing.py`.
- [x] 1.4 Verify invariance under a random rational invertible $T$ (poles, eigenvector and residue maps), reusing `derive_alp_inverse_propagator_congruence.py`.

## 2. Limits and diagnostics

- [ ] 2.1 Verify the weak-mixing limit symbolically against the first-order results of Alda:2026nrd; report the scaling of the residual with $f_\pi/f_a$.
- [ ] 2.2 Implement stability diagnostics ($K>0$, $M^2\ge0$, real nonnegative poles) and show analytically ($\det K=1$) and numerically that $K>0$ at LO; test the unstable-vacuum branch.
- [ ] 2.3 Implement near-degeneracy detection, the two-level formula check against `finite-alp-meson-mixing.md`, and the degenerate-subspace output.
- [ ] 2.4 Implement pole tracking along a scan in $m_{a,0}^2$ and in $f_\pi/f_a$; check avoided crossings and label stability.

## 3. Records and integration

- [ ] 3.1 Write the memory record(s) (derived/numerical, with inputs linked and the scripts' minimal invocations verified), including the explicit statement of the LO truncation limits and of what is deferred (NLO, loops, resonances, complex poles).
- [ ] 3.2 Update every affected `INDEX.md`; run all new and related scripts together; check links.
- [ ] 3.3 Check each alp-sjp.2 acceptance criterion against the evidence for the LO part; add a note to the Beads issue; create a follow-up bead for the NLO, loop and resonance extension and link it to alp-sjp.2.
