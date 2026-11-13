# Tasks

Scripts are written by the `coder` agent. Each memory record has the front matter required by AGENTS.md, links its inputs, and goes into its `INDEX.md` and every parent index. Track with `bd` (alp-sjp.6, mass-shift part only). Do not close alp-sjp.6. Do not commit unless asked. Results are analytical; exact rational checks against the solver verify them and are not results. No expansion in $f_\pi/f_a$ is part of this change.

## 1. Exact relations

- [x] 1.1 Script and record: charged-meson zero shift (link the decoupling record, add the $m_q$-dependent unmixed masses), exact $K^0$ pole, trace sum rule with explicit $K^{-1}$ and in the sources, product rule $\det M^2$. Verify symbolically in the 6x6 neutral block and with exact rationals against `solve_alp_meson_spectrum.py`. Record in `memory/alp/shared/`.

## 2. Finite-mixing poles

- [x] 2.1 Illustration, not a result of this change: the $O(\epsilon^2)$ shift of a generic meson pole, shown as the leading behaviour of the exact equation of 2.2. Record and script exist and are kept.
- [x] 2.2 Script and record: exact pole equations of the neutral mesons and the ALP (Schur form, canonical frame in invariants, poles as eigenvalues of the canonical mass matrix), with the invariance under constant axial redefinitions.
- [x] 2.3 Script and record: exact interlacing bounds and the arrow-matrix spectrum.
- [x] 2.4 Script and record: isolated meson pole in closed form with its exact error bound, and isolated ALP pole in closed form with its fixed-point iteration (both labelled as approximations under the far-apart assumption).
- [x] 2.5 Script and record: explicit exact pole equations for $\pi^0$, $\eta$, $\eta'$ at diagonal $\chi_0$ (the secular polynomial and the per-meson residues), general $m_q$ kept as in 2.2. No isospin-limit rotation.
- [x] 2.6 Script and record: $K^0$ doublet, exact pole and its CP-conserving and CP-violating reductions; consistency with the exact pole of 1.1.

## 3. Checks

- [x] 3.1 Script: the exact trace and product sum rules of 1.1 checked against the arrow-matrix eigenvalues and against `solve_alp_meson_spectrum.py` with exact rationals.

## 4. Limits and closing

- [x] 4.1 Record: validity domain of the isolated approximations (the spacing condition), level-repulsion sign rule, exact two-level formula for near degeneracy (link the existing record), the $f_a\to\infty$ limit returning the unmixed masses, and the regime statement: small mixing keeps no shift at order $f_\pi/f_a$, large mixing uses the exact poles.
- [x] 4.2 Update every affected `INDEX.md`; run all new scripts and tests together; check links and front matter.
- [x] 4.3 Add a note to alp-sjp.6 listing what is done (LO analytical pole equations and their approximations) and what remains (decay constants, residues, mixing observables, other multiplets, NLO/loops, uncertainties); keep it open. Report changed files and a proposed commit message.
