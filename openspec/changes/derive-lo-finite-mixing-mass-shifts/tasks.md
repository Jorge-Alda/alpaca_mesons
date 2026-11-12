# Tasks

Scripts are written by the `coder` agent. Each memory record has the front matter required by AGENTS.md, links its inputs, and goes into its `INDEX.md` and every parent index. Track with `bd` (alp-sjp.6, mass-shift part only). Do not close alp-sjp.6. Do not commit unless asked. Results are analytical; exact rational checks against the solver verify them and are not results.

## 1. Exact relations

- [x] 1.1 Script and record: charged-meson zero shift (link the decoupling record, add the $m_q$-dependent unmixed masses), exact $K^0$ pole, trace sum rule with explicit $K^{-1}$ and in the sources, product rule $\det M^2$. Verify symbolically in the 6x6 neutral block and with exact rationals against `solve_alp_meson_spectrum.py`. Record in `memory/alp/shared/`.

## 2. Meson shifts

- [ ] 2.1 Script and record: $O(\epsilon^2)$ shift of a generic meson pole $j$ from the Schur complement, written with $D_{ja}(m_{j,0}^2)$, then in the invariants $c_q^{AP}$, $c_{GA}$, $\chi_0$, $m_0^2$, $m_{a,0}^2$. Proof of invariance under a symbolic constant chiral redefinition.
- [ ] 2.2 Script and record: explicit isospin-limit forms for $\pi^0$, $\eta$, $\eta'$ (closed $\eta$-$\eta'$ rotation) and general-$m_q$ form through $O$; check the $m_u\to m_d$ limit.
- [ ] 2.3 Script and record: $K^0$ doublet shift via $|D_{6a}|^2+|D_{7a}|^2$, CP-conserving and aligned-phase reductions, CP-violating case; consistency with the exact pole of 1.1.
- [ ] 2.4 Script and record: $O(\epsilon^4)$ meson shifts by one more iteration; verify the coefficient against the weak-mixing expansion and by residual exponents on exact rational sources.

## 3. ALP shift

- [ ] 3.1 Script and record: $\delta m_a^2$ at $O(\epsilon^2)$ and $O(\epsilon^4)$, split into kinetic, anomaly and exchange parts, plus the additive $\operatorname{tr}S_2$ term; no-meson and $c_G=0$ limits.
- [ ] 3.2 Script and record: pole-aligned frame: exact implicit equation $D'_{aa}(m_a^2)=0$, equality with the unaligned expansion, and the same for a meson pole; check that the sum of all shifts equals the trace sum rule of 1.1 at $O(\epsilon^2)$.

## 4. Limits and closing

- [ ] 4.1 Record: validity domain, level-repulsion sign rule, exact two-level formula for near degeneracy (link the existing record), $\epsilon\to0$ and weak-mixing checks against Alda:2026nrd.
- [ ] 4.2 Update every affected `INDEX.md`; run all new scripts and tests together; check links and front matter.
- [ ] 4.3 Add a note to alp-sjp.6 listing what is done (LO analytical mass shifts) and what remains (decay constants, residues, mixing observables, other multiplets, NLO/loops, uncertainties); keep it open. Report changed files and a proposed commit message.
