# Tasks

Scripts are written by the `coder` agent; each memory record follows the front-matter rules of AGENTS.md and is added to `memory/alp/shared/INDEX.md` (and parent indexes) with its script. Track with `bd` (alp-sjp.1).

## 1. Hadronic transformation of sources and anomaly

- [x] 1.1 Derive the transformed $r_\mu,\ell_\mu,\chi$ and the $\widehat U$ representation for general $g_{R,L}(x)$ with a script (extending `derive_alp_invariant_coefficient_functions.py` or a new one); verify that $\partial_\mu\widehat U=W_R^\dagger(D_\mu U)W_L$ and the composition law of two redefinitions hold symbolically through order 3, with noncommuting $U(2)$ matrices.
- [x] 1.2 Derive the anomalous shift of $\Theta(x)$ and of the meson-level topological source for a general singlet $\theta(x)$; verify the Jacobian shift and $X_{TC}$ invariance in a script, including the $4\pi/\alpha_s$ conversion to Alda:2026nrd's normalization.
- [x] 1.3 Write the memory record for 1.1–1.2 (derived, linked to the partonic record and `alp-qcd-topological-cp.md`); verify the invocation example in the record runs as written.

## 2. Nonlinear coordinate map and transformed vertices

- [x] 2.1 Derive the all-orders meson map $U'=g_RUg_L^\dagger$, the compensator $h$ and the transformed $u_\mu,\Gamma_\mu,\chi_\pm$ with a script; verify invariance of $\langle u_\mu u^\mu\rangle$, $\langle\chi_+\rangle$ symbolically.
- [ ] 2.2 Expand the transformed Lagrangian to the meson order needed for $aPP$ and $aPPP$ vertices; verify with a script that the original and transformed amputated vertices differ by exactly the pole-vector and contact terms, and that dropping the contact terms fails the check.
- [ ] 2.3 Write the memory record for 2.1–2.2 and cross-link `alp-ccwz-arbitrary-order-expansion.md`; verify the script reproduces the stated first-order limit $c_q^{VS},c_q^{AP}$.

## 3. Combined inverse propagator and pole alignment

- [ ] 3.1 Derive $D'(s)=T^{\mathsf T}D(s)T$ and the invariance of poles, null vectors, residues and the LSZ amplitude for generic $T$; verify with exact rational data for $N_f=3$ in a script.
- [ ] 3.2 Derive the single-pole aligned coordinates, their regularity condition and the exceptional-pole alternative; verify the identity $m_k^2K'_{i\hat a}-M^{2\prime}_{i\hat a}=0$ and agreement with `derive_alp_finite_mixing.py`.
- [ ] 3.3 Derive the multi-pole alignment conditions for shift-type maps and for general $T$; verify with a two-pole counterexample and the pencil diagonalization in a script.
- [ ] 3.4 Show with a script that the linear alignment shift of the finite-mixing record is the constant-generator case of $g_{R,L}(x)$, with the generator-to-$w_i^{(a)}$ relation stated.
- [ ] 3.5 Write the memory record for 3.1–3.4, stating the cases in which alignment fails; cross-link `finite-alp-meson-mixing.md`.

## 4. Integration

- [ ] 4.1 Run all new scripts together with the existing finite-mixing and invariant-coefficient scripts and verify they pass; update every affected `INDEX.md` and run the repository's index and link checks.
- [ ] 4.2 Check each acceptance criterion of alp-sjp.1 against the records and scripts, add a note to the Beads issue with the evidence, and state remaining limits for alp-sjp.2, .3, .5 and .6.

## Workflow follow-up

- Propose a commit message after the files are modified; commit only on request.
- Close alp-sjp.1 in Beads only after 4.2 is confirmed.
