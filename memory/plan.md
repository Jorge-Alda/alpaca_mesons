# Plan for the research project

## 1. Contact interactions from the lowest order chiral Lagrangian

- $a\to 3P$ from flavour-conserving interactions: DONE
- $a\to 2P$ from flavour-violating interactions: DONE

## 1b. Field-redefinition invariants of the ALP sources

- First-order dressed invariants $\widehat S_1$, $\widehat P_1$, $\widehat\Theta_1$, $\widehat\Theta_{\gamma,1}$: DONE
- Second-order dressed invariants $\widehat S_2$, $\widehat P_2$, $\widehat\Theta_2$, $\widehat\Theta_{\gamma,2}$ (derivative Taylor convention, compact commutator forms): memory DONE, latex draft MISSING
- Next step: use $\widehat S_2$, $\widehat P_2$ in the LO meson vertices at order $a^2/f_a^2$.

## 1c. Radial mode in the model of DiLuzio:2026poh

- $\sigma\to P_iP_j$, polar and Cartesian, general sources: memory DONE
- $a\to P_iP_jP_k$ with $\sigma$ exchange, general sources, polar coordinates ([record](alp/three-pseudoscalar/alp-three-pseudoscalar-radial-exchange.md)): memory DONE, latex draft DONE (`sec:murayama-three-general` and following subsections of `latex/murayama_matching.tex`)
- Quartic vertex in the basis without ALP–meson mixing ([record](alp/three-pseudoscalar/murayama-unmixed-basis-quartic-vertex.md)), parity rule removing cubic vertices ([record](alp/three-pseudoscalar/murayama-cubic-vertex-parity-rule.md)), only-up assembly with $\sigma$ exchange ([record](alp/three-pseudoscalar/murayama-only-up-three-pseudoscalar-amplitude.md)), contact cancellation for three distinct neutral daughters ([record](alp/three-pseudoscalar/murayama-distinct-neutral-contact-cancellation.md)): memory DONE, latex DONE, draft equations tested against the direct calculation.
- $\sigma$ exchange in `scripts/implementation/only_up_decays.py` (explicit `radial_mass_squared` argument): DONE.
- Vertex and $\sigma$-coupling formulas (draft and the `murayama-*` records) rewritten in invariant form: traces of the original $\widehat S_{rs}$, $\widehat P_{rs}$ with $\mathbb M^{-1}_{m^2}(\widehat P_{01})$, $N(m^2)$ and pole masses; no $w$ vectors, shifted source coefficients or auxiliary node symbols. DONE, tested equation by equation.
- Pending decision (user): $w_k^{(r)}$ still appears in the draft's pole-direction, chiral-shift and eq:murayama-T definitions, and $\Omega_r=\sum_kw_k^{(r)}\lambda_k$ in [radial-decay-polar-cartesian-general-sources.md](alp/shared/radial-decay-polar-cartesian-general-sources.md) and its script.
- Next step: numerical widths of the only-up model with $\sigma$ exchange as a function of $m_\sigma$ (plots in `plots/murayama_3P.pdf` still omit it); treat the $\sigma$ width when $m_\sigma<m_a-m_r$ (currently raises `ValueError`).

## 1d. Validity of the LO theory at finite $f_\pi/f_a$

- Removal of ALP derivative couplings and on-shell frame equality ([record](alp/shared/alp-dressed-frame-derivative-coupling-removal.md)): memory DONE, latex MISSING
- Tree-level S-wave unitarity ([record](alp/shared/alp-lo-partial-wave-unitarity.md)): memory DONE, latex MISSING
- One-loop power counting with ALP lines ([record](alp/shared/alp-loop-power-counting.md)) and beyond-LO singlet terms ([record](alp/shared/alp-singlet-derivative-terms-energy-growth.md)): memory DONE, latex MISSING
- Loop-to-tree crossing at finite mixing (OpenSpec change `estimate-loop-tree-crossing-finite-mixing`): task 1.1, comparison of four candidate observables against tests (a)-(d) with exact leading powers of $\epsilon$ and frame behaviour ([record](alp/shared/alp-loop-observable-candidates-comparison.md)): memory DONE. Open: the pure-ALP tadpole of the amplitude is exactly a shift of $\widehat P_1$, which weakens design D3; the non-absorbability of the mixing-induced pieces is a conjecture.
- Task 1.2, selection of the observable ([record](alp/shared/alp-loop-observable-selection-a-to-3pi0.md)): memory DONE. Selected: coefficients of $\bar B(s_{ij};m_A,m_P)$ ($P=\pi^0,\eta,\eta'$) and $\bar B(s_{ij};m_A,m_A)$ in $A\to3\pi^0$. They are finite and $\mu$-independent, and no operator with at most four derivatives reproduces them. All $s_{ij}$-independent parts, including the $A_0(m_A^2)$ logarithm, are absorbed. The loop does not follow the tree: $(A,\eta)$ survives when $\langle\lambda_3\widehat P_1\rangle=0$.
- Task 1.2 was revised with the scope of the change: the class used is an on-shell amplitude with the ALP at first order in $1/f_a$ at tree level; the bubble selection above is not used by the estimate.
- Task 2.1, exact tree of $A\to P_{n_1}P_{n_2}P_{n_3}$ at finite mixing in $\mathbb M^{-1}_s(\mathcal C)$ and the invariants ([record](alp/shared/alp-finite-mixing-three-meson-tree.md)): memory DONE, latex MISSING. Pending user approval: the display shorthand $\Phi_\ell=\sum_kv^{(\ell)}_k\lambda_k$, the assumption $\langle\widehat S_1\rangle=0$, and the daughter sign convention (continuity with the unmixed vector).
- Task 3.1, power counting of the dressed-frame vertices in the invariants and loop classes of the amplitude ([record](alp/shared/alp-dressed-vertex-power-counting-invariants.md)): memory DONE, latex MISSING. Each ALP end costs $\epsilon/\sqrt{N(m_A^2)}$; classes at relative $\epsilon^2$: closed $A$ line, $A$ line in bubbles/triangles (cubic $\widehat S_1$ vertex only with a $K^0$-like or $\pi^\pm K^\mp$ partner), meson loops with a mixing insertion, meson one-point functions from a closed $A$ line; $\epsilon^4$: two $A$ lines. Corrected 2026-11-22: $\widehat S_1=-i[c_q^{VS},m_q]$ has zero diagonal (only $d$-$s$ entries); the $A$ one-point from a closed meson line vanishes at this order; the meson one-point keeps a $\widehat S_1$ term through $Y_{ds}(\widehat S_1)_{sd}$ (CP odd). 30 tests pass. Representative-dependent entries of `alp-loop-power-counting.md` marked superseded. Open (user): treatment of one-point functions (kept or removed by vacuum re-alignment at one loop). Inconsistency left: the tree record's (task 2.1) test inputs still use a diagonal $\widehat S_1$; its scripts should be corrected the same way.
- Next step (this change): task 3.2, naturalness assumption and choice of the loop-factor mass per class (design D3), with every symbol defined or flagged.
- Loop-to-tree crossing with the radial mode $\sigma$ (only-up and general sources, regime of interest $\kappa_u\in[0.1,1]$): IN PROGRESS 2026-11-25, OpenSpec change `phi-dependent-loop-tree-crossing-murayama`: groups 1 to 5 done, 6.1 record done. Scope fixed by the user: $\widehat\Theta_1=0$, meson-only loops excluded, full numerator for loops with an internal ALP-like or radial line. Crossing record: [murayama-radial-crossing-implicit](alp/shared/murayama-radial-crossing-implicit.md).
  - Starting point: the [crossing record](alp/shared/alp-loop-crossing-exact-tree-condition.md) is an implicit relation $R(\epsilon)=R_{\max}$ valid for any $\epsilon$, with the exact tree denominator. Its one approximation is the truncation of the loop numerator at the leading mixing-dependent order. The leading-order inversions are secondary and not the result.
  - Main task: remove that truncation. Numerator = complete one-loop amplitude (meson-only loops included) with exact pole projectors, as coefficients of scalar one-loop functions; denominator = exact $\sigma$-inclusive tree; $R_{\max}$ symbolic. The fixed Dalitz point (fixed independent invariants) of the existing record is retained. The renormalization-scale ($\mu$) dependence must be shown to cancel against analytically specified counterterms; the result is $\mu$-independent.
  - Two steps: (1) general result with dressed sources as functions of the complex radial-phase field $\phi$ ($\phi$-valued invariants), of which the $n_q$ case of DiLuzio:2026poh ([coefficient record](alp/shared/alp-murayama-radial-phase-mass-coefficients.md)) is a special case; (2) specialization to the up-only model, $f_a=\sqrt2m_u/\kappa_u$, with $\kappa_u$ and $\lambda$ analytical parameters.
  - Tasks: (T1) $\phi$-valued radial mass and ALP-like pole equation, with the $n_q$ and up-only reductions and the analytic comparison of $m_\sigma$ with $m_a$ (coder); (T2) $\phi$-valued one-loop vertex inventory, polar and Cartesian (coder); (T3) bibliography for scalar one-loop functions and the U(3) and $\phi$-sector counterterm bases (bibliography-manager); (T4) general $\phi$-valued one-loop amplitude with counterterms, $\mu$ cancellation, polar = Cartesian check, $\epsilon\to0$ limit, before any specialization (coder); (T5) implicit crossing without the numerator truncation, general then up-only (orchestrator with brainstormer review, coder for scripts); (T6) Cartesian magnitude cross-check (coder). T1–T3 can start.
  - Caveats only (no calculation): ALP-like/$\eta'$ near-degeneracy region; validity of LO chiral perturbation theory at $m_a\sim1$–$2.3$ GeV.
  - Counterterms and $\mu$: handled by the naturalness assumption, with no absorption analysis, as in the ALP-only record ([record](alp/shared/murayama-radial-naturalness-assumption.md)).
- Next step: coupled-channel S-wave analysis that includes the constant ALP entries together with the growing $PP\to PP$ entries, exact in the mixing. Also decide how to express the $O(x^2)$ dressed sources, which depend on $c_q^A$ and $c_q^P$ separately (user decision). Then compute the full one-loop ALP-meson self-energies including the $aP$ bubbles.

## 2. Interactions with Vector mesons

### 2.1 Even-parity sector

- Resonance ChiPT Lagrangian: memory DONE, latex draft MISSING
- Vector-photon mixing: memory DONE, latex draft MISSING
- Exchange contributions to $a \to 3P$: memory DONE, latex draft MISSING

### 2.2 Odd-parity sector

- Complete Lagrangian: IN PROGRESS
- Fit to the parameters: IN PROGRESS; use the [focused preliminary fit](rcht/odd-parity/odd-preliminary-fit.md) with one fixed mass/data setup and separate source hypotheses, then address material input sensitivities before further optimization.
- Calculation of $a \to 2V$, $a \to \gamma V$, $a \to VPP$, $a\to \gamma PP$ and mixing contribution to $a\to 2\gamma$: MISSING

## 3. Decays into 4 particles

- Efficient integration of phase space for 4 particles final states: FUTURE WORK
- Contributions from 5-pseudoscalar WZW terms: FUTURE WORK
- Contributions from $a \to 2V \to 4P$

## 4. Interactions with axial-vector mesons

FUTURE WORK
