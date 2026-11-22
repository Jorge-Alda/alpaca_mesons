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
- Next step (this change): task 1.2, select the observable or state that none passes, after fixing the reading of test (c) (absorption by a local counterterm polynomial in the Lagrangian parameters) and checking the conjectured non-absorbable pieces; consult the brainstormer on whether D3 must be revised.
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
