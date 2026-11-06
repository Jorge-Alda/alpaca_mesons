# Tasks

## 1. Goodness of fit at accepted minima

- [ ] 1.1 (coder) Add an assessment script that loads the four accepted endpoints from `joint_vector_endpoint_comparison.json`, refuses non-accepted endpoints, and writes total/per-block $\chi^2$, dof counting, nominal $p$-values, full standardized residuals and ranked largest pulls at $(4,128,128)$; verify with a focused unittest including the refusal case.
- [ ] 1.2 Run it and write `data/vector_sector/joint_vector_fit_adequacy.json`; verify the total $\chi^2$ matches the comparison ledger to $10^{-8}$.
- [ ] 1.3 Record per-block results, pull rankings and the BaBar residual energy pattern per channel in a memory record linked from `memory/rcht/odd-parity/INDEX.md`; verify index links.

## 2. Block tensions

- [ ] 2.1 (coder) Add block-exclusion refits warm-started from `basin_probes/even_shared_widths`, reporting optimizer status, remaining-block $\chi^2$, excluded-block pulls and near-singular columns; verify with a test on a reduced block.
- [ ] 2.2 Run all six exclusions and store the ledger; verify each entry is independently interpretable.

## 3. Sensitivity variations

- [ ] 3.1 (coder) Add likelihood options for the fixed variation list in `design.md`, with defaults reproducing the baseline; verify bitwise-equal residuals at defaults by unittest.
- [ ] 3.2 Run each variation (with basin-switch reruns where scaled distance exceeds 1), high-order re-evaluation, and acceptance gates; store all outcomes including failures.
- [ ] 3.3 Record the variation table ($\chi^2$, dof, $p$, block $\chi^2$, parameter shifts in baseline local errors) in memory; verify against the ledger.

## 4. Verdict

- [ ] 4.1 Write the adequacy verdict with the stated criterion; if inadequate, name the responsible blocks/rows and the needed extension and create a follow-up bead under alp-uzh; verify the bead exists and is linked.
- [ ] 4.2 Run the focused test suite and confirm historical fit artifacts are unchanged; update alp-uzh.6 and alp-uzh.7 accordingly.
