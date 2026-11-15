# Tasks

## 1. Literature and conventions

- [x] 1.1 Locate and register in `bibliography.bib` (INSPIRE BibTeX, via `bibliography-manager`) the sources for the QCD short-distance behavior of pseudoscalar TFFs (Wilson coefficient and its $\alpha_s$ correction, light-cone expansion, singlet axial current anomalous dimension, flavor decoupling); write extract records; verify every cited key exists in `bibliography.bib` and each extract is indexed.
- [x] 1.2 Record the conventions shared with Alda:2026nrd and Messerli:2025rnv (sources $U_a,Y_a$, tensor $R$, metric, $\epsilon^{0123}$) in a short memory record; verify it links the existing hard-current and tree-TFF records.

## 2. Source-resolved matrix elements and scheme

- [x] 2.1 Derive $\langle0|\bar q_f\gamma_\sigma\gamma_5q_f|a(p)\rangle$ for $f=u,d,s$ at LO ChPT in terms of $U_a,Y_a$ and the contact coefficient, split into isovector, octet, and singlet parts; verify with a SymPy script (reproducing the $4/9,1/9,1/9$ charge weights) and a paired memory record, and that no meson decay constant is used as an ALP input.
- [x] 2.2 Derive the role of the singlet current, anomaly, and gluonic operator, and the scheme/scale and renormalization-group evolution of the matrix elements at the order fixed in the literature; verify with a script checking the RG solution against its differential equation and a record stating scheme, scale, and remaining named inputs.
- [x] 2.3 Derive the heavy-flavor ($c,b$) threshold matching and the below-threshold heavy-quark contribution; verify the continuity of the coefficient across thresholds symbolically and give the active-flavor table for $2m_e<m_a<2m_\tau$.

## 3. Asymptotic TFF in each regime

- [x] 3.1 Complete the nearly symmetric limit: assemble the source-linear coefficient from tasks 2.1–2.3 and the existing local OPE; verify against `scripts/derive_alp_hard_current_ope.py --test` and Bose/Ward constraints.
- [x] 3.2 Derive the fixed-ratio doubly virtual limit and the singly virtual limit by the light-cone expansion, with the leading-twist structure and the declared distribution-amplitude assumption; verify that it reduces to task 3.1 at equal virtualities and that the assumption enters as a named parameter.
- [ ] 3.3 Derive the leading $m_a^2/Q^2$ corrections and the timelike/Euclidean continuation statement for the hard kernel; verify the $m_a\to0$ and large-$Q^2$ limits by script and state the supported virtuality domain per regime.
- [ ] 3.4 Cross-check the method (not the coefficient) against the known $\pi^0$ asymptotic structure by running the same derivation with pion sources; verify agreement and that the ALP result is not built from it.

## 4. Matching to the low-energy form

- [ ] 4.1 Define the low-energy benchmark forms for matching tests (tree-level poles from `alp-doubly-virtual-tff.md` and a symbolic dispersive form) independent of the vector fit; verify each is labelled conditional and linked to its record.
- [ ] 4.2 Derive the matching prescriptions (sharp join, interpolation, effective pole) preserving $F(0,0)$ and discontinuities and avoiding double counting; verify with a script that tests normalization, continuity, and absence of double-counted strength.
- [ ] 4.3 Quantify the matching-scale variation and the spread among alternative completions on the benchmark forms; verify the variations are reported separately from hadronic and scheme uncertainty.

## 5. Records and closure

- [ ] 5.1 Write the paired memory records for tasks 2–4 under `memory/alp/lepton-pair/` with required front matter, and update the containing, parent, and `memory/README.md` indexes; verify links resolve.
- [ ] 5.2 Update the status in `alp-isoscalar-tff-input-gates.md` and `alp-lepton-pair-calculation-plan.md`; list unresolved inputs (singly virtual shape, fit-dependent matching) as blockers on `alp-6os.1`; verify all scripts run and `--test` passes.
