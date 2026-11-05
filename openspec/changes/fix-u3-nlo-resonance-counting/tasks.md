# Tasks

## 1. Setup and inputs

- [x] 1.1 Claim beads issue alp-gii.1 (`bd update alp-gii.1 --claim`) and verify that `bd show alp-gii.1` reports it in progress
- [x] 1.2 Ask the user to approve a notation for residual local constants (design Open Questions) and verify that the approved symbol is recorded in the design and used nowhere before approval
- [x] 1.3 Through `bibliography-manager`, check that every source cited by the design (Kaiser:1998ds, Kaiser:2000gs, Herrera-Siklody:1996tqr, Gasser:1983yg, Ecker:1988te, Ecker:1989yg, Rosell:2004mn, Pich:2010sm, Alda:2026nrd) is in `bibliography.bib` with an INSPIRE record. Add any further record consulted for the $N_c$ scaling of resonance couplings. Verify with `grep` on each key
- [x] 1.4 Create `memory/nlo/counting/INDEX.md` with front matter. Add a topic-level entry to `memory/nlo/INDEX.md`, and verify that all relative links resolve

## 2. Field content and power counting

- [x] 2.1 Write `memory/nlo/counting/field-content.md`: the nonet with $\eta_0$, the ALP sources of Alda:2026nrd, and the $S$, $P$, $V$, $A$ nonets (spin one as antisymmetric tensors), with chiral transformations and links to their LO Lagrangians. Mark missing $S$/$P$/$A$ LO records as inputs from Ecker:1988te. Verify that the front matter has all six required fields
- [x] 2.2 Through `coder`, create `scripts/derive_nlo_u3_resonance_counting.py` with an order table that assigns $\delta$ orders to every operator in `memory/nlo/lagrangian.md`, to the resonance couplings and masses, and to the low-energy limit of single-resonance exchange. Verify that `.venv/bin/python scripts/derive_nlo_u3_resonance_counting.py --orders` prints the table and asserts $L_5,L_8,\Lambda_{1,2}\to\delta$, $L_4,L_6,L_7,L_{18},L_{25}\to\delta^2$, and exchange limit $\to\delta$
- [x] 2.3 Write `memory/nlo/counting/power-counting.md`, which explains step by step the counting rules and the output of 2.2, including unexpanded resonance propagators, and labels assumptions. Verify that every order in the record matches the script output

## 3. Loop classes

- [x] 3.1 Through `coder`, extend the script with a classification of one-loop topologies (Goldstone-only with LO vertices; with an ALP-source vertex; with an internal resonance line; with a resonance vertex; Goldstone-loop dressing of a resonance propagator), giving each an order and an included/excluded flag. Verify that `--loops` prints the classification and asserts that only the Goldstone-only classes are included
- [x] 3.2 Write `memory/nlo/counting/loop-classes.md`, which explains 3.1, cites Rosell:2004mn and Pich:2010sm for the excluded resonance loops, and gives the width-prescription rule with a link to `memory/rcht/foundations/vector-width-lineshape-prescriptions.md`. Verify that the links resolve and the classification matches the script

## 4. Subtraction scheme and counterterm basis

- [x] 4.1 Through `coder`, add the $\lambda(\mu_\chi)$ scheme, the $U(3)$ $\Gamma_i$ table ($L_4$–$L_8$, $L_{18}$, $L_{25}$, $L_{47}$; $\Lambda_{1,2}$ without a pole) and the one-loop running. Include a check that $L_i^r(\mu_2)-L_i^r(\mu_1)=-\Gamma_i\log(\mu_2/\mu_1)/16\pi^2$ and the $Z_A$ labels for $\Lambda_1$, $\Lambda_2$, $2L_5+3L_{18}$ and $2L_8-3L_{25}$. Verify that `--scheme` runs and that every assertion passes for $\mu_\chi\in\{0.5,0.77,1.0\}\,$GeV
- [x] 4.2 Write `memory/nlo/counting/subtraction-scheme-and-scales.md`: the default $\mu_\chi=M_\rho$, the variation range, and the separation of the QCD scale. Translate between the conventions of Kaiser:1998ds and Herrera-Siklody:1996tqr. Verify that the coefficients agree with `memory/nlo/one-loop-effects-and-renormalized-lecs.md` and with the script
- [x] 4.3 Write `memory/nlo/counting/counterterm-basis.md`, a table with one row per constant: $\delta$ order, $\Gamma_i$, $Z_A$ behaviour, contributing resonance multiplets, and the truncation under which a zero $\Gamma_i$ holds. Include a marked row for further ALP-source operators. Verify that every row is reproduced by `--scheme --orders`

## 5. Resonance exchange and residual constants

- [x] 5.1 Write `memory/nlo/counting/resonance-exchange-and-residual-constants.md`. Use the approved notation from 1.2 to define residual constants, the separation scale, the $\mu_\chi$ independence of exchange, and the labelled saturation assumption. Refer the matching to alp-gii.3. Verify consistency with the explicit-resonance rule in `memory/nlo/resonance-saturation-and-dynamical-masses.md`
- [x] 5.2 Add one-line cross-links from `memory/nlo/one-loop-effects-and-renormalized-lecs.md` and `memory/nlo/resonance-saturation-and-dynamical-masses.md` to the new records, without altering their content. Update their `updated` field and verify the links

## 6. Validity range and truncation

- [ ] 6.1 Through `coder`, add the truncation prescription $\delta_{\rm est}(s)=\max(1/N_c,s/\Lambda_\chi^2,M_K^2/\Lambda_\chi^2)$, with the excluded resonance-loop size added in quadrature. Verify that `--truncation 0.5 1.0` prints $\delta_{\rm est}$ at $\sqrt s=0.5$ and $1.0\,$GeV with all inputs stated in GeV
- [ ] 6.2 Write `memory/nlo/counting/validity-range-and-truncation.md` covering the range up to about $1\,$GeV, the failure near resonance poles and above $1\,$GeV, and the inputs of the prescription with their sources. Verify that the numbers match 6.1

## 7. Integration check

- [ ] 7.1 Run `.venv/bin/python scripts/derive_nlo_u3_resonance_counting.py` (all checks) and verify exit code 0. Check that every `memory/nlo/counting/*.md` record has valid front matter and status/provenance values, and that `memory/README.md` → `memory/nlo/INDEX.md` → `memory/nlo/counting/INDEX.md` links resolve
- [ ] 7.2 Check each acceptance criterion of alp-gii.1 against the records (loop classes, counterterms, scale/scheme, singlet treatment, validity range). Then close the issue with `bd close alp-gii.1`, and verify with `bd ready` that alp-gii.2, alp-gii.3 and alp-gii.8 are no longer blocked by it

## Workflow follow-up

- Propose a commit message for the changed files. Do not commit without explicit approval.
- Archive the change with `/openspec-archive-change` after review.
