# Tasks

## 1. Literature, conventions and process inventory

- [ ] 1.1 Delegate to `bibliography-manager`: add INSPIRE BibTeX entries and `memory/extracts/<key>/` digests for the axial-RChT references (Ecker:1988te, Jamin:2000wn, Roig:2013baa, GomezDumm:2010ra and any others it finds); verify each key exists in `bibliography.bib` and has an extract INDEX entry
- [ ] 1.2 Create `memory/rcht/axial-even-parity/` with `INDEX.md`, link it from `memory/rcht/INDEX.md`; verify the indexes follow the AGENTS.md navigation contract
- [ ] 1.3 Write the axial inventory record (axial external states and axial exchange for every existing ALP process, with exclusions and reasons); verify every process in `memory/alp/INDEX.md` is classified
- [ ] 1.4 Write the field/flavor/phase convention record (tensor $A_{\mu\nu}$, nonet matrix, ideal basis, Levi-Civita/metric), including the dictionary to Alda:2026nrd; verify the matrix and normalization are reproduced by a script test (`coder`)

## 2. Quadratic sector

- [ ] 2.1 `coder`: script deriving kinetic and mass terms, neutral isoscalar mixing angle and generator rotation for the axial nonet; verify ideal limit and run with the documented invocation
- [ ] 2.2 Write the paired memory record with assumptions on omitted mass operators; verify front matter, `$` delimiters and index entry
- [ ] 2.3 Write the power-counting and width-treatment record (narrow-width criteria per state, exchange propagator prescription); verify it states validity limits and cites sources

## 3. Operator basis

- [ ] 3.1 `coder`: script enumerating even-parity one-axial-field operators from $u_\mu$, $f_\pm^{\mu\nu}$, $\chi_\pm$, $\nabla_\mu$ and reducing by Schouten/Bianchi/EOM; verify the independent operator count and that it contains the literature operators
- [ ] 3.2 Write the operator-basis record with normalizations, completeness argument and the independent-versus-fixed coefficient table; verify each coefficient is marked and any new notation is listed for user approval
- [ ] 3.3 Ask the user to approve new coefficient notation before the vertex records use it; verify approval is recorded

## 4. Vertices and equivalences

- [ ] 4.1 `coder`: script expanding the Lagrangian to extract $A$–$P$, $A$–$V$, $AVP$, $APP$, $A\gamma P$, $AV\gamma$ and ALP-source vertices and removing quadratic mixing; verify symbolic diagonality check passes
- [ ] 4.2 `coder`: script deriving tensor↔Proca field redefinitions and induced contacts; verify amplitude equality test on a representative process
- [ ] 4.3 `coder`: Ward-identity, crossing and ideal-limit tests for the source vertices; verify all tests pass
- [ ] 4.4 Write derivation records for 4.1–4.3 with step-by-step explanation, input links and verification section; verify each record names its script and invocation

## 5. Handoff and integration

- [ ] 5.1 Write the coefficient/contact handoff table mapping each independent coefficient and contact to vertices and inventory channels; verify every inventory channel is covered or excluded
- [ ] 5.2 Update `memory/plan.md` §4, `memory/rcht/INDEX.md` and parent indexes; run all scripts once and verify every link and index entry resolves
- [ ] 5.3 Update bead alp-ct8.1 with results and propose a commit message; verify `bd show alp-ct8.1` reflects completion state (close only if authorized)
