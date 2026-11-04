# Tasks

## 1. Establish the reproducible audit inputs

- [x] 1.1 Inventory the global-fit and conditional $G_V$ inputs, their source paths and hashes, row definitions, mass roles, width-sharing rules, and covariance assumptions; verify the inventory against the stored JSON metadata and fitting code.
- [x] 1.2 Ask the coder agent to implement a deterministic audit script under `scripts/` that reads the existing fit and measurement artifacts without overwriting them, and emits machine-readable block and comparison results; verify it runs from the documented invocation.

## 2. Quantify the shared-block diagnosis

- [x] 2.1 Use the audit script to reproduce the stored joint block contributions and selected shared-even residual diagnostic; verify the reported values agree with the baseline JSON to numerical precision and remain explicitly diagnostic subsets.
- [x] 2.2 Compute the charged/neutral $K^*$ and correlated $\phi\to K\bar K$ pair minima under their documented covariance, and evaluate only source-backed mass, width, saturation, or covariance alternatives; verify each result records its fixed/free inputs and likelihood comparison.
- [x] 2.3 Record shared-width, duplicate-measurement, $K^*\to K\pi$ saturation, and phase-space versus tensor-mass conventions, including any proposed alternative that cannot be evaluated due to missing source information; verify no conditional $G_V$ $\chi^2$ is added to the joint total.

## 3. Document and preserve the result

- [x] 3.1 Add one atomic numerical memory record with input links, source keys, assumptions, equations or calculation steps, output paths, reproduction command, and limitations; verify its front matter and cited inputs follow the repository memory protocol.
- [x] 3.2 Update `memory/rcht/odd-parity/INDEX.md` and every affected parent index with relative links and one-line topic descriptions; verify all new links resolve.
- [x] 3.3 Run the documented audit command from a clean stored-input state and verify the output hashes and summary match the memory record; report any baseline input-hash mismatch rather than silently updating the baseline.
