# HEP-Ph Research Protocol

Use this repository as a reproducible research record. Load only the protocol section and memory files needed for the current task.

## Evidence Record

- Add every consulted source to `bibliography.bib` using its INSPIRE-HEP BibTeX record before relying on it.
- Record each learned concept, definition, equation, or result in a small `memory/*.md` file. Cite its BibTeX key.
- Give every memory file YAML front matter with `title`, `status`, `created`, `updated`, `bibtex_keys`, and `provenance`.
- Valid `status` values are `literature`, `derived`, `numerical`, `conjectured`, and `mixed`.
- Valid `provenance` values are `literature`, `analytical-derivation`, `numerical-computation`, `assumption`, and `conjecture`.
- Label assumptions and approximations explicitly; do not turn them into established results.

## Reproducibility

- Put analytical and numerical work in executable Python scripts under `scripts/`, using `scipy` where appropriate.
- Pair every derivation script with a step-by-step explanation in `memory/`.
- In each derivation record, identify every input and link it to the relevant memory file or cite the stated assumption.
- Keep scripts deterministic, state units and conventions, and include a minimal invocation example.

## Source Handling

- Store downloaded source material under `memory/sources/<bibtex-key>/`; this path is intentionally ignored by Git.
- For arXiv work, download source archives rather than rendered PDFs when available.
- Store modular paper digests in `memory/extracts/<bibtex-key>/`.

## Communication

- Keep terminal output to modified paths and necessary follow-up questions.
- Respond concisely. Do not praise the user.
