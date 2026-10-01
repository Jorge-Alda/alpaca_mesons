# HEP-Ph Research Protocol

Use this repository as a reproducible research record. Load only the protocol section and memory files needed for the current task.

## Hierarchical Memory and Retrieval

- Organize research notes under the hierarchical topic structure rooted at `memory/`. New records must live in the most specific existing topic subfolder; create a new nested topic folder only when no existing branch is appropriate.
- Begin every memory retrieval by reading `memory/README.md`, then crawl the relevant `INDEX.md` files from the root toward the target topic. Load only the indexed records needed for the task; do not scan or bulk-load unrelated memory files, extracts, or source archives.
- Treat every `INDEX.md` as a navigation contract. Each directory containing memory records must have an index linking to its immediate files and subfolders, with a one-line description for every entry. Keep links relative and valid after moves.
- When creating, splitting, renaming, or moving a memory record, update the containing index and every affected parent index in the same change. Preserve the record's YAML front matter, equations, citations, evidence labels, assumptions, provenance, and substantive content.
- Keep records atomic: one concept, derivation, calculation, result, or decision per file. Split mixed records when doing so does not lose context; carry the relevant front matter and references into each resulting record and link related records through their indexes or text.
- Literature extracts remain under `memory/extracts/<bibtex-key>/` and downloaded material under `memory/sources/<bibtex-key>/`; retrieve them through their indexes only when the selected memory record requires source-level detail.

## Evidence Record

- Add every consulted source to `bibliography.bib` using its INSPIRE-HEP BibTeX record before relying on it.
- Record each learned concept, definition, equation, or result in a small `memory/*.md` file. Cite its BibTeX key.
- Give every memory file YAML front matter with `title`, `status`, `created`, `updated`, `bibtex_keys`, and `provenance`.
- Valid `status` values are `literature`, `derived`, `numerical`, `conjectured`, and `mixed`.
- Valid `provenance` values are `literature`, `analytical-derivation`, `numerical-computation`, `assumption`, and `conjecture`.
- Label assumptions and approximations explicitly; do not turn them into established results.
- In Markdown files, always use `$...$` and `$$...$$` as LaTeX delimiters.

## Reproducibility

- Put analytical and numerical work in executable Python scripts under `scripts/`, using `scipy` where appropriate.
- When both analytical and numerical derivations are available, prefer the analytical derivation.
- Do not treat an analytical derivation as a check: delegate its calculation to an executable Python script.
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
- Use the same conventions and definitions as Alda:2026bar
