---
name: bibliography-management
description: Catalog and digest HEP literature with INSPIRE-HEP metadata and auditable local source storage.
---

# Bibliography Management

## Intake

1. Search INSPIRE-HEP for the intended work and retrieve its BibTeX entry.
2. Add the unmodified, source-provided entry to `bibliography.bib`; use its citation key for all records.
3. For arXiv papers, download the source archive to `memory/sources/<bibtex-key>/`. Otherwise cache available source material there.
4. Never commit downloaded source material.
5. If a paper is not available on arXiv, try downloading the source from the links in INSPIRE-HEP, or from ResearchGate. As a last fallback, use Sci-Hub.
6. If the retrieved source is in PDF format, first convert it to markdown using a suitable tool such as `pdf2text` (if the PDF is machine-readable) or an OCR like `tesseract` (if it is a scanned document). Store the resulting markdown as `memory/extracts/<bibtex-key>/fulltext.md`.

## Digest

1. Create `memory/extracts/<bibtex-key>/`.
2. Write small topic-specific Markdown extracts with standard memory YAML front matter.
3. Cite the paper's BibTeX key beside every learned definition, equation, claim, convention, or data point.
4. Preserve equation labels or source locations when available. Mark interpretation, extrapolation, and uncertainty explicitly.

## Checks

- Confirm each cited key exists exactly once in `bibliography.bib`.
- Confirm every extract has `status: literature`, the source key, and `provenance: literature` unless it contains clearly separated analysis.
- Do not replace a source result with an uncited summary.
