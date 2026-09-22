---
title: Research memory index
status: mixed
created: 2026-09-22
updated: 2026-09-22
bibtex_keys: []
provenance: assumption
---

# Research Memory

Each file records one small, auditable unit of research knowledge. Use lower-case kebab-case names, for example `renormalization-conventions.md`.

```yaml
---
title: Short descriptive title
status: literature # literature | derived | numerical | conjectured | mixed
created: YYYY-MM-DD
updated: YYYY-MM-DD
bibtex_keys: [Author:Year]
provenance: literature # literature | analytical-derivation | numerical-computation | assumption | conjecture
---
```

For derivations, include sections for assumptions, inputs and provenance, method, result, and verification. Link each input to its memory record.
