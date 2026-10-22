---
title: Research memory index
status: mixed
created: 2026-09-22
updated: 2026-10-19
bibtex_keys: []
provenance: assumption
---

# Research Memory

Each record is one auditable concept, derivation, calculation, or decision. Records are grouped below by topic; the records themselves retain their complete original content and front matter. Use lower-case kebab-case names, for example `renormalization-conventions.md`.

## Topic indexes

- [Plan](plan.md) - Long-term project planning.
- [ALP](alp/INDEX.md) — ALP derivations and decay calculations organized by pseudoscalar, photon, and vector final states, with related shared conventions and audits.
- [Kaon states](kaon/INDEX.md) — neutral-kaon projections and related selection-rule diagnostics.
- [K pi rescattering](kpi/INDEX.md) — production kernels, Omnes functions, and isospin breaking.
- [RChT](rcht/INDEX.md) — foundations, odd-sector conventions and operators, normalized kernels, literature prescription matching and flat inputs, focused original-operator preliminary fit, conditional calibration/profiles, finite ALP target/null audits and matching/completion plans.
- [Literature extracts](extracts/INDEX.md) — source-backed modular paper digests, including Alda:2026nrd's commutator/source conventions and the Bauer photon-operator normalization and tensor convention.

The ignored [source archive](sources/INDEX.md) contains downloaded source material rather than research records.

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
