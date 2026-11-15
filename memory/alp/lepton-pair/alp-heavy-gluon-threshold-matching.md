---
title: Leading heavy-flavor matching of the ALP gluonic coupling
status: derived
created: 2026-11-15
updated: 2026-11-15
bibtex_keys: [Alda:2026nrd, Messerli:2025rnv]
provenance: analytical-derivation
---

# Scope and inputs

Use the conventions and natural units of the [high-energy conventions](alp-high-energy-tff-conventions.md). The input is the invariant $c_{GP}=c_G+\sum_f c_{ff}^P$ and the heavy-flavor triangle's zero-momentum limit, as recorded in the [triangle calculation](alp-heavy-quark-triangle.md) following Alda:2026nrd. The matching is at leading order in $\alpha_s$ and $1/f_a$, for external gluon momenta and $m_a$ small compared with the heavy-quark mass $m_h$. It does not establish local matching in the charmonium region.

# Derivation

The heavy-quark triangle's pseudoscalar term approaches the constant proportional to $c_{hh}^P$, while its derivative term vanishes, as $m_h$ tends to infinity. Equivalently, matching the invariant $c_{GP}$ above and below the threshold gives
$$
c_G^{(n_f-1)}+\sum_{f\ne h}c_{ff}^P=c_G^{(n_f)}+\sum_f c_{ff}^P,
\qquad c_G^{(n_f-1)}=c_G^{(n_f)}+c_{hh}^P.
$$
For $c_{GA}=c_G+\sum_f c_{ff}^A$ this implies
$$
c_{GA}^{(n_f-1)}=c_{GA}^{(n_f)}-c_{hh}^A+c_{hh}^P.
$$
This is a source matching relation, not a derivation of the renormalized singlet-current matching coefficient. Its corrections include powers of external momenta over $m_h$ and QCD corrections. When $m_a$ approaches or exceeds $2m_c$, local charm matching is not justified by this expansion.

# Verification and remaining limitations

The [script](../../../scripts/derive_alp_symmetric_tff.py) solves the matching equation, checks invariance under the stated heavy-flavor rotation, and checks the heavy triangle's limits and the two Ward contractions and Bose symmetry of the leading tensor. Run:

```text
.venv/bin/python scripts/derive_alp_symmetric_tff.py
.venv/bin/python scripts/derive_alp_symmetric_tff.py --test
```

The script also prints a conditional sum of the contact term, heavy triangle, and light axial-current term from the [hard-current OPE](alp-hard-current-ope.md). This sum is not established as the complete coefficient: the overlap between the explicit hard light-quark contribution and the renormalized axial matrix element must be fixed, and the [singlet-current record](alp-singlet-axial-current-evolution.md) leaves the gluonic Wilson coefficient and threshold matching open. The [heavy-flavor record](alp-heavy-flavor-thresholds.md) states the threshold-continuity assumption. No numerical uncertainty is assigned to these omissions.
