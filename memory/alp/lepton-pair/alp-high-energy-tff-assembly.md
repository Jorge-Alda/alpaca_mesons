---
title: Assembly of the ALP two-photon form factor from quark flavors and virtuality regimes
status: mixed
created: 2026-11-15
updated: 2026-11-15
bibtex_keys: [Alda:2026nrd, Messerli:2025rnv]
provenance: analytical-derivation
---

# Scope

This record states how the pieces derived in the records of the high-energy completion add up to $F_{a\gamma^*\gamma^*}(q_1^2,q_2^2)$, which pieces exist, and which do not yet.
It contains no new derivation and writes an expression only where it is complete and verified in the record that supplies it. Conventions: [shared conventions](alp-high-energy-tff-conventions.md);
$F_{a\gamma^*\gamma^*}$ is the scalar of $i\mathcal M=ie^2RF_{a\gamma^*\gamma^*}$ [@Messerli:2025rnv].

# Decomposition by quark flavor

At leading order in $\alpha_s$ and order $1/f_a$,
$$
F_{a\gamma^*\gamma^*}(q_1^2,q_2^2)=F^{(u,d,s)}(q_1^2,q_2^2)+F_c(q_1^2,q_2^2)+F_b(q_1^2,q_2^2).
$$

- $F^{(u,d,s)}$ is the part of the form factor that comes from the $u,d,s$ quarks, from the photon operator with coefficient $c_\gamma$ (counted once, here), and from the gluonic coupling. It is described by the
  hadronic representation at low virtualities (the [tree TFF](alp-doubly-virtual-tff.md) and the other low-energy records indexed in [INDEX.md](INDEX.md)) and by the short-distance expansion at high virtualities.
- $F_h$, $h=c,b$, is the contribution of the heavy flavor: the one-loop triangle with its ALP couplings,
  $$
  F_h=\frac{2N_cQ_h^2}{4\pi^2f_a}\Bigl[c_{hh}^P+\bigl(c_{hh}^A-c_{hh}^P\bigr)B_4\Bigr],
  $$
  with $B_4=B_4(4m_h^2/m_a^2,\,4m_h^2/q_1^2,\,4m_h^2/q_2^2)$ of the [function record](alp-loop-function-b4.md); see the [triangle record](alp-heavy-quark-triangle.md).
  $F_h$ is a contribution to $F_{a\gamma^*\gamma^*}$ and is added to $F^{(u,d,s)}$, not substituted for part of it.

**No double counting.** The hadronic representation of $F^{(u,d,s)}$ has only $u,d,s$ degrees of freedom and the photon contact term enters it once, through $c_\gamma$ (the coefficient in the Lagrangian). The loop of $h$ is a separate
diagram and appears only in $F_h$. The derivative-free axial matrix element of $h$ is not used for $h=c,b$, because the local expansion with that operator alone does not reproduce $F_h$
([triangle record](alp-heavy-quark-triangle.md)).

**Where additivity is not complete.**
1. The heavy flavors with a pseudoscalar coupling change the effective ALP–gluon coupling seen by the three-flavor theory. The three-flavor sources $\bm{\mathcal C}$ and $c_{GA}$ of the
   [axial-current record](alp-axial-current-matrix-element.md) must then be evaluated with that effective gluonic coupling, not with the gluonic coupling of the full theory. The relation between the two
   is not derived in these records.
2. Corrections of relative order $\alpha_s$ to $F_h$ and mixed diagrams with gluons connecting the heavy-quark loop to the light sector are not included.
3. For $m_a>2m_c$ the free-charm triangle has a $c\bar c$ cut and no resonances; the treatment of charm in the charmonium region is an open decision.

# Pieces by virtuality regime

| Regime | Piece | Supplied by | Status |
|---|---|---|---|
| Low virtualities | $F^{(u,d,s)}$, hadronic | tree, isovector, mixed-vector and channel-input records of [INDEX.md](INDEX.md) | under construction; the vector-fit calibration is not converged |
| Low virtualities | $F_h$ | triangle record | derived: $F_h\to\frac{2N_cQ_h^2}{4\pi^2f_a}c_{hh}^P$ plus corrections suppressed by $1/m_h^2$, with the first correction $-\frac{2N_cQ_h^2}{4\pi^2f_a}(c_{hh}^A-c_{hh}^P)\frac{m_a^2+q_1^2+q_2^2}{12m_h^2}$ |
| High virtualities | operator contribution of $u,d,s$ with its matrix element, singlet evolution, and $\alpha_s$ coefficient | [axial-current record](alp-axial-current-matrix-element.md), [singlet evolution record](alp-singlet-axial-current-evolution.md) | derived for nearly symmetric virtualities, at the order and with the named inputs stated there |
| High virtualities | term independent of $q^2$ from the derivative couplings of $u,d,s$ | massless limit of the triangle (per massless quark, $\frac{2N_cQ_q^2}{4\pi^2f_a}c_{qq}^A$), anchored to the repository's $C_{a\gamma\gamma}$ | the per-quark limit is derived; its assembly with $c_\gamma$ and the operator contribution, and the Bose and Ward checks, are task 3.1 and are not done |
| High virtualities | $F_h$ | triangle record | derived: $F_h\to\frac{2N_cQ_h^2}{4\pi^2f_a}c_{hh}^A$ with corrections of relative order $\frac{m_h^2}{|q^2|}\ln\frac{|q^2|}{m_h^2}$ (equal virtualities) or $\frac{m_h^2}{|q^2|}\ln^2\frac{|q^2|}{m_h^2}$ (one real photon) |
| One virtuality large, the other fixed | $F^{(u,d,s)}$ | light-cone expansion | not derived (task 3.2) |
| Between the hadronic and hard regimes | join of $F^{(u,d,s)}$ | matching (tasks 4.1 to 4.3) | not done |

No expression for the total $F_{a\gamma^*\gamma^*}$ in the hard region is written here, because the part of $F^{(u,d,s)}$ listed as task 3.1 is not yet assembled.

# Invariance

Under an ALP-dependent axial rotation, $c_q^A\to c_q^A+\kappa$ shifts the explicit couplings $c_\gamma$ and $c_{hh}^A,c_{hh}^P$ together; the combinations $c_{hh}^A-c_{hh}^P$ and
$c_\gamma+2N_c\langle Q^2c_q^A\rangle+2N_c\sum_hQ_h^2c_{hh}^A$ are unchanged. For the heavy flavors this is verified in the script of the triangle record; for the light flavors it is the
statement of the [LSZ record](../../rcht/odd-parity/odd-alp-source-lsz-assembly.md). The total form factor must be invariant, which is the check to be applied when 3.1 is assembled.
