---
title: KK pi vector-current coefficient rank audit
status: mixed
created: 2026-10-29
updated: 2026-10-29
bibtex_keys: [GomezDumm:2010ra, Alda:2026nrd, Ruiz-Femenia:2003jdx]
provenance: numerical-computation
---

# Scope and inputs

This audit asks whether the vector-current form factors for
$\tau^-\to K^+K^-\pi^-\nu_\tau$ and
$\tau^-\to K^-K^0\pi^0\nu_\tau$ cover the VJP/VVP directions and the
remaining VPPP combinations identified in the
[ALP vector-channel overview](alp-vector-channel-coefficient-overview.md).
The channel formulas are the $F_3$ expressions in Eqs. (30) and (35) and
the $C^R,C^{RR},D^R,E^R$ polynomials in Appendix A of [@GomezDumm:2010ra].
The source extract records the literal five-operator VPPP basis and its
normalization [here](../../extracts/GomezDumm:2010ra/vppp-basis-and-omega-three-pion.md).

Assume tree level, exact strong isospin, ideal $\omega$--$\phi$ mixing,
nonzero fixed even-parity couplings, and independently variable odd
coefficients. Do not impose the short-distance relations of
[@GomezDumm:2010ra] as data. The canonical Levi-Civita/coupling conversion
reverses all odd coefficients together, so it leaves these ranks invariant.
The calculation treats zero-width propagators as rational functions for the
coefficient audit. It is not a fit to measured $\tau$ or $e^+e^-$ data.

## Reproduction

Run:

```text
.venv/bin/python scripts/audit_kkpi_coefficient_rank.py
```

The [script](../../../scripts/audit_kkpi_coefficient_rank.py) builds exact
SymPy coefficient matrices. It keeps meson masses symbolic for the stripped
sector maps and uses rational generic fixed inputs as a finite exact
certificate for full rational-amplitude ranks. Polynomial coefficient rows
represent variation over kinematics; they are not statistically independent
observations. The [transcript](../../../data/vector_sector/kkpi_coefficient_rank.txt)
stores the printed matrices' ranks and nullspaces. Five focused unit tests
are in [the paired test file](../../../scripts/test_audit_kkpi_coefficient_rank.py).

# Sector projections

Write the VPPP coordinates as $(g_1,g_2,g_3,g_4,g_5)$. The local charged-mode
polynomial $D^R$ spans
$$
\left\{g_1-g_3,\quad g_2,\quad
(m_K^2+m_\pi^2)g_4+m_K^2g_5\right\},
$$
with rank three when its energy dependence is retained. At a single fixed
$Q^2$, it has rank two. The neutral-mode polynomial $E^R$ spans only
$$
g_1-g_3+2g_2,
$$
so its rank is one. Combining the two $KK\pi$ local projections remains rank
three. The two isoscalar $3\pi$ combinations
$$
A_{3\pi}=g_1-g_2-g_3,\qquad B_{3\pi}=2g_4+g_5
$$
raise the charged or combined VPPP rank to four when
$m_K^2\ne m_\pi^2$; with the neutral mode alone they raise it only to three.
The charged mode plus these isoscalar combinations therefore spans the full
four-dimensional VPPP space in the overview. In the exact $SU(3)$ mass limit,
the charged-plus-isoscalar rank falls to three, so the separation of $g_4$
from $g_5$ specifically relies on mass breaking.

The two VJP numerators in the channels span four combinations:
$$
c_1-c_2+c_5,\quad
-(c_1-c_2-c_5+2c_6),\quad
c_1+c_2+8c_3-c_5,\quad c_4.
$$
The VVP numerator pair spans $d_3$ and $d_1+8d_2-d_3$. Thus the $KK\pi$
form factors contain the preliminary fit's missing VJP direction
$c_1-c_2+c_5$, while the preliminary fit already supplies the additional
independent VJP and VVP combinations needed to separate their sector maps.
The $KK\pi$ VVP projection alone has rank two and does not resolve all three
overview VVP directions.

# Full form-factor rank and fit combination

Ranks below are for coefficient functions across kinematics, with the full
$F_3$ pole structure retained. “Preliminary span” means the seven measured
coefficient rows in the original low-energy fit. “Isoscalar span” means the
two rows $A_{3\pi},B_{3\pi}$, included as independent information for this
rank comparison; the repository's proposed isoscalar fit has not yet
produced a numerical result.

| Included rows | $K^+K^-\pi^-$ | $K^-K^0\pi^0$ | Both modes |
|---|---:|---:|---:|
| Local VPPP only | 3 | 1 | 3 |
| Local VPPP plus isoscalar span | 4 | 3 | 4 |
| Full $F_3$ form factor | 8 | 4 | 8 |
| Full $F_3$ plus preliminary span | 11 | 9 | 11 |
| Full $F_3$ plus preliminary and isoscalar spans | 12 | 11 | 12 |

The twelve is the maximum relevant coefficient rank, namely the overview's
five VJP, three VVP, and four VPPP combinations. The charged mode together
with the preliminary and isoscalar spans reaches this rank; adding the neutral
mode does not increase it further. The neutral mode remains one direction
short of full rank even after those rows, because it cannot separate the two
mass-source VPPP coefficients. Without the preliminary information, the full
$KK\pi$ form factors have rank eight, below the sum of the separated VJP,
VVP, and VPPP sector ranks: one VJP/VPPP displacement compensates between the
sectors in both full amplitudes. The script gives an explicit symbolic
compensation. The preliminary rows remove that degeneracy.

These are amplitude ranks, not information ranks. BaBar Tables 2 and 3 are
now transcribed as inclusive charge-resolved spectra and fitted jointly in
[babar-kkpi-unified-vector-coefficient-fit.md](babar-kkpi-unified-vector-coefficient-fit.md).
The full electromagnetic $KK\pi$ form factors retain VPPP sensitivity, and
the integrated-rate coefficient Jacobian has rank ten in the chosen
combination basis. The four $g_i$ combinations remain weakly constrained:
their profile minima include zero and their outer profile points have
optimizer/bound failures. The $K^*K$ isobar-only Table 7 fit still has zero
local VPPP sensitivity.
