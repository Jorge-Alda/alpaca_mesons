---
title: ArXiv Alda:2026nrd neutral generators in the amplitude kernels
status: numerical
created: 2026-10-13
updated: 2026-10-13
bibtex_keys: [Alda:2026nrd]
provenance: numerical-computation
---

# Alda:2026nrd inputs for neutral amplitudes

## Inputs and convention

Use the [arXiv-v1 source record](../../extracts/Alda:2026nrd/neutral-meson-generators.md),
including its mass matrix, rounded fitted central parameters and displayed
generator phases. The [fixed selected-input JSON](../../../data/vector_sector/alda_2026nrd_arxiv_v1.json)
pins the archive hash and source equation/table labels. No local manuscript
value is an input to this implementation.

The [helper](../../../scripts/implementation/alda_2026nrd_neutral_inputs.py) reconstructs
the full real three-by-three neutral mass matrix, diagonalizes it, orders the
states by eigenmass and fixes their signs against the published generator
entries. It includes the pion/isospin admixtures. It then applies the source
field-convention map

$$
T_P=\frac{F}{\sqrt2 f_\pi}\lambda_P,\qquad
\chi=2\operatorname{diag}(B_0m_u,B_0m_d,B_0m_s).
$$

This is a conversion of the same canonical fields, not an independent
rescaling of physical states. At $F=f_\pi=0.093$ GeV the reconstructed
$T_P$ matrices are orthonormal; rounded printed generator entries are not
renormalized or treated as exact eigenvectors. The anomalous singlet mass
enters the neutral mass matrix but is absent from the quark spurion.

## Explicit kinematic choice

Call `alda_2026nrd_neutral_inputs(mass_scheme='tree')` or
`alda_2026nrd_neutral_inputs(mass_scheme='physical')`. Each pi0/eta/etap entry supplies
`TP`, `chi`, `F`, `mp`, `tree_mass` and `mass_scheme`.

The tree scheme uses the fitted matrix eigenvalues. The physical scheme uses
the on-shell masses in Alda:2026nrd's arXiv input table while retaining its generators
and quark-spurion inputs. In particular the LO eta eigenmass is about
0.49410 GeV, whereas its physical kinematic mass is 0.547862 GeV.
The latter substitution is an explicitly phenomenological kinematic choice;
it neither changes the eigenvectors nor represents an NLO refit.
No mass-spurion trace is inferred from either eta mass.

## Verification and uncertainty status

Run [the executable verification](../../../scripts/verify_alda_2026nrd_neutral_inputs.py):

```text
.venv/bin/python scripts/verify_alda_2026nrd_neutral_inputs.py
```

It verifies the reconstructed matrices against the arXiv's displayed rounded
diagonals, orthonormality, the eta tree eigenmass and the explicit $F/f_\pi$
convention conversion. Precision beyond the rounded central parameter inputs
is not claimed. Posterior correlations and chiral truncation effects have not
been propagated through this helper. Stage-4 ranks use fixed central inputs.

The stage-3 vector block uses an isospin-averaged spurion and omits rho–isoscalar
mixing, whereas these pseudoscalar/odd-vertex inputs retain Alda:2026nrd's small isospin
breaking. This is a declared partial-isospin approximation, not a fully
isospin-resolved mass/source action. Charged/neutral vector scenarios remain
separate. A complete treatment must revisit the even-vector block and its
calibration rather than silently adding rho mixing to the existing angles.
