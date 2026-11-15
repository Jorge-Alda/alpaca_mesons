---
title: RunDec package for running and decoupling of alpha_s and quark masses
status: literature
created: 2026-11-15
updated: 2026-11-15
bibtex_keys: [Chetyrkin:2000yt]
provenance: literature
---

# Source and scope

Read from the abstract and procedure list of `memory/sources/Chetyrkin:2000yt/src/rundec.tex`
(arXiv hep-ph/0004189) [@Chetyrkin:2000yt]. The formulas were not extracted individually.

# Statements

- The paper collects formulas for the strong coupling and quark masses at different
  scales and numbers of active flavors with QCD corrections to three and
  sometimes four loops, and provides a Mathematica package.
- The package includes decoupling routines (`DecAsUpMS`, `DecAsDownMS`, and mass analogs;
  usage examples at the lines near 2090-2130 of the TeX file) for stepping the
  number of active flavors up or down in $\overline{\rm MS}$.
- Use: reference implementation for numerical cross-checks of the coupling
  matching in [Chetyrkin:1997sg](../Chetyrkin:1997sg/alphas-threshold-matching.md) and
  [Chetyrkin:1997un](../Chetyrkin:1997un/alphas-decoupling-constant.md); no
  axial-current or TFF content.

# Flags

- No numerical values are recorded here. Any numerical use requires extracting the
  explicit equations first.
