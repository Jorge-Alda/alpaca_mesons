---
title: Numerical and phase input status for the isoscalar vector three-pion decay
status: mixed
created: 2026-10-25
updated: 2026-10-25
bibtex_keys: [Alda:2026nrd, Niecknig:2012sj, Hoferichter:2023bjm, Hoferichter:2025lcz]
provenance: literature
---

# Scope

This record checks whether the dispersive $V_i\to3\pi$ solution needed in the
strict ALP pole residue is available as a complex numerical input, with
correlated variations and a phase matched to the ALP production amplitude.
The required amplitude convention and pole residue are defined in the
[vector partial-wave record](alp-three-pion-vector-partial-wave.md).

# Inputs and provenance

- The source convention and required relative phase are from the
  [vector partial-wave record](alp-three-pion-vector-partial-wave.md), using
  the no-mixing convention of Alda:2026nrd [@Alda:2026nrd].
- The Khuri--Treiman decay equation, numerical solution strategy, and
  $\omega,\phi\to3\pi$ fits are in Niecknig:2012sj [@Niecknig:2012sj]. Its
  arXiv source archive is listed in
  [the source index](../../sources/Niecknig:2012sj/INDEX.md).
- The later current-amplitude fits are documented in the indexed records for
  Hoferichter:2023bjm and Hoferichter:2025lcz
  [@Hoferichter:2023bjm; @Hoferichter:2025lcz]. They concern the
  electromagnetic-current matrix element and do not supply an ALP-side
  $V_i\to3\pi$ phase match.

# Result

Niecknig:2012sj supplies the coupled single-variable dispersive equation and
its numerical solution as plots, together with the amplitude convention
$$
\mathcal M(V_i\to3\pi)=i\epsilon_{\mu\nu\alpha\beta}
n^\mu p_+^\nu p_-^\alpha p_0^\beta\mathcal F(s,t,u),
\qquad \mathcal F(s,t,u)=\mathcal F(s)+\mathcal F(t)+\mathcal F(u),
$$
and normalized shape functions such as $\mathcal F_a$ for a unit subtraction
constant [@Niecknig:2012sj, Eqs. (1), (7), (14), Figs. 7--8]. The physical
normalization is fixed from the decay width or a Dalitz fit. The paper
compares Bern and Madrid $\pi\pi$ P-wave phases, once- and twice-subtracted
representations, and reports upper/lower amplitude bands; it does not provide
a machine-readable complex $F_i^{\rm on}(t,u)$ grid or a joint set of
correlated parameter replicas [@Niecknig:2012sj, Secs. 4.4--5.2]. The
archived arXiv source contains the TeX and figures, but no numerical grid
files. The later three-pion current-fit sources likewise do not fill that
gap for the vector decay function [@Hoferichter:2023bjm;
@Hoferichter:2025lcz].

The phase in the displayed decay-amplitude convention fixes the kinematic
tensor convention, while width and Dalitz data determine an overall modulus
and shape. They do not determine the absolute complex phase relative to the
source-resolved ALP production residue. The complex subtraction parameter
used in some Dalitz fits is an internal parameter of that decay solution; it
does not by itself fix the relative phase $\zeta_i$ in the ALP pole
factorization. That phase still requires both amplitudes to be matched in a
common LSZ and S-matrix convention, as stated in the linked partial-wave
record [@Niecknig:2012sj; @Alda:2026nrd].

# Reproducible route and remaining requirements

A numerical decay grid can be reproduced by solving the published
Khuri--Treiman equation with a declared $\pi\pi$ P-wave phase, Omnès
subtraction prescription, and analytic continuation of the crossed-channel
angular average. Its scale can then be fixed using the selected measured
decay width or a Dalitz likelihood. To propagate uncertainty, preserve the
joint variations of phase input, subtraction prescription, fitted complex
subtraction parameters, and width or data normalization; the published
upper/lower curves are sensitivity bounds, not a covariance matrix and not
random replicas. A new numerical implementation must therefore declare its
own input ensemble and cannot label independent quoted marginal errors as
correlated replicas.

The remaining gates for using this solution in the ALP spectral calculation
are: (1) choose the phase-shift and subtraction prescription and provide its
numerical input; (2) obtain the experimental normalization inputs and their
uncertainties, or select the published decay-width normalization; (3)
generate the complex Dalitz-grid solution and documented correlated
variations; and (4) match its global phase to the ALP production residue in
the same external-state and LSZ convention. No ALP production amplitude or
relative phase follows from the $V_i$ decay fit alone.
