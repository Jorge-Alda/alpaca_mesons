---
title: First physics draft of the odd-parity vector and ALP sector
status: mixed
created: 2026-10-14
updated: 2026-10-15
bibtex_keys: [Alda:2026nrd, Ecker:1988te, Ecker:1989yg, Ruiz-Femenia:2003jdx, Graf:2020yxt, Kaiser:2000gs, Bauer:2020jbp]
provenance: assumption
---

# Scope and provenance

The user requested a first physics-only draft covering odd-intrinsic-parity
vector mesons and ALPs, operators, Hilbert series, assumptions and limitations.
The draft is integrated into the existing [main.tex](../../../latex/main.tex)
through [vectors.tex](../../../latex/vectors.tex), following the user's explicit
manuscript location and modular structure. Its [README](../../../latex/README.md) identifies the
source records and executable derivations supporting its equations. The draft
reuses established repository results and introduces no new calculation.
The existing title, author, introduction and pseudoscalar material are preserved.

# Scientific selection

The scientific organization was reviewed by the brainstormer agent before
drafting: conventions and anomaly, the original seven VJP/four VVP operators,
the bounded covariant extension, its finite Hilbert polynomial, source/LSZ
assembly, and power-counting limitations. Draft prose was delegated to the
drafter agent. Sources already have inspected repository bibliography entries.

Inputs are [RChT foundations](../foundations/rcht-basics.md),
[the original operator extract](../../extracts/Ruiz-Femenia:2003jdx/vector-odd-parity.md),
[the bounded action count](odd-higher-covariant-operator-count.md),
[the Hilbert polynomial](odd-bounded-hilbert-series.md),
[the parent filtration](odd-hilbert-power-counting.md),
[hard derivative counting](odd-resonance-derivative-counting.md),
[source/LSZ assembly](odd-alp-source-lsz-assembly.md), and
[operator relevance](odd-operator-alp-relevance.md).

Assumptions and derived results remain distinguished. The count of 228 applies
only to the declared seed quotient; it is not an unrestricted U(3) EFT count.
Raw character multiplicities are independently derived, while the quotient
inherits exact relation certificates. Formal parent counts do not establish
observable ranks or controlled convergence. Independent-source differentiation
precedes physical source restriction, and fixed topological contacts remain.
The photon-operator normalization follows Bauer's Eqs. (39)/(54); it is not
attributed to Alda's WZW-excluding action. The meson contribution is a separate
leading-nonet WZW/source derivation. No fit implementation is included.

# Verification

The integrated draft builds with `latexmk -pdf -interaction=nonstopmode
-halt-on-error main.tex` from `latex/`, using the existing local `utphys` style.
The six-page PDF after the 2026-10-15 narrative revision has no final LaTeX
warnings, unresolved citations/references or overfull boxes. The earlier
Hilbert/table and ALP-source pages and the revised prose/comment page were
rendered and visually inspected. Both repository and manuscript
`git diff --check` pass.
These checks validate presentation, not a fresh scientific derivation. The
corresponding scripts are linked in the manuscript README and were not rerun
solely to prepare this text. The obsolete separate draft was removed after the
integrated manuscript compiled successfully.

# Manuscript conventions

Following the user's correction, operator families are written directly with
$\chi_+$ and $i\chi_-$ rather than auxiliary $X,Y$ definitions. The ALP source
section uses the original $\mathbb M_a^{-1}(C)$, $c_q^A$, $c_q^P$, $m_q$ and
$m_a^2$ notation without new reduced source-weight symbols. Flavour traces use
angular brackets. The extension is motivated through independent sources,
singlet and OZI structures, flavour breaking and allowed derivative corrections;
its bounded completeness and short-distance limitations remain explicit.

The subsequent convention audit removed the $\Phi_R,F$ comparison from the
manuscript. Every displayed interaction uses the declared metric $(+---)$ and
$\epsilon^{0123}=+1$, hence $\epsilon_{0123}=-1$, directly, following
[@Bauer:2020jbp, immediately after Eq. (1)]. The photon sign audit uses the derivative form

$$
\mathcal L_{a\gamma\gamma}^{\mathrm{contact}}
=\frac{e^2}{8\pi^2f_a}C_{a\gamma\gamma}a\,
\epsilon_{\mu\nu\rho\sigma}\partial^\mu A^\nu\partial^\rho A^\sigma.
$$

Its sign was audited by the coder against the existing
[signed phase derivation](../../alp/photon/alp-photon-odd-action-phase-matching.md)
and the [canonical conversion transcript](../../../data/vector_sector/levi_civita_convention.md).
The recorded derivative-tensor contraction has two-photon factor $-2$ in the
canonical tensor, yielding the negative scalar contact coefficient
$-C_{a\gamma\gamma}/(4\pi^2f_a)$. The manuscript's positive dual coefficient
and its odd couplings preserve the source action under the conversion in
[the convention audit](levi-civita-convention.md). The old reference-backend
transcript instead used the opposite tensor and a positive scalar coefficient
$C_{a\gamma\gamma}/(4\pi^2f_a)$. This is an application of an existing
executable result together with the executed conversion, not a new anomaly matching.
The photon-operator normalization is supplied by [@Bauer:2020jbp, Eqs. (39)/(54)];
independent local/anomaly matching of the enlarged RChT action remains open.
Angular brackets and bold
matrix notation are applied consistently in the manuscript. The brainstormer
review also scoped the $d_4$ canonical two-boson zero and retained only the
declared single-odd tree $a\to3P$ zero without suggesting that loops alone
remove its external-momentum constraint.

# Planning

The user's prose preference is to name the physical object or operation rather
than use vague process vocabulary. In the manuscript, neutral-meson
contributions, the solution of the mass-operator equation and the field
redefinition that removes bilinear ALP--meson mixing replace loose uses of
``profile''. Assumptions and missing calculations are stated explicitly rather
than summarized as ``conditional''. The Hilbert discussion names the covariant
reduction and exact three-flavour trace calculation that establish the quotient,
without implying that the raw character count independently proves it.

This drafting milestone records the currently established physics, rather than
changing the ongoing calibration plan. The 2026-10-15 narrative revision was
delegated to the drafter and reviewed with the brainstormer. Its physics argument
now proceeds from the independent-source problem to the enlarged operator space,
the useful nested truncations, ALP source insertion, and the physical scope of the
resulting amplitudes. The basis and ALP subsections end by stating what they
establish and which next subsection addresses the remaining issue.

The user's requested contiguous `\CGPT{...}` comments distinguish documentation
of completed calculations in a future appendix from actual missing work. The
[neutral cubic maps](odd-power-counted-channel-maps.md), LO source/LSZ assembly,
Hilbert quotient and parent filtration are established calculations, not marked
as missing derivations. The [finite neutral target audit](odd-alp-measured-null-audit.md)
already includes canonical vector factors and the signed photon response.
Independent hadronic local/anomaly matching, charged and higher-field ALP
targets, full physical-domain/width completion and quantitative convergence
remain separate unfinished goals. The next writing step is to collect the
existing explicit representatives, reduction certificates and kernel/LSZ
details in the future appendix while resolving those physics gaps before
claiming complete decay predictors. The overarching
physical-predictor completion gates remain in the existing
[next-steps map](next-steps.md). Numerical ALP predictions or uncertainty bands
must wait for those gates and cannot be inferred from the draft's operator
counts. Fit implementation remains deferred by the user's instruction.
