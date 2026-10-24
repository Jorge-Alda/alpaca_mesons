---
title: Recommended route to a reliable constraint on the ALP TFF null direction
status: mixed
created: 2026-10-24
updated: 2026-10-24
bibtex_keys: [Alda:2026nrd, Knecht:2001qf, BaBar:2009rrj, BaBar:2018zpn, Adlarson:2016hpp, NA60:2016nad]
provenance: assumption
---

# Purpose and present status

The fit includes three classes of finite-virtuality observables: vector Dalitz
spectra, singly virtual pion TFF data, and doubly virtual $\eta'$ TFF data.
Together they probe one direction in the preliminary four-dimensional null
space, but the current original-eleven-coefficient tree model describes them
poorly when combined. In particular, the existing solution predicts the
BaBar double-tag $\eta'$ TFF too large by factors of roughly 13--42 before
refitting; after refitting, the low-energy and NA60 blocks become much worse.
The result is recorded in the
[three-class fit note](doubly-virtual-etap-joint-fit.md). Therefore the
direction is algebraically visible but not yet reliably measured.

The objective is to constrain the ALP doubly virtual TFF variation with
data through shared odd-parity coefficients, while keeping the physical
meson source matrices used in the data fit distinct from the ALP no-mixing
source. The ALP-specific projection is documented in the
[rho-pole null record](../../alp/lepton-pair/alp-isovector-tff-fit-null-projection.md).

# What a short-distance-compatible TFF means

The resonance TFF must have the QCD leading powers and coefficients in each
hard Euclidean limit. For one large virtuality and one real photon,
$$
F_P(-Q^2,0)=\frac{C_P^{\mathrm{single}}}{Q^2}
+O(Q^{-4})
\quad (Q^2\to\infty),
$$
up to the logarithmic evolution required at the chosen accuracy. For two
large virtualities at fixed positive ratio,
$$
F_P(-Q_1^2,-Q_2^2)=
\frac{C_P^{\mathrm{double}}(Q_1^2/Q_2^2)}{Q_1^2+Q_2^2}
+\text{power-suppressed terms},
$$
with a coefficient derived for the relevant physical meson and factorization
convention. These equations specify leading powers; they are not permission
to impose an arbitrary coefficient or transfer the pion or $\eta'$ coefficient
to an ALP. The pion short-distance input is indexed in
[Knecht:2001qf](../../extracts/Knecht:2001qf/pion-tff-short-distance.md).
The source-specific ALP coefficient still needs its own QCD matching.

The currently retained tree expression must be checked along several paths:
the singly virtual limit $Q_1^2\to\infty$, $Q_2^2=0$; the symmetric
double-virtual path $Q_1^2=Q_2^2\to\infty$; and fixed unequal ratios
$Q_1^2/Q_2^2$. A valid completion must remove any constant or growing term
where QCD predicts inverse-power falloff. It must keep Bose symmetry, both
photon placements, electromagnetic Ward identities, the measured real-photon
normalization, and the low-energy poles/cuts already represented.

The practical prescription is to derive the high-energy terms from the
source-resolved QCD matrix elements or operator-product expansion, then join
them to the low-energy resonance/dispersive representation at a declared
matching scale. State what is included on each side so rho-region strength is
not counted twice. Vary the matching scale and the allowed completion
parameters to estimate the matching uncertainty. A fitted extra pole or
polynomial with no matching derivation is only a diagnostic parameterization,
not a short-distance completion. Evaluate the hard limits for the pion,
physical $\eta'$, and ALP source separately.

# What a better vector Dalitz description means

The A2 and NA60 observable is the finite-timelike-virtuality spectrum for
$\omega\to\pi^0\ell^+\ell^-$. Improving it requires a predictor for the
measured bin observable, not just a different fitted slope. The next model
should:

1. Include the full direct VJP and VVP-to-photon conversion contributions,
   their interference, and the physical flavor/vector responses with the
   same signs and normalization used in the real-photon widths.
2. Replace the fixed zero-width internal $\rho$ pole by a finite-width
   propagator with the appropriate energy dependence and analytic structure
   over the measured dilepton-mass range. Derive or source the line shape and
   vary its documented parameters; do not tune a width to absorb odd-coupling
   mismatches.
3. Fold the prediction over each experimental bin and its phase-space and
   acceptance definition where available. Retain the bin-center result as a
   comparison, and keep A2 and NA60 normalization/systematic nuisances
   separate. Use published covariance if obtainable; otherwise label the
   diagonal approximation and test its effect.
4. Check normalization at $q^2=0$, the threshold behavior, the real-photon
   limit, residuals versus $q^2$, and stability when fitting A2 and NA60
   separately or withholding part of one spectrum. A single integrated
   $\chi^2$ or nonzero null derivative does not establish a satisfactory
   shape model.

These improvements test the Dalitz description itself. If dispersive pion
rescattering or additional channels are added, specify their matching to the
resonance terms and avoid counting the same cut or pole contribution twice.

# Recommended work sequence and decision test

1. **Freeze the present baseline.** Retain the three-class output and its
   block residuals as a reference. Do not use its small local coefficient
   errors as calibrated uncertainties.
2. **Complete the meson TFF high-energy behavior.** Derive and verify the
   single- and double-virtual asymptotic coefficients, build one matched
   predictor for all TFF rows, and report each data block's residuals and
   sensitivity to matching choices. Keep the current 11-parameter fit and
   any extension as separately labelled models.
3. **Improve the Dalitz predictor.** Implement the finite-width, bin-folded
   description and establish whether it describes both spectra with the
   same coefficients and nuisance treatment.
4. **Refit jointly only after the component descriptions pass their own
   checks.** Profile the even-sector, line-shape, normalization and matching
   nuisances supported by data. Report coefficient rank, singular values,
   profile intervals, and per-block pulls; test changes under reasonable
   fit-range and covariance assumptions.
5. **Project the accepted coefficient region onto the ALP target.** Use the
   ALP no-mixing sources to evaluate the spread in
   $\delta F_{a\gamma^*\gamma^*}(s_1,s_2)$ over the kinematics relevant to
   the lepton loop. A useful constraint requires both acceptable data
   description and a stable, finite target profile after nuisance
   propagation. Report any remaining null-sensitive component rather than
   choosing its zero representative.

The neutral meson TFF classes are expected to probe only the same one
old-null component found so far. The other three preliminary null directions
remain unconstrained by these rows and require observables with different
operator projections or independent theory input. This is a fit-planning
recommendation, not a claim that the required improved predictors or a
successful joint calibration have already been obtained.

# Inputs and provenance

- The present fit and its conditional failures are in
  [doubly-virtual-etap-joint-fit.md](doubly-virtual-etap-joint-fit.md).
- The original retained neutral TFF and its hard-limit caveat are in
  [the source-resolved TFF note](../../alp/lepton-pair/alp-doubly-virtual-tff.md).
- The ALP old-null projection is in
  [the rho-pole TFF note](../../alp/lepton-pair/alp-isovector-tff-fit-null-projection.md).
- The vector Dalitz baseline is in
  [the conditional A2/NA60 fit](vector-dalitz-conditional-fit.md).
- The experimental TFF blocks are indexed in the [CLEO extract](../../extracts/CLEO:1997fho/tff-tables-and-errors.md),
  [BaBar single-virtual extract](../../extracts/BaBar:2009rrj/pi0-tff-table.md),
  and [BaBar double-virtual extract](../../extracts/BaBar:2018zpn/doubly-virtual-eta-prime-tff.md).
