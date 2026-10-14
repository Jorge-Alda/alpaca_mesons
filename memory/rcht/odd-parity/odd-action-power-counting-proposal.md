---
title: Proposed momentum, quark-mass and large-Nc hierarchy for the odd vector action
status: mixed
created: 2026-10-14
updated: 2026-10-14
bibtex_keys: [Kaiser:2000gs, Alda:2026nrd, Ruiz-Femenia:2003jdx]
provenance: assumption
---

# Scope and evidence status

Saved at the user's request after discussion of restricting the amplitude
combinations. This is a proposed truncation and its limitations, not an
implemented reduction or a demonstrated uniform power counting for all ALP
kinematics. Saving the proposal does not establish its numerical adequacy.

Inputs: the [bounded action count](odd-higher-covariant-operator-count.md),
[higher two-boson combinations](odd-higher-two-boson-combinations.md),
[operator relevance and mismatch assessment](odd-operator-alp-relevance.md),
[nonet/OZI discussion](u3-nonet-odd-parity-extension-options.md), and
[low-energy calibration plan](odd-breaking-low-energy-calibration-plan.md).
Conventions remain those of Alda:2026nrd, canonical tensor vectors, the fixed
derivative ceiling and independent ALP derivative/mass sources. No vector EOM
is introduced. No independent topological coefficient functions are added.

Kaiser--Leutwyler supplies the systematic simultaneous momentum, quark-mass
and large-$N_c$ expansion in Goldstone-sector $U(3)$ ChPT [@Kaiser:2000gs;
https://arxiv.org/html/hep-ph/0007101]. Applying this organization to explicit
resonances below is a model assumption, not a result established by that
Goldstone-sector construction. This record preserves a proposal rather than
performing a new analytical calculation.

# Proposed bookkeeping

In natural units define the dimensionless parameters

$$
\epsilon_p=\frac{p_{\rm soft}^2}{\Lambda^2},\qquad
\epsilon_\chi=\frac{B_0m_q}{\Lambda^2},\qquad
\epsilon_N=\frac1{N_c}.
$$

Here $p_{\rm soft}$ denotes soft momenta or virtuality variations in the stated
prediction domain. The quark-mass symbol denotes the relevant flavour
insertion, not a replacement by a physical pseudoscalar mass. $\Lambda$ is
the declared breakdown scale of the chosen resonance model; no numerical
value or automatic identification with $M_V$ has been selected. The existing
powers of $M_V$ in operator normalization do not by themselves establish a
breakdown scale or a numerical hierarchy of coefficients.

Assign the three parameters the same nominal correction order only as an
explicit truncation assumption. They need not be numerically equal at
$N_c=3$. Count additional factors relative to the lowest-order operator in
its VJP or VVP sector; leading mass/source structures are retained rather
than discarded merely because they contain $\chi_-$. Canonical field/source
normalizations must be included before comparing physical amplitudes.

| Level | Retained structures |
|---|---|
| Leading | Lowest-order single-trace derivative and mass/source operators in the declared space |
| First corrections | Leading double traces; single-trace $\chi_+$ dressings; single-trace two-derivative $\chi_-$ extensions |
| Deferred mixed corrections | Double-trace $\chi_+$ dressings and double-trace two-derivative mass/source extensions |

Under the equal-order bookkeeping, an added $\chi_+$ insertion plus an extra
trace carries $\epsilon_\chi\epsilon_N$; two additional soft derivatives plus
an extra trace carry $\epsilon_p\epsilon_N$. Such terms lie beyond the proposed
first-correction truncation. Other retained seed blocks, including charged
$f_-$ structures, must be classified consistently at the operator level.
This table is not a claim that the original eleven exhaust every leading
covariant or independent-source structure.

# Limitations and required safeguards

**Mass dimension.** Increasing dimension by two within a fixed VJP or VVP
sector suggests an additional inverse-scale factor. It does not establish
numerical suppression without a momentum/mass expansion and coefficient
assumptions. Dimension-six VVP and dimension-seven VJP operators cannot be
ranked against each other solely by their engineering dimensions.

**Resonance momenta.** Derivatives acting on an on-shell vector probe momenta
of order $M_V$, and are not automatically soft. Keep resonance propagators
explicit rather than expanding across their poles. Separate soft virtuality
variations from resonance-mass factors before assigning derivative order.
For $a\to VV$, ordinary soft-momentum counting can fail; a resonance-region
prescription or a restricted domain is needed before the proposed hierarchy
is advertised as controlled there. The same caution applies to radiative
calibration involving on-shell vectors.

**Large-$N_c$ and singlets.** One additional flavour trace generically has
nominal relative $1/N_c$ OZI suppression, separately from momentum/spurion
counting [@Kaiser:2000gs]. This does not prove uniform channel-by-channel
numerical suppression. A nominally suppressed structure may be the first
nonzero contribution in a channel. Singlet fields themselves are retained;
trace counting does not impose universal singlet/octet couplings or remove
the leading anomaly/topological mass sector.

**Tree accuracy.** Keeping subleading tree operators while deferring loops
and self-energies remains the selected phenomenological tree truncation. It
must not be called a complete next-order calculation in a resonance
$1/N_c$ expansion. No loop/self-energy computation is required merely to
save or explore this hierarchy; the user's tree-level accuracy choice stands.

**Sources and contacts.** Keep independent derivative and pseudoscalar ALP
weights, and all contacts tied to retained parent operators. Do not delete
source-sensitive terms by applying an on-shell meson-mass substitution before
source differentiation.

**No implicit priors.** Power counting organizes omissions and sensitivity
studies. It does not set retained coefficients to fixed values, equate them,
or introduce a Gaussian naturalness/OZI penalty into the likelihood.

# Relevance to the matching objectives

The proposed first-correction truncation retains single-trace flavour-breaking
derivative structures that could address the source-flavour mismatch, plus
leading double traces for independent singlet sensitivity. Their sufficiency
is unestablished; omitted mixed flavour-breaking/OZI corrections may have to
be promoted if the enlarged matching equations or low-energy adequacy require
them. Failure would be reported as a limitation of the truncation.

This hierarchy does not promise to repair the chiral common-pole
correlator/TFF mismatch. Quark-spurion dressings vanish in the strict chiral
limit and cannot be used to claim that repair. Low-energy counting is not a
justification for extrapolating the retained action to arbitrarily hard
momenta or forcing exact asymptotic constraints.

# Implementation and outstanding count

Apply the hierarchy to the parent operator columns before projecting or
reducing amplitudes. The saved independent combinations can mix single- and
double-trace parents; their representative labels are not valid filters for
an entire combination. Preserve the original full basis and compute the
retained-column image and its exact rank separately.

The current 51 higher neutral combinations are an intermediate stripped
vertex rank. No power-counting-filtered count has yet been computed. Next
classify retained parents by sector, additional spurion/derivative order and
trace topology, resolve the treatment of derivatives on resonances, then
recompute the joint leading/higher physical-channel maps. Validate the
hierarchy by comparing nested truncations, correction sizes and channel
adequacy within a declared domain before propagating an ALP uncertainty.
