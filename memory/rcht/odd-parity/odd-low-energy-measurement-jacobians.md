---
title: Fixed-scenario low-energy measured-width Jacobians and conditional ranks
status: mixed
created: 2026-10-14
updated: 2026-10-14
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx, ParticleDataGroup:2024cfk]
provenance: numerical-computation
---

# Result and evidence status

All 64 retained mass/source scenarios have the same **conditional local
measurement ranks**, after including the original even-sector and shared
total-width observations once:

| Parent filtration | Leading | First corrections | Full bounded space |
|---|---:|---:|---:|
| All additional derivatives soft | 7 | 12 | 12 |
| Hard-vector derivative parents promoted | 12 | 12 | 12 |

The original eleven-coupling subspace has rank seven. The parent-action
dimensions for the first row are $(12,101,228)$; for the promoted row they
are $(35,136,228)$. Thus the full bounded action has 216 locally unmeasured
action directions in this snapshot, beyond its raw presentation redundancies.
These are fixed-scenario measured-width Jacobians, distinct from the earlier
arbitrary-source family and stripped-vertex ranks.

The numerical information is conditional on a diagonal covariance for the
PDG summaries. Their complete covariance is unavailable; some summaries
share underlying measurements. Broad-rho NWA and constrained-fit K-star
entries must not be promoted into precision independent data by this rank
calculation. Three separate controls quantify their role:

| Measurement control | Rows | Soft ranks | Promoted ranks |
|---|---:|---|---|
| Exclude broad-rho summaries | 9 | $(7,9,9)$ | $(9,9,9)$ |
| Exclude constrained-fit K-star summaries | 10 | $(7,10,10)$ | $(10,10,10)$ |
| Exclude both sets | 7 | $(6,7,7)$ | $(7,7,7)$ |

These control ranks hold in every retained scenario. Removing rows does not
supply their missing covariance or establish model adequacy.

# Inputs and assumptions

Use the source-backed [PDG 2024 snapshot](../../extracts/ParticleDataGroup:2024cfk/odd-parity-low-energy-measurements.md),
its [data JSON](../../../data/odd_parity/odd_low_energy_measurements.json),
the [bounded numerical kernels](odd-bounded-sm-numerical-kernels.md),
[Alda input conversion](alda-neutral-amplitude-inputs.md),
[stage-3 fixed sectors](../foundations/vector-stage3-calibration.md),
and the [selected calibration plan](odd-breaking-low-energy-calibration-plan.md).
Sources are already catalogued in `bibliography.bib`; the snapshot pins
seven official PDG source PDF hashes and printed pages.

The twelve rows are nine radiative branching fractions and three two-photon
partial widths. The pion uses PDG's cited direct PrimEx partial width with
the two quoted errors combined in quadrature; its mean-life average is not
added as a second datum. Seven shared resonance total widths are represented
once. The inclusive $\eta'\to\rho^0\gamma$ entry is deferred because the
current predictor lacks its nonresonant contribution. Hard CLEO/BaBar,
near-origin slopes without an accepted snapshot, and three-pion data are
excluded. Scale factors are not applied to quoted errors again.

Natural units are GeV, $F=f_\pi=0.093$ GeV, operator scale $M_V=0.770$ GeV,
$\alpha$ and the electron mass are the fixed stage-3 inputs. Both mass
branches, four charged/neutral input choices, two phi electronic-width
summaries and four signed source hypotheses remain separate: 64 scenarios.
The near-ideal preference is a tie-breaker, not a likelihood term.

The LO neutral eigenmasses and charged masses
$m_{P,ij}^2=(\chi_i+\chi_j)/2$ enter the local vertex invariants. This is
necessary for the same bounded source-preserving Goldstone-EOM quotient.
Physical Alda input-table masses enter the phase-space factors separately
as a stated phenomenological kinematic replacement. It is not an NLO
source/mass refit. A physical-mass replacement inside the local polynomial
would instead require its nonzero EOM remainders; it cannot silently be
called the same 228-dimensional action map.

The even mass/source block retains exact isospin, while odd spurions and
neutral generators retain Alda isospin breaking. This partial-isospin
truncation and broad-vector NWA limits remain. Alda central inputs are fixed:
posterior covariance and chiral/model truncation effects are not propagated.
Both charged/neutral rho vertices share the selected rho mass in each
scenario, and both K-star vertices share its selected K-star mass. Their
measured total-width nuisances remain separate rather than being averaged.

# Calculation and reproduction

Run

```text
.venv/bin/python scripts/verify_odd_bounded_amplitudes.py
.venv/bin/python scripts/derive_odd_low_energy_jacobians.py
```

The numerical implementation evaluates all 241 normalized raw parent
gradients, including charged $f_-$ terms, their composite derivative,
the signed full spectral response and both photon placements. The exact
13 inherited action relations determine the parent quotient image of each
retained-column set. Filters act on parent columns before amplitude reduction.

Write a real tree-level channel amplitude as
$G_i=G_{0i}+A_{ij}g_j$, with the inherited WZW offset only for the neutral
two-photon form factor. Its phase-space factor is $k_i$, including the
electromagnetic convention appropriate to that amplitude. The executable
SymPy derivation differentiates

$$
y_i=\frac{k_iG_i^2}{\Gamma_i},\qquad
\frac{\partial y_i}{\partial g_j}
=\frac{2k_iG_i}{\Gamma_i}A_{ij},\qquad
\frac{\partial y_i}{\partial\Gamma_i}
=-\frac{k_iG_i^2}{\Gamma_i^2}.
$$

For a direct partial width omit $\Gamma_i$. Each radiative gradient includes
its explicit $e$; the two-photon scalar strips $e^2$ and uses the established
identical-photon factor in $k_i$. Amplitudes have dimension GeV$^{-1}$ and
$k_i$ has dimension GeV$^3$.

The common reference point is the explicit deterministic original-coupling
seed recorded in the artifact, with new couplings zero. It satisfies no OPE
constraint and is neither a fit nor a coefficient prior. All reference
amplitudes are nonzero. Consequently each measured-width row differs from
its corresponding amplitude row by a nonzero scalar, so their local ranks
agree. At zero odd couplings the radiative-width Jacobian vanishes: that
singular point does not establish a global blind direction. Width-only data
also leave discrete amplitude-sign ambiguities that a local Jacobian does
not resolve.

The runner whitens with the stated summary errors. It constructs nuisance
columns for the declared mass-action coordinates, the permitted source
coefficients and standardized shared total widths. Their original four
mass, three electronic-width and seven total-width rows are appended once;
stage-3 fitted parameter errors are not added as priors. Mass matrices come
from the existing executable SymPy action derivation. Nonlinear even-sector
response derivatives use deterministic central differences at two step sizes.

The nuisance constraint block has full local column rank in every scenario.
Thus profiling those *locally constrained* nuisance columns preserves the
odd rank. This does not mean nuisance uncertainty is zero: the projected
singular values still encode its effect. The marginal rank is computed as
$\operatorname{rank}[J_{\rm nuisance},J_{\rm odd}]
-\operatorname{rank}J_{\rm nuisance}$ and the artifact also exports the
orthogonal nuisance projection.

Numerical rank uses nonzero column normalization and relative SVD tolerance
$10^{-10}$, retaining measurement whitening. Unnormalized whitened singular
values are also saved. Their scale depends on the declared raw coefficient
basis and arbitrary reference point; no fitted covariance or naturalness
information is inferred. Pivot-column reconstructions export the measured
combination maps. The raw null generators are defined explicitly by those
maps; inherited action-relation directions must then be divided out to
obtain the physical unmeasured action space.

# Verification and implications

The [JSON artifact](../../../data/vector_sector/odd_low_energy_jacobians.json)
contains each separate measured and nuisance Jacobian, retained parent IDs,
combination maps, whitened/projection singular values, all controls, the
reference point, input hashes and limitations.
The [SymPy transcript](../../../data/vector_sector/odd_low_energy_jacobians.md)
records the observable derivatives and nested rank results.

All 64 runs recover the held signed spectral sectors. Every local gradient
annihilates the 13 action relations; the largest normalized residual is
$2.49\times10^{-16}$. Width derivatives agree with independent simultaneous
finite displacements of the eleven original coefficients. The largest
relative change on halving the nonlinear nuisance step is
$6.30\times10^{-8}$. The independent kernel verification additionally
passes 42 legacy-amplitude checks, 723 charged raw-column checks and 4338
photon-parent Bose checks. All retained-column reconstructions pass.

At soft first corrections the twelve summary observables already span
twelve independent local directions. The full action adds no measured
direction here. Hard derivative promotion reaches that rank already at
its leading level. These saturated models can interpolate this snapshot
locally; an improved in-sample residual alone cannot validate their power
counting or numerical hierarchy. The narrow/nonconstrained control likewise
saturates at seven directions after soft first corrections.

The subsequent [conditional calibration](odd-low-energy-calibration.md)
now executes 1792 joint scenario/model/control fits with local observable
uncertainties and held-out null tests for the unsaturated models. No
preferred truncation or ALP uncertainty band is established. Its next
gate is sign/nuisance profiling and target-null projection, retaining
all controls and source alternatives. Do not fit all
217 new action coefficients, conceal null directions with priors, or assume
the 216 unmeasured directions are irrelevant for ALP channels. Their target
projection requires the separate ALP sensitivity audit before propagation.
