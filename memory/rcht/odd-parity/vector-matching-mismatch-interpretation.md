---
title: Interpretation and possible remedies for the two vector matching mismatches
status: mixed
created: 2026-10-13
updated: 2026-10-13
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx, Roig:2013baa, Knecht:2001qf]
provenance: assumption
---

# Scope of the explanation

Saved at the user's request on 2026-10-13. The results below are established
in executable records; their interpretation and proposed remedies are
explicitly distinguished from demonstrated completions. This is not a new
calculation or a physical fit. Conventions are tree-level tensor RChT,
natural units and the positive-anomaly/positive-strange phi phases.

## First mismatch: different hard limits of a restricted action

The [chiral-octet correlator calculation](vector-chiral-octet-matching.md)
takes all three invariants hard, including the pseudoscalar-density momentum.
The on-shell single-virtual TFF instead keeps one photon real and the
pseudoscalar on shell; the symmetric double-virtual TFF keeps the
pseudoscalar on shell and makes both photons hard. They are independent
kinematic constraints on the same assumed resonance action.

In its massless common-pole control, correlator matching requires

$$
d_3=-\frac{N_CM_V^2}{64\pi^2F_V^2}+\frac{F^2}{8F_V^2}.
$$

Together with the other matched combinations it leaves
$\lim F_{\pi^0}(-Q^2,0)=F/(3M_V^2)$, violating single-virtual
$1/Q^2$ falloff. Changing $d_3$ to remove the $F^2/(8F_V^2)$ term,
while enforcing both TFF constant cancellations, instead gives
$\lim Q^2F_{\pi^0}(-Q^2,-Q^2)=0$, whereas the symmetric pion OPE
requires $2F/3$ in this phase convention. Thus vanishing alone is not
sufficient matching. These are exact results for the current common-pole
kernel, not a theorem excluding every theory with a single vector multiplet.

The literature likewise distinguishes vector-only correlator relations
from the conditions after adding pseudoscalar resonances
[@Roig:2013baa; see the
[extract](../../extracts/Roig:2013baa/tff-correlator-compatibility.md)].
The nonzero symmetric coefficient and distinction between LMD and LMD+V
are documented in [the pion target extract](../../extracts/Knecht:2001qf/pion-tff-short-distance.md)
[@Knecht:2001qf].

## Second mismatch: source breaking without interaction completion

The [physical-channel matching record](vector-physical-tff-matching.md)
defines $W,S,R_i,T_i,q$ explicitly. The response has
$W(-Q^2)=S/Q^2+O(Q^{-4})$, so its leading hard flavor structure is
the source moment $S$, independent of pole masses. The WZW term instead
has the electromagnetic flavor structure $q(Q,Q)$.

With the present universal odd coefficients, cancellation of all three
physical neutral-channel double-virtual constants, after removing
single-virtual growth, requires $S=\Lambda Q$. The declared source action
has $S=F_VQ+z_V\{X,Q\}+y_V\operatorname{tr}(XQ)I$; for
$x_l\ne x_s$ the proportionality equations force $z_V=y_V=0$.
The retained source-breaking models therefore cannot meet these necessary
conditions without changing the action beyond masses and source residues.

All 48 fixed-central source-breaking scenarios fail those cancellations.
The 16 universal controls pass the necessary zero-coefficient equations
but do not establish leading-coefficient matching and already fail the
held leptonic adequacy diagnostic. These are separate model/input hypotheses,
not 64 independent experimental tests. The reported source residual is
an algebraic flavor mismatch, not a decay-width systematic or chi-squared.

## Plausible reasons and their evidence status

- **Incomplete odd interaction breaking:** the even action was enlarged while
  odd vertices were held at the common single-trace form. This restriction is
  explicit, and the missing flavor freedom is established by the $S\propto Q$
  condition. A particular sufficient operator completion has not been derived.
- **Insufficient resonance/momentum structure:** the common-pole control
  cannot satisfy the displayed hard requirements simultaneously. Additional
  structures are a possible remedy, but their sufficiency for the full
  finite-spurion nonet has not been shown.
- **Extrapolating a finite low-energy truncation to arbitrarily hard momenta:**
  higher-order contributions need not remain small in that extrapolation.
  Failure of exact asymptotic cancellation does not quantify the error of a
  low-energy decay calculation. This is an EFT interpretation, not a numerical
  estimate of the omitted terms.

Changing a mixing branch or pole masses alone cannot remove the leading
source obstruction, since those masses enter at the next spectral moment.
The exact flavor trace argument also cannot be evaded by relabeling
pseudoscalar eigenvectors within the same complete diagonal basis. The
existing Ward, basis, common-limit and direct-kernel expansion checks support
an algebraic origin; they do not test every omitted physical effect.

## Remedies and limitations

| Route | Intended role | Unestablished work |
|---|---|---|
| Odd flavor-breaking, singlet and single/double-trace interactions | Supply additional channel-dependent flavor tensors and calibrate low-energy deviations | Independent operator reduction, normalized amplitudes, identifiability and adequacy |
| Pseudoscalar resonances | Modify full correlator matching and its compatibility with TFF falloff | Implement their source/interactions and verify leading coefficients in this action; literature compatibility conditions cannot be copied without the extra terms |
| Additional vector multiplet | Add poles and numerator freedom; heavier states may participate in ultraviolet flavor cancellations | Derive and match the extended action, rather than inserting a second pole phenomenologically |
| Restricted low-energy calibration | Test radiative decays and real-photon normalizations/slopes without imposing the incompatible hard prescription | Define the dataset/domain beforehand, propagate input errors and report remaining unidentified combinations |

Roig:2013baa provides a pseudoscalar-resonance compatibility example;
Knecht:2001qf discusses the two-vector LMD+V ansatz. Neither establishes
a complete solution for our finite-spurion $U(3)$ action. Adding mass-type
double traces alone must not be advertised as canceling double-virtual
constants: their numerator has no hard derivative, so two internal-vector
responses make them subleading in that limit.

For the user-selected tree-level ALP accuracy, these results do not by
themselves require loops or momentum-dependent mixing. On 2026-10-13 the
user chose the [combined odd-extension/low-energy plan](odd-breaking-low-energy-calibration-plan.md).
The hard-TFF/complete-OPE objective is retained as a distinct later extension,
not a prerequisite for a demonstrably adequate low-energy calibration.
