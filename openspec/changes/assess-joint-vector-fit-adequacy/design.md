# Design

## Context

See `proposal.md`. Four accepted minima exist (`joint_vector_endpoint_comparison.json`), with common-order $\chi^2\approx225.73$ for 75 residuals and 32 coordinates. The previous non-converged snapshot already showed the largest contributions from the shared total widths (7 rows) and the 45 BaBar $KK\pi$ bins. The baseline assumptions (negative single-trace source, near-ideal inverse-mass branch, fixed excited poles with 2.5 GeV regulator, no $K_2^*$, diagonal PDG covariance apart from the $\phi\to K\bar K$ correlation $-0.72$) are recorded in `babar-kkpi-global-vector-sector-fit-baseline.md` and `vector-stage5-joint-fit-specification.md`.

## Goals / Non-Goals

**Goals:**
- Diagnose where the $\chi^2$ comes from and whether a justified variation removes it.
- Keep every variation reproducible and switchable, with exact baseline recovery.

**Non-Goals:**
- New operators, loops, or new data sets. These belong to the follow-up issue if the verdict requires them.
- Profile intervals and identifiability (alp-uzh.7).
- Inflating errors to reach $\chi^2/\mathrm{dof}\approx1$ without a physical or experimental justification.

## Decisions

- **Assess all four accepted minima, refit from the best.** Goodness-of-fit numbers are reported for all four; block exclusions and variations are warm-started from `basin_probes/even_shared_widths`. Because the four minima differ by $<10^{-3}$ in $\chi^2$, multi-starting each variation is unnecessary; a variation whose refit moves far (scaled distance above 1) is re-run from the other three accepted endpoints to check for a basin switch.
- **Quadrature.** Refits use $(2,16,16)$ with four-bin blocks and `max_nfev=5000`; each final endpoint is re-evaluated at $(4,128,128)$, and only the high-order $\chi^2$ is reported. Variations failing the certified tolerances are reported as not qualified.
- **Variations as explicit likelihood options.** Each variation is an option of the likelihood class with a default that reproduces the baseline; a regression test checks bitwise-equal residuals at the default. Alternatives such as a separate likelihood copy per variation were rejected to avoid divergence.
- **Variation list (fixed before running):** BaBar upper cut 1.70, 1.60 GeV; excited-pole masses and widths $\pm1\sigma$ (PDG), regulator 2.0 and 3.0 GeV; tensor products fixed to zero; source hypotheses universal, single-trace positive, single-plus-double-trace; second inverse-mass branch; $\phi\to K\bar K$ correlation set to zero; alternative $\phi$ electronic-width summary; 5% fractional systematic on BaBar bins added in quadrature (assumption, labelled). Window and covariance variations change the degrees of freedom; it is recorded per variation.
- **Deferred variation: single-plus-double-trace source (decided 2026-11-07 by the user).** The double-trace source term $y\,\mathrm{tr}(\chi_+Q)\mathbb{1}$ has not been derived for the inclusive $KK\pi$ current, and the likelihood has no $y$ coordinate. Deriving it would add a new operator contribution, which is a non-goal of this change. It is therefore not run here. It is recorded as a candidate model extension in the follow-up bead, and the adequacy verdict is stated as conditional on its absence.
- **Excited-pole and tensor-pole variations (decided 2026-11-07 after consulting the brainstormer).** Each of the eight fixed excited poles and the $a_2(1320)$ has its mass and its width shifted on its own by $+1\sigma$ and by $-1\sigma$, giving 36 refits. The uncertainties are the PDG 2024 values in `memory/rcht/odd-parity/pdg2024-excited-vector-and-tensor-poles.md`. The $a_2$ width uses the OUR ESTIMATE uncertainty of 5 MeV. The asymmetric $\phi(2170)$ width uncertainty is applied as quoted, without symmetrising. The shifts are applied around the hard-coded baseline values, which stay unchanged. A linear $\Delta\chi^2$ screen at the best endpoint is used only to order the refits and to cross-check the refit results, never to drop any of them. These refits show how much the fit depends on the fixed pole inputs; they do not replace the baseline in the verdict.
- **Adequacy criterion.** Nominal $p>0.01$ for the total, with no single block having $p<10^{-3}$. This threshold is an assumption and stated as such; the nominal $p$-value is not corrected for the bounded and nonlinear parameters.
- **Verdict ownership.** The assessment script outputs numbers; the verdict and the proposed physics extension are written in the memory record, and the follow-up bead is created with `bd create --parent alp-uzh`.

## Risks / Trade-offs

- [A variation needs physics not yet implemented, e.g. a $K_2^*$ contribution] → Not added here; recorded as a candidate extension for the follow-up issue.
- [Block exclusion of shared-width rows leaves some coordinates unconstrained] → Flag near-singular Jacobian columns and report those refits as conditional.
- [Many refits are expensive] → Bounded-block path, one start per variation; timings recorded.
- [Threshold choice drives the verdict] → Report all $p$-values so the reader can apply another threshold.

## Migration Plan

Add options and driver, run, store new ledgers; existing outputs untouched. Rollback removes the new options and files.
