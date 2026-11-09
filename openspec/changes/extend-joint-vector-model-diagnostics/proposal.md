# Proposal

## Why

The joint vector-sector fit is inadequate: $\chi^2=225.7$ for 43 degrees of freedom ($p=1.3\times10^{-26}$) at all four accepted minima (alp-uzh.6). Three blocks fail $p<10^{-3}$: the shared total widths ($\Gamma_{K^{*+}}$ pull $-8.3$, $\Gamma_{K^{*0}}$ pull $-4.4$), $\phi\to K\bar K$, and the odd block ($\rho^0\to\eta\gamma$). The verdict holds without the double-trace source term and under all 47 variations. The identifiability study (alp-uzh.7) cannot start until a model meets the adequacy criterion or a documented reason shows that it cannot. The conjectured causes are recorded in `memory/rcht/odd-parity/joint-vector-fit-inadequacy-candidate-causes.md`; none has been tested.

## What Changes

- Run the cheap diagnostics listed in the candidate-causes record before adding any operator:
  - recompute the on-shell $\Gamma_{K^*}$ at the best point and compare it with the tree-level RChT width at physical masses;
  - compute the irreducible $\chi^2$ of the two $K^*$ width rows for any common width;
  - rotate the whitened $\phi\to K\bar K$ residuals back to the $K^+K^-$ and $K^0\bar K^0$ channels;
  - refit without the $\eta'$ radiative rows;
  - refit without the mass block, to see which mass rows move.
- Add switchable model extensions to the joint likelihood, each recovering the baseline residuals exactly when off, tried in the priority order of the issue:
  1. separate charged and neutral $K^*$ masses and widths, and a revised $K\pi$-saturation / shared-width treatment;
  2. SU(3) breaking in the vector–pseudoscalar–pseudoscalar coupling (a free $K^*K\pi$ to $\rho\pi\pi$ ratio);
  3. the double-trace source term $y\,\mathrm{tr}(\chi_+Q)\mathbb{1}$, previously deferred (2026-11-07), now derived for the vector mass spectrum and the $KK\pi$ current;
  4. isospin breaking in the even sector ($\rho$–$\omega$ mixing, $\phi\to K^+K^-$ versus $K^0\bar K^0$);
  5. the $\eta'$ and odd-block tension.
- Refit each extension with the acceptance gates of the baseline, and re-apply the adequacy criterion of `assess-joint-vector-fit-adequacy`. Stop adding extensions once the criterion is met.
- Record the outcome: an adequate extended model, or a documented reason why the model cannot meet the criterion.

## Capabilities

### New Capabilities

- `joint-vector-model-extension`: Diagnostics and switchable extensions of the joint vector-sector likelihood aimed at the inadequacy found in alp-uzh.6, ending in an adequacy result or a documented reason for failure.

### Modified Capabilities

None.

## Impact

- New option flags in `scripts/fit_babar_kkpi_global_vector_sector.py` and a new driver script, with focused tests. All script changes are delegated to the `coder` agent.
- New analytical derivations for any extension that adds an operator (double-trace source, SU(3)-breaking VPP coupling, isospin-breaking terms), each with a derivation script and a memory record.
- New ledgers under `data/vector_sector/` and indexed memory records under `memory/rcht/odd-parity/`.
- Baseline likelihood, accepted endpoints and the alp-uzh.6 ledgers are unchanged.
- Closes alp-uzh.10; unblocks alp-uzh.7 only if the adequacy criterion is met.
