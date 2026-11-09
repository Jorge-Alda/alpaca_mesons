---
title: Candidate causes of the joint vector-fit inadequacy and diagnostics to separate them
status: conjectured
created: 2026-11-09
updated: 2026-11-09
bibtex_keys: [Alda:2026nrd]
provenance: conjecture
---

# Candidate causes of the inadequacy

This record lists conjectures that follow the [adequacy verdict](joint-vector-fit-adequacy-verdict.md). None of them has been tested. They are ordered by priority for the follow-up issue.

1. **$K^*$ width prediction and isospin treatment.** This candidate needs no new operator.
   - The fit uses one shared $K^*$ width coordinate with $K\pi$ saturation, and counts the charged and neutral $K^*$ width measurements once each. See the [baseline record](babar-kkpi-global-vector-sector-fit-baseline.md).
   - **Diagnostics:**
     - Recompute the on-shell $\Gamma_{K^*}$ at the best point outside the fit, and compare it with the tree-level RChT width at physical masses.
     - Compute the irreducible $\chi^2$ of the two width rows for any common width.
     - Refit with separate charged and neutral $K^*$ widths and masses.
   - $K\pi$ saturation alone is an assumption in force. It is not expected to produce a shift of several MeV, but this expectation is not sourced.
2. **SU(3) breaking in the vector–pseudoscalar–pseudoscalar coupling.**
   - **Motivation:** $\phi\to K\bar K$ and $\Gamma_\phi$ pull in the direction opposite to the $K^*$ widths.
   - **Diagnostic:** free one $K^*K\pi$ to $\rho\pi\pi$ coupling ratio, and check whether the shared-width and $\phi\to K\bar K$ blocks improve together.
3. **The deferred double-trace source term $y\,\mathrm{tr}(\chi_+Q)\mathbb{1}$.**
   - This term would move the $\omega,\phi$ masses independently of $\rho,K^*$, provided the source enters the vector mass spectrum. That is to be confirmed.
   - **Diagnostic:** in the refit without the mass block, see which mass rows move.
4. **Isospin breaking in the even sector beyond $K^*$.** This covers $\rho$–$\omega$ mixing and the difference between $K^+K^-$ and $K^0\bar K^0$ in $\phi$ decays.
   - **Diagnostic:** rotate the whitened $\phi\to K\bar K$ residuals back to the two channels. A residual mainly in their ratio favours this candidate; a residual mainly in their sum favours candidate 2.
5. **$\eta'$ and odd-sector tension.**
   - **Diagnostic:** refit without the $\eta'$ radiative rows, and check whether both $\Gamma_{\eta'}$ and the odd block recover.
   - $\rho^0\to\eta\gamma$ would remain as a separate tension.
