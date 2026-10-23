---
title: ALP dilepton two-photon cut and width normalization
status: derived
created: 2026-10-23
updated: 2026-10-23
bibtex_keys: [Messerli:2025rnv]
provenance: analytical-derivation
---

# Assumptions

Use the project's unnormalized scalar $F_{a\gamma^*\gamma^*}$ from the [tree-level doubly virtual TFF record](alp-doubly-virtual-tff.md), with tensor convention $i\mathcal M=i e^2 R F$. Define $F_{a\gamma\gamma}=F_{a\gamma^*\gamma^*}(0,0)$ and normalize $\widetilde F=F/F_{a\gamma\gamma}$. Work at leading QED order and $s=m_a^2>4m_\ell^2$. The only discontinuity evaluated here is the two-real-photon cut. The loop convention and normalized branching-ratio formula are taken from the [Messerli:2025rnv source extract](../../extracts/Messerli:2025rnv/dilepton-loop-formalism.md). No source profile or numerical ALP couplings are selected.

# Derivation

The executable [derivation](../../../scripts/derive_alp_lepton_pair_normalization.py) evaluates the pseudoscalar lepton spin trace, the two-photon polarization sum, and two-body phase space. It uses the coupling identification $g_{\ell}=2\alpha^2m_\ell F_{a\gamma\gamma}\mathcal A_\ell$ implied by the reduced loop normalization. With $\sigma_\ell=\sqrt{1-4m_\ell^2/s}$, the results are
$$
\sum_{\mathrm{spins}}|\bar u i\gamma_5 v|^2=2s,
\qquad \sum_{\mathrm{pol.}}|R|^2=\frac{s^2}{2},
$$
$$
\Gamma(a\to\gamma\gamma)=\frac{\pi\alpha^2}{4}|F_{a\gamma\gamma}|^2s^{3/2},
\qquad
\Gamma(a\to\ell^+\ell^-)_{\mathrm{TFF}}=
\frac{\alpha^4m_\ell^2|F_{a\gamma\gamma}|^2\sqrt{s}}{2\pi}
\sigma_\ell|\mathcal A_\ell(s)|^2.
$$
Consequently,
$$
\frac{\Gamma(a\to\ell^+\ell^-)_{\mathrm{TFF}}}{\Gamma(a\to\gamma\gamma)}=
2\sigma_\ell\left(\frac{\alpha}{\pi}\right)^2\frac{m_\ell^2}{s}|\mathcal A_\ell(s)|^2.
$$
Cutting both photon propagators in the reduced loop and integrating their relative angle gives
$$
\operatorname{Im}_{\gamma\gamma}\mathcal A_\ell(s)=
\frac{\pi}{2\sigma_\ell}\log\!\left(\frac{1-\sigma_\ell}{1+\sigma_\ell}\right),
\qquad \lim_{\sigma_\ell\to0^+}\operatorname{Im}_{\gamma\gamma}\mathcal A_\ell=-\pi.
$$
The reduced cut is negative in the stated convention. In the massless-lepton limit the two-photon contribution to the width vanishes as $m_\ell^2$, as required by helicity suppression.

# Verification and limits

Run `.venv/bin/python scripts/derive_alp_lepton_pair_normalization.py` for the exact expressions. The script contains the angular integrand and primitive used to obtain the cut. This fixes the width normalization and absorptive part, but not the full amplitude or a numerical width: its dispersive real part requires the source-specific doubly virtual ALP TFF over all virtualities, its high-energy matching, and any consistently matched local ALP-lepton term. Other open TFF cuts can also contribute to the imaginary part; therefore the squared cut alone is not asserted as a lower bound.
