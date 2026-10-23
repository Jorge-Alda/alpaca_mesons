---
title: Dilepton loop normalization and two-photon cut
status: literature
created: 2026-10-23
updated: 2026-10-23
bibtex_keys: [Messerli:2025rnv]
provenance: literature
---

# Reduced loop amplitude

Messerli:2025rnv defines the unnormalized transition form factor by the signed tensor in its Eq. (9), with $\epsilon^{0123}=+1$. The normalized TFF is $\widetilde F_{P\gamma^*\gamma^*}=F_{P\gamma^*\gamma^*}/F_{P\gamma\gamma}$, where $F_{P\gamma\gamma}=F_{P\gamma^*\gamma^*}(0,0)$. With $q^2=M_P^2$, the reduced two-photon amplitude is [@Messerli:2025rnv, Eqs. (6)--(9)]
$$
\mathcal A_\ell(q^2)=\frac{2i}{\pi^2q^2}\int d^4k\,
\frac{q^2k^2-(q\cdot k)^2}{k^2(q-k)^2[(p-k)^2-m_\ell^2]}
\widetilde F_{P\gamma^*\gamma^*}(k^2,(q-k)^2).
$$
The normalized rate is
$$
\frac{\mathrm{Br}[P\to\ell^+\ell^-]}{\mathrm{Br}[P\to\gamma\gamma]}=
2\sigma_\ell(q^2)\left(\frac{\alpha}{\pi}\right)^2
\frac{m_\ell^2}{q^2}|\mathcal A_\ell(q^2)|^2,
\qquad \sigma_\ell(q^2)=\sqrt{1-\frac{4m_\ell^2}{q^2}}.
$$

# Two-real-photon discontinuity

The cut through the two real photons is independent of the off-shell TFF after this normalization:
$$
\operatorname{Im}_{\gamma\gamma}\mathcal A_\ell(q^2)=
\frac{\pi}{2\sigma_\ell(q^2)}\log\!\left(\frac{1-\sigma_\ell(q^2)}{1+\sigma_\ell(q^2)}\right).
$$
This cut gives the unitarity lower bound only if no other imaginary contributions are present. Messerli:2025rnv notes additional hadronic cuts for the heavier $\eta^{(\prime)}$ channels; the same caution applies to a source-dependent ALP TFF. [@Messerli:2025rnv, Eqs. (13)--(15)]
