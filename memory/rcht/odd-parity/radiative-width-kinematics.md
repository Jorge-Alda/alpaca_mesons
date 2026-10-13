---
title: Normalized kinematics and widths for vector radiative decays
status: derived
created: 2026-10-04
updated: 2026-10-13
bibtex_keys: [Ruiz-Femenia:2003jdx, Ecker:1989yg, ParticleDataGroup:2024cfk]
provenance: analytical-derivation
---

# Kinematics and widths for $V\leftrightarrow P\gamma$

The normalized action-derived amplitude, including direct VJP terms, linear-$f_-$ contributions, and the complete neutral-vector exchange, is given in [radiative-sm-amplitude.md](radiative-sm-amplitude.md). The executable derivation contracts the published tensor operators with their Fourier phases and verifies the canonical $\omega\to\pi^0\gamma$ amplitude of Ruiz-Femenia et al. [@Ruiz-Femenia:2003jdx, Eq. (27)]. No coefficient is calibrated to data.

## Amplitude coefficient

For $i\mathcal M=i\epsilon_{\alpha\beta\rho\sigma}\epsilon_V^\alpha\epsilon_\gamma^\beta q^\rho p_P^\sigma G$, $G$ has units GeV$^{-1}$. The direct and exchange contributions are

$$
G_{\rm direct}=\frac{2\sqrt2e}{FM_{\rm ext}M_V}
\left[F_Q\big((c_1+c_2-c_5)m_P^2+(-c_1+c_2+c_5-2c_6)M_{\rm ext}^2\big)+4c_3F_H-4c_4F_K\right],
$$

$$
G_{\rm exch}=-\frac{4eF_V}{FM_{\rm ext}}\sum_{\alpha=\rho^0,\omega,\phi}
\frac{w_\alpha}{M_\alpha^2}
\left[F_{Q,\alpha}\big((d_1-d_3)m_P^2+d_3M_{\rm ext}^2\big)+4d_2F_{H,\alpha}\right],
\qquad G=G_{\rm direct}+G_{\rm exch}.
$$

$M_V$ normalizes the operators, $M_{\rm ext}$ is the physical external-vector mass, and $M_\alpha$ is the internal neutral-vector pole mass. The same external vector mass is used in both decay orientations.

## Widths

With the initial-state spin average included,

$$
\Gamma(V\to P\gamma)=\frac{|G|^2(M_{\rm ext}^2-m_P^2)^3}{96\pi M_{\rm ext}^3},\qquad
\Gamma(P\to V\gamma)=\frac{|G|^2(m_P^2-M_{\rm ext}^2)^3}{32\pi m_P^3}.
$$

The different denominators account for the spin average over an initial vector in the first orientation and an initial pseudoscalar in the second. Closed channels are excluded.

## Reproduction and scope

```text
.venv/bin/python scripts/derive_radiative_sm_amplitude.py
.venv/bin/python -m unittest scripts.implementation.test_radiative_sm_verification
```

The implementation is [`scripts/implementation/radiative_sm.py`](../../../scripts/implementation/radiative_sm.py). It requires explicit flavor matrices and the mass spurion $\chi=2B_0m_q$. Physical $\eta,\eta'$ inputs and a matched action for split internal masses remain unspecified. The source-specific $\omega\to\pi^0\gamma$ predictor remains available in `radiative_odd_parity.py`.
