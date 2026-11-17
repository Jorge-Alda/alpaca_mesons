---
title: Only-up angular-mode inputs for finite-mixing decays
status: mixed
created: 2026-11-17
updated: 2026-11-17
bibtex_keys: [DiLuzio:2026poh, Alda:2026nrd, ParticleDataGroup:2024cfk]
provenance: assumption
---

# Only-up angular-mode inputs

The user specifies $n_u=1$, $n_d=n_s=0$, identifies the angular mode with the ALP, and sets $f_\varphi=f_a$. With the radial field frozen, the phase sources and their sign in the conventions of Alda:2026nrd are given in the [phase-source record](alp-murayama-phase-mass-sources.md). The relation $m_u=\kappa_u f_\varphi/\sqrt2$ follows from DiLuzio:2026poh, so $f_a=\sqrt2m_u/\kappa_u$ for positive $\kappa_u$.

The specified inputs are $m_{\pi^\pm}=0.13957061$ GeV, $f_\pi=0.0923$ GeV, $m_u=0.00220$ GeV, $m_d=0.00469$ GeV and $m_s=0.09274$ GeV, with $B_0=m_{\pi^\pm}^2/(m_u+m_d)$. The user specifies $m_{\eta'}^2=m_1^2+\frac23B_0(m_u+m_d+m_s)$, fixing $m_1^2$ from the central value $m_{\eta'}=957.78\pm0.06$ MeV in the [2025 PDG listing](../../extracts/ParticleDataGroup:2024cfk/eta-prime-mass-2025-update.md) [@ParticleDataGroup:2024cfk]. This prescription sets the singlet diagonal entry to the measured mass squared; it is not in general a physical eigenvalue after mixing. See the [quadratic matrix extract](../../extracts/DiLuzio:2026poh/angular-mode-couplings.md).

The central inputs give $B_0=2.8272794159$ GeV and $m_1^2=0.7295546296$ GeV$^2$. The input and decay adapters are [`only_up_inputs.py`](../../../scripts/implementation/only_up_inputs.py) and [`only_up_decays.py`](../../../scripts/implementation/only_up_decays.py). They use LO angular-mode sources and omit radial exchange, loops and resonances. The existing finite-mixing routines have pole-coordinate singularities near exceptional poles.

## Verification

The [adapter checks](../../../scripts/implementation/test_only_up_decays.py) compare the four neutral eigenvalues and signed ALP leg at $\kappa_u=0.01$ to the matrix of DiLuzio:2026poh, with the meson coordinates reversed when translating conventions. They check the source derivatives, photon contractions, common WZW orientation and the closed three-pion threshold. At $\kappa_u=0.3$, the $\pi^0\pi^+\pi^-$ channel is open and the adapter width agrees with the underlying implementation, giving $0.0003113357$ GeV. This is a strong-mixing LO example and does not establish higher-order accuracy.

On 2026-11-17 the following command passed 26 checks, including 21 existing WZW checks:

```sh
.venv/bin/python -m unittest scripts.implementation.test_only_up_decays scripts.implementation.test_alp_photon_pp_wzw_finite_mixing
```

Three existing three-pseudoscalar checks also passed after mapping their deleted original module name to the present finite implementation. Their wrapper comparison checks integration of the existing implementation, not an independent derivation of its amplitude.

## WZW photon widths

The photon width functions in `only_up_decays.py` use physical pole masses and the WZW coefficients. These are WZW contributions to the widths. Vector-meson contributions interfere with WZW contributions at the amplitude level; the full width requires the coherent sum of amplitudes.

The electromagnetic input is fixed to $\alpha_{\rm em}=1/137$, with $e=\sqrt{4\pi\alpha_{\rm em}}$. For the coefficient $C_{a\gamma\gamma}$ defined in the [finite-mixing photon record](../photon/alp-lo-finite-mixing-photon-coefficient.md), the two-photon function evaluates
$$
\Gamma_{\gamma\gamma}^{\mathrm{WZW}}=\frac{\alpha_{\rm em}^2m_a^3}{64\pi^3f_a^2}|C_{a\gamma\gamma}|^2.
$$
For the momentum-independent WZW coefficient $\mathcal A_{\mathrm{WZW}}$ and the phase-space conventions of the [photon plus two pseudoscalars width record](../photon-two-pseudoscalar/alp-photon-pp-width.md), integration over the second Dalitz invariant gives
$$
\Gamma_{\gamma PP}^{\mathrm{WZW}}=\frac{|\mathcal A_{\mathrm{WZW}}|^2}{6144\pi^3m_a^3S_{PP}}\int_{(m_1+m_2)^2}^{m_a^2}ds\,\frac{(m_a^2-s)^3\lambda(s,m_1^2,m_2^2)^{3/2}}{s^2}.
$$
Here $S_{PP}=2$ for identical daughters and $1$ otherwise; the integral is zero at or below threshold. The remaining integration uses adaptive quadrature. Neutral kaon pole masses are degenerate for these diagonal sources.

The adapter checks compare the one-dimensional integral with an independent two-dimensional Dalitz integral and cover threshold, identical daughters, daughter exchange and two-photon normalization. On 2026-11-17 the command above passed 27 checks. At $\kappa_u=0.3$, the WZW widths are $6.14279\times10^{-8}$ GeV for $\gamma\gamma$ and $3.93070\times10^{-8}$ GeV for $\gamma\pi^+\pi^-$. These numerical examples use the LO model at strong mixing.
