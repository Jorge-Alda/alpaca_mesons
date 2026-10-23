---
title: Plan for fitting vector-to-pseudoscalar dilepton form factors
status: mixed
created: 2026-10-23
updated: 2026-10-23
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx, Adlarson:2016hpp, NA60:2016nad]
provenance: assumption
---

# Plan for $V\to P\ell^+\ell^-$ differential form factors

## Goal and scope

Use the measured vector Dalitz spectra to constrain the finite-virtuality odd-parity $VP\gamma^*$ vertex in the preliminary-plus-$G_V$ fit. The first fit block is $\omega\to\pi^0e^+e^-$ from Adlarson:2016hpp and $\omega\to\pi^0\mu^+\mu^-$ from NA60:2016nad. Their electron and muon channels probe the same hadronic form factor over complementary $q^2=m_{\ell\ell}^2$ ranges. The associated $\eta\to\ell^+\ell^-\gamma$ measurements in these papers are not vector Dalitz data and must not be combined as if they measured the same form factor. Source-backed summaries are in the [A2 extract](../../extracts/Adlarson:2016hpp/omega-pi-electron-spectrum.md) and [NA60 extract](../../extracts/NA60:2016nad/omega-pi-muon-spectrum.md).

The first pass should fit $\omega\to\pi^0\ell^+\ell^-$. Extend to other $VP$ channels only after the corresponding experimental spectra, flavor factors, and mixing conventions have been catalogued. Do not use a fitted pole slope as a substitute for the binned spectrum.

## Amplitude and differential rate

For $V(p,\epsilon_V)\to P(k)\ell^+(p_+)\ell^-(p_-)$, set $q=p_++p_-=p-k$ and $q^2=m_{\ell\ell}^2$. Define the transition form factor by the tensor structure

$$
\mathcal M(V\to P\gamma^*)= e\,F_{VP}(q^2)\,\epsilon_{\mu\nu\alpha\beta}\,
\epsilon_V^\mu p^\nu q^\alpha \epsilon_\gamma^{*\beta},
$$

with an explicit convention note mapping this scalar to the repository's signed canonical epsilon convention. The measured distribution is proportional to

$$
\frac{d\Gamma(V\to P\ell^+\ell^-)}{dq^2}
\propto \frac{\alpha}{q^2}\,\beta_\ell(q^2)
\left(1+\frac{2m_\ell^2}{q^2}\right)
\lambda^{3/2}(m_V^2,m_P^2,q^2)\,|F_{VP}(q^2)|^2,
\qquad \beta_\ell(q^2)=\sqrt{1-\frac{4m_\ell^2}{q^2}},
$$

where $\lambda(x,y,z)=x^2+y^2+z^2-2xy-2xz-2yz$. Fix the overall rate normalization against the measured $V\to P\gamma$ width using the same tensor and phase convention, and fit the $q^2$ dependence. Before implementation, verify the experiment's exact form-factor normalization and bin observable from its paper; include the normalization factor in the executable data adapter rather than relying on the proportionality above.

## Theory input

Use the retained odd-parity VJP and VVP Lagrangian and the even vector-photon response from the preliminary-plus-$G_V$ setup. At tree level, the virtual photon couples through both the direct VJP vertex and vector exchange from VVP followed by photon-vector mixing. Evaluate the complete channel amplitude at finite $q^2$, retaining the declared vector propagator prescription and flavor/mixing inputs. The starting point and sign convention are the repository's [neutral TFF derivation](radiative-tff-amplitude.md), [radiative SM amplitudes](radiative-sm-amplitude.md), and [fit Lagrangian and null basis](preliminary-fit-lagrangian-and-null-directions.md). Translate their all-incoming momenta to the physical decay kinematics explicitly.

Use the no-mixing ALP basis only for the ALP source if that sector is later propagated; it does not remove physical $\omega$--$\phi$ mixing from the SM vector propagator or external vector state. Keep the ALP basis convention and the vector mixing prescription as separate inputs. The first fit should follow the same vector-sector prescription as the preliminary-plus-$G_V$ fit and not silently replace it with ideal mixing or a single-pole slope model.

## Data and likelihood

1. Transcribe each published differential spectrum with its bin edges, measured quantity, statistical and systematic uncertainties, and covariance where reported. Preserve the experiment's normalization convention and any radiative-width input.
2. Convert each bin into the same predicted observable used by the experiment, including bin integration and any stated acceptance or resolution treatment. If only bin-center values are provided, document that approximation and compare it with bin integration.
3. Fit the A2 electron and NA60 muon spectra jointly with the existing low-energy odd-sector data, counting shared external widths once and keeping experiment-specific normalization nuisances explicit. Use source-level covariance; if it is unavailable, state the diagonal approximation and keep that block's rank diagnostic conditional.
4. Retain the current fitted even-sector quantities and $G_V$ according to the existing preliminary fit's treatment. Profile only nuisance parameters that the published likelihood supports; report fixed-nuisance and profiled results separately.
5. Check residuals and stability under the fit range and bin choices. In particular, do not interpret a poor high-$q^2$ description as evidence for a shifted odd coupling before checking whether the retained tree-level resonance model describes the spectrum.

## Null-direction and identifiability check

The radiative widths use $q^2=0$ and cannot by themselves determine the finite-virtuality deformation associated with $\delta c_1=\delta c_5$. The vector Dalitz spectra are promising because they sample nonzero $q^2$ in the same $VP\gamma^*$ amplitude. For every accepted fit point, compute the Jacobian of all old and new observables with respect to the original coefficient vector, project it onto the preliminary-fit null basis, and report singular values and the resolved linear combinations. Repeat with normalization and shape nuisances projected out. This is the test of whether the data actually lift the $c_1=c_5$ null; a nonzero amplitude derivative alone is not a precision measurement.

Do not claim separate determination of $c_1$ and $c_5$ from this channel unless the complete projected Jacobian supports it. The known common shift $\delta c_2=\delta c_5=\delta c_6$ is invisible to the neutral $P\gamma^*\gamma^*$ form factors, and the current radiative baseline also leaves it unprobed. Explicitly test whether the $VP\gamma^*$ rows lift that direction in the full channel implementation. If not, preserve it as a null coordinate and identify an additional observable before quoting individual $c_i$ values. Also carry through the $c_7$ and $d_4$ directions; their insensitivity in the current preliminary fit is a measured-channel statement, not a general operator identity.

## Completion criteria

The fit block is ready for inclusion in the global analysis when it has source-linked bin tables, documented covariance and normalization, a validated finite-$q^2$ amplitude in the canonical convention, a reproducible joint likelihood, residual plots, nuisance profiles, and the projected-null Jacobian. The main output should be the posterior or profile for data-constrained combinations and their correlations, not a preferred zero-null representative masquerading as separately measured couplings.

The conditional first-pass calculation is recorded in [the fit result](vector-dalitz-conditional-fit.md). Its old-null Jacobian adds one shape-sensitive direction, but the fixed tree-level model gives a poor joint fit and is not accepted as a coupling calibration. The result uses the existing analytical inputs [Alda:2026nrd, Ruiz-Femenia:2003jdx] and the published vector Dalitz measurements [Adlarson:2016hpp, NA60:2016nad].
