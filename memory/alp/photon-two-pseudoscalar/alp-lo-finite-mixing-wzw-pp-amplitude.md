---
title: LO WZW amplitude for a -> gamma P P at finite ALP-meson mixing in dressed invariants
status: derived
created: 2026-11-14
updated: 2026-11-14
bibtex_keys: [Alda:2026nrd, Wess:1971yu, Bauer:2020jbp]
provenance: analytical-derivation
---

# LO WZW amplitude for $a\to\gamma P_1P_2$ at finite ALP-meson mixing in dressed invariants

## Inputs

- Gauged $\gamma3P$ WZW term of the [WZW $\gamma3\pi$ record](../../rcht/odd-parity/wzw-gamma3pi-contribution.md), with the common orientation $\mathrm{orientation}=-1$ fixed by the signed neutral-pion matching in [`scripts/derive_wzw_external_currents.py`](../../../scripts/derive_wzw_external_currents.py) and recorded in the [Levi-Civita convention record](../../rcht/odd-parity/levi-civita-convention.md) [@Wess:1971yu].
- Direct axial-source insertion $\Phi\to\Phi-2\epsilon c_q^Aa$, $\epsilon=f_\pi/f_a$, from the [WZW record](wzw-alp-photon-pp-amplitude.md) [@Alda:2026nrd].
- The dressed invariants and the tadpole source $\mathcal C=-B_0\{m_q,c_q^{AP}\}-\frac{m_0^2}{N_f}c_{GA}\mathbf 1$ of the [photon coefficient record](../photon/alp-lo-finite-mixing-photon-coefficient.md), and the mass operator $\mathbb M_s(X)=\{B_0m_q,X\}+\frac{m_0^2}{N_f}\operatorname{tr}(X)\mathbf 1-sX$.
- Physical pole masses $m_a$ (ALP), $m_1$, $m_2$ (daughters) of the neutral system, from the [spectrum solver](../../../scripts/solve_alp_meson_spectrum.py).

## Assumptions

- LO, tree level, stable simple poles; no ALP–meson degeneracy among the external poles.
- Sources commute with $Q$.
- Only the leading WZW sector. Vector-meson exchange, loops, NLO and widths are excluded.
- No exchange diagram: the $\gamma PP'$ three-point vertex vanishes by momentum conservation in this sector.
- Parity-even photon bremsstrahlung from $a\to P^+P^-$ (minimal coupling) has a different Lorentz structure and is not included.

## Identities used (not assumptions)

- Terms with two or three ALP-component legs contribute nothing: the Levi-Civita contraction gives $\epsilon^{\mu\nu\rho\sigma}\partial_\mu a\,\partial_\nu a=0$, so each such term is zero identically.
- Each external leg enters the alternating trace only through its flavour matrix, the meson component of its pole eigenvector plus its ALP component times $-2\epsilon c_q^A$. In the canonical frame this matrix is the meson-component matrix of the pole alone, so the $c_q^A$ dependence cancels leg by leg.

## Leg matrices

For a neutral pole at physical mass $m$ with nonzero ALP component, the leg is the flavour matrix
$$
\frac{2\epsilon\,\mathbb M_{m^2}^{-1}(\mathcal C)}{\sqrt{1+2\epsilon^2\operatorname{tr}\big[\big(\mathbb M_{m^2}^{-1}(\mathcal C)\big)^2\big]}},
$$
with the normalization written out explicitly, so no further symbol is needed. For the ALP leg $m=m_a$. For a charged daughter (no ALP admixture) the leg is $\lambda^\dagger$, the outgoing generator. In the limit $\epsilon\to0$ a neutral daughter leg tends to its unmixed generator, with residuals of order $\epsilon^2$.

## Amplitude

For both neutral daughters, the amplitude is
$$
\boxed{\begin{aligned}
\mathcal M(a\to\gamma P_1P_2)&=\epsilon_{\mu\nu\rho\sigma}\varepsilon_\gamma^\mu p_1^\nu p_2^\rho k^\sigma\;
\frac{N_c e}{3\pi^2f_a^3}\,
\prod_{m\in\{m_a,m_1,m_2\}}\frac{1}{\sqrt{1+2\epsilon^2\operatorname{tr}\big[\big(\mathbb M_{m^2}^{-1}(\mathcal C)\big)^2\big]}}\\
&\quad\times\sum_{\sigma\in S_3}\operatorname{sgn}(\sigma)\,
\operatorname{tr}\Big(Q\,\mathbb M^{-1}_{m_{\sigma(1)}^2}(\mathcal C)\,\mathbb M^{-1}_{m_{\sigma(2)}^2}(\mathcal C)\,\mathbb M^{-1}_{m_{\sigma(3)}^2}(\mathcal C)\Big)
\end{aligned}}
$$
Here $(m_{\sigma(1)},m_{\sigma(2)},m_{\sigma(3)})$ is the list $(m_a,m_1,m_2)$ reordered by the permutation $\sigma\in S_3$, and $\operatorname{sgn}(\sigma)$ is its sign. The overall sign is $+$, with $\mathrm{orientation}=-1$ already inserted. The prefactor $\frac{N_c e}{3\pi^2f_a^3}$ is $\frac{N_c e}{24\pi^2f_\pi^3}(2\epsilon)^3$, from the factor $2\epsilon$ on each neutral leg.

For a charged daughter, the leg $\mathbb M^{-1}_{m^2}(\mathcal C)$ is replaced by $\lambda^\dagger$ of that daughter, and that daughter contributes no factor $2\epsilon/\sqrt{\cdots}$ and no normalization factor. In general the prefactor is $\frac{N_c e}{24\pi^2 f_\pi^3}$ times $2\epsilon/\sqrt{\cdots}$ for each neutral leg.

**Weak-mixing limit.** With the ALP leg and charged daughters, and neutral daughters as unmixed generators, the amplitude reduces to
$$
\mathcal M\to\epsilon_{\mu\nu\rho\sigma}\varepsilon_\gamma^\mu p_1^\nu p_2^\rho k^\sigma\;
\frac{N_c e}{12\pi^2 f_a f_\pi^2}\Big[\operatorname{tr}\big(Q\,\mathbb M_{m_a^2}^{-1}(\mathcal C)\,\lambda_1^\dagger\lambda_2^\dagger\big)-\operatorname{tr}\big(Q\,\mathbb M_{m_a^2}^{-1}(\mathcal C)\,\lambda_2^\dagger\lambda_1^\dagger\big)-\operatorname{tr}\big(Q\lambda_1^\dagger\,\mathbb M_{m_a^2}^{-1}(\mathcal C)\,\lambda_2^\dagger\big)+\operatorname{tr}\big(Q\lambda_1^\dagger\lambda_2^\dagger\,\mathbb M_{m_a^2}^{-1}(\mathcal C)\big)+\operatorname{tr}\big(Q\lambda_2^\dagger\,\mathbb M_{m_a^2}^{-1}(\mathcal C)\,\lambda_1^\dagger\big)-\operatorname{tr}\big(Q\lambda_2^\dagger\lambda_1^\dagger\,\mathbb M_{m_a^2}^{-1}(\mathcal C)\big)\Big],
$$
which is the no-mixing table of the [WZW record](wzw-alp-photon-pp-amplitude.md).

**Invariance.** $\mathbb M_s^{-1}(\mathcal C)$ is built from dressed invariants and $m_q$ only, and the leg normalization is a function of the same object. So the full amplitude is unchanged by the redefinitions of the Q-commuting class, with the same caveat as the [photon record](../photon/alp-lo-finite-mixing-photon-coefficient.md).

## Checks

From [`scripts/derive_alp_lo_photon_channels.py`](../../../scripts/derive_alp_lo_photon_channels.py), tested in [`scripts/test_derive_alp_lo_photon_channels.py`](../../../scripts/test_derive_alp_lo_photon_channels.py) (18 tests):

1. Weak-mixing limit: the WZW table of the [WZW record](wzw-alp-photon-pp-amplitude.md) is reproduced, including the charged-channel coefficients and the overall prefactor. The $\pi^+\pi^-$ coefficient tends to $-2(\mathbb M^{-1})_{uu}$, with residual of order $\epsilon^2$.
2. Leg matrices from the solver agree with the invariant leg formula to $9\times10^{-14}$ at the physical poles.
3. Ward identity: contraction with $\varepsilon_\gamma\to k$ vanishes for arbitrary non-commuting leg matrices, so the result holds for pole-projected daughters.
4. Bose symmetry: the flavour factor is antisymmetric in the two daughter legs and vanishes for identical daughters.
5. Invariance under constant $\kappa$ (exact) and under random nonlinear Q-commuting redefinitions (amplitudes unchanged to $10^{-9}$ with the solver).

## Reproducibility

```sh
.venv/bin/python scripts/derive_alp_lo_photon_channels.py
.venv/bin/python -m unittest discover -s scripts -p test_derive_alp_lo_photon_channels.py
```

Both ran on 2026-11-14 and passed. The [existing WZW script](../../../scripts/derive_wzw_alp_photon_two_pseudoscalar.py) still passes its seven tests.

## Validity and open points

- **Zero-ALP neutral poles.** For a neutral pole with zero ALP component (the exact $K^0$ pole), $\mathbb M_{m^2}^{-1}$ is singular and the leg formula above does not apply. Its leg is that pole's own unit generator in the $K^0$ eigenspace. This case is not implemented.
- **Real neutral external fields.** For real neutral external states the six-term trace is purely imaginary, which is a phase convention of the external state, not a physical effect.
- **Orientation of the photon coupling.** The sign of the $\gamma\gamma$ coefficient is not tied to $\mathrm{orientation}=-1$ in any record; it is not checked here.
- **Scope.** LO only. Resonance parts of the [existing $a\to\gamma PP$ record](alp-photon-pp-amplitude.md) must be added coherently and are not included.
