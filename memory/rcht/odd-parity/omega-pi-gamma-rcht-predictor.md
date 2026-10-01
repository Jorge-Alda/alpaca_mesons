---
title: On-shell omega to neutral-pion gamma RChT predictor
status: mixed
created: 2026-09-30
updated: 2026-09-30
bibtex_keys: [Ruiz-Femenia:2003jdx]
provenance: analytical-derivation
---

# Scope and conventions

This record documents the first fit-compatible physical calibration predictor in the repository. It is only the on-shell $\omega\to\pi^0\gamma$ result printed by Ruiz--Femenia et al.; it is literature-formula algebra, not a fresh operator derivation [@Ruiz-Femenia:2003jdx]. The formula inputs are transcribed in the [on-shell radiative extract](../../extracts/Ruiz-Femenia:2003jdx/on-shell-radiative-amplitudes.md), while the five OPE inputs are recorded in the [odd-vector extract](../../extracts/Ruiz-Femenia:2003jdx/vector-odd-parity.md).

Natural units are used. $F,F_V,M_V,M_\omega,m_\pi$ and the caller-supplied $\Gamma_\omega^{\rm total}$ are in GeV, the scalar coefficient $G_{\omega\pi\gamma}$ is in $\mathrm{GeV}^{-1}$, the partial width is in GeV, and the branching fraction is dimensionless. $M_V$ is the common resonance-operator normalization mass and internal-vector mass appearing in the published approximation; it is not replaced by the external on-shell $M_\omega$.

## Literature evidence

With

$$
i\mathcal M=i\epsilon_{\alpha\beta\rho\sigma}
\epsilon_\omega^\alpha\epsilon_\gamma^\beta q^\rho k^\sigma
G_{\omega\pi\gamma},
$$

the uneliminated Eq. (27) is [@Ruiz-Femenia:2003jdx]

$$
G_{\rm direct}=\frac{2\sqrt2e}{M_\omega M_VF}
\left[(c_2-c_1+c_5-2c_6)M_\omega^2
+(c_1+c_2+8c_3-c_5)m_\pi^2\right],
$$

$$
G_\rho=-\frac{4eF_V}{M_V^2M_\omega F}
\left[d_3M_\omega^2+(d_1+8d_2-d_3)m_\pi^2\right],
\qquad G_{\omega\pi\gamma}=G_{\rm direct}+G_\rho.
$$

Thus $c_4,c_7,d_4$ are accepted as members of the full coupling records but are naturally absent from this channel. Applying all five published OPE equations gives Eq. (28),

$$
G_{\omega\pi\gamma}^{\rm OPE}=\frac e{F_V}\left[
\frac{N_C}{8\pi^2}\frac{M_\omega}{F}
-\frac F2\frac{M_\omega}{M_V^2}
\left(1+\frac{m_\pi^2}{M_\omega^2}\right)\right].
$$

The generic scalar-coefficient width normalization is

$$
\Gamma(V\to P\gamma)=\frac{|G|^2(M_V^{\rm ext\,2}-m_P^2)^3}
{96\pi M_V^{\rm ext\,3}}.
$$

After $F_V=\sqrt2F$ and $\alpha=e^2/(4\pi)$, this becomes the printed Eq. (29),

$$
\Gamma(\omega\to\pi\gamma)=\frac{\alpha}{192}M_\omega
\left(1-\frac{m_\pi^2}{M_\omega^2}\right)^3
\left[\frac{N_C}{4\pi^2}\frac{M_\omega^2}{F^2}
-\frac{M_\omega^2}{M_V^2}\left(1+\frac{m_\pi^2}{M_\omega^2}\right)\right]^2.
$$

## Derived evidence and method

[`derive_omega_pi_gamma_rcht.py`](../../../scripts/derive_omega_pi_gamma_rcht.py) performs the exact SymPy calculation in these steps:

1. construct the two published Eq. (27) mass combinations and their relative sign;
2. solve the five OPE equations for $c_3,c_2,c_6,d_1,d_3$ without assigning the remaining free couplings;
3. substitute them into the direct and $\rho$ terms separately and simplify their sum;
4. assert exact equality to Eq. (28);
5. insert Eq. (28) into the generic width, then apply $F_V=\sqrt2F$ and $e^2=4\pi\alpha$;
6. assert exact equality to Eq. (29) and verify that $c_4,c_7,d_4$ never occur.

Invocation:

```text
.venv/bin/python scripts/derive_omega_pi_gamma_rcht.py
```

## Numerical implementation evidence

[`radiative_odd_parity.py`](../../../scripts/implementation/radiative_odd_parity.py) provides frozen, validated `OmegaPiGammaInputs`, separate signed direct and $\rho$ coefficients, their full uneliminated sum, the constrained Eq. (28) coefficient, the generic width, and one-element runtime predictors `rcht.omega_pi0_gamma.width` and `rcht.omega_pi0_gamma.branching_fraction`. Predictor parameters are `c1` through `c7` and `d1` through `d4`. The width configuration has exactly `e`, `F`, `FV`, `MV`, `Momega`, and `mpi`; the branching configuration additionally requires positive `GammaOmegaTotal` and returns

$$
\mathcal B(\omega\to\pi^0\gamma)=\frac{\Gamma(\omega\to\pi^0\gamma)}{\Gamma_\omega^{\rm total}}.
$$

The width predictor is unchanged. Both predictors use the uneliminated coupling so they can calibrate fitted coefficient combinations. The branching predictor is connected to the [PDG datum](../../extracts/ParticleDataGroup:2024cfk/omega-pi0-gamma-datum.md), but no fit is run.

The implementation fixes the external antisymmetric-tensor residue to the published tree-level canonical value $Z_\omega=1$ exactly. It does not expose a residue nuisance parameter.

Tests independently transcribe Eq. (27), build a point satisfying all five OPE relations and compare it with Eq. (28), reproduce Eq. (29), check linear coupling and quadratic-width scaling, preserve the direct/$\rho$ interference sign, verify both registry predictions and immutability, compare the branching fraction directly with partial width divided by total width, and reject missing, unexpected, nonfinite, nonpositive-total-width, and closed-channel inputs:

```text
.venv/bin/python -m unittest scripts.implementation.test_radiative_odd_parity
.venv/bin/python -m unittest discover -s scripts/implementation -t .
```

## Explicit assumptions and limits

- **Literature assumptions:** ideal $\omega$--$\phi$ mixing, one common $M_V$, tree-level narrow-width on-shell kinematics, and canonical $Z_\omega=1$ [@Ruiz-Femenia:2003jdx].
- **Implementation scope:** the single-resonance antisymmetric-tensor RChT formula is used; HLS is not used.
- **No silent generalization:** there is no SU(3) or channel substitution. The source explicitly warns that $K^*\to K\gamma$ and $\phi\to\eta\gamma$ require flavor-resolved short-distance matching [@Ruiz-Femenia:2003jdx]. Other $V\to P\gamma$ predictors remain blocked until that matching is implemented.
- **Data boundary:** this predictor module itself introduces no experimental value or likelihood; the separate validated loader connects the PDG datum. No physical fit has yet been run because the first-release problem remains rank deficient.
