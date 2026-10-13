---
title: Numerical implementation of charge-resolved K pi S waves in ALP to K Kbar pi
status: numerical
created: 2026-09-23
updated: 2026-09-23
bibtex_keys:
  - Alda:2026nrd
  - Omnes:1958hv
  - ParticleDataGroup:2024cfk
provenance: numerical-computation
---

# Numerical $K\pi$ $S$ Waves in $a\to K\bar K\pi$

## Scope and inputs

The implementation is [`scripts/implementation/kpi_s_wave_decay.py`](../../scripts/implementation/kpi_s_wave_decay.py). It covers the four charge-conserving flavor modes $K^+\pi^0K^-$, $K^0\pi^0\bar K^0$, $K^+\pi^-\bar K^0$, and $K^0\pi^+K^-$, as well as their physical-neutral-kaon projections $\pi^0K_SK_S$, $\pi^0K_LK_L$, $\pi^0K_SK_L$, $\pi^+K^-K_{S,L}$, and $\pi^-K^+K_{S,L}$. The projections include indirect CP violation through the configurable $\epsilon_K$ convention documented in [the PDG extract](../extracts/ParticleDataGroup:2024cfk/neutral-kaon-mixing.md) [@ParticleDataGroup:2024cfk]. The internal ordering is $(K,\pi,\bar K)$, with $s=(p_K+p_\pi)^2$ and $t=(p_\pi+p_{\bar K})^2$.

The script evaluates the flavor traces of the leading-order $a\Phi^3$ interaction with diagonal hatted couplings, nondegenerate fitted $B_0m_q$, physical neutral generators, and physical charged-neutral masses in the conventions of `Alda:2026nrd` [@Alda:2026nrd]. All energies are in GeV and widths are returned in GeV.

## Band decomposition and final-state interactions

After eliminating the $K\bar K$ invariant, the contact amplitude is affine in $s$ and $t$. The implementation assigns the $s$ and $t$ slopes to their corresponding $K\pi$ bands and splits the invariant-independent term equally between them. This equal split is an explicit convention: it is invariant under exchange of the two bands and reproduces the full leading-order contact amplitude exactly when the final-state-interaction factors equal one, but only the sum of constants is fixed for non-charge-conjugate bands.

A caller may supply a complex scalar function `fsi(invariant, pair)`. This realizes

$$
\mathcal M(s,t)=P_s(s)F_s(s)+P_t(t)F_t(t).
$$

The default is $F_s=F_t=1$. The normalized scalar implementation in [the numerical Omnes record](kpi-omnes-numerical.md) can be supplied directly or dispatched by physical pair [@Omnes:1958hv].

For $K_S$ and $K_L$, the default $F_s=F_t=1$ is the direct coherent projection of the LO flavor contact amplitude for the selected `epsilon_k`. Passing a nontrivial scalar function keyed by a physical pair such as `(KS, pi0)` is an additional effective physical-basis dressing assumption. It does not implement coupled $K^0$--$\bar K^0$ propagation or derive the physical-pair function from flavor-basis scattering. Repeated $K_SK_S$ and $K_LK_L$ states include the $2!$ phase-space symmetry factor.

### Example with the Pelaez--Rodas phase

The following example uses the central constrained $I=1/2$ elastic phase of Pelaez and Rodas [@Pelaez:2016tgi] for both charge-conjugate $K\pi$ bands in $a\to K^+K^-\pi^0$:

```python
import numpy as np

from scripts.implementation import (
    OmnesFunction,
    PairOmnesFSI,
    PelaezRodas2016ElasticPhase,
    TabulatedOmnes,
    s_wave_decay_width_fast,
)

phase = PelaezRodas2016ElasticPhase(isospin="1/2", fit="CFD")

omega = OmnesFunction(
    phase=phase,
    threshold=phase.threshold,
    phase_input_max=phase.maximum,
    subtraction_point=0.0,
    tail="constant",
    # Explicit assumption: freeze the elastic phase above the K-eta threshold.
    asymptotic_phase=phase(phase.maximum),
)

alp_mass = 1.2
kpi_min = (0.493677 + 0.1349768) ** 2
kpi_max = (alp_mass - 0.493677) ** 2
omega_fast = TabulatedOmnes(
    omega,
    minimum=kpi_min,
    maximum=kpi_max,
    relative_tolerance=1e-6,
)

fsi = PairOmnesFSI(
    {
        ("K+", "pi0"): omega_fast,
        ("K-", "pi0"): omega_fast,
    }
)

c_q_A = np.diag([1.0, 1.0, 0.0])
c_q_P = np.zeros((3, 3))

width = s_wave_decay_width_fast(
    alp_mass=alp_mass,
    c_q_A=c_q_A,
    c_q_P=c_q_P,
    channel=("K+", "K-", "pi0"),
    f_a=1000.0,
    fsi=fsi,
    relative_tolerance=1e-6,
)
print(width)  # GeV
```

Here `phase.threshold` and `phase.maximum` use the isospin-limit masses and elastic domain of the cited fit, while `kpi_min` and `kpi_max` cover the physical invariant range of each charged band. `TabulatedOmnes` validates its interpolation against omitted direct evaluations before the width is integrated. The constant continuation is not supplied by Pelaez and Rodas; it is an explicit high-energy assumption required because the Omnes integral extends to infinity. This scalar example also omits the $I=3/2$ source component, physical charge-channel mixing, and inelastic coupled-channel dynamics. It should therefore be used only when those approximations are appropriate.

## Assumptions and omissions

The implemented dressing is a charge-resolved diagonal independent-isobar approximation. Physical masses retain split $K\pi$ thresholds, and the flavor traces retain source-level isospin breaking. The implementation does not solve the coupled charge-basis Omnes matrix, and therefore does not reproduce charge-exchange rescattering or the exact two-threshold cusp. It also omits dynamical $m_d-m_u$ corrections to $K\pi$ scattering, electromagnetic corrections, and crossed-channel three-body rescattering. These omissions must not be interpreted as a complete isospin-breaking calculation.

## Verification

Run

```bash
.venv/bin/python -m unittest scripts.implementation.test_kpi_s_wave_decay
.venv/bin/python -m unittest scripts.implementation.test_omnes scripts.implementation.test_kpi_phases
```

The deterministic tests verify exact recovery of the full flavor-trace contact amplitude for unity dressing, agreement with the analytical isospin-limit $K^+\pi^0K^-$ kernel, all four flavor modes and seven physical-neutral-kaon modes, charged-neutral threshold splitting, threshold and coupling validation, quadratic coupling scaling, constant-FSI scaling, one-neutral-kaon projection, two-neutral-kaon interference, and identical-particle factors.
