---
title: Bitwise-identical reuse of repeated calculations in the isospin-breaking VPP likelihood
status: numerical
created: 2026-11-10
updated: 2026-11-10
bibtex_keys: [Alda:2026nrd, Qin:2024ulb]
provenance: numerical-computation
---

# Bitwise-identical reuse of repeated calculations in the isospin-breaking $VPP$ likelihood

Module: [`scripts/joint_vector_isospin_breaking_vpp_reuse.py`](../../../scripts/joint_vector_isospin_breaking_vpp_reuse.py) (class `IsospinBreakingVppLikelihoodWithReuse`); tests: [`scripts/test_joint_vector_isospin_breaking_vpp_reuse.py`](../../../scripts/test_joint_vector_isospin_breaking_vpp_reuse.py). The driver's `FitLikelihood` of [joint-vector-isospin-breaking-vpp-restart-option.md](joint-vector-isospin-breaking-vpp-restart-option.md) now inherits from this class. The original `scan_point`/`zero_control` ledger was produced with `IsospinBreakingVppLikelihood` itself.

```bash
.venv/bin/python -m unittest scripts.test_joint_vector_isospin_breaking_vpp_reuse
.venv/bin/python scripts/joint_vector_isospin_breaking_vpp_reuse.py   # prints the bitwise comparison at orders (1,3,4)
```

## Scope and assumptions

- The physics is unchanged. The option of [joint-vector-isospin-breaking-vpp-likelihood-option.md](joint-vector-isospin-breaking-vpp-likelihood-option.md), in `joint_vector_isospin_breaking_vpp_option.py`, is not modified.
- The class removes three repeated calculations identified in [joint-vector-isospin-breaking-vpp-evaluation-cost.md](joint-vector-isospin-breaking-vpp-evaluation-cost.md).
- Equality is tested, not assumed. The tests compare IEEE-754 bit patterns of the residual and the Jacobian against the parent class.
- Units and conventions are those of the option record [@Alda:2026nrd].

## Method

1. **Zero-excitation running propagator.** With $r\neq0$, `mode_scaled` follows each $K^*$ call, and the $s$-channel $\rho$ call, of `running_propagator(name, x, beta, ...)` with `running_propagator(name, x, (0,0), ...)` at the same invariant $x$.

   In the reference propagator, `fit_babar_kkpi_g_coefficients.propagator`, the excited-multiplet terms of Qin:2024ulb [@Qin:2024ulb] enter as $\beta_1 D_1+\beta_2 D_2$. The ground pole is computed once on the invariant and broadcast over the excitation choices. With $\beta=0$, both choices are therefore the ground pole from the first call, provided $D_1, D_2$ are finite.

   The second call is thus $\mathrm{BW}+(D_0-D_0)$. Here $D_0$ is that ground pole and $\mathrm{BW}=1/(m^2-x-im\Gamma(x))$ is the rescaled running Breit-Wigner of the first call. `GroundPropagatorMemory` keeps $\mathrm{BW}$ and $D_0$, and returns that expression when four things agree: the name, the invariant object, the reference-parameter object, and the mass and shared width. Any other call evaluates the running-propagator expression written out operation by operation.
2. **Duplicate reference residual.** The option's `jacobian` calls the inherited Jacobian, whose first step is `residual(values)`. It then calls `residual(values)` again as the reference for the $r_1,r_2,r_3$ forward differences. The class keeps the first result and uses it as that reference. The forward-difference steps and arithmetic are the same as before.
3. **Vertex factors.** $T(r)/T(0)$ of the $K^*$ and $\phi$ modes, and of the $s$-channel modes, were recomputed for every energy-bin block. The class keeps them per key, where the key is the exact bytes of $r$ and of the flavor matrices. It hands out a fresh copy each time, because the option updates the dictionary in place. At most 8 keys are kept.

## Result (numerical)

The figures below are at the restart point and $(2,16,16)$; see the [cost record](joint-vector-isospin-breaking-vpp-evaluation-cost.md).

| | option class | reuse class | change |
|---|---|---|---|
| residual, restart $r$ | 0.0725 s | 0.0454 s | $-37\%$ |
| Jacobian, restart $r$ | 1.518 s | 0.939 s | $-38\%$ |
| Jacobian, $r=0$ | 1.030 s | 0.906 s | $-12\%$ |

The approximate saving per Jacobian at the restart $r$ splits as follows:

- item 1: $\approx0.47$ s;
- item 3: $\approx0.08$ s;
- item 2: one residual, $\approx0.045$ s with the reuse class (0.073 s with the option class).

The remaining extra cost over the baseline (ratio 1.22) is mainly the three forward-difference residuals for $r_1,r_2,r_3$. These cannot be restricted to one channel, because $r$ enters every block.

## Verification

`scripts.test_joint_vector_isospin_breaking_vpp_reuse` has 6 tests, run on 2026-11-10, all passing:

- bit patterns of the residual and the $75\times35$ Jacobian equal those of `IsospinBreakingVppLikelihood` at $(2,16,16)$ with four-bin blocks, at three points: the restart endpoint, the same point with $r=0$, and the same point with $r$ at the scan point $(-0.350,1.855,0.729)$ GeV$^{-2}$;
- the driver's `FitLikelihood` inherits from the class;
- the module globals that are patched temporarily are restored after a Jacobian;
- `GroundPropagatorMemory` agrees bit for bit with `running_propagator` for $K^*$ and $\rho$, for both the full and the zero-excitation call, with one evaluation and one reuse;
- no reuse happens when the invariant object, the shared width or the name differs;
- remembered factors equal the parent values and are copies.

The benchmark also confirmed bitwise equality at $(2,16,16)$ for the restart point and $r=0$. The restart dry run printed the same output, apart from the elapsed time, before and after the driver switched to this class. Equality at points other than these three was not tested. The argument above covers finite propagators only.
