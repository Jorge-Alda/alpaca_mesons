---
title: Diagnostic record for the pi0 K_L K_S vector-width slowdown
status: mixed
created: 2026-09-25
updated: 2026-09-30
bibtex_keys: [Alda:2026nrd, Ecker:1989yg, ParticleDataGroup:2024cfk]
provenance: analytical-derivation
---

# Diagnostic for the $\pi^0K_LK_S$ Vector-Width Slowdown

## Scope and evidence status

This record diagnoses the slowdown and records the implemented correction. The executable calculation is [`scripts/diagnose_pi0_kl_ks_vector_width.py`](../../scripts/diagnose_pi0_kl_ks_vector_width.py). For direct file invocation it adds the repository root to `sys.path`, then imports through the `scripts.implementation` package namespace; it does not use legacy top-level `implementation.*` imports and never launches the nested production width integral.

The evidence status is:

- **Exact analytical result under the stated assumptions:** the two allowed ordered flavor projections cancel.
- **Numerical-computation result before correction:** with the default complex $\epsilon_K$, the vector builders left a floating-point residual at machine cancellation scale; with $\epsilon_K=0$, the evaluated components were exactly zero in the representative calculation.
- **Diagnosed numerical root cause:** `epsabs=0` asked nested adaptive quadrature to satisfy a relative-only criterion on roundoff noise around an analytically zero vector integrand.
- **Production status:** corrected on 2026-09-26. Vector contributions with a common spectator, pair, pole mass, and pole width now have their flavor/order coefficients aggregated before kinematics. A group whose aggregate is zero at the existing scale-relative roundoff criterion is removed. Empty vector-pole and matching models return exact `0.0j` before Dalitz or pole evaluation.

## Conventions, units, assumptions, and inputs

The neutral-kaon convention and bra-level conjugation are those in [the neutral-kaon projection record](neutral-kaon-amplitude-projection.md):

$$
|K_S\rangle=p|K^0\rangle+q|\bar K^0\rangle,
\qquad
|K_L\rangle=p|K^0\rangle-q|\bar K^0\rangle
$$

with $p=(1+\epsilon_K)/\sqrt{2(1+|\epsilon_K|^2)}$ and $q=(1-\epsilon_K)/\sqrt{2(1+|\epsilon_K|^2)}$ [@ParticleDataGroup:2024cfk]. The default value is the one recorded in [the numerical-width record](../alp/three-meson/alp-three-meson-widths-numerical.md). Natural units are used. Masses and decay constants are in GeV, Dalitz invariants are in $\mathrm{GeV}^2$, invariant amplitudes are dimensionless, and the bounded inner integral of $|\mathcal M|^2ds_{23}$ is in $\mathrm{GeV}^2$.

The representative inputs are explicit assumptions, not a phenomenological benchmark:

$$
m_a=1.6\ \mathrm{GeV},\qquad
f_a=100\ \mathrm{GeV},\qquad
c_q^A=\operatorname{diag}(0.8,-0.3,0.15),\qquad
c_q^P=0.
$$

The ChPT conventions and fitted pseudoscalar inputs follow [the numerical-width record](../alp/three-meson/alp-three-meson-widths-numerical.md) and its Alda:2026nrd source [@Alda:2026nrd]. The pole and matching amplitudes use the assumptions listed in [the vector-exchange record](../alp/three-meson/a3p-vector-exchange-amplitude.md): current constant-width pole inputs, ideal vector mixing, and $G_V=f_\pi/\sqrt2$. The local term required by tensor-to-Proca matching is kept separate from the pole [@Ecker:1989yg]. No approximation to $\epsilon_K$ is made. The only flavor assumption used for the exact zero is that the relevant strangeness-zero amplitude is Bose/order symmetric,

$$
\mathcal A_{0\bar0}=\mathcal A_{\bar00},
$$

with the equal neutral-kaon mass used by the implementation. Components with $K^0K^0$ or $\bar K^0\bar K^0$ are excluded by the strangeness-conserving flavor structure for the diagonal representative coupling. The implementation additionally aggregates only terms with identical numerical kernels: the spectator and pair labels and the complete `VectorPoleParameters` value must agree. Thus independently overridden $K^{*0}$ and $\bar K^{*0}$ masses or widths are not silently treated as degenerate. The roundoff classification uses $|\sum_i c_i|\leq64\epsilon_{\rm mach}\sum_i|c_i|$; this is an explicit numerical representation rule for an analytically established cancellation, not a fixed physical-amplitude cutoff.

## Stepwise analytical derivation

1. The final-state bra coefficients on the $K_L$ leg are $(p^*,-q^*)$ and those on the $K_S$ leg are $(p^*,q^*)$.
2. The $K^0\bar K^0$ ordering therefore carries $p^*q^*$, while the $\bar K^0K^0$ ordering carries $-q^*p^*$.
3. SymPy constructs the projection rather than transcribing its result:

   $$
   \mathcal A(\pi^0K_LK_S)
   =p^*q^*\mathcal A_{0\bar0}
   -q^*p^*\mathcal A_{\bar00}.
   $$

4. Substitution of $\mathcal A_{0\bar0}=\mathcal A_{\bar00}=\mathcal A_{\rm sym}$ and symbolic simplification gives exactly

   $$
   \boxed{\mathcal A(\pi^0K_LK_S)=0}.
   $$

This cancellation is independent of the numerical value or phase of $\epsilon_K$: the same product $p^*q^*$ multiplies the two opposite signs. It is not an $\epsilon_K$ suppression approximation.

## Numerical method and outputs

At the deterministic interior point

$$
s_{12}=0.8892238342761362\ \mathrm{GeV}^2,
\qquad
s_{23}=1.5368660756870405\ \mathrm{GeV}^2,
$$

the script builds the current ChPT and corrected vector-pole and vector-matching callables for both $\pi^0K_LK_S$ and $\pi^0K^0\bar K^0$. It prints every retained vector exchange with spectator, pair, flavor factor, and term. It also prints the two explicitly projected flavor orderings, sector component sums, sums of absolute term magnitudes, and cancellation residuals.

For the noncancelled flavor channel, the representative coherent component sum is approximately

$$
\mathcal M(\pi^0K^0\bar K^0)
=2.0004\times10^{-2}-1.3553\times10^{-2}i,
$$

whereas both corrected physical vector builders retain zero contributions and return exactly zero. Regression mocks prove that the pair-invariant construction and pole denominator are not evaluated; scalar type and finiteness validation is preserved. The noncancelled flavor channel still retains two contributions in each vector model and gives the same representative nonzero values. The current direct ChPT physical-state construction and the explicitly projected ChPT flavor orderings also return exact zero at this point. The exact printed values are authoritative because they are generated from the current working tree.

## Bounded quadrature diagnosis

The script evaluates only one inner $s_{23}$ integral at fixed $s_{12}$, with `limit=24` and `epsrel=1e-6`. This bound avoids launching the nested production width. Before correction, the default complex $\epsilon_K$ and `epsabs=0` exhausted all 24 subintervals on an integral of order $10^{-35}\,\mathrm{GeV}^2$. With the current package implementation, the physical-channel integrand is exactly zero and QUADPACK converges on its initial 21-point panel, both with `epsabs=0` and with the stated diagnostic absolute tolerance. The same holds for the explicit $\epsilon_K=0$ case. The noncancelled flavor-channel integral used to define the diagnostic scale is $1.8222312564801699\times10^{-3}\,\mathrm{GeV}^2$ in the direct invocation.

The diagnostic absolute tolerance is chosen as $10^{-12}$ times the corresponding noncancelled flavor-channel integral scale. This factor is a stated diagnostic precision assumption, not a universal production tolerance. A production tolerance should instead be tied to the requested absolute width accuracy and a channel/component scale.

These checks exclude threshold closure, absent underlying flavor couplings, and a genuinely tiny physical vector effect: phase space is open, the definite-flavor channel and its vector terms remain nonzero, and only the physical projection vanishes. Finite constant pole widths also keep the representative noncancelled integrand away from a real-axis pole singularity.

## Implemented correction and verification

[`vector_meson_exchange.py`](../../scripts/implementation/vector_meson_exchange.py) now groups contributions by their complete kinematic kernel identity and aggregates flavor/order coefficients before evaluating invariants, numerators, local matching functions, or propagators. Exactly cancelled groups are removed while nonzero groups retain their original labeled contribution metadata. This is channel-independent; there is no $\pi^0K_LK_S$ name check. The pole and matching callables short-circuit empty contribution tuples to exact zero. Regression tests in [`test_vector_meson_exchange.py`](../../scripts/implementation/test_vector_meson_exchange.py) verify exact zero and patch the Dalitz helper to fail if kinematics are reached; they also verify that $\pi^0K^0\bar K^0$ remains nonzero.

No production `epsabs` policy was changed. A nonzero scale-aware `epsabs` remains a general integration-robustness recommendation, but it is not a substitute for removing this established exact vector cancellation.

Modified paths are:

- `scripts/implementation/vector_meson_exchange.py`;
- `scripts/implementation/test_vector_meson_exchange.py`;
- `scripts/diagnose_pi0_kl_ks_vector_width.py`;
- `memory/pi0-kl-ks-vector-width-slowdown-diagnostic.md`.

Run:

```bash
.venv/bin/python scripts/diagnose_pi0_kl_ks_vector_width.py
git diff --check
```

The diagnostic asserts the exact SymPy cancellation, empty physical vector models and pre-kinematic exact-zero calls for default complex $\epsilon_K$, retained nonzero definite-flavor models, exact numerical zero for $\epsilon_K=0$, and one-panel convergence of the bounded integral. Package-form test and diagnostic results are generated with the repository `.venv`; the direct script invocation uses only the repository-root path insertion documented above. The last recorded full run, on 2026-09-28, had all 74 implementation tests and the direct diagnostic passing; this task changes documentation only and does not relabel that run as a 2026-09-30 validation.
