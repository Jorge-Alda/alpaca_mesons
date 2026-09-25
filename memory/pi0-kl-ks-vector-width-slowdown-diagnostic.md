---
title: Diagnostic record for the pi0 K_L K_S vector-width slowdown
status: mixed
created: 2026-09-25
updated: 2026-09-25
bibtex_keys: [Alda:2026bar, Ecker:1989yg, ParticleDataGroup:2024cfk]
provenance: analytical-derivation
---

# Diagnostic for the $\pi^0K_LK_S$ Vector-Width Slowdown

## Scope and evidence status

This record diagnoses the established slowdown without changing the production implementation or tests. The executable calculation is [`scripts/diagnose_pi0_kl_ks_vector_width.py`](../scripts/diagnose_pi0_kl_ks_vector_width.py). It imports the current implementation as top-level modules from `scripts/implementation/`, matching that implementation's present import arrangement, and never launches the nested production width integral.

The evidence status is:

- **Exact analytical result under the stated assumptions:** the two allowed ordered flavor projections cancel.
- **Numerical-computation result:** with the default complex $\epsilon_K$, the current vector builders leave only a floating-point residual at machine cancellation scale; with $\epsilon_K=0$, the evaluated components are exactly zero in this representative calculation.
- **Diagnosed numerical root cause:** `epsabs=0` asks nested adaptive quadrature to satisfy a relative-only criterion on roundoff noise around an analytically zero integrand.
- **Production status:** no production fix has yet been applied.

## Conventions, units, assumptions, and inputs

The neutral-kaon convention and bra-level conjugation are those in [the neutral-kaon projection record](neutral-kaon-amplitude-projection.md):

$$
|K_S\rangle=p|K^0\rangle+q|\bar K^0\rangle,
\qquad
|K_L\rangle=p|K^0\rangle-q|\bar K^0\rangle
$$

with $p=(1+\epsilon_K)/\sqrt{2(1+|\epsilon_K|^2)}$ and $q=(1-\epsilon_K)/\sqrt{2(1+|\epsilon_K|^2)}$ [@ParticleDataGroup:2024cfk]. The default value is the one recorded in [the numerical-width record](alp-three-meson-widths-numerical.md). Natural units are used. Masses and decay constants are in GeV, Dalitz invariants are in $\mathrm{GeV}^2$, invariant amplitudes are dimensionless, and the bounded inner integral of $|\mathcal M|^2ds_{23}$ is in $\mathrm{GeV}^2$.

The representative inputs are explicit assumptions, not a phenomenological benchmark:

$$
m_a=1.6\ \mathrm{GeV},\qquad
f_a=100\ \mathrm{GeV},\qquad
c_q^A=\operatorname{diag}(0.8,-0.3,0.15),\qquad
c_q^P=0.
$$

The ChPT conventions and fitted pseudoscalar inputs follow [the numerical-width record](alp-three-meson-widths-numerical.md) and its Alda source [@Alda:2026bar]. The pole and matching amplitudes use the assumptions listed in [the vector-exchange record](a3p-vector-exchange-amplitude.md): current constant-width pole inputs, ideal vector mixing, and $G_V=f_\pi/\sqrt2$. The local term required by tensor-to-Proca matching is kept separate from the pole [@Ecker:1989yg]. No approximation to $\epsilon_K$ is made. The only flavor assumption used for the exact zero is that the relevant strangeness-zero amplitude is Bose/order symmetric,

$$
\mathcal A_{0\bar0}=\mathcal A_{\bar00},
$$

with the equal neutral-kaon mass used by the implementation. Components with $K^0K^0$ or $\bar K^0\bar K^0$ are excluded by the strangeness-conserving flavor structure for the diagonal representative coupling.

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

the script builds the current ChPT, vector-pole, and vector-matching callables for both $\pi^0K_LK_S$ and $\pi^0K^0\bar K^0$. It prints every retained vector exchange with spectator, pair, flavor factor, and term. It also prints the two explicitly projected flavor orderings, sector component sums, sums of absolute term magnitudes, and cancellation residuals.

For the noncancelled flavor channel, the representative coherent component sum is approximately

$$
\mathcal M(\pi^0K^0\bar K^0)
=2.0004\times10^{-2}-1.3553\times10^{-2}i,
$$

whereas the physical $\pi^0K_LK_S$ sum is at the $10^{-18}$ level despite individual retained projected/vector terms at ordinary nonzero scales. The script asserts that the physical residual is within $512\epsilon_{\rm mach}$ times the absolute-term sum and that the flavor amplitude is nonzero. The exact printed values are authoritative because they are generated from the current working tree.

## Bounded quadrature diagnosis

The script evaluates only one inner $s_{23}$ integral at fixed $s_{12}$, with `limit=24` and `epsrel=1e-6`. This bound prevents the minutes-long nested production calculation. For the default complex $\epsilon_K$ and `epsabs=0`, QUADPACK exhausts all 24 subintervals while trying to resolve an integral of order $10^{-35}\,\mathrm{GeV}^2$. The same residual converges immediately when supplied a nonzero diagnostic absolute tolerance. In the explicit $\epsilon_K=0$ case, the current arithmetic evaluates the amplitude and integrand as exact numerical zero, so even `epsabs=0` terminates after the initial panel.

The diagnostic absolute tolerance is chosen as $10^{-12}$ times the corresponding noncancelled flavor-channel integral scale. This factor is a stated diagnostic precision assumption, not a universal production tolerance. A production tolerance should instead be tied to the requested absolute width accuracy and a channel/component scale.

These checks exclude threshold closure, absent couplings or flavor traces, a genuinely tiny physical vector effect, and a need for more subdivisions as explanations: phase space is open, the flavor channel and retained vector terms are nonzero, the exact projection is zero, and subdivision only attempts to resolve floating noise. Finite constant pole widths also keep this representative integrand away from a real-axis pole singularity.

## Recommended correction and verification

The preferred correction is to aggregate the two flavor orderings before numerical evaluation and return an exact zero when their symbolic/flavor coefficient cancels. This avoids manufacturing a residual that quadrature can mistake for signal. Independently, adaptive integration should use a nonzero, scale-aware `epsabs` in addition to `epsrel`; this is a robustness requirement, not a substitute for the exact flavor cancellation. No such production correction is made in this diagnostic record.

Run:

```bash
.venv/bin/python scripts/diagnose_pi0_kl_ks_vector_width.py
git diff --check
```

The script asserts the exact SymPy cancellation, machine-scale numerical cancellation for default complex $\epsilon_K$, exact numerical zero for $\epsilon_K=0$ at the representative point, bounded failure of relative-only quadrature on the residual, and convergence with a nonzero scale-aware absolute tolerance.
