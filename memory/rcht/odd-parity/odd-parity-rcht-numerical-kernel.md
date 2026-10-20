---
title: Non-HLS odd-parity RChT numerical kernel
status: numerical
created: 2026-09-28
updated: 2026-09-30
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx, Ecker:1989yg]
provenance: numerical-computation
---

# Non-HLS Odd-Parity RChT Numerical Kernel

## Evidence, scope, and assumptions

[`scripts/implementation/odd_parity_rcht.py`](../../../scripts/implementation/odd_parity_rcht.py) implements the operator-level proposal derived in the [two-source ALP RChT record](two-source-alp-rcht-odd-parity-extension.md), the [direct VJP insertion record](odd-parity-coupling-prescriptions-and-direct-alp-vgamma.md), and the [exact photon--vector diagonalization](../foundations/rcht-photon-vector-diagonalization.md). The $VJP/VVP$ basis and five short-distance relations are literature inputs [@Ruiz-Femenia:2003jdx]; the ALP spurion insertion follows Alda:2026nrd's conventions [@Alda:2026nrd]. The numerical organization and validation policy are implemented results, not additional physics claims.

Natural units are used. $e,c_i,d_i,\mathbb M_a^{-1}(C),c_q^{AP}$ and generators are dimensionless; $f_a,M_V,F_V,F$ have mass dimension one, while `b0_mq` has mass dimension two. `b0_mq` denotes $B_0m_q$, so no second factor of $B_0$ is inserted. Flavor traces retain the repository normalization. Source generators may be non-Hermitian (for charged fields), while $\mathbb M_a^{-1}(C),c_q^{AP},B_0m_q$ are finite Hermitian $3\times3$ matrices.

Assumptions are explicit: tree level, one resonance multiplet, the Goldstone vacuum, $f_-=0$, and electromagnetic-compatible ALP spurions. The eleven $c_i,d_i$ are finite real inputs. No on-shell reduction, integration by parts, fitted coupling value, width prescription, or uncertainty is inferred. **HLS is intentionally excluded by user instruction.** Ideal $\omega$--$\phi$ mixing is used only by the labeled `ideal_vector_generators()` convenience function and by the pre-existing vector-exchange implementation.

## Stepwise method and outputs

1. `ALPFlavorSpurions` copies and freezes the three Hermitian inputs. For source matrices $T_A,T_B$, `odd_parity_flavor_kernels` computes without assuming source Hermiticity
   $$
   D=\operatorname{tr}(\{T_A,T_B\}\mathbb M_a^{-1}(C)),\qquad
   M=\operatorname{tr}(\{T_A,T_B\}
   \{c_q^{AP}+\mathbb M_a^{-1}(C),B_0m_q\}).
   $$
2. The seven unreduced VJP coefficients are
   $$
   \frac1{f_a}(4eD,4eD,-8eM,0,4eD,4eD,4eD),
   $$
   and weighting multiplies component $i$ by $c_i/M_V$. The four distinct VVP Lorentz structures remain separate as
   $$
   \frac1{f_a}(-2D,4M,-2D,-2D),
   $$
   with componentwise $d_i$ weighting. In particular, $d_2$ remains when $\mathbb M_a^{-1}(C)=0$ but $c_q^{AP}\ne0$.
3. `anomaly_seed` requires the caller to state whether sources are identical and evaluates
   $$
   c_{AB}+sN_c\operatorname{tr}(\{T_A,T_B\}
   (c_q^A+\mathbb M_a^{-1}(C))),\qquad
   s=1\ \text{(identical)},\quad s=2\ \text{(distinct)}.
   $$
   `photon_anomaly_seed` evaluates $c_\gamma+2N_c\operatorname{tr}(Q^2(c_q^A+\mathbb M_a^{-1}(C)))$ directly.
4. `ope_linear_system` returns read-only $A,b$ for $A(c_1,\ldots,c_7,d_1,\ldots,d_4)^T=b$, using exactly the five Ruiz--Femenia relations recorded in the [coupling prescription](odd-parity-coupling-prescriptions-and-direct-alp-vgamma.md). `ope_residuals` returns $A\theta-b$. `OPEPrior` accepts only a symmetric positive-definite $5\times5$ covariance, and `ope_prior_chi2` evaluates the quadratic form by a SciPy Cholesky triangular solve. The user-authorized `fixed_ope_prior` builder implements the explicitly assumed Gaussian scale model in [the fit-infrastructure record](correlated-odd-parity-fit-infrastructure.md), evaluated once from caller-supplied reference `OPEParameters`; it is not parameter-dependent during fitting.
5. For real $\boldsymbol\epsilon$ and a real orthogonal rotation $O$, `photon_vector_transform` computes exactly
   $$
   Z=1+\boldsymbol\epsilon^T\boldsymbol\epsilon,\quad
   K_V^{1/2}=I-\frac{\boldsymbol\epsilon\boldsymbol\epsilon^T}{\sqrt Z(\sqrt Z+1)},
   $$
   $$
   H=K_V^{1/2}O,\qquad \boldsymbol\rho=O^T\boldsymbol\epsilon/\sqrt Z,
   \qquad
   S=\begin{pmatrix}Z^{-1/2}&-\boldsymbol\rho^T\\Z^{-1/2}\boldsymbol\epsilon&H\end{pmatrix}.
   $$
   `transform_anomalous_bilinear` requires a real symmetric $\Lambda_{\rm old}$ and returns $S^T\Lambda_{\rm old}S$ without silently symmetrizing.

There is **no automatic conversion from unordered $C_{AB}$ coefficients to $\Lambda$**. Diagonal multiplicities, source couplings, and tensor-to-Proca normalization must be supplied before calling the congruence, as explained in the [two-source record](two-source-alp-rcht-odd-parity-extension.md).

## Verification and invocation

Run exactly:

```text
.venv/bin/python -m unittest scripts.implementation.test_odd_parity_rcht
.venv/bin/python -m unittest scripts.implementation.test_vector_meson_exchange
.venv/bin/python -m unittest discover -s scripts/implementation -t .
git diff --check
```

The deterministic tests compare every coefficient directly, check seeded unitary flavor covariance and axial reparametrization invariance, expose the OPE null directions, compare the Cholesky result with a manual diagonal-covariance result, verify congruence and quadratic-form invariance, test zero-mixing blocks, and check ideal-generator normalization and charge conjugation. Validation tests cover malformed, nonfinite, non-Hermitian, nonpositive, nonreal, nonorthogonal, nonsymmetric, and non-SPD inputs. Existing vector-exchange tests verify that replacing its duplicated private generators with this public source causes no behavior change.

## Decisions, implementation tasks, and blocker

The [implementation decision record](odd-parity-rcht-implementation-decisions.md) fixes the first-prediction baseline: leading-OZI single-trace $U(3)$ continuation, common $c_i,d_i$, ideal $\omega$--$\phi$ mixing and physical masses as declared approximations, narrow on-shell external vectors, and the Ecker antisymmetric-tensor external-state normalization. Independent double traces, singlet-specific operators, and precision line shapes are upgrades triggered by the target dataset, not global blockers. The generic background-source anomaly helper retains caller-supplied $c_{AB}$, but the pure-RChT baseline uses WZW/counterterm normalization plus $c_i,d_i$ and does not require a free hadronic $c_{AB}$.

Alda:2026nrd Eq. (III.37) requires

$$
|m_a^2-m_i^2|\gg\frac{f_\pi}{f_a}
\left|\operatorname{tr}(C\lambda_i^\dagger)\right|
\approx\frac{f_\pi}{f_a}m_i^2.
$$

The excluded real-mass neighborhoods scale as $f_\pi/f_a$ and are generally negligible in measure or scan resolution for $f_a\gg f_\pi$, but Eq. (III.37) is a pointwise validity condition. A finite caller-supplied `m_a_inverse_C` is necessary but not sufficient. The calling scan must evaluate the exact overlap at every point and mask/flag a violation or replace the perturbative treatment by a coupled mixing calculation. Matrix magnitude alone is not the criterion. For unstable states, closeness of the real mass parameters is distinct from closeness of the complex poles; a retained near-pole point requires a coupled complex-propagator analysis. This operator kernel is not passed the pseudoscalar spectrum and therefore cannot perform that pointwise validation itself [@Alda:2026nrd].

> **The remaining blocker is validated physical data ingestion and observable predictors for the calibration channels. The fixed-prior generic fit layer is ready, but no physical $c_i,d_i$ fit has been run.**

Physical on-shell reduction, external-tensor normalization in amplitudes, width models for internal vectors, Ward-identity tests, and optional spectral conversion of unstable external-vector channels are implementation tasks. This kernel itself continues to supply operator-level coefficients rather than physical widths or fitted couplings.

For the current distinction between implemented ALP insertions and missing
complete ALP amplitudes, use the 2026-10-13
[target audit](alp-target-amplitude-audit.md). New projected null checks do not
promote this operator-level kernel into a normalized physical predictor.
