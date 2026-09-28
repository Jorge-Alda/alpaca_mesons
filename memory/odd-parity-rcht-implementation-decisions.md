---
title: Decision-ready baseline for odd-parity RChT implementation
status: mixed
created: 2026-09-28
updated: 2026-09-28
bibtex_keys: [Alda:2026bar, Bai:2025fvl, Fujiwara:1984mp, Ruiz-Femenia:2003jdx, Roig:2013baa, Kaiser:2000gs, Herrera-Siklody:1996tqr, Ecker:1988te, Ecker:1989yg, Rosell:2004mn, Pich:2010sm, Hanhart:2012wi, GomezDumm:2013sib, Jamin:2000wn, ParticleDataGroup:2024cfk]
provenance: assumption
---

# Odd-Parity RChT Implementation Decisions

## Scope and evidence

This record turns the evidence and options in [the $U(3)$ nonet record](u3-nonet-odd-parity-extension-options.md) and [the vector line-shape record](vector-width-lineshape-prescriptions.md) into an implementation baseline. Literature fixes the Ruiz--Femenia operator basis, OZI hierarchy, Ecker antisymmetric-tensor convention, and line-shape constraints; the selected truncations below are declared recommendations, not additional literature results [@Ruiz-Femenia:2003jdx; @Kaiser:2000gs; @Ecker:1989yg]. **HLS remains excluded.** Natural units and Alda notation are used throughout [@Alda:2026bar].

> **Recommended baseline**
>
> 1. Use a defined leading-OZI, single-trace $U(3)$ nonet model based on the seven $VJP$ and four $VVP$ Ruiz--Femenia operators, with common $c_i,d_i$ and all independent double traces set explicitly to zero. Common singlet/octet couplings are an ansatz, not a complete $U(3)$ result.
> 2. Use $T_0=\mathbf1/\sqrt3$, $\operatorname{tr}(T_0^2)=1$, for the vector singlet and $\lambda_0=\sqrt{2/3}\,\mathbf1$, $\operatorname{tr}(\lambda_0^2)=2$, for Alda's pseudoscalar singlet.
> 3. Use ideal $\omega$--$\phi$ mixing and physical vector masses only as declared first-pass approximations.
> 4. Treat a vector as external only in a justified narrow-width approximation (NWA), sufficiently far from daughter thresholds. Then use its LSZ matrix element and no Breit--Wigner factor. Broad $\rho/K^*$ states or channels near daughter threshold require an inclusive spectral daughter observable.
> 5. For internal vectors use physical thresholds and energy-dependent P-wave widths for $\rho$ and $K^*$; use a narrow constant $\omega$ width unless the $3\pi$ line shape is resolved; for $\phi$ near its pole resolve charged and neutral $K\bar K$ widths and include $3\pi$ at precision level.
> 6. In pure antisymmetric-tensor RChT, use WZW with the declared Bardeen-type counterterm convention to fix the anomalous variation, and use $c_i,d_i$ plus any explicitly retained invariant local odd-parity operators for resonance interactions. Do not infer those invariant terms from the anomaly or add a free hadronic $c_{AB}$ by default.
> 7. In the NWA encode the Ecker external-state LSZ matrix element with pole mass and residue, using tree-level canonical $Z_V=1$ only at the stated baseline accuracy.
> 8. Do not invent numerical fit defaults. Dataset selection and covariance/theory-error modeling remain the genuine numerical decision.

## Nonet and OZI decision

“Singlet/double trace” is not one prescription. The baseline is a defined leading-single-quark-loop/OZI model: continue the existing single traces to the full nonet, use the two normalized singlet generators above, take common $c_i,d_i$ for octet and singlet directions, and set every independent double-trace coefficient to zero. Generic extra flavor traces correspond to additional quark-loop topologies and OZI-suppressed coefficient counting, with one nominal $1/N_c$ suppression per extra trace [@Kaiser:2000gs]. This trace counting is separate from the large-$N_c$ counting of powers and derivatives of the anomalous source $X_{TC}$. Common singlet/octet $c_i,d_i$ are therefore a model ansatz and operator truncation, not a symmetry theorem. Ideal neutral-vector mixing and physical masses are likewise phenomenological baseline choices. This model is sufficient for first $aV\gamma$ and $aVV$ predictions and is not a global architecture blocker.

Upgrade when fitting $\eta'$, $\omega/\phi$, or other OZI-sensitive data at subleading $1/N_c$ accuracy. Concrete additions include the breaking operator $i\epsilon_{\mu\nu\rho\sigma}\langle V^{\mu\nu}V^{\rho\sigma}\rangle\langle\chi_-\rangle$ alongside its single-trace counterpart, and trace-factorized derivative structures such as $\epsilon_{\mu\nu\rho\sigma}\langle V^{\mu\nu}f_+^{\rho\alpha}\rangle\langle\nabla_\alpha u^\sigma\rangle$ or analogous $VVP$ factorizations. A complete $U(3)$ upgrade must additionally enumerate parity-allowed coefficient functions of $X_{TC}$ and operators with derivatives of $X_{TC}$, then reduce the candidate set using integration by parts, equations of motion, and $U(3)$ trace identities. These examples are upgrade templates, not a claim of a complete independent basis [@Roig:2013baa; @Herrera-Siklody:1996tqr].

The source conventions must not be identified by name alone. Herrera-Siklody et al. define

$$
X_{\rm HS}=\langle\log U\rangle+\hat\theta,
\qquad \hat\theta=i\theta_{\rm HS},
$$

whereas Alda Eq. (III.14), with $U=u^2$ on a common logarithm branch, gives

$$
X_{TC}=2i\left(\langle\log u\rangle+i c_G\frac{a}{f_a}\right)
=i\langle\log U\rangle-2c_G\frac{a}{f_a}.
$$

Consequently the exact algebraic relation is

$$
X_{TC}=iX_{\rm HS}+\theta_{\rm HS}-2c_G\frac{a}{f_a}.
$$

Only under the explicit source identification $\theta_{\rm HS}=2c_Ga/f_a$ does this reduce to $X_{TC}=iX_{\rm HS}$, equivalently $X_{\rm HS}=-iX_{TC}$. Because the two papers' underlying source normalizations and signs have not been matched in this repository, uses of Herrera-Siklody source functions are translated by an explicitly declared convention rather than a claimed universal equality [@Herrera-Siklody:1996tqr; @Alda:2026bar]. $X_{TC}$-dependent coefficients can contribute to dynamical ALP amplitudes after expanding $X_{TC}$ in the ALP field. Their omission here is an explicit large-$N_c$/operator truncation, not a consequence of topological observables being absent.

## Width decision by use case

An external vector receives no Breit--Wigner factor. Reported two-body $a\to V\gamma$ and $a\to VV$ predictions may use on-shell vectors only when the state is narrow enough and the kinematics are sufficiently far from its daughter threshold for an NWA to be controlled. A broad $\rho$ or $K^*$, or any vector produced near daughter threshold, should instead be reported through the inclusive daughter channel with a normalized spectral integral. Widths are therefore not blockers for the odd operator basis, but they can be mandatory for the physical observable even when a formal on-shell two-body amplitude exists.

For an internal $\rho$ or $K^*$, the baseline is an energy-dependent P-wave width with physical $\pi\pi$ or $K\pi$ thresholds. For an internal $\omega$, a narrow constant width is acceptable unless the calculation resolves the $3\pi$ line shape. For an internal $\phi$ near its pole, use separate charged- and neutral-$K\bar K$ energy-dependent widths and add $3\pi$ when precision requires it. A constant-width Breit--Wigner is only a compatibility or local-fit alternative; Dyson-resummed or dispersive treatments define precision systematics rather than the baseline [@Rosell:2004mn; @Pich:2010sm; @Hanhart:2012wi; @GomezDumm:2013sib].

After any dressing, verify $k_\mu\mathcal M^\mu=0$ for every external photon and for the complete amplitude, not only an isolated graph. Also verify the low-energy/anomaly normalization and that contact, vector--photon-mixing, and exchange terms use one consistent dressing prescription.

## Anomaly coefficient and external tensor

The numerical $c_{AB}$ is not needed by the selected pure-RChT baseline. WZW, together with the chosen Bardeen-type counterterm convention, fixes the anomalous variation and current assignment; it does not fix invariant local odd-parity operators [@Bai:2025fvl; @Fujiwara:1984mp]. The $c_i,d_i$ and any explicitly retained local invariant terms remain independent homogeneous couplings. The generic background-source helper retains $c_{AB}$ because that convention can be useful, but it requires explicit caller input and must not define hadronic resonance couplings by default.

For a stable vector, or an unstable vector treated as external within a controlled NWA, the Ecker kinetic convention gives the LSZ matrix element

$$
\langle0|V_{\mu\nu}(0)|V(p,\lambda)\rangle
=i\frac{\sqrt{Z_V}}{M_{\rm pole}}
\left(p_\mu\varepsilon_\nu^{(\lambda)}-p_\nu\varepsilon_\mu^{(\lambda)}\right),
$$

up to a consistently fixed one-particle phase [@Ecker:1988te; @Ecker:1989yg]. Tree-level canonical normalization sets $Z_V=1$; beyond tree level the pole residue and generally complex pole must be defined by the chosen propagator scheme. Encoding this factor and deriving physical amplitudes are implementation tasks, not requests for an arbitrary user prescription. Broad vectors are not promoted to exact external asymptotic states by this formula.

## Near-degeneracy domain

Alda Eq. (III.37) requires

$$
|m_a^2-m_i^2|\gg\frac{f_\pi}{f_a}
\left|\operatorname{tr}(C\lambda_i^\dagger)\right|
\approx\frac{f_\pi}{f_a}m_i^2.
$$

The approximate right-hand side shows that the excluded bands scale as $f_\pi/f_a$ and are negligible in generic broad scans for $f_a\gg f_\pi$. The actual pointwise test must nevertheless use the exact overlap $|\operatorname{tr}(C\lambda_i^\dagger)|$. A finite $\mathbb M_a^{-1}(C)$ is necessary but not sufficient: every sampled point that violates the inequality is locally outside the perturbative-mixing domain and must be masked or warned, or treated with a coupled mixing calculation. This local exception is not a global architecture blocker, and matrix magnitude alone is not the validity test. For stable or narrow pseudoscalars, the replacement is exact coupled diagonalization of the relevant quadratic system. For unstable pseudoscalars, validity is controlled by separation of the complex poles relative to the off-diagonal mixing, and the replacement is a coupled complex propagator/pole analysis rather than an ordinary Hermitian exact diagonalization [@Alda:2026bar].

## Calibration decision

The minimal calibration set is:

- $V\to P\gamma$ radiative widths, constraining direct-plus-exchange $VJP/VVP$ combinations;
- $P\to\gamma\gamma$ widths and transition-form-factor normalizations/slopes, constraining WZW/counterterm normalization and photon--vector contributions;
- $e^+e^-\to P\gamma$ or $VP$ line shapes, constraining timelike combinations, phases, and internal-vector prescriptions;
- $\omega/\phi\to3\pi$, constraining odd three-meson amplitudes, direct terms, and the $\omega/\phi$ line-shape sector.

Use PDG correlations when available. Otherwise document an independent-Gaussian experimental approximation and add a separately stated theory covariance for truncation, line-shape, and symmetry-breaking uncertainty [@ParticleDataGroup:2024cfk]. No central values, covariance defaults, or priors are selected here.

## Blocker classification

**Genuine blocker for numerical predictions:** the calibration dataset, its experimental covariance, and the theory-covariance/prior model for the identifiable $c_i,d_i$ combinations and any promoted nuisance parameters.

**Implementation tasks, not physics-choice blockers:** on-shell amplitude reduction; the fixed external-tensor normalization; internal-vector line-shape code by the use cases above; spectral daughter distributions; and complete-amplitude Ward-identity tests.

**Nonblockers at baseline:** independent double traces and singlet-specific couplings within the declared model; nonideal $\omega$--$\phi$ mixing; a free hadronic $c_{AB}$; precision Dyson/dispersive line shapes; and generic construction of $\mathbb M_a^{-1}(C)$ away from the local Alda degeneracy domain. Pointwise validity masking/warnings and coupled treatment of any intentionally retained violating points remain required implementation work.
