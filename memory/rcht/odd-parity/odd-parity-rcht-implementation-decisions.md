---
title: Decision-ready baseline for odd-parity RChT implementation
status: mixed
created: 2026-09-28
updated: 2026-09-30
bibtex_keys: [Alda:2026nrd, Bai:2025fvl, Fujiwara:1984mp, Ruiz-Femenia:2003jdx, Roig:2013baa, Kaiser:2000gs, Herrera-Siklody:1996tqr, Ecker:1988te, Ecker:1989yg, Rosell:2004mn, Pich:2010sm, Hanhart:2012wi, GomezDumm:2013sib, Jamin:2000wn, ParticleDataGroup:2024cfk]
provenance: assumption
---

# Odd-Parity RChT Implementation Decisions

## Scope and evidence

This record turns the evidence and options in [the $U(3)$ nonet record](u3-nonet-odd-parity-extension-options.md) and [the vector line-shape record](../foundations/vector-width-lineshape-prescriptions.md) into an implementation baseline. Literature fixes the Ruiz--Femenia operator basis, OZI hierarchy, Ecker antisymmetric-tensor convention, and line-shape constraints; the selected truncations below are declared recommendations, not additional literature results [@Ruiz-Femenia:2003jdx; @Kaiser:2000gs; @Ecker:1989yg]. **HLS remains excluded.** Natural units and Alda notation are used throughout [@Alda:2026nrd].

> **Recommended baseline**
>
> 1. Use a defined leading-OZI, single-trace $U(3)$ nonet model based on the seven $VJP$ and four $VVP$ Ruiz--Femenia operators, with common $c_i,d_i$ and all independent double traces set explicitly to zero. This is trace-topology bookkeeping, not a proof that every omitted multi-trace coefficient has the same suppression. Common singlet/octet couplings are an ansatz, not a symmetry consequence or complete $U(3)$ result.
> 2. Use $T_0=\mathbf1/\sqrt3$, $\operatorname{tr}(T_0^2)=1$, for the vector singlet and $\lambda_0=\sqrt{2/3}\,\mathbf1$, $\operatorname{tr}(\lambda_0^2)=2$, for Alda's pseudoscalar singlet.
> 3. Use ideal $\omega$--$\phi$ mixing and physical vector masses only as declared first-pass approximations.
> 4. Treat a vector as external only in a justified narrow-width approximation (NWA), sufficiently far from daughter thresholds. Then use its LSZ matrix element and no Breit--Wigner factor. Broad $\rho/K^*$ states or channels near daughter threshold require an inclusive spectral daughter observable.
> 5. For internal vectors use physical thresholds and energy-dependent P-wave widths for $\rho$ and $K^*$; use a narrow constant $\omega$ width unless the $3\pi$ line shape is resolved; for $\phi$ near its pole resolve charged and neutral $K\bar K$ widths and include $3\pi$ at precision level.
> 6. In pure antisymmetric-tensor RChT, use WZW with the declared external-source/current assignment and Bardeen-type counterterm convention to reproduce the corresponding anomalous variation. This does not uniquely fix finite invariant local terms or resonance matching. Use $c_i,d_i$ plus any explicitly retained invariant local odd-parity operators for resonance interactions; do not infer those terms from the anomaly or add a free hadronic $c_{AB}$ by default.
> 7. In the NWA encode the Ecker external-state matrix element with a declared real on-shell mass and residue prescription, using the tree-level canonical limit $Z_V=1$ only at the stated baseline accuracy; do not present an unstable complex pole as an exact LSZ state.
> 8. Use the user-authorized fixed OPE theory prior below. This is an explicit assumption, not a literature result. Dataset ingestion and physical observable predictors remain missing.

## Nonet and OZI decision

“Singlet/double trace” is not one prescription. The baseline is a defined leading-single-quark-loop/OZI model: continue the existing single traces to the full nonet, use the two normalized singlet generators above, take common $c_i,d_i$ for octet and singlet directions, and set every independent double-trace coefficient to zero. In the standard large-$N_c$ generating-functional organization, an additional flavor trace generically represents an additional quark-loop topology and motivates one nominal $1/N_c$ suppression [@Kaiser:2000gs]. That rule is only coefficient/topology bookkeeping: field normalizations, singlet insertions, anomalous-source dependence, and channel-specific leading terms can modify the net scaling, so trace number alone does not prove a numerical suppression for every operator. It is also separate from the large-$N_c$ counting of powers and derivatives of Alda's anomalous source $X_{TC}$. Common singlet/octet $c_i,d_i$ are therefore a model ansatz and operator truncation, not a symmetry theorem. Ideal neutral-vector mixing and physical masses are likewise phenomenological baseline choices. This model is sufficient for first $aV\gamma$ and $aVV$ predictions and is not a global architecture blocker.

Upgrade when fitting $\eta'$, $\omega/\phi$, or other OZI-sensitive data at subleading $1/N_c$ accuracy. Concrete additions include the breaking operator $i\epsilon_{\mu\nu\rho\sigma}\langle V^{\mu\nu}V^{\rho\sigma}\rangle\langle\chi_-\rangle$ alongside its single-trace counterpart, and trace-factorized derivative structures such as $\epsilon_{\mu\nu\rho\sigma}\langle V^{\mu\nu}f_+^{\rho\alpha}\rangle\langle\nabla_\alpha u^\sigma\rangle$ or analogous $VVP$ factorizations. A complete $U(3)$ upgrade must additionally enumerate parity-allowed coefficient functions of $X_{TC}$ and operators with derivatives of $X_{TC}$, then reduce the candidate set using integration by parts, equations of motion, and $U(3)$ trace identities. These examples are upgrade templates, not a claim of a complete independent basis [@Roig:2013baa; @Herrera-Siklody:1996tqr].

The source conventions must not be identified by name alone. Denote the Herrera-Siklody et al. invariant by

$$
X_{\rm HS}=\langle\log U\rangle+\hat\theta,
$$

while Alda Eq. (III.14), with $U=u^2$ on a common logarithm branch, defines the distinct convention

$$
X_{TC}=2i\left(\langle\log u\rangle+i c_G\frac{a}{f_a}\right).
$$

The recorded conventions do not derive the normalization and sign map between $\hat\theta$ and Alda's ALP source $c_Ga/f_a$. Therefore no exact factor relation between $X_{\rm HS}$ and $X_{TC}$ is asserted here. The common-branch identity $\langle\log U\rangle=2\langle\log u\rangle$ matches only the mesonic logarithms; it does not supply the missing source map. Importing a Herrera-Siklody coefficient function consequently requires a separately stated and derived source identification [@Herrera-Siklody:1996tqr; @Alda:2026nrd]. Coefficient functions of Alda's $X_{TC}$ can contribute to dynamical-ALP amplitudes after expanding $X_{TC}$ in the ALP field. Their omission here is an explicit large-$N_c$/operator truncation, not a consequence of topological observables being absent.

## Width decision by use case

An external vector receives no Breit--Wigner factor. Reported two-body $a\to V\gamma$ and $a\to VV$ predictions may use on-shell vectors only when the state is narrow enough and the kinematics are sufficiently far from its daughter threshold for an NWA to be controlled. A broad $\rho$ or $K^*$, or any vector produced near daughter threshold, should instead be reported through the inclusive daughter channel with a normalized spectral integral. Widths are therefore not blockers for the odd operator basis, but they can be mandatory for the physical observable even when a formal on-shell two-body amplitude exists.

For an internal $\rho$ or $K^*$, the baseline is an energy-dependent P-wave width with physical $\pi\pi$ or $K\pi$ thresholds. For an internal $\omega$, a narrow constant width is acceptable unless the calculation resolves the $3\pi$ line shape. For an internal $\phi$ near its pole, use separate charged- and neutral-$K\bar K$ energy-dependent widths and add $3\pi$ when precision requires it. A constant-width Breit--Wigner is only a compatibility or local-fit alternative; Dyson-resummed or dispersive treatments define precision systematics rather than the baseline [@Rosell:2004mn; @Pich:2010sm; @Hanhart:2012wi; @GomezDumm:2013sib].

After any dressing, verify $k_\mu\mathcal M^\mu=0$ for every external photon and for the complete amplitude, not only an isolated graph. Also verify the low-energy/anomaly normalization and that contact, vector--photon-mixing, and exchange terms use one consistent dressing prescription.

## Anomaly coefficient and external tensor

The numerical $c_{AB}$ is not needed by the selected pure-RChT baseline. Once the external sources, current assignment, regulator, and Bardeen-type counterterm convention are fixed, WZW reproduces the corresponding anomalous variation. This statement does not fix finite invariant local odd-parity counterterms, does not by itself determine resonance matching, and is not a convention-independent assignment of the anomaly among currents [@Bai:2025fvl; @Fujiwara:1984mp]. The $c_i,d_i$ and any explicitly retained local invariant terms remain independent homogeneous couplings. The generic background-source helper retains $c_{AB}$ because that convention can be useful, but it requires explicit caller input and must not define hadronic resonance couplings by default.

For a stable vector, the Ecker kinetic convention gives the LSZ matrix element

$$
\langle0|V_{\mu\nu}(0)|V(p,\lambda)\rangle
=i\frac{\sqrt{Z_V}}{M_V}
\left(p_\mu\varepsilon_\nu^{(\lambda)}-p_\nu\varepsilon_\mu^{(\lambda)}\right),
$$

up to a consistently fixed one-particle phase, where $Z_V$ is the residue assigned to the stable transverse pole in the chosen field convention [@Ecker:1988te; @Ecker:1989yg]. Tree-level canonical normalization has $Z_V=1$ and $M_V$ equal to the tree-level mass parameter. An unstable vector has a complex pole and no exact one-particle LSZ state. Calling it “external” means an NWA projection onto a real on-shell state, with the residue treatment and accuracy declared; the displayed matrix element is then an NWA prescription, not an exact matrix element at the complex pole. Broad vectors and near-daughter-threshold configurations require the inclusive daughter channel instead.

## Near-degeneracy domain

Alda Eq. (III.37) requires

$$
|m_a^2-m_i^2|\gg\frac{f_\pi}{f_a}
\left|\operatorname{tr}(C\lambda_i^\dagger)\right|
\approx\frac{f_\pi}{f_a}m_i^2.
$$

The approximate right-hand side shows that the excluded real-mass neighborhoods scale as $f_\pi/f_a$ and are generally negligible in measure or scan resolution when $f_a\gg f_\pi$. This is not permission to retain a violating point. The actual pointwise test must use the exact overlap $|\operatorname{tr}(C\lambda_i^\dagger)|$: every sampled point that violates the inequality is locally outside the perturbative-mixing domain and must be masked/flagged or treated with a coupled mixing calculation. A finite $\mathbb M_a^{-1}(C)$ is necessary but not sufficient, and its matrix magnitude is not the validity criterion. Equation (III.37) is a real tree-level mass-separation condition. For stable or sufficiently narrow pseudoscalars, a retained violating point requires coupled diagonalization of the relevant quadratic system. For unstable states, real values of $m_a^2$ and $m_i^2$ can appear close without the complex poles being close, or conversely; validity is controlled by complex-pole separation relative to the off-diagonal mixing and requires a coupled complex propagator/pole analysis rather than ordinary Hermitian diagonalization [@Alda:2026nrd].

## Calibration decision

The minimal calibration set is:

- $V\to P\gamma$ radiative widths, constraining direct-plus-exchange $VJP/VVP$ combinations;
- $P\to\gamma\gamma$ widths and transition-form-factor normalizations/slopes, testing the chosen WZW/current/counterterm convention and constraining invariant local and photon--vector contributions;
- $e^+e^-\to P\gamma$ or $VP$ line shapes, constraining timelike combinations, phases, and internal-vector prescriptions;
- $\omega/\phi\to3\pi$, constraining odd three-meson amplitudes, direct terms, and the $\omega/\phi$ line-shape sector.

Use PDG correlations when available. Otherwise document an independent-Gaussian experimental approximation and add separately stated uncertainty models for line-shape and symmetry-breaking effects [@ParticleDataGroup:2024cfk]. No physical central values or fitted couplings are selected here.

### Selected fixed OPE theory prior

The user authorizes option 2 with a fixed Gaussian covariance for the five literature OPE residual definitions. This statistical model is an **assumption, not literature**. Evaluate the structural scales and covariance once at the reference point

$$
N_c=3,\qquad M_V=0.770\,\mathrm{GeV},\qquad F=0.092\,\mathrm{GeV},
\qquad F_V=\sqrt2F,
$$

using baseline truncation fraction $\tau=1/N_c$ and equicorrelation $\rho=0$. Do not update this covariance as $c_i,d_i$ vary. Required sensitivity cases are $\tau=0.2,0.5$ and $\rho=0.5$. The exact formula and numerical baseline widths are recorded in [the fixed-prior fit-infrastructure record](correlated-odd-parity-fit-infrastructure.md).

## Blocker classification

**Genuine blocker for a physical fit:** validated data ingestion and observable amplitudes/predictors for the selected channels. The generic fixed-covariance minimizer and baseline OPE prior now exist, but no $c_i,d_i$ fit is executable or has been run.

**Implementation tasks, not physics-choice blockers:** on-shell amplitude reduction; the fixed external-tensor normalization; internal-vector line-shape code by the use cases above; spectral daughter distributions; and complete-amplitude Ward-identity tests.

**Nonblockers at baseline:** independent double traces and singlet-specific couplings within the declared model; nonideal $\omega$--$\phi$ mixing; a free hadronic $c_{AB}$; precision Dyson/dispersive line shapes; and generic construction of $\mathbb M_a^{-1}(C)$ away from the local Alda degeneracy domain. Pointwise validity masking/warnings and coupled treatment of any intentionally retained violating points remain required implementation work.
