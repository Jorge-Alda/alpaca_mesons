---
title: Odd-parity coupling prescriptions and direct ALP-vector-photon VJP insertions
status: mixed
created: 2026-09-28
updated: 2026-09-28
bibtex_keys: [Alda:2026nrd, Bai:2025fvl, Ruiz-Femenia:2003jdx, Ecker:1989yg, Fujiwara:1984mp]
provenance: analytical-derivation
---

# Odd-Parity Coupling Prescriptions and Direct $aV\gamma$ Insertions

## Scope and evidence labels

This record answers two questions: how to determine the odd-parity $c_i,d_i$, and which direct $aV\gamma$ vertices follow from the exact seven-operator $VJP$ basis. The symbolic derivation is [`scripts/derive_direct_alp_vgamma_vjp.py`](../../../scripts/derive_direct_alp_vgamma_vjp.py). Inputs are the [Ruiz--Femenia operator extract](../../extracts/Ruiz-Femenia:2003jdx/vector-odd-parity.md), the [two-source ALP derivation](two-source-alp-rcht-odd-parity-extension.md), the [framework comparison](../foundations/vector-meson-framework-comparison.md), and the [background-field/RChT matching extract](../../extracts/Bai:2025fvl/no-axial-vector-and-rcht-matching.md). [@Ruiz-Femenia:2003jdx; @Bai:2025fvl]

Labels below mean: **literature** for statements established in the cited sources, **derived assessment** for comparisons inferred from those facts, and **recommendation** for the proposed implementation strategy.

## Part 1: prescriptions for $c_i,d_i$

The three prescriptions do not determine the same object. Prescription A constrains combinations in a truncated RChT basis, B estimates combinations visible in chosen observables, and C fixes a background-field/HLS realization only after selecting otherwise free homogeneous terms.

| Criterion | A: QCD short-distance-constrained RChT | B: correlated phenomenological fit | C: HLS/background-field WZW plus selected homogeneous action |
|---|---|---|---|
| Predictivity | **Pro:** gives analytic relations among couplings. **Con:** only combinations are fixed; $c_4,c_7,d_4$ are untouched in the cited VVP matching. | **Pro:** data determine combinations entering measured amplitudes. **Con:** prediction depends on calibration channels and covariance treatment. | **Pro:** a selected homogeneous action can be compact and highly predictive. **Con:** anomaly matching does not fix that selection. |
| Assumptions | Single-resonance, leading-$1/N_c$ and selected OPE/form-factor asymptotics. | Dataset, truncation, line shapes, errors, correlations, and theory-error model must be specified. | Hidden-gauge/background-source matching, plus chosen homogeneous coefficients; common VMD, $a=2$, ideal-mixing choices are additional assumptions. |
| Parameter identifiability | Rank-deficient: e.g. $d_1+8d_2$, not $d_1,d_2$ separately, is constrained. | A global multi-channel fit can expose flat directions and correlations but cannot identify absent structures. | Parameters set by a homogeneous-action choice are assumed rather than identified; freeing them restores degeneracies. |
| $SU(3)$ breaking | The cited constraints are stated to survive beyond the chiral limit, but they do not constitute a complete broken-$SU(3)$ resonance theory. | Breaking and singlet mixing can be fitted if sufficiently diverse channels are included, at the cost of parameters. | Requires explicit breaking operators and nonideal mixing; neither follows from WZW matching. |
| Widths and phases | Leading large-$N_c$ result is narrow and tree level; widths are extra dynamical input. | Can constrain width/self-energy prescriptions and phases, but these must preserve Ward identities and low-energy matching. | WZW matching supplies no propagator or width prescription; phenomenological insertions add model dependence. |
| Field-redefinition invariance | Holds for complete amplitudes only when all induced local terms are retained; a finite representation-dependent truncation can differ. | Amplitude-level fits can be invariant if the transformed complete basis is used; fitted individual $c_i,d_i$ are basis dependent. | A source-to-propagating-field identification is not a literal field identity; matched local terms are required. |
| Double counting | Adding local odd-parity terms already saturated by retained vectors can double count and spoil asymptotics. | Contact and exchange pieces must be fitted in one declared basis, not counted twice. | A matched HLS/WZW result is an alternative representation/cross-check, not an additive extra contribution to RChT. |
| ALP flavor extrapolation | The operator basis accepts general ALP spurions, but QCD constraints on hadronic combinations do not determine every flavor-changing ALP amplitude. | Standard hadronic data calibrate universal resonance couplings, while ALP spurion traces and poorly sampled breaking directions remain extrapolations. | Covariant source insertion organizes flavor, but a selected homogeneous action and idealized flavor assignments can hide uncertainty in general Alda flavor directions. |

**Literature-established facts.** The cited OPE matching gives

$$
4c_3+c_1=0,\qquad c_1-c_2+c_5=0,\qquad
c_5-c_6=\frac{N_C}{64\pi^2}\frac{M_V}{\sqrt2F_V},
$$

$$
d_1+8d_2=-\frac{N_C}{64\pi^2}\frac{M_V^2}{F_V^2}
+\frac{F^2}{4F_V^2},\qquad
d_3=-\frac{N_C}{64\pi^2}\frac{M_V^2}{F_V^2}
+\frac{F^2}{8F_V^2},
$$

and does not fix all coefficients. Strict VMD fails in the cited $\omega\to3\pi$ test. Gauged WZW requires its counterterm, while homogeneous anomalous terms remain free. [@Ruiz-Femenia:2003jdx; @Bai:2025fvl; @Fujiwara:1984mp]

**Derived assessment.** None of A--C alone simultaneously supplies QCD asymptotics, identifiable uncertainties, realistic symmetry breaking and widths, and representation-independent complete amplitudes.

**Recommendation.** Use the full RChT odd sector together with WZW and its required counterterm; impose QCD short-distance relations as uncertainty-bearing constraints or priors rather than exact numerical identities; perform a correlated multi-channel fit including a consistent width and $SU(3)$-breaking treatment; and match an HLS/background-field realization to the same poles and low-energy terms as a cross-check. The matched HLS result is **not additive physics**. Compare complete amplitudes, not direct/exchange pieces or individual basis coefficients.

## Part 2: exact $VJP$ basis and source insertions

### Conventions, units, and assumptions

Natural units are used. The fields $a$, $f_a$, $B_0$, $m_q$, and $M_V$ have mass dimension one; $e$, $c_i$, $d_i$, $\mathbb M_a^{-1}(C)$, $c_q^{AP}$, $Q$, and flavor generators are dimensionless. We use $\epsilon_{0123}=+1$, $V=T_VV$, and

$$
\mathcal L_{VJP}=\sum_{i=1}^{7}\frac{c_i}{M_V}\mathcal O^i_{VJP}.
$$

The derivation keeps one physical ALP at $O(1/f_a)$, uses the Goldstone vacuum and the standard electromagnetic source, and performs no integration by parts, equation-of-motion reduction, short-distance substitution, or on-shell reduction except the explicit $c_7$ check below. Flavor matrices are fully general $3\times3$ matrices in the symbolic calculation; no matrix commutativity is assumed.

$$
u_\mu\big|_{a_{\rm physical}}=-2\frac{\partial_\mu a}{f_a}\mathbb M_a^{-1}(C),
\qquad
\chi_-\big|_{a_{\rm physical}}=-4iB_0\frac{a}{f_a}
\{c_q^{AP}+\mathbb M_a^{-1}(C),m_q\},
\qquad
f_+^{\mu\nu}=-2eQF^{\mu\nu}.
$$

For the standard electromagnetic-compatible Goldstone vacuum, $f_-^{\mu\nu}=0$ because the charge matrix commutes with the ALP spurions. Alda permits flavor-changing entries only among equal-charge quarks, so this commutation is the relevant framework assumption.

### Published seven-operator basis

The exact basis is [@Ruiz-Femenia:2003jdx]

$$
\begin{aligned}
\mathcal O^1_{VJP}&=\epsilon_{\mu\nu\rho\sigma}\langle\{V^{\mu\nu},f_+^{\rho\alpha}\}\nabla_\alpha u^\sigma\rangle,\\
\mathcal O^2_{VJP}&=\epsilon_{\mu\nu\rho\sigma}\langle\{V^{\mu\alpha},f_+^{\rho\sigma}\}\nabla_\alpha u^\nu\rangle,\\
\mathcal O^3_{VJP}&=i\epsilon_{\mu\nu\rho\sigma}\langle\{V^{\mu\nu},f_+^{\rho\sigma}\}\chi_-\rangle,\\
\mathcal O^4_{VJP}&=i\epsilon_{\mu\nu\rho\sigma}\langle V^{\mu\nu}[f_-^{\rho\sigma},\chi_+]\rangle,\\
\mathcal O^5_{VJP}&=\epsilon_{\mu\nu\rho\sigma}\langle\{\nabla_\alpha V^{\mu\nu},f_+^{\rho\alpha}\}u^\sigma\rangle,\\
\mathcal O^6_{VJP}&=\epsilon_{\mu\nu\rho\sigma}\langle\{\nabla_\alpha V^{\mu\alpha},f_+^{\rho\sigma}\}u^\nu\rangle,\\
\mathcal O^7_{VJP}&=\epsilon_{\mu\nu\rho\sigma}\langle\{\nabla^\sigma V^{\mu\nu},f_+^{\rho\alpha}\}u_\alpha\rangle.
\end{aligned}
$$

### Stepwise SymPy derivation

1. Substitute $f_+=-2eQF$ and $u_\mu=-2(\partial_\mu a)\mathbb M_a^{-1}(C)/f_a$. SymPy expands the general component matrices and obtains the common product $(-2e)(-2/f_a)=+4e/f_a$ with flavor trace $\operatorname{tr}(\{T_V,Q\}\mathbb M_a^{-1}(C))$.
2. Retaining every Lorentz contraction separately gives

$$
\begin{aligned}
\mathcal O^1_{aV\gamma}&=\frac{4e}{f_a}\operatorname{tr}(\{T_V,Q\}\mathbb M_a^{-1}(C))
\epsilon_{\mu\nu\rho\sigma}V^{\mu\nu}F^{\rho\alpha}\partial_\alpha\partial^\sigma a,\\
\mathcal O^2_{aV\gamma}&=\frac{4e}{f_a}\operatorname{tr}(\{T_V,Q\}\mathbb M_a^{-1}(C))
\epsilon_{\mu\nu\rho\sigma}V^{\mu\alpha}F^{\rho\sigma}\partial_\alpha\partial^\nu a,\\
\mathcal O^5_{aV\gamma}&=\frac{4e}{f_a}\operatorname{tr}(\{T_V,Q\}\mathbb M_a^{-1}(C))
\epsilon_{\mu\nu\rho\sigma}(\partial_\alpha V^{\mu\nu})F^{\rho\alpha}\partial^\sigma a,\\
\mathcal O^6_{aV\gamma}&=\frac{4e}{f_a}\operatorname{tr}(\{T_V,Q\}\mathbb M_a^{-1}(C))
\epsilon_{\mu\nu\rho\sigma}(\partial_\alpha V^{\mu\alpha})F^{\rho\sigma}\partial^\nu a,\\
\mathcal O^7_{aV\gamma}&=\frac{4e}{f_a}\operatorname{tr}(\{T_V,Q\}\mathbb M_a^{-1}(C))
\epsilon_{\mu\nu\rho\sigma}(\partial^\sigma V^{\mu\nu})F^{\rho\alpha}\partial_\alpha a.
\end{aligned}
$$

3. For $\mathcal O^3$, the scalar factors give

$$
i(-2e)\left(-\frac{4iB_0}{f_a}\right)=-\frac{8eB_0}{f_a},
$$

because $i^2=-1$. Thus

$$
\mathcal O^3_{aV\gamma}=-\frac{8eB_0}{f_a}
\operatorname{tr}\!\left(\{T_V,Q\}
\{c_q^{AP}+\mathbb M_a^{-1}(C),m_q\}\right)
\epsilon_{\mu\nu\rho\sigma}aV^{\mu\nu}F^{\rho\sigma}.
$$

4. $\mathcal O^4_{aV\gamma}=0$ under the stated $f_-=0$ framework assumption.

These are operator-level source insertions. Multiplication by $c_i/M_V$ gives their contribution to $\mathcal L_{VJP}$.

### The $c_7$ on-shell check

$\mathcal O^7$ is an off-shell $aV\gamma$ vertex. For an external on-shell antisymmetric-tensor vector,

$$
V^{\mu\nu}=\frac{p_V^\mu\varepsilon_V^\nu-p_V^\nu\varepsilon_V^\mu}{M_V},
$$

the script evaluates all four free-$\rho$ components of

$$
\epsilon_{\mu\nu\rho\sigma}p_V^\sigma V^{\mu\nu}
$$

as zero. The result follows from antisymmetry and repeated factors of $p_V$; no transversality or mass-shell scalar equation is needed. Therefore the **direct** $c_7$ term vanishes for this external on-shell tensor vector. This does not remove $c_7$ from general off-shell Green functions and does not automatically remove all induced or local terms after field redefinitions.

### Indirect terms and physical interpretation

Every $VVP$ coefficient $d_i$ can induce $aV\gamma$ through the standard $F_V$ photon--vector mixing. In particular, $d_2$ contributes through the physical $\chi_-$ insertion and must not be dropped. The division into direct $VJP$ and indirect $VVP$ contributions is representation and field-basis dependent; only the complete amplitude, including WZW/counterterm and induced local terms, is physical.

## Verification

Run

```text
.venv/bin/python scripts/derive_direct_alp_vgamma_vjp.py
```

The script derives the derivative and mass flavor factors with general component matrices, preserves the five derivative Lorentz contractions, verifies $\mathcal O^4=0$ at $f_-=0$, and explicitly contracts the on-shell $\mathcal O^7$ tensor. All assertions must print as passed.
