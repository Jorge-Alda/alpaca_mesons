---
title: Stage-4 nonideal vector response in radiative and two-photon amplitudes
status: derived
created: 2026-10-13
updated: 2026-10-15
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx, Ecker:1988te, Ecker:1989yg]
provenance: analytical-derivation
---

Canonical tensor and odd-coupling conventions follow [the Levi-Civita audit](levi-civita-convention.md). The linked original scripts and saved artifacts use the explicitly identified reference backend; apply that record's signed coefficient and tensor adapter when interpreting them in the canonical convention.

# Action-level amplitude integration

## Inputs and assumptions

Use the canonical fields and actual Lorentz contractions in [radiative SM amplitude](radiative-sm-amplitude.md) and [neutral TFF amplitude](radiative-tff-amplitude.md), the physical generators in [nonideal mass action](../foundations/vector-nonideal-mass-action.md), the [mass extension](../foundations/vector-mass-extension.md), and the signed source residues derived in [leptonic widths](../foundations/vector-leptonic-width.md). Calibration branches are recorded in [stage 3](../foundations/vector-stage3-calibration.md). Pseudoscalar coefficients are independent explicit inputs; their Alda provenance is the arXiv version indexed in [neutral meson generators](../../extracts/Alda:2026nrd/neutral-meson-generators.md), not the local manuscript.

Assume tree level, one vector nonet, canonical tensor kinetic terms, real vacuum spurions, the declared mass/source invariants and the existing single-trace odd action. Natural units, masses in GeV, metric $(+---)$ and $\epsilon^{0123}=+1, \epsilon_{0123}=-1$ apply. This extension does not assert completeness of the flavor-breaking odd action or establish its short-distance matching. No loop, kinetic-breaking or momentum-dependent self-energy prescription is introduced.

Both mass branches remain in stage 3. The user's selection rule is to prefer the smaller departure from ideal mixing when all other considerations are equal; it supplies no likelihood factor or sign determination. Source-sign and charged/neutral input scenarios remain independent choices.

For the numerical calibration application, the even-vector mass/source block
retains stage 3's exact-isospin projection, while the odd mass-spurion and
neutral-pseudoscalar inputs retain Alda's full fitted isospin breaking.
This explicitly declared partial-isospin truncation omits rho–isoscalar mixing;
it does not establish a complete isospin-resolved vector action.

## Executable derivation

Run [the script](../../../scripts/derive_vector_stage4_amplitudes.py):

```text
.venv/bin/python scripts/derive_vector_stage4_amplitudes.py
```

1. Insert $V_i=r_iF_\gamma$ into the tensor quadratic action and factor out $F_\gamma^2/4$. SymPy differentiates $\sum_i[(M_i^2-s)r_i^2-2\sqrt2eR_ir_i]$, solves for the response, and assembles

$$
V/F_\gamma=\sqrt2eW(s),\qquad W(s)=\sum_i\frac{R_iT_i}{M_i^2-s}.
$$

Here $R_i$ is signed and has mass dimension one. The source-action trace calculation gives

$$
R_i=F_V\operatorname{tr}(T_iQ)+z_V\operatorname{tr}(T_i\{X,Q\})
+y_V\operatorname{tr}(T_i)\operatorname{tr}(XQ),\qquad
X=\chi_+|_{\Phi=0}=2\chi.
$$

For the fixed positive-strange convention $T_\omega=\cos\theta_V T_{\omega_I}-\sin\theta_V T_{\phi_I}$ and $T_\phi=\sin\theta_V T_{\omega_I}+\cos\theta_V T_{\phi_I}$. The same rotation acts on the residues. SymPy verifies that the physical diagonal response equals the ideal-coordinate inverse mass-matrix response by multiplying by $K(s)=U^T\operatorname{diag}(M_\omega^2-s,M_\phi^2-s)U$. This is basis covariance; holding a diagonal split mass matrix fixed while rotating only generators is not that transformation.

2. Reuse the actual Fourier-index VJP/VVP contractions, including both VVP assignments, to generate generic flavor-pair polynomials. For any matrices $A,B$, with the external-leg conjugations of the baseline, define

$$
q(A,B)=\operatorname{tr}(\{A,B\}T_P),\qquad
h(A,B)=\operatorname{tr}(\{A,B\}\{T_P,\chi\}),
$$

$$
C(A,B;s,t)=q(A,B)[(c_1+c_2-c_5)m_P^2+(-c_1+c_2+c_5-2c_6)s+(c_1-c_2+c_5)t]+4c_3h(A,B),
$$

$$
H(A,B;s,t)=q(A,B)[(d_1-d_3)m_P^2+d_3(s+t)]+4d_2h(A,B).
$$

These are contraction outputs, not replacements for the Lorentz calculation. The script keeps physical neutral $T_P=\operatorname{diag}(t_{Pu},t_{Pd},t_{Ps})$ independent of vector generators; no eta/eta-prime normalization is inferred from a vector rotation.

3. Insert the response into each internal vector slot. The explicit photon building block remains $f_+=-2eQF_\gamma$: its $Q$ is not replaced by $W$.

$$
F_P(s_1,s_2)=\frac{\sqrt2N_Cq(Q,Q)}{8\pi^2F}
+\frac4{FM_V}[C(W(s_1),Q;s_1,s_2)+C(W(s_2),Q;s_2,s_1)]
-\frac{4\sqrt2}{F}H(W(s_1),W(s_2);s_1,s_2).
$$

$M_V$ remains the operator-normalization scale. Every term has dimension GeV$^{-1}$; $W$ has dimension GeV$^{-1}$. The WZW anomaly coefficient is unchanged by the declared even-vector extensions.

4. For a real-photon radiative decay, use the physical external vector matrix $T_V$ and its canonical wavefunction factor $1/M_{\rm ext}$; this factor is never replaced by an internal pole denominator. With $F_K=\operatorname{tr}(T_V[[T_P,Q],\chi])$,

$$
C_V=q(T_V,Q)[(c_1+c_2-c_5)m_P^2+(-c_1+c_2+c_5-2c_6)M_{\rm ext}^2]+4c_3h(T_V,Q)-4c_4F_K,
$$

$$
G=\frac{2\sqrt2e}{FM_{\rm ext}M_V}C_V
-\frac{4e}{FM_{\rm ext}}H(T_V,W(0);M_{\rm ext}^2,0).
$$

Crossing retains the external vector mass and the existing two width formulas. $c_4$ is absent from diagonal neutral $T_P$ channels but remains in charged radiative channels through the baseline commutator.

5. Expand $u=1+izT_P/(\sqrt2F)$ and its conjugate, with $z$ a scalar field bookkeeper. SymPy differentiates $\chi_+=u^\dagger\chi u^\dagger+u\chi u$ and $f_+=-e(uQu^\dagger+u^\dagger Qu)$ and obtains zero coefficients at first SM-pseudoscalar order. Consequently the declared even source and mass extensions add no direct one-SM-$P$ vertex at this tree order; they enter these amplitudes through physical masses, generators and source responses. This statement uses a real vacuum spurion and is not an ALP/external-source assertion or a completeness claim for odd breaking operators.

## Verification and limitations

The executable verifies both photon Ward identities with momentum substituted for polarization in the unreduced indexed contractions, Bose symmetry with virtualities and flavor responses exchanged, rotation covariance, and the universal-source/common-pole reduction $W(s)=F_VQ/(M_V^2-s)$. It compares the generalized TFF with the complete prior action-derived amplitude in that limit. Thus the published normalization checks of that baseline carry over in its limiting domain.

The old five OPE constraints and old rank statements were derived for the unbroken action. They are not automatically constraints or ranks of this enlarged action. Applying them can be a separately labelled diagnostic, but renewed high-energy matching is required before claiming matching of split masses and source residues. Physical meson inputs, interference-sensitive residue signs, interaction-breaking sufficiency and omitted self-energies remain explicit limits; no fit or precision branching-fraction result is produced here.

## Reusable implementation and observations

[`vector_sector.py`](../../../scripts/implementation/vector_sector.py) supplies
the canonical neutral sector, trace-derived source residues and $W(s)$.
`VectorSector.from_stage3` requires a declared rho/K-star scenario, phi-width
choice and source hypothesis. The default `branch='near_ideal'` implements
the user-selected tie-breaker; `branch='other'` retains the second solution.
The archived stage-3 input checksum and selected angle are validated.

[`vector_sector_amplitudes.py`](../../../scripts/implementation/vector_sector_amplitudes.py)
implements the formulas above, both crossed radiative widths, the real
two-photon width and the $Q^2|F|$ predictor. Its `VectorAmplitudeInputs` takes
explicit $T_P,\chi,F,M_V,m_P$ and a sector. $M_V$ is the operator normalization
scale, not automatically an internal or external pole mass.

The existing CLEO and BaBar observation factories in
[`odd_parity_data.py`](../../../scripts/implementation/odd_parity_data.py) accept
these new inputs as well as the old common-pole `TFFInputs`.
`load_vector_omega_pi0_gamma_observation` connects the same checksum-pinned
omega branching datum to the new sector. Their values, covariance policy,
asymmetric-error conversion and normalization nuisance are unchanged.
Mapping the helper's `etap` to the CLEO table's `eta_prime` is explicit.
No fit is run, and data connection does not validate the high-Q-squared model.

Use [the verification runner](../../../scripts/verify_vector_stage4.py) to run
focused new tests and existing amplitude/data regressions without triggering
the package's unrelated eager imports:

```text
.venv/bin/python scripts/verify_alda_neutral_inputs.py
.venv/bin/python scripts/verify_vector_stage4.py
```

Tests exercise independent inverse-matrix/spectral response agreement,
source normalization, phase and degenerate-basis invariance, both branches,
old-amplitude limits, crossing, published pion width, pole/invalid inputs,
and unchanged likelihood/finite predictions in the actual data factories.
