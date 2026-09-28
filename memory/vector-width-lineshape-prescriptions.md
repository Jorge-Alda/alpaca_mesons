---
title: Width and line-shape prescriptions for odd-parity antisymmetric-tensor RChT amplitudes
status: mixed
created: 2026-09-28
updated: 2026-09-28
bibtex_keys: [Ecker:1988te, Ecker:1989yg, Ruiz-Femenia:2003jdx, Rosell:2004mn, Pich:2010sm, Hanhart:2012wi, GomezDumm:2013sib, Jamin:2000wn, Omnes:1958hv, Watson:1952ji, ParticleDataGroup:2024cfk]
provenance: literature
---

# Width and line-shape choices

## Scope and conventions

This is a prescription comparison for tree-level odd-intrinsic-parity amplitudes with antisymmetric-tensor vectors. The odd vertex is kept tree level; a line shape is inserted only when a vector is an internal line or when a production form factor is being modeled. The antisymmetric-tensor representation does not alter the physical two-body cut, and a field-representation change requires all induced local terms. [@Ecker:1989yg; @Ruiz-Femenia:2003jdx]

For a vector $X$ decaying to two pseudoscalars $a,b$, define

$$
q_{ab}(s)=\frac{\sqrt{\lambda(s,m_a^2,m_b^2)}}{2\sqrt{s}},
\qquad \lambda(x,y,z)=x^2+y^2+z^2-2xy-2xz-2yz,
$$

and the P-wave width model

$$
\Gamma_{X\to ab}(s)=\Gamma_{X\to ab}(M_X^2)\frac{M_X^2}{s}
\left[\frac{q_{ab}(s)}{q_{ab}(M_X^2)}\right]^3
\theta\!\left(s-(m_a+m_b)^2\right).
$$

Thus the width turns on as $(s-s_{th})^{3/2}$ for a two-pseudoscalar P wave. A three-pion channel has a three-body phase-space onset and cannot be represented by the same single $q^3$ factor. The $K^*$ channel is the $K\pi$ P wave; the rho channel is the pi pi P wave. These channel assignments are consistent with resonance chiral partial-wave treatments. [@Jamin:2000wn; @ParticleDataGroup:2024cfk]

The useful dressed denominator is

$$
D_X(s)=\frac{1}{M_{X,0}^2-s-\Sigma_X(s)},
\qquad \operatorname{Im}\Sigma_X(s+i0)=+\sqrt{s}\,\Gamma_X(s),
$$

above an open cut, so that the denominator contains $-i\sqrt{s}\,\Gamma_X(s)$ for the stated $+i0$ convention. A constant-width Breit--Wigner replaces the self-energy by an imaginary constant, while a Dyson construction retains a renormalized complex $\Sigma_X(s)$ and its real part. [@Rosell:2004mn; @Pich:2010sm]

## Prescription comparison

| Prescription | Threshold and analyticity | Unitarity and Ward-identity risk | Fit complexity | External versus internal vector |
|---|---|---|---|---|
| Zero width / narrow state, $D=1/(M^2-s)$ | No absorptive cut; meromorphic tree-level pole. Correct as leading-$1/N_c$ RChT, not near an open broad resonance. [@Ecker:1988te] | No loop unitarity; no width-induced gauge issue. | Minimal. | Use for an external vector only under a controlled NWA, away from daughter thresholds; also usable for low-resolution/off-resonance exchange. |
| Constant-width BW, $D=1/(M^2-s-iM\Gamma)$ | Nonzero width below threshold and no physical branch point; analyticity is only phenomenological. | Can spoil low-energy limits, crossing, and gauge cancellations if only one internal line is dressed. The photon Ward identity is safest when the photon appears through $F_{\mu\nu}$, but vector-photon mixing and contact terms must not be dressed inconsistently. | One width per state; easy fit. | Avoid for precision internal exchange; harmless as a declared narrow-window fit ansatz. Never attach it to an external asymptotic state. |
| Energy-dependent width, $-i\sqrt{s}\Gamma(s)$ | Correct leading threshold power for specified channels, but the real dispersive part and coupled-channel cuts are omitted unless added separately. | Better elastic behavior than constant width, but not exact unitarity; the same gauge-consistent set of self-energy, mixing, and contact pieces is still required. | Moderate: masses, normalizations, thresholds, channel fractions, and possible Blatt--Weisskopf or form-factor parameters. | Practical internal-line baseline when the invariant range crosses one known two-body threshold and the amplitude is not used as a precision form factor. |
| Dyson-resummed RChT propagator | Has the loop cuts and a dispersive real part after renormalization; channels and counterterms define the analytic structure. [@Rosell:2004mn; @Pich:2010sm] | Best field-theory option within the specified truncation, but resumming a propagator alone is not a unitary amplitude. Vertices, mixing, subtraction terms, and local odd-parity terms must be treated consistently. | High: loop functions, subtraction scheme, correlated couplings, and coupled channels. | Use for precision internal rho and, with a dedicated calculation, K*; not justified for an external vector merely to assign it a width. |
| Dispersive / Omnes, $F(s)=P(s)\Omega(s)$ | Explicit right-hand cut and exact elastic phase; analyticity is controlled by subtractions and phase continuation. [@Omnes:1958hv; @Watson:1952ji; @Hanhart:2012wi] | Elastic unitarity is strongest; inelasticity requires a coupled-channel Omnes matrix. Electromagnetic Ward identities remain a property of the full tensor amplitude and its subtraction/contact terms, not of $\Omega$ alone. | Highest for a general multi-channel fit; phase input and continuation are required. [@GomezDumm:2013sib] | Best for external-current form factors or internal exchange when data/scattering phases exist; not a universal replacement for an RChT propagator in every three-body amplitude. |

## State/channel assessment

| State and dominant channel | Threshold behavior | Where prescriptions differ materially |
|---|---|---|
| rho: pi pi P wave | $(s-4m_{pi}^2)^{3/2}$; broad resonance and strong elastic phase motion. | Use a running two-body $\pi\pi$ P-wave width at baseline. Differences are material across the rho peak and in precision spectra; use dispersive or loop-resummed treatment for precision. A final-state rho generally calls for the inclusive $\pi\pi$ daughter observable rather than an external-state NWA. [@Hanhart:2012wi; @GomezDumm:2013sib; @Rosell:2004mn] |
| $K^*$: $K pi$ P wave | $(s-(m_K+m_pi)^2)^{3/2}$ with unequal-mass momentum; nearby inelastic thresholds can matter above the elastic window. | Use a running two-body $K\pi$ P-wave width at baseline. It is material in $K\pi$ invariant-mass scans and near threshold; a final-state $K^*$ requires a spectral daughter treatment when the NWA is not controlled. [@Jamin:2000wn] |
| omega: $3pi$ (often through rho pi) | Three-body threshold at $3m_pi$ with a nontrivial subchannel integral; narrow pole. | A constant width is acceptable as a declared local/narrow-window model. Use an explicit $3\pi$ self-energy or spectral function when resolving the line shape or daughter threshold; do not model it with a two-body $q^3$ width. [@ParticleDataGroup:2024cfk; @Ruiz-Femenia:2003jdx] |
| phi: $K\bar K$ P wave, plus $3pi$ | $K\bar K$ onset is P-wave and close to the physical phi pole; charged and neutral thresholds are split. | Use separate charged- and neutral-$K\bar K$ running P-wave contributions plus $3\pi$ when the pole region is resolved. A constant width is only a local narrow-window alternative. [@ParticleDataGroup:2024cfk] |

## Ward identities and external vectors

For a photon represented through $F_{\mu\nu}$, each tree-level odd vertex is transverse in the photon momentum. A width insertion does not by itself destroy that algebraic transversality. The risk arises when vector-photon mixing, longitudinal/contact terms, or several diagrams are dressed differently: then the cancellations required by the complete current can fail. This is a derived implementation warning, not a claim that every energy-dependent width violates a Ward identity. The safe rule is to test $k_\mu\mathcal M^\mu=0$ after dressing the complete amplitude, and to preserve the low-energy normalization and local terms. [@Ecker:1989yg; @Rosell:2004mn]

A stable vector is an asymptotic state. An unstable vector may be treated as external only within a controlled NWA, when it is narrow enough for the requested accuracy and the production point is sufficiently far from daughter thresholds. In that regime do not replace its external LSZ matrix element by a Breit--Wigner factor. For broad $\rho/K^*$ states or near daughter threshold, the physical prediction is an inclusive daughter-channel observable with a spectral integral. An internal vector carries a variable $s$ and may need a line shape whenever the integration samples a threshold, pole, or rapidly varying phase.

## Recommendation and systematics

**Practical baseline:** use a zero-width external-vector treatment only under the NWA conditions above. For internal rho and $K^*$ use running two-body P-wave widths with physical thresholds. For omega use a constant width only as a local/narrow-window approximation, otherwise an explicit $3\pi$ line shape. For phi use separate charged/neutral $K\bar K$ running widths plus $3\pi$ when resolving the pole. State the real-part convention, keep the odd-parity numerator and all local/contact terms consistent, and verify the photon Ward identity and soft/low-energy normalization for the full amplitude.

**Systematic alternatives:** (1) narrow-state versus baseline comparison; (2) constant-width BW only as a deliberately local fit ansatz; (3) replace the baseline denominator by the renormalized Dyson-resummed RChT result for rho and K* when loop data and counterterms are available; (4) replace the relevant production form factor by a dispersive/Omnes representation, varying subtraction constants, phase continuation, and inelastic treatment. For phi, vary charged and neutral $K\bar K$ thresholds separately. For omega, vary the three-pion spectral model rather than pretending that its width has universal two-body P-wave kinematics.

The baseline is a modeling choice, not a unitary completion. A claim of exact unitarity should be reserved for an explicitly matched dispersive or coupled-channel construction; a claim of a Dyson-resummed RChT result requires the loop, subtraction, mixing, and vertex ingredients to be specified together. [@Rosell:2004mn; @Pich:2010sm; @Hanhart:2012wi; @GomezDumm:2013sib]
