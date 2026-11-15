---
title: Finite ALP mass in the hard photon propagator
status: derived
created: 2026-11-15
updated: 2026-11-15
bibtex_keys: [Alda:2026nrd, Hoferichter:2020lap]
provenance: analytical-derivation
---

# Inputs and assumptions

Use natural units and the tensor, source, phase and metric conventions in [the conventions record](alp-high-energy-tff-conventions.md). The normalized flavor distribution amplitude and its $a_2^f$ ansatz, the source-linear prefactor, renormalization scheme and threshold inputs are those of [the light-cone result](alp-light-cone-high-energy-tff.md). The free hard quark propagator of [Hoferichter:2020lap](../../extracts/Hoferichter:2020lap/pseudoscalar-light-cone-tff.md) is used at tree level. The calculation below isolates its mass dependence; it is not the complete correction at relative order $m_a^2/Q^2$. Higher-twist two- and multi-particle matrix elements also enter at that order and are not fixed by the leading-twist distribution amplitude.

# Derivation

For $Q_i^2=-q_i^2>0$ the scalar denominator for either photon ordering, after the symmetric distribution is used, is
$$
u Q_1^2+(1-u)Q_2^2+u(1-u)m_a^2.
$$
Expanding the propagator in $m_a^2$ gives
$$
\frac{1}{uQ_1^2+(1-u)Q_2^2}-\frac{u(1-u)m_a^2}{[uQ_1^2+(1-u)Q_2^2]^2}+O(m_a^4).
$$
Insert this expression in the flavor sum and prefactor in the light-cone record. For its declared Gegenbauer ansatz, the script derives at equal virtualities and, separately, for one real photon,
$$
\frac1{Q^2}-\frac{m_a^2}{Q^4}\frac{7-3a_2^f}{35},\qquad
\frac{3(1+a_2^f)}{Q^2}-\frac{(2+3a_2^f)m_a^2}{Q^4}.
$$
The second term in each expression quantifies the propagator contribution omitted in task 3.2, not the full power correction. The scalar projection retains $p^2=m_a^2$; the apparent $1/m_a^2$ cancels against $p_\mu\langle0|\bar q\gamma^\mu\gamma_5q|a\rangle$, as in the light-cone record. The finite physical matrix element rather than its projection ratio should be used when taking $m_a\to0$.

# Continuation and domains

For complex virtualities the tree denominator is $-u q_1^2-(1-u)q_2^2+u(1-u)m_a^2-i0$. Thus the Euclidean expression continues with the Feynman prescription; real timelike values can cross a propagator pole and acquire an imaginary part. The continuation must be made before integration, retaining $i0$; an absolute value or a Euclidean series is not a prescription for crossing a cut. The hard coefficient alone does not determine timelike hadronic resonance or threshold effects.

For fixed nonzero spacelike ratio, require both $Q_i^2\gg m_a^2$ and hadronic scales, with a fixed ratio and a hard scale well above the active quark masses for their massless treatment. Nearly equal virtualities satisfy the same conditions and additionally $|q\cdot p|\ll|q^2|$. The one-real-photon endpoint is conditional on leading-twist endpoint factorization; the displayed integral is finite for the assumed distribution, but the power expansion and higher-twist terms are not uniformly controlled there. For a fixed nonzero soft photon virtuality, the leading-power singly virtual limit is conditional on factorization and its finite-soft-virtuality power corrections are not given. If $Q^2$ is not parametrically larger than $m_a^2$, neither the truncated mass series nor an unquantified omission is justified. Timelike poles, hadronic cuts and heavy-flavor thresholds require their own continuation and inputs as in the [symmetric result](alp-symmetric-high-energy-tff.md).

# Verification

Run `.venv/bin/python scripts/derive_alp_finite_mass_hard_tff.py` and the same command with `--test`. SymPy expands the denominator and integrates the stated ansatz. Unit tests check the massless limit, large-$Q^2$ limit, equal and singly virtual integrals and photon exchange. These algebraic checks do not determine the omitted higher-twist matrix elements.
