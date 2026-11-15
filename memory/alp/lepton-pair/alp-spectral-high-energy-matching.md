---
title: Conditional spectral matching of the ALP transition form factor
status: derived
created: 2026-11-15
updated: 2026-11-15
bibtex_keys: [Alda:2026nrd, Messerli:2025rnv, Hoferichter:2018kwz]
provenance: analytical-derivation
---

# Inputs and assumptions

Use natural units and the conventions and source-linear light-current coefficient from [the light-cone calculation](alp-light-cone-high-energy-tff.md), the conditional [low-energy benchmarks](alp-low-energy-matching-benchmarks.md), and the spectral construction of [Messerli:2025rnv](../../extracts/Messerli:2025rnv/dispersive-tff-asymptotic-matching.md) and [Hoferichter:2018kwz](../../extracts/Hoferichter:2018kwz/pion-tff-dispersive-asymptotic-representation.md). The result below assumes a common asymptotic distribution amplitude $6u(1-u)$ for the relevant light-flavor current matrix elements, tree-level hard matching and a fixed renormalization scale. The flavor trace is over the three light quarks; $Q$ is their charge matrix, and $\langle\cdots\rangle$ denotes the flavor trace. The source matrices are kept general and are not assumed to commute. The coefficient is obtained from ALP current matrix elements, not a meson coefficient. The source-dependent contact and massive-heavy-quark contributions of [the symmetric assembly](alp-symmetric-high-energy-tff.md) require separate treatment and must not be counted twice. Non-asymptotic distribution amplitudes and their evolution are not reproduced by this spectral ansatz.

# Spectral contribution

The ALP current coefficient multiplying the normalized convolution is
$$
\frac{4f_\pi^2}{f_a}\left\langle Q^2\mathbb M_a^{-1}(\bm{\mathcal C})\right\rangle_\mu,
$$
with $\bm{\mathcal C}=-B_0\{m_q,c_q^{AP}\}-(m_0^2/N_f)c_{GA}\mathbf1$ and $\mathbb M_a$ as in [the axial-current record](alp-axial-current-matrix-element.md). When the LO chiral approximation is unavailable, replace the coefficient by $2p_\mu\langle0|\bar q\gamma^\mu\gamma_5Q^2q|a(p)\rangle_\mu/(i m_a^2)$. Multiplying it by three gives the coefficient of the source-inspired continuum integral
$$
\frac{12f_\pi^2}{f_a}\left\langle Q^2\mathbb M_a^{-1}(\bm{\mathcal C})\right\rangle_\mu
\int_{s_m}^\infty dx\,
\frac{q_1^2q_2^2}{(x-q_1^2-i0)^2(x-q_2^2-i0)^2}.
$$
This term vanishes at $q_1^2=0$ or $q_2^2=0$. At equal spacelike virtualities, the integral is $Q^4/[3(s_m+Q^2)^3]$, so the expression tends to the ALP coefficient divided by $Q^2$. For two hard spacelike virtualities at fixed nonzero ratio, the finite lower limit changes only subleading terms; the zero-threshold integral reproduces the asymptotic distribution-amplitude convolution. The tree-level statement does not supply the unequal-virtuality NLO Wilson coefficient or nonperturbative threshold effects. Keeping $s_m>0$ defines its analytic continuation with high-energy spectral support; the added high-energy discontinuity is not a pre-existing low-energy cut.

# Normalization and effective pole

For a conditional low-energy form with $F^{\rm low}(0,0)$ and singly virtual $1/Q^2$ coefficient, the spectral term vanishes when either photon is real. A symmetric effective pole with residue equal to the missing real-photon normalization and mass squared equal to the missing singly virtual coefficient divided by that residue enforces both conditions when the residue is nonzero. The mass squared must be positive for a physical stable pole. A vanishing residue or an incompatible sign requires a different spectral model. The pole changes the spectral strength and is a model for omitted contributions, not an exact determination. It is suppressed as $1/Q^4$ for equal hard virtualities and hence does not alter the above leading symmetric coefficient.

The explicit low-energy cut is preserved only if the added spectral term and effective pole have no support on its assigned low-energy spectral region. The continuum and pole describe new model strength; assigning them disjoint support from an actual low-energy density and subtracting any pole or continuum already contained in it are prerequisites, not consequences of the algebra. The tree benchmark has stable poles and the dispersive benchmark is only a component; they cannot simply be added without the overlap prescription discussed in [the mixed-channel overlap record](alp-tff-rho-isoscalar-overlap.md). The omission of mixed low/high spectral regions is a model assumption, not an established physical zero. Numerical matching remains conditional on the low-energy inputs.

# Verification and alternatives

Run `.venv/bin/python scripts/derive_alp_spectral_matching.py` and append `--test`. Ten symbolic checks verify the origin and the effective-pole normalization using the explicit two-vector tree benchmark and the symbolic double-spectral integral as inputs, photon exchange for the pole and continuum integrands, equal hard limit and its matrix-source coefficient, a fixed unequal spacelike ratio, the singly virtual pole constraint, spacelike continuity, variation of the continuum lower bound, and zero overlap of an illustrative pair of disjoint spectral supports. For the dispersive example the check treats its integral as a symbolic linear input; it does not evaluate the unknown density. The disjoint-support example does not establish the physical spectral density or its overlap. A sharp join or smooth interpolation restricted to Euclidean momenta can provide a comparison but is not by itself an analytic prescription preserving the physical cuts; see [the analyticity audit](alp-matching-analyticity-audit.md). The size of their differences is not evaluated here.
