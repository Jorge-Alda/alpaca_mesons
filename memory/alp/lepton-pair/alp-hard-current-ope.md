---
title: Leading local hard-current expansion with explicit ALP axial matrix elements
status: derived
created: 2026-10-25
updated: 2026-10-25
bibtex_keys: [Alda:2026nrd, Messerli:2025rnv]
provenance: analytical-derivation
---

# Scope, inputs, and assumptions

This is a direct free-quark derivation of the leading parity-odd local term
in the electromagnetic-current product. It advances the short-distance
part of the [isoscalar TFF plan](alp-isoscalar-tff-plan.md), without assigning
source-specific ALP matrix elements. The scalar convention and source
boundaries are those of the [tree TFF](alp-doubly-virtual-tff.md) and
[input-requirement record](alp-isoscalar-tff-input-gates.md)
[@Alda:2026nrd; @Messerli:2025rnv]. The bibliography entries were inspected.
No physical-meson asymptotic coefficient is used as an ALP input.

Use natural units, metric $(+---)$, $\epsilon^{0123}=+1$, and
$\gamma_5=i\gamma^0\gamma^1\gamma^2\gamma^3$. The electromagnetic
current excludes $e$:
$$
J_{\rm em}^{\mu}=\sum_f Q_f\bar q_f\gamma^\mu q_f.
$$
The calculation assumes a massless free hard-quark propagator and leading
order in $\alpha_s$. The momenta obey
$$
p=q_1+q_2,\qquad q=\frac{q_1-q_2}{2},\qquad
q_1=q+\frac p2,\quad q_2=-q+\frac p2.
$$
Take $-q^2$ large with fixed soft $p$ and quark masses negligible compared
with the hard scale. In this limit the two virtualities become nearly equal
and deep Euclidean. This is not a derivation for arbitrary unequal fixed
ratios or for a singly virtual photon.

# Executable Wick and Clifford derivation

The paired [script](../../../scripts/derive_alp_hard_current_ope.py) forms
both single contractions of the symmetrically placed currents,
$$
i\int d^4x\,e^{iq\cdot x}
T\{J_{\rm em}^{\mu}(x/2)J_{\rm em}^{\nu}(-x/2)\}.
$$
Use $S_F(k)=i\slashed k/(k^2+i0)$, expand the uncontracted quark fields
about the origin, and retain their local, derivative-free term. The two
contractions contain $S_F(q)$ and $S_F(-q)=-S_F(q)$. The Fourier prefactor
$i$ therefore gives
$$
-\frac1{q^2+i0}\sum_f Q_f^2\bar q_f
\bigl(\gamma^\mu\slashed q\gamma^\nu
      -\gamma^\nu\slashed q\gamma^\mu\bigr)q_f.
$$
There is no extra factor $N_c$: color is already summed in the remaining
quark bilinear. The disconnected vacuum contraction does not contribute to
the ALP matrix element considered here.

SymPy constructs the standard Dirac matrices and projects this difference
onto the axial Clifford basis using traces. It obtains
$$
\operatorname{tr}(\gamma^0\gamma^1\gamma^2\gamma^3\gamma_5)=-4i,
$$
$$
\gamma^\mu\slashed q\gamma^\nu
-\gamma^\nu\slashed q\gamma^\mu
=-2i\epsilon^{\mu\nu\rho\sigma}q_\rho\gamma_\sigma\gamma_5.
$$
Consequently the leading parity-odd local operator is
$$
\boxed{
i\int d^4x\,e^{iq\cdot x}
T\{J_{\rm em}^{\mu}(x/2)J_{\rm em}^{\nu}(-x/2)\}
\;\supset\;
\frac{2i}{q^2+i0}\epsilon^{\mu\nu\rho\sigma}q_\rho
\sum_f Q_f^2\bar q_f\gamma_\sigma\gamma_5q_f .
}
$$
The displayed sign follows the explicit propagator, Fourier prefactor,
gamma-five trace, and lowered momentum convention above. It is not an
independent phase choice.

For the three light flavors the axial combination is explicitly
$$
\frac49\bar u\gamma_\sigma\gamma_5u
+\frac19\bar d\gamma_\sigma\gamma_5d
+\frac19\bar s\gamma_\sigma\gamma_5s.
$$
Additional active flavors require their own matrix elements and threshold
matching. The full unbounded loop cannot be described by silently keeping
only three massless flavors at every hard scale.

# ALP matrix element and scalar extraction

Insert the vacuum-to-ALP matrix element of the displayed local operator.
In the existing current-tensor convention its leading contribution is
$$
T^{\mu\nu}\;\supset\;
\frac{2i}{q^2+i0}\epsilon^{\mu\nu\rho\sigma}q_\rho
\sum_fQ_f^2\langle0|\bar q_f\gamma_\sigma\gamma_5q_f|a(p)\rangle.
$$
Lorentz covariance makes each axial-current matrix element parallel to
$p_\sigma$. The script evaluates the nondegenerate ALP rest frame, contracts
with the lowered two-photon epsilon tensor, and divides by its scalar norm.
It never divides by a tensor. The resulting scalar, in the convention
$T^{\mu\nu}=\epsilon^{\mu\nu\rho\sigma}q_{1\rho}q_{2\sigma}
F_{a\gamma^*\gamma^*}$, is
$$
\boxed{
F_{a\gamma^*\gamma^*}\;\supset\;
\frac{2i}{(q^2+i0)m_a^2}
\sum_f Q_f^2p^\sigma
\langle0|\bar q_f\gamma_\sigma\gamma_5q_f|a(p)\rangle .
}
$$
This expression applies at $m_a>0$ in the stated hard limit; it is the
covariant form of the rest-frame result
$2i\langle0|\bar q\gamma^0\gamma_5q|a\rangle/(m_a q^2)$ for one
unit-charge flavor. It has dimension $\mathrm{GeV}^{-1}$ with relativistic
one-particle normalization. The apparent $m_a^{-2}$ is only the scalar
projection of a longitudinal matrix element, not a claim of a massless-ALP
singularity. No new decay-constant notation or numerical matrix element is
assigned.

The hard-momentum Ward contractions vanish before taking the ALP matrix
element. After using its longitudinal form, the leading tensor also obeys
both physical photon Ward contractions. These statements do not establish
the full operator Ward identities at nonzero soft momentum: derivative
operators, omitted here, are part of that expansion. Photon exchange
$\mu\leftrightarrow\nu$, $q\to-q$ leaves the leading term unchanged.

# Remaining matching and physical limits

The independent no-mixing $U_a,Y_a$ sources specify the supplied low-energy
description. They do not, by themselves, determine the renormalized QCD
matrix elements in the boxed expression. Their matching must specify each
ALP source, the axial-current normalization, scale and scheme, active-flavor
thresholds, and singlet mixing where relevant. The direct local photon
contribution is separate and has not been absorbed into this quark-current
coefficient. Higher-order Wilson coefficients and higher-dimensional
operators are not computed here.

The local expansion controls nearly symmetric hard photons at fixed soft
$p$. It does not determine the singly virtual limit or all unequal-ratio
hard limits; resummed derivative/light-cone information and its source
matrix elements may be required there. It also supplies no timelike
continuation through hadronic cuts. Thus the user-specified
$2m_e<m_a<2m_\tau$ with the full two-photon loop still requires the
remaining matching and continuation described in the plan. A high-energy
ALP TFF has not been obtained by expanding the retained resonance poles.

# Reproducibility and checks

Run
```text
.venv/bin/python scripts/derive_alp_hard_current_ope.py
.venv/bin/python scripts/derive_alp_hard_current_ope.py --test
```
Both invocations completed. Four unit tests verify the Clifford algebra and
gamma-five orientation, reconstruct every ordered gamma difference from the
trace-derived axial coefficients, check the Wick/Fourier sign and hard Ward
and Bose relations, and verify scalar reconstruction, physical photon Ward
identities after the ALP matrix element, and light-flavor weights. The
calculation uses standard SymPy matrix multiplication, traces,
Levi-Civita tensors, simplification, and LaTeX output.
