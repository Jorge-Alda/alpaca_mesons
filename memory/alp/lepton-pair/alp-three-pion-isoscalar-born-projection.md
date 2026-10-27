---
title: Exact-isospin I=0 projection of the ALP three-pion Born current
status: derived
created: 2026-10-27
updated: 2026-10-27
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Scope and inputs

This derives the charged component of the $I=0$ projection of the full-charge
LO current for $a(p)+\gamma^*(q)\to\pi^+(p_+)\pi^-(p_-)\pi^0(p_0)$. The
starting Born current, including its source traces and photon coupling, is
given in the [LO current record](alp-three-pion-current-born-amplitude.md).
Use the no-mixing diagonal ALP sources of Alda:2026nrd
[@Alda:2026nrd], exact isospin in the pion triplet, $T_0=\lambda_3$,
common on-shell pion masses, and natural units. No isoscalar replacement of
the photon charge matrix is made.

The unique $I=0$ tensor in $1\otimes1\otimes1$ is fully antisymmetric. For
the charged component define the projector by
$$
P_0J^\mu=\frac16\sum_{\sigma\in S_3}\operatorname{sgn}(\sigma)
J^\mu(p_{\sigma(+)},p_{\sigma(-)},p_{\sigma(0)}),
$$
with even permutations carrying plus signs and odd permutations minus
signs. This is the charged component of the projector; it is not the matrix
element with a separately normalized $I=0$ three-pion state.

# Executable derivation

Run `.venv/bin/python scripts/derive_alp_three_pion_isospin_projection.py`.
The script obtains the source vertex from the existing Born derivation,
applies all six signed permutations, and contracts the resulting current
with $q_\mu$. It keeps independent pion and photon vector components until
the Ward contraction.

For the neutral pion leg $r$, the source vertex produced by the existing
traces is
$$
\mathcal V(r)=\frac{-B_0\langle X_PT_0\rangle
-2m_a^2\langle c_q^AT_0\rangle
+6(p\cdot p_r)\langle c_q^AT_0\rangle}{3f_af_\pi},
$$
where $X_P=\{c_q^P,m_q\}$ as in the [LO current
record](alp-three-pion-current-born-amplitude.md), and $T_0=\lambda_3$.
The six-term antisymmetrization cancels the constant pseudoscalar-source
part. With $D_r=q^2-2q\cdot p_r$, the executed result is
$$
\begin{aligned}
P_0J^\mu={}&\frac{4e\langle c_q^AT_0\rangle}{3f_af_\pi}\Bigg[
 (p\cdot p_0-p\cdot p_-)
 \frac{(2p_+-q)^\mu}{D_+}
 +(p\cdot p_+-p\cdot p_0)
 \frac{(2p_--q)^\mu}{D_-}\\
&\hspace{38mm} +(p\cdot p_--p\cdot p_+)
 \frac{(2p_0-q)^\mu}{D_0} \Bigg].
\end{aligned}
$$
Thus the projected charged component depends only on the axial-source trace
$\langle c_q^AT_0\rangle$. It vanishes when the light-quark part of
$c_q^A$ is isoscalar. It is transverse: using on-shell pions,
$q\cdot(2p_r-q)=-D_r$, SymPy gives $q_\mu P_0J^\mu=0$.

# Limits and verification

The paired derivation script and six focused tests pass. The tests check the
Ward identity, charged-pion exchange antisymmetry, cancellation of the
pseudoscalar-source constants, the isospin-symmetric axial-source zero, and
agreement between the matrix-trace and expanded source expressions.
The derivation is an isospin projection only: it does not perform the
$J=1$ angular projection needed for the spin-one unitarity integral, add
rescattering, or determine resonance and regular terms. It therefore does
not complete the physical isoscalar $a\gamma^*\to3\pi$ amplitude.
