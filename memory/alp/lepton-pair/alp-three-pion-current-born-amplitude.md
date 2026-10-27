---
title: LO ALP electromagnetic-current amplitude to three pions
status: derived
created: 2026-10-25
updated: 2026-10-27
bibtex_keys: [Alda:2026nrd, Hoferichter:2025lcz]
provenance: analytical-derivation
---

# Scope and conventions

This record derives the gauged leading-order chiral contribution to the
crossed current matrix element
$$
\langle\pi^+(p_+)\pi^-(p_-)\pi^0(p_0)|J^\mu_{\rm em}|a(p)\rangle,
$$
and hence to $a(p)+\gamma^*(q)\to\pi^+\pi^-\pi^0$. It is not the same
matrix element as the anomalous current amplitude
$$
\langle\pi^+\pi^-\pi^0|J^\mu_{\rm em}|0\rangle
=-\epsilon^{\mu\nu\alpha\beta}p_{+\nu}p_{-\alpha}p_{0\beta}
\mathcal F(s,t,u;q^2),
$$
used in the three-pion Khuri--Treiman treatment [@Hoferichter:2025lcz]. The
ALP matrix element here has ordinary-parity vector structures and depends on
the ALP and photon momenta in addition to the three-pion Dalitz variables.
The Khuri--Treiman scalar cannot be reused as this production amplitude.

Use natural units, metric $(+---)$, incoming $p$ and $q$, outgoing pion
momenta, and $p+q=p_++p_-+p_0$. The photon amplitude includes $e$; divide by
$e$ for the current matrix element. The diagonal ALP source matrices
$c_q^A,c_q^P$ are inputs already in the no-mixing basis of Alda:2026nrd
[@Alda:2026nrd]. Write
$$
X_P=\{c_q^P,m_q\},\qquad
T_+=\frac{\lambda_1+i\lambda_2}{\sqrt2},\qquad
T_-=T_+^\dagger,
$$
where $\lambda_1,\lambda_2$ are Gell-Mann matrices normalized by
$\operatorname{tr}(\lambda_i\lambda_j)=2\delta_{ij}$. Denote the physical
neutral-pion generator by $T_0$. The source profile is not fixed to a
benchmark. Expand
$$
\Phi=\pi^+T_++\pi^-T_-+\pi^0T_0.
$$

# LO derivation

Start from the one-ALP, three-Goldstone terms in the LO chiral Lagrangian,
already derived in the [cubic-interaction record](../three-pseudoscalar/alp-phi3-lo-chipt-derivation.md):
$$
\mathcal L_{a\Phi^3}=
\frac{\partial_\mu a}{6f_af_\pi}
\left\langle c_q^A[\Phi,[\Phi,\partial^\mu\Phi]]\right\rangle
-\frac{B_0a}{6f_af_\pi}
\left\langle\Phi^3\{c_q^P,m_q\}\right\rangle.
$$
The coefficients of $\partial\pi^+$, $\partial\pi^-$ and $\partial\pi^0$
in the kinetic trace are the following matrix traces:
$$
\begin{aligned}
K_+&=\left\langle c_q^A\left(
[T_0,[T_-,T_+]]+[T_-,[T_0,T_+]]
\right)\right\rangle,\\
K_-&=\left\langle c_q^A\left(
[T_0,[T_+,T_-]]+[T_+,[T_0,T_-]]
\right)\right\rangle,\\
K_0&=\left\langle c_q^A\left(
[T_+,[T_-,T_0]]+[T_-,[T_+,T_0]]
\right)\right\rangle.
\end{aligned}
$$
These are the full coefficients obtained by expanding
$\langle c_q^A[\Phi,[\Phi,\partial^\mu\Phi]]\rangle$. For diagonal
$c_q^A$ and $T_0=\lambda_3$, the explicit SymPy trace evaluation gives
$$
K_+=K_-=-4\langle c_q^AT_0\rangle,
\qquad K_0=-2K_+,
\qquad M_{+-0}=2\langle X_PT_0\rangle.
$$
The coefficient of
$\pi^+\pi^-\pi^0$ in the mass-source trace is
$$
M_{+-0}=\sum_{\sigma\in S_3}
\operatorname{tr}\!\left[
X_P T_{\sigma(1)}T_{\sigma(2)}T_{\sigma(3)}
\right],
\qquad (T_1,T_2,T_3)=(T_+,T_-,T_0).
$$

For either charged-pion pole graph, shift the momentum entering that charged
leg of the ALP vertex by $p_\pm\mapsto p_\pm-q$. Use the incoming-photon
identity
$$
p\cdot p_-=m_a^2+p\cdot q-p\cdot p_+-p\cdot p_0
$$
to express the two shifted source vertices in common invariants. In trace
notation they are
$$
\mathcal V_+=\mathcal V_-\equiv\mathcal V_{a\pi^+\pi^-\pi^0}
=\frac{1}{6f_af_\pi}
\left[K_+(m_a^2-3p\cdot p_0)-B_0M_{+-0}\right].
$$

Gauging the charged derivatives produces a contact-current term
proportional to $K_+-K_-$. Before applying the evaluated trace identities,
the current is
$$
J^\mu_{a\to3\pi,\rm Born}=e\left[
\mathcal V_+\frac{(2p_+-q)^\mu}{D_+}
-\mathcal V_-\frac{(2p_--q)^\mu}{D_-}
-\frac{K_+-K_-}{6f_af_\pi}p^\mu
\right],
\qquad D_\pm=q^2-2q\cdot p_\pm+i0.
$$
The evaluated identities $K_+=K_-$ and $\mathcal V_+=\mathcal V_-$ then
reduce this to
$$
J^\mu_{a\to3\pi,\rm Born}=e\,\mathcal V_{a\pi^+\pi^-\pi^0}
\left[\frac{(2p_+-q)^\mu}{D_+}
-\frac{(2p_--q)^\mu}{D_-}\right].
$$
On the external pion mass shell,
$q\cdot(2p_\pm-q)=-D_\pm$, so
$$
q_\mu J^\mu_{a\to3\pi,\rm Born}
=-e\mathcal V_{a\pi^+\pi^-\pi^0}
+e\mathcal V_{a\pi^+\pi^-\pi^0}=0.
$$
This Ward identity follows after independently constructing both shifted
vertices and the contact term; the zero contact coefficient is a result, not
an omitted diagram assumption.

# Isoscalar-current limit and interpretation

On the light-quark subspace, the isoscalar electromagnetic generator is
$Q^{(0)}_\ell=\frac16\mathbf1_2$. Since
$[Q^{(0)}_\ell,T_\pm]=0$, it gives no minimal gauged LO pion current. The full
electromagnetic generator also has an isovector component; the Born
expression above is not by itself an isoscalar spectral projection.
Separately, for isoscalar source matrices and the ideal isovector neutral
generator $T_0=\lambda_3$, the nonradiative source trace vanishes. Neither
zero implies that the full three-pion contribution to the ALP TFF vanishes.

The first isoscalar resonance contribution can arise from an ALP--vector--
virtual-photon vertex followed by a vector-to-three-pion amplitude. It has two
epsilon tensors before contraction and is not contained in this LO Born
result. The existing [on-shell $a\to V\gamma$ record](../../rcht/odd-parity/alp-vgamma-amplitude.md)
does not give the required off-shell vertex. The available dispersive
$V\to3\pi$ framework supplies a subamplitude, but the full ALP production
amplitude still needs its kinematic projection and source-dependent
normalization. No resonance-dressed or Khuri--Treiman solution is claimed
here.

# Reproducibility

Run:
```text
.venv/bin/python scripts/derive_alp_virtual_photon_three_pion_born.py
```
The executed SymPy derivation prints the matrix-trace source formula and
reports the contact coefficient, shifted-vertex difference, and Ward
contraction in LaTeX. It returns zero for each of those quantities. No
numerical ALP point, rescattering input,
resonance term, or dispersive subtraction is selected.
