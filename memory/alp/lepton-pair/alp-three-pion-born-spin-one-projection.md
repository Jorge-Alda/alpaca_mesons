---
title: Spin-one orientation projection of the isoscalar ALP three-pion Born current
status: derived
created: 2026-10-27
updated: 2026-10-27
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Scope and conventions

This derives the normalized total-orientation average of the exact-isospin
$I=0$ LO Born current in the spin-one channel selected by the vacuum
three-pion current. It uses the trace-based [Born projection](alp-three-pion-isoscalar-born-projection.md)
and the geometric current and Lorentz conventions in the
[three-pion discontinuity record](alp-three-pion-discontinuity.md), with the
no-mixing source convention of Alda:2026nrd [@Alda:2026nrd]. Let
$$
k=p+q=p_++p_-+p_0,\qquad x=k^2,
$$
and work in the $k$ rest frame. The average is normalized to one over
$SO(3)$ at fixed Dalitz invariants. This is a total three-pion orientation
projection; it is distinct from the pairwise $\pi\pi$ P-wave projection.
Explicitly, for a function $G$ of the orientation,
$$
\langle G\rangle_\Omega=\frac{1}{8\pi^2}
\int_0^{2\pi}d\alpha\int_0^\pi\sin\beta\,d\beta
\int_0^{2\pi}d\gamma\,G(\alpha,\beta,\gamma).
$$

The result is the overlap selected by the geometric vacuum current. It is
not a complete helicity decomposition of the Born amplitude and adds no
rescattering, direct $V3\pi$ term, or strict vector pole.

# Born input and angular variables

Set
$$
C=\frac{4e\langle c_q^AT_0\rangle}{3f_af_\pi},\qquad
D_r=q^2-2q^0E_r+2\mathbf q\cdot\mathbf p_r+i0,
$$
and use the cyclic triples $(r,s,t)=(+ ,0,-),(-,+,0),(0,-,+)$. The Born
current is
$$
J^\mu_{\mathrm{Born},I=0}
=C\sum_{(r,s,t)}(p\cdot p_s-p\cdot p_t)
\frac{(2p_r-q)^\mu}{D_r}.
$$
Here $p=k-q$, so in this frame
$$
p\cdot(p_s-p_t)=(\sqrt{x}-q^0)(E_s-E_t)
+\mathbf q\cdot(\mathbf p_s-\mathbf p_t).
$$
For each pion define
$$
a_r=q^2-2q^0E_r,\qquad b_r=2|\mathbf q|\,|\mathbf p_r|,
\qquad
I_r=\int_{-1}^{1}\frac{1-z^2}{a_r+b_rz+i0}\,dz.
$$
The analytic expression below is first obtained for $a_r>b_r>0$. Elsewhere
continue from the original denominator with its Feynman prescription; the
logarithms must use that same continuation.

# Orientation projection

Write the spatial geometric current as
$$
\mathbf J_{3\pi}=-\sqrt{x}\,(\mathbf p_+\times\mathbf p_-),
$$
which follows from $\epsilon^{0123}=+1$ and the lowered pion momenta in the
[discontinuity convention](alp-three-pion-discontinuity.md). SymPy performs
the azimuthal orientation integral and obtains the factor
$|\mathbf q|^2(1-z^2)/2$. The three cyclic tensor contractions have equal
ratio $2/\sqrt{x}$. The resulting overlap and tensor average are
$$
\left\langle \mathbf J_{3\pi}\cdot
  (\mathbf J_{\mathrm{Born},I=0}\times\mathbf q)\right\rangle_\Omega
=\frac{C|\mathbf q|^2|\mathbf J_{3\pi}|^2}{\sqrt{x}}
\sum_{r=+,-,0} I_r,
$$
$$
\left\langle (J_{3\pi})_i
(J_{\mathrm{Born},I=0})_j\right\rangle_\Omega
=\frac{C|\mathbf J_{3\pi}|^2}{2\sqrt{x}}
\left(\sum_{r=+,-,0} I_r\right)\epsilon_{ijk}q_k.
$$
Thus the current component reconstructed in this selected spin-one channel is
$$
\left(\mathbf J_{\mathrm{Born},I=0}\right)_{J=1}
=-\frac{3C}{2\sqrt{x}}
\left(\sum_{r=+,-,0} I_r\right)
(\mathbf J_{3\pi}\times\mathbf q),
\qquad
\left(J^0_{\mathrm{Born},I=0}\right)_{J=1}=0.
$$
The component is linear in $\langle c_q^AT_0\rangle$ and vanishes for an
isoscalar light-quark axial source.

# Closed angular integral

Polynomial division and integration give
$$
I_r=\frac{2a_r}{b_r^2}
+\frac{b_r^2-a_r^2}{b_r^3}
\left[\log(a_r+b_r)-\log(a_r-b_r)\right],
\qquad b_r\ne0,
$$
with the logarithms analytically continued from the original $D_r+i0$.
SymPy gives the regular limit
$$
\lim_{b_r\to0}I_r=\frac{4}{3a_r}.
$$
The initial derivation domain excludes angular pion-pole crossings; there the
prescribed continuation, rather than independent endpoint logarithms,
defines the integral.

# Reproducibility and limits

Run the [derivation script](../../../scripts/derive_alp_three_pion_born_spin_one.py)
and its [focused test](../../../scripts/test_alp_three_pion_born_spin_one.py):
```text
.venv/bin/python scripts/derive_alp_three_pion_born_spin_one.py
.venv/bin/python -m unittest discover -s scripts -p test_alp_three_pion_born_spin_one.py
```
The derivation script executes. Five focused tests pass, including the
analytic primitive, the $b_r\to0$ limit, cyclic orientation signs, an
independent four-vector epsilon-tensor sign calculation, and direct-current
angular quadrature at nonsingular kinematics. The projected time component
vanishes. The result determines only the component selected by the geometric
vacuum current. A full amplitude, its other helicity structures, rescattering,
and the spectral input remain open.
