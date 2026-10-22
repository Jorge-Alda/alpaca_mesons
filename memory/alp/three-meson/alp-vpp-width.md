---
title: Invariant width for an ALP to a vector and two pseudoscalars
status: derived
created: 2026-10-21
updated: 2026-10-22
bibtex_keys: [Alda:2026nrd, Ecker:1988te, Ecker:1989yg, Ruiz-Femenia:2003jdx]
provenance: analytical-derivation
---

# Scope and inputs

This record reduces the strong tree amplitude in [the $a\to VPP$ amplitude
record](alp-vpp-amplitude.md) to a two-dimensional invariant phase-space
integral. Write $m_V=m_\omega$ or $m_\phi$ for the physical on-shell external
vector mass, so $k^2=m_V^2$. Keep $M_V$ for the distinct vector-mass parameter
in the RChT Lagrangian and in the amplitude coefficients. Use $m_V$, not
$M_V$, in kinematics, the polarization projector, and phase-space boundaries.
Internal exchanges
are coherent complex functions, with any line shape specified in those
functions. Natural units, metric $(+---)$, and the amplitude phases of the
parent record are used. Pseudoscalar amplitudes and generators follow the
Alda:2026nrd convention in the parent record: $f_\pi=93\,\mathrm{MeV}$,
$\Phi=\sum_i\lambda_i\phi_i$, and
$\operatorname{tr}(\lambda_i\lambda_j)=2\delta_{ij}$. Thus the width
formula here uses the already Alda-normalized coefficients from that record;
there is no extra generator rescaling in phase space or in the polarization
sum. The even scalar-source term is restricted to
$|\Delta S|=1$: in the $(u,d,s)$ basis retain only $ds,sd$ entries of
$c_q^{VS}$. No new literature input is used.

The executable derivation is
[`scripts/derive_alp_vector_two_pseudoscalar_width.py`](../../../scripts/derive_alp_vector_two_pseudoscalar_width.py):

```text
.venv/bin/python scripts/derive_alp_vector_two_pseudoscalar_width.py
```

# Invariant amplitude and polarization sum

Let $k$ be the vector momentum and $p_1,p_2$ the pseudoscalar momenta. Collect
the parent-record coefficients into
$$
\mathcal M=\varepsilon_V^\mu J_\mu,\qquad
J_\mu=F_{\rm odd}(s,t)\epsilon_{\mu\nu\rho\sigma}
p_1^\nu p_2^\rho k^\sigma+F_1(s,t)p_{1\mu}+F_2(s,t)p_{2\mu}.
$$
Here $F_{\rm odd}$ is the full coefficient of the odd structure, and $F_1,F_2$
are the coefficients of $\varepsilon_V\cdot p_1$ and
$\varepsilon_V\cdot p_2$ after inserting the restricted $c_q^{VS}$ into the
two even flavour traces. All exchange assignments contributing to each
coefficient are summed before squaring. $u$ is eliminated using
$$
s+t+u=m_a^2+m_V^2+m_1^2+m_2^2.
$$

The polarization projector is $-g^{\mu\nu}+k^\mu k^\nu/m_V^2$. The odd
vector is orthogonal to $k,p_1,p_2$, so odd-even interference vanishes
pointwise, while even-even interference remains. Define
$$
q_{ij}=p_i\cdot p_j-\frac{(k\cdot p_i)(k\cdot p_j)}{m_V^2},
\qquad
\Delta_G=\det\begin{pmatrix}
 m_V^2&k\cdot p_1&k\cdot p_2\\
k\cdot p_1&m_1^2&p_1\cdot p_2\\
k\cdot p_2&p_1\cdot p_2&m_2^2
\end{pmatrix}.
$$
Then
$$
\sum_{\rm pol}|\mathcal M|^2=
\Delta_G|F_{\rm odd}|^2-q_{11}|F_1|^2-q_{22}|F_2|^2
-2q_{12}\operatorname{Re}(F_1F_2^*).
$$
The Gram determinant is nonnegative in physical phase space. The signs in the
even terms follow from $q_{ij}$ as defined above, equivalently the transverse
polarization contraction $-p_i\cdot p_j+(k\cdot p_i)(k\cdot p_j)/m_V^2$.

# Two-variable width

For $s=(p_1+p_2)^2$, $t=(k+p_2)^2$, the width is
$$
\Gamma=\frac{1}{256\pi^3m_a^3 S_{PP}}
\int_{(m_1+m_2)^2}^{(m_a-m_V)^2}ds
\int_{t_-(s)}^{t_+(s)}dt\;\sum_{\rm pol}|\mathcal M(s,t)|^2,
$$
where $S_{PP}=2!$ for identical pseudoscalars and $S_{PP}=1$ otherwise.
The boundaries, with $\lambda(x,y,z)=x^2+y^2+z^2-2xy-2xz-2yz$, are
$$
t_\pm(s)=m_V^2+m_2^2+
\frac{(m_a^2-m_V^2-s)(s+m_2^2-m_1^2)
\pm\sqrt{\lambda(m_a^2,m_V^2,s)}\sqrt{\lambda(s,m_1^2,m_2^2)}}{2s}.
$$
Use the displayed pair as integration limits, ordered at each physical $s$ as
the smaller and larger value. The integral vanishes
below $m_a<m_V+m_1+m_2$. For identical pseudoscalars the amplitude must already
be Bose symmetrized; $S_{PP}$ then supplies the phase-space factor $1/2!$.

# Verification and limits

The executable script expresses the Gram contractions, Mandelstam elimination,
cosine-to-$t$ Jacobian, and phase-space normalization symbolically. Its
normalization follows $d\Phi_3=ds\,dt/(128\pi^3m_a^2)$, hence the displayed
width density after multiplication by $1/(2m_a)$. No spin average is made for
the spin-zero initial state. The result retains the parent amplitude's tree
truncation; it does not prescribe internal resonance line shapes or include
rescattering, loops, or electromagnetic conversions.
