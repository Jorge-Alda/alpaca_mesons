---
title: Mandelstam representation and integration reduction for ALP to four pseudoscalars
status: derived
created: 2026-10-23
updated: 2026-10-23
bibtex_keys: [Alda:2026nrd, Ecker:1988te, Ecker:1989yg, ParticleDataGroup:2024cfk, Ruiz-Femenia:2003jdx]
provenance: analytical-derivation
---

# Scope and inputs

This record gives the invariant phase space and width for the retained
tree-level axial/pseudoscalar-source amplitude in the
[four-pseudoscalar source derivation](alp-four-pseudoscalar-axial-pseudoscalar.md).
The amplitude model, no-mixing source convention, pseudoscalar and vector
normalizations, and exclusions are those stated there. In particular,
electromagnetic interactions, loops, and interactions outside the retained
vector-exchange fit truncation are not included. The propagators use the
constant-width line shape from the linked $a\to3P$ exchange record, with
PDG total-width values; this is a phenomenological implementation choice,
not a claim that the PDG specifies a universal Breit--Wigner model
[@ParticleDataGroup:2024cfk]. We use natural units, metric $(+---)$, and
$\epsilon^{0123}=+1$.

The symbolic derivation is implemented in
[`scripts/derive_alp_four_pseudoscalar.py`](../../../scripts/derive_alp_four_pseudoscalar.py):

```text
.venv/bin/python scripts/derive_alp_four_pseudoscalar.py --phase-space
```

# Mandelstam form of the amplitude and width

Write the coherent amplitude as
$$
\mathcal M=\epsilon(p_1,p_2,p_3,p_4)F(s_{ij}),\qquad
s_{ij}=(p_i+p_j)^2,
$$
where $F$ is the full sum of the $aVV$ and axial $aVP$ contributions in the
source derivation. Momentum conservation fixes
$$
\sum_{i<j}s_{ij}=m_a^2+2\sum_{i=1}^4m_i^2.
$$
Thus five of the six pair invariants are independent. Let
$$
G_{ij}=p_i\cdot p_j=\begin{cases}
m_i^2,&i=j,\\
\dfrac{s_{ij}-m_i^2-m_j^2}{2},&i\ne j.
\end{cases}
$$
Then
$$
|\mathcal M|^2=-\det G\,|F(s_{ij})|^2.
$$
This removes all explicit four-vectors and Levi-Civita tensors from the
squared amplitude.

For distinguishable labeled particles, the Lorentz-invariant width is
$$
\Gamma=\frac{1}{16384\pi^6m_a^3S}\int_{\mathcal D}
ds_{12}\,ds_{34}\,dz_1\,dz_2\,du\,
\sqrt{-\det G}\,|F(s_{ij})|^2,
$$
where $u=s_{13}$, $z_1=s_{13}+s_{14}$, $z_2=s_{13}+s_{23}$, and
$S=\prod_h n_h!$ is the identical-particle factor for each species $h$
with multiplicity $n_h$. The sixth pair invariant is eliminated by the
sum rule. This is the measure after the two orientations corresponding to
the azimuthal angle have been included. For four identical particles,
$S=4!$.

The outer pair masses range over
$$
(m_1+m_2)^2\le s_{12}\le(m_a-m_3-m_4)^2,
\qquad
(m_3+m_4)^2\le s_{34}\le(m_a-\sqrt{s_{12}})^2.
$$
At fixed $s_{12},s_{34}$, the two summed invariants range independently over
$$
z_1^{\rm c}-r_1\le z_1\le z_1^{\rm c}+r_1,\qquad
z_2^{\rm c}-r_2\le z_2\le z_2^{\rm c}+r_2,
$$
with the centers and half-ranges
$$
\begin{aligned}
z_1^{\rm c}&=2m_1^2+m_3^2+m_4^2+
\frac{(s_{12}+m_1^2-m_2^2)(m_a^2-s_{12}-s_{34})}{2s_{12}},\\
r_1&=\frac{\sqrt{\lambda(m_a^2,s_{12},s_{34})\,
\lambda(s_{12},m_1^2,m_2^2)}}{2s_{12}},\\
z_2^{\rm c}&=m_1^2+m_2^2+2m_3^2+
\frac{(s_{34}+m_3^2-m_4^2)(m_a^2-s_{12}-s_{34})}{2s_{34}},\\
r_2&=\frac{\sqrt{\lambda(m_a^2,s_{12},s_{34})\,
\lambda(s_{34},m_3^2,m_4^2)}}{2s_{34}},
\end{aligned}
$$
where $\lambda(x,y,z)=x^2+y^2+z^2-2xy-2xz-2yz$. At fixed five
coordinates, physical momenta require the Gram matrix to have the physical
signature and $-\det G\ge0$.

# Exact inner invariant and reduction in integration dimension

At fixed $s=s_{12}$, $t=s_{34}$, $z_1$, and $z_2$, the determinant is
quadratic in $u=s_{13}$. If
$$
-\det G=A u^2+B u+C,
$$
then the executable derivation obtains $A=-\lambda(m_a^2,s,t)/16$. Define
$$
c=-\frac{B}{2A},\qquad
r^2=\frac{B^2-4AC}{4A^2}.
$$
The physical interval is $c-r\le u\le c+r$, and
$$
-\det G=\frac{\lambda(m_a^2,s,t)}{16}
\bigl[r^2-(u-c)^2\bigr].
$$
Consequently, with $u=c+r\cos\theta$, the inner part of the width integral
is proportional to
$$
\frac{\sqrt{\lambda(m_a^2,s,t)}}{4}r^2
\int_0^\pi d\theta\,\sin^2\theta\,
|F(c+r\cos\theta)|^2.
$$

This proves a conditional exact reduction. If $F$ is independent of $u$,
the last integral is $\pi\sqrt{\lambda(m_a^2,s,t)}r^2|F|^2/8$. If
$|F(c+r\cos\theta)|^2$ is polynomial in $\cos\theta$, the moments
$$
\int_0^\pi\sin^2\theta\cos^{2n}\theta\,d\theta
=\frac{\pi\,C_n}{2^{2n+1}},\qquad
\int_0^\pi\sin^2\theta\cos^{2n+1}\theta\,d\theta=0,
$$
with $C_n=\binom{2n}{n}/(n+1)$, give the integral algebraically.

For a rational $F$ with constant pole parameters and no pole on the interval,
partial fractions reduce the same integral to polynomial moments and
$$
J(z)=\int_{-r}^{r}\frac{\sqrt{r^2-x^2}}{x-z}\,dx
=\pi\bigl(\sqrt{z^2-r^2}-z\bigr),
$$
where the square root is chosen to approach $z$ at infinity. A pole of order
$k$ follows by differentiating $J(z)$ $k-1$ times and dividing by
$(k-1)!$. This applies to the full coherent rational tree amplitude,
including crossed terms, with the constant-width propagators specified in
the linked source record. It reduces the numerical integration from five
dimensions to four. It is an analytic method in principle; the resulting
partial-fraction expressions can be ill-conditioned near coincident poles
and need stable limiting formulas in a production integrator. The finite
widths move vector poles away from the real integration domain, removing the
zero-width on-shell divergence. Energy-dependent widths and nonrational form
factors require a separate inner-integration analysis.

# Special two-dimensional case

For a single $aVV$ pairing, $F=F(s_{12},s_{34})$ depends only on the two
pair masses, so all three relative angles can be integrated. Its contribution
is
$$
\frac{d\Gamma}{ds_{12}\,ds_{34}}=
\frac{\lambda(m_a^2,s_{12},s_{34})^{3/2}
\lambda(s_{12},m_1^2,m_2^2)^{3/2}
\lambda(s_{34},m_3^2,m_4^2)^{3/2}}
{1179648\pi^5m_a^3s_{12}^2s_{34}^2S}
|F(s_{12},s_{34})|^2.
$$
This is not the generic answer for the complete amplitude: the three
crossed $aVV$ pairings interfere, and the axial $aVP$ terms contain
three-particle invariants $(p_a-p_i)^2$. Those terms generally restore
dependence on all five independent Mandelstam variables. Bose symmetry does
not remove that dependence.

# Inputs and provenance

The amplitude and vector-tree truncation are taken from the linked source
derivation and fit record. The line shape and PDG width inputs follow the
[$a\to3P$ vector-exchange record](../three-pseudoscalar/a3p-vector-exchange-amplitude.md)
and the [light-vector input record](../../extracts/ParticleDataGroup:2024cfk/light-vector-pole-parameters.md)
[@ParticleDataGroup:2024cfk]. The five-invariant Gram determinant,
phase-space normalization, physical coordinate bounds, and conditional
analytic integration are derived in the cited SymPy script. The four-
dimensional reduction assumes the stated constant widths and does not apply
unchanged to energy-dependent widths.
