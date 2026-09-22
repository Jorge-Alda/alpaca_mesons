---
title: ALP to three-meson amplitudes and decay widths
status: derived
created: 2026-09-22
updated: 2026-09-22
bibtex_keys: [Alda:2026bar]
provenance: analytical-derivation
---

# ALP to Three-Meson Amplitudes and Decay Widths

## Inputs and assumptions

The interaction terms are the three generic-generator component Lagrangians derived in [the generic-generator expansion](alp-three-meson-generic-generators.md). Their matrix-level origin and literature provenance are recorded in [the leading-order chiral derivation](alp-phi3-lo-chipt-derivation.md) and [the Appendix A extract](extracts/Alda:2026bar/appendix-a-expansion.md) [@Alda:2026bar]. The calculation is implemented in [`scripts/derive_alp_three_meson_amplitudes.py`](../scripts/derive_alp_three_meson_amplitudes.py).

Use natural units, the mostly-minus metric, and

$$
a(p)\longrightarrow \phi_1(p_1)\phi_2(p_2)\phi_3(p_3),
\qquad
s_{rs}=(p_r+p_s)^2.
$$

All particles are on shell and $p=p_1+p_2+p_3$. The amplitude convention is $S=1+iT$, so that the vertex is $i\mathcal M$. Couplings and flavour traces may in general be complex; squared amplitudes therefore mean absolute squares. The calculation includes only the direct leading-order interactions in the input record. It does not add external-leg mixing or higher-order chiral corrections.

Each width below is understood to vanish unless $m_a\ge m_1+m_2+m_3$.

The useful kinematic identities are

$$
p\cdot p_1=\frac{m_a^2+m_1^2-s_{23}}{2},\qquad
p\cdot p_2=\frac{m_a^2+m_2^2-s_{13}}{2},\qquad
p\cdot p_3=\frac{m_a^2+m_3^2-s_{12}}{2},
$$

and

$$
s_{12}+s_{13}+s_{23}=m_a^2+m_1^2+m_2^2+m_3^2.
$$

With the Fourier convention implicit in $S=1+iT$, a term $C\,\partial_\mu a\,\partial^\mu\phi_r$ contributes $C\,p\cdot p_r$ to the decay amplitude.

## Functional differentiation

Replacing a field by the sum of its labeled external legs and differentiating once with respect to every leg gives the vertex multiplicities. Define in addition

$$
T_{iii}^{P}=\left\langle\lambda_i^3P\right\rangle.
$$

For three identical fields, the coefficient of $x_1x_2x_3$ in $(x_1+x_2+x_3)^3$ is $3!=6$. For two identical fields, the coefficients extracted from the kinetic and symmetry-breaking terms are respectively

$$
2p\cdot p_3-p\cdot p_1-p\cdot p_2,
\qquad 2.
$$

No extra vertex multiplicity occurs for three different fields because the ordered-generator tensors already contain all field orderings inherited from the Lagrangian.

## Three identical mesons

For $a(p)\to\phi_i(p_1)\phi_i(p_2)\phi_i(p_3)$, the kinetic interaction vanishes and

$$
\boxed{
\mathcal M_{iii}=6M T_{iii}^{P}.
}
$$

Here and below $K=1/(6f_af_\pi)$ and $M=-B_0/(6f_af_\pi)$ as in the input record.

## Two identical and one different meson

Label $a(p)\to\phi_i(p_1)\phi_i(p_2)\phi_j(p_3)$, with $i\ne j$. Functional differentiation gives

$$
\mathcal M_{iij}
=K T_{iij}^{A}
\left(2p\cdot p_3-p\cdot p_1-p\cdot p_2\right)
+2M T_{iij}^{P}.
$$

Using $m_1=m_2=m_i$, $m_3=m_j$, and the Mandelstam sum rule yields

$$
\boxed{
\mathcal M_{iij}(s_{12})
=\frac{K}{2}T_{iij}^{A}
\left(m_a^2+3m_j^2-3s_{12}\right)
+2M T_{iij}^{P}.
}
$$

The result is symmetric under $p_1\leftrightarrow p_2$, as required.

## Three different mesons

Label $a(p)\to\phi_i(p_1)\phi_j(p_2)\phi_k(p_3)$, with $i,j,k$ pairwise distinct. The amplitude is

$$
\mathcal M_{ijk}
=K\left(
T_{ij;k}^{A}p\cdot p_3
+T_{ik;j}^{A}p\cdot p_2
+T_{jk;i}^{A}p\cdot p_1
\right)
+M T_{ijk}^{P}.
$$

Therefore

$$
\boxed{
\begin{aligned}
\mathcal M_{ijk}(s_{12},s_{13},s_{23})
={}&\frac{K}{2}\Big[
T_{ij;k}^{A}(m_a^2+m_k^2-s_{12})\\
&+T_{ik;j}^{A}(m_a^2+m_j^2-s_{13})
+T_{jk;i}^{A}(m_a^2+m_i^2-s_{23})
\Big]
+M T_{ijk}^{P}.
\end{aligned}
}
$$

Only two Mandelstam variables are independent. The Jacobi identity gives

$$
T_{ij;k}^{A}+T_{ik;j}^{A}+T_{jk;i}^{A}=0,
$$

which may be used together with the Mandelstam sum rule to rewrite the amplitude, but the symmetric three-variable form keeps the differentiated-field origin manifest.

## Mandelstam decay widths

For labeled masses $(m_1,m_2,m_3)$, define the Kallen function

$$
\lambda(x,y,z)=x^2+y^2+z^2-2xy-2xz-2yz.
$$

The physical Dalitz region can be parametrized by

$$
(m_1+m_2)^2\le s_{12}\le(m_a-m_3)^2,
$$

and

$$
\begin{aligned}
s_{23}^{\pm}(s_{12})={}&m_2^2+m_3^2
+\frac{(s_{12}+m_2^2-m_1^2)(m_a^2-s_{12}-m_3^2)}{2s_{12}}\\
&\pm\frac{\sqrt{\lambda(s_{12},m_1^2,m_2^2)
\lambda(m_a^2,s_{12},m_3^2)}}{2s_{12}}.
\end{aligned}
$$

For a spin-zero ALP and spin-zero mesons,

$$
\boxed{
\frac{d^2\Gamma}{ds_{12}\,ds_{23}}
=\frac{1}{S\,256\pi^3m_a^3}
\left|\mathcal M(s_{12},s_{23})\right|^2,
}
$$

where $S=\prod_r n_r!$ is the final-state symmetry factor. Hence

$$
\boxed{
\Gamma
=\frac{1}{S\,256\pi^3m_a^3}
\int_{(m_1+m_2)^2}^{(m_a-m_3)^2}ds_{12}
\int_{s_{23}^-}^{s_{23}^+}ds_{23}
\left|\mathcal M(s_{12},s_{23})\right|^2.
}
$$

The three cases use

$$
\begin{array}{c|c|c|c}
\text{channel}&(m_1,m_2,m_3)&S&\mathcal M\\ \hline
iii&(m_i,m_i,m_i)&3!&\mathcal M_{iii}\\
iij&(m_i,m_i,m_j)&2!&\mathcal M_{iij}(s_{12})\\
ijk&(m_i,m_j,m_k)&1&\mathcal M_{ijk}(s_{12},s_{13},s_{23})
\end{array}
$$

with $s_{13}=m_a^2+m_1^2+m_2^2+m_3^2-s_{12}-s_{23}$ in the double integral. For the first two cases the $s_{23}$ integration is immediate:

$$
\Gamma_{iii}
=\frac{1}{3!\,256\pi^3m_a^3}
\int_{4m_i^2}^{(m_a-m_i)^2}ds_{12}\,
\left[s_{23}^+(s_{12})-s_{23}^-(s_{12})\right]
\left|6M T_{iii}^{P}\right|^2,
$$

$$
\Gamma_{iij}
=\frac{1}{2!\,256\pi^3m_a^3}
\int_{4m_i^2}^{(m_a-m_j)^2}ds_{12}\,
\left[s_{23}^+(s_{12})-s_{23}^-(s_{12})\right]
\left|\frac{K}{2}T_{iij}^{A}(m_a^2+3m_j^2-3s_{12})
+2M T_{iij}^{P}\right|^2,
$$

where the Dalitz boundaries in each equation use that channel's masses.

## Verification

Run

```bash
.venv/bin/python scripts/derive_alp_three_meson_amplitudes.py
```

The script prints the field-polynomial inputs, extracted amplitudes, Mandelstam substitutions, phase-space boundaries, and widths in LaTeX. Assertions verify the factorial coefficients, the reduction of $\mathcal M_{iij}$ to a function of $s_{12}$, the Mandelstam sum rule, and the Dalitz-boundary width.
