---
title: One-ALP three-pseudoscalar interactions at LO in chiPT
status: derived
created: 2026-09-22
updated: 2026-09-22
bibtex_keys: [Alda:2026bar]
provenance: analytical-derivation
---

# One-ALP Three-Pseudoscalar Interactions at LO

## Scope and inputs

The calculation retains terms linear in the ALP and cubic in the pseudoscalar matrix $\Phi$ from

$$
\mathcal L_{\mathrm{LO}}=\frac{f_\pi^2}{4}
\left\langle u_\mu u^\mu+\chi_+\right\rangle.
$$

The definitions $\phi=\Phi/(2f_\pi)$, $\mathcal C_\phi[X]=[X,\phi]$, $\mathcal A_\phi[X]=\{X,\phi\}$, and the all-orders expressions for $u_\mu$ and $\chi_+$ are taken from [the Appendix A extract](../../extracts/Alda:2026bar/appendix-a-expansion.md) [@Alda:2026bar]. The executable derivation is [`scripts/derive_alp_phi3_lo_chipt.py`](../../../scripts/derive_alp_phi3_lo_chipt.py).

Assumptions and truncations:

- Only $O(1/f_a)$ terms are retained.
- The calculation is at leading chiral order and outside the Wess-Zumino-Witten sector.
- No flavour commutativity is assumed.
- Lorentz indices are contracted as displayed; flavour brackets denote a trace.

## Kinetic term

Separate the vielbein into its ALP-independent mesonic part and its one-ALP part. The terms needed at cubic meson order are

$$
u_\mu^{(0)}=-2\partial_\mu\phi
+\frac{1}{3}\mathcal C_\phi^2[\partial_\mu\phi]+O(\phi^5),
$$

$$
u_\mu^{(a,A)}=2\frac{\partial_\mu a}{f_a}c_q^A
-\frac{\partial_\mu a}{f_a}\mathcal C_\phi^2[c_q^A]+O(\phi^4),
$$

$$
u_\mu^{(a,V)}=-2i\frac{\partial_\mu a}{f_a}
\left(\mathcal C_\phi[c_q^V]-\frac{1}{6}\mathcal C_\phi^3[c_q^V]+\cdots\right).
$$

The vector series is odd in $\Phi$, while $u_\mu^{(0)}$ is also odd. Their products therefore contain only even powers of $\Phi$ and do not contribute at cubic order. The two axial cross terms give

$$
\left.\mathcal L_{\mathrm{kin}}\right|_{a\Phi^3}
=f_\pi^2\frac{\partial_\mu a}{f_a}
\left(
\left\langle \partial^\mu\phi\,\mathcal C_\phi^2[c_q^A]\right\rangle
+\frac{1}{3}\left\langle c_q^A\mathcal C_\phi^2[\partial^\mu\phi]\right\rangle
\right).
$$

Using trace cyclicity in the Appendix identity

$$
\left\langle X\mathcal C_\phi^2[Y]\right\rangle
=\left\langle \mathcal C_\phi^2[X]Y\right\rangle,
$$

the two traces coincide. Restoring $\Phi=2f_\pi\phi$ gives

$$
\boxed{
\left.\mathcal L_{\mathrm{kin}}\right|_{a\Phi^3}
=\frac{\partial_\mu a}{6f_af_\pi}
\left\langle c_q^A[\Phi,[\Phi,\partial^\mu\Phi]]\right\rangle .
}
$$

The double commutator is unchanged if the opposite convention for $\mathcal C_\phi$ is used.

## Symmetry-breaking term

Write

$$
\chi_S=2B_0\left(m_q+i\frac{a}{f_a}[c_q^S,m_q]\right),
\qquad
\chi_P=2B_0\frac{a}{f_a}\{c_q^P,m_q\}.
$$

At one-ALP order, the scalar source appears in $2\cos(\mathcal A_\phi)[\chi_S]$ and therefore multiplies only even powers of $\Phi$. It gives no $a\Phi^3$ term. The pseudoscalar source contributes through

$$
2\sin(\mathcal A_\phi)[\chi_P]
=2\mathcal A_\phi[\chi_P]
-\frac{1}{3}\mathcal A_\phi^3[\chi_P]+O(\phi^5).
$$

Consequently,

$$
\left.\mathcal L_{\chi}\right|_{a\Phi^3}
=-\frac{B_0a}{48f_af_\pi}
\left\langle
\{\Phi,\{\Phi,\{\Phi,\{c_q^P,m_q\}\}\}\}
\right\rangle.
$$

Repeated trace cyclicity yields

$$
\left\langle\mathcal A_\Phi^3[X]\right\rangle
=8\left\langle\Phi^3X\right\rangle,
$$

and hence

$$
\boxed{
\left.\mathcal L_{\chi}\right|_{a\Phi^3}
=-\frac{B_0a}{6f_af_\pi}
\left\langle\Phi^3\{c_q^P,m_q\}\right\rangle .
}
$$

## No-mixing parametrization

For the one-ALP prescription in the Appendix extract,

$$
\widehat c_q^A=-\mathbb M_a^{-1}(C),
\qquad
\widehat c_q^P=-c_q^{AP}-\mathbb M_a^{-1}(C).
$$

The direct cubic interactions therefore become

$$
\boxed{
\left.\mathcal L_{\mathrm{kin}}\right|_{a\Phi^3}
=-\frac{\partial_\mu a}{6f_af_\pi}
\left\langle\mathbb M_a^{-1}(C)
[\Phi,[\Phi,\partial^\mu\Phi]]\right\rangle,
}
$$

$$
\boxed{
\left.\mathcal L_{\chi}\right|_{a\Phi^3}
=\frac{B_0a}{6f_af_\pi}
\left\langle\Phi^3
\{c_q^{AP}+\mathbb M_a^{-1}(C),m_q\}\right\rangle.
}
$$
