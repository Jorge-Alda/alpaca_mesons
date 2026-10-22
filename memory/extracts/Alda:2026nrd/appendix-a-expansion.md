---
title: Appendix A CCWZ expansion formulas
status: literature
created: 2026-09-22
updated: 2026-10-22
bibtex_keys: [Alda:2026nrd]
provenance: literature
---

# Appendix A CCWZ Expansion Formulas

Define the dimensionless pseudoscalar matrix

$$
\phi = \frac{\Phi}{2f_\pi},
$$

and the commutator and anticommutator maps

$$
\mathcal C_\phi[A]=[\phi,A],\qquad
\mathcal A_\phi[A]=\{\phi,A\}.
$$

These conventions follow Appendix A, Eqs. (A.1)-(A.4), of `Alda:2026nrd` [@Alda:2026nrd].

The ordering is verified directly in the archived arXiv-v1
[`main_arXiv.tex`](../../sources/Alda:2026nrd/main_arXiv.tex), lines 2042--2046.
The explicit linear vielbein at line 2109 contains
$-i(\partial_\mu a)/(f_\pi f_a)[\bm\Phi,\bm{c_q^V}]$.
This fixes the sign of odd powers of the commutator map. Calculations that
previously combined the reversed map with the displayed sine formulas require
an individual audit; changing this definition does not establish that every
downstream derivation has been rechecked.

At linear order in the ALP, the all-orders CCWZ vielbein is

$$
u_\mu =
-2\frac{\sin\mathcal C_\phi}{\mathcal C_\phi}[\partial_\mu\phi]
+2\frac{\partial_\mu a}{f_a}\cos\mathcal C_\phi[c_q^A]
-2i\frac{\partial_\mu a}{f_a}\sin\mathcal C_\phi[c_q^V].
$$

This is Eq. (A.13) of `Alda:2026nrd` [@Alda:2026nrd].

Writing $\chi=\chi_S+i\chi_P$, the symmetry-breaking building block is

$$
\chi_+=2\cos\mathcal A_\phi[\chi_S]
+2\sin\mathcal A_\phi[\chi_P],
$$

where

$$
\chi_S=2B_0\left(m_q+i\frac{a}{f_a}[c_q^S,m_q]\right),
\qquad
\chi_P=2B_0\frac{a}{f_a}\{c_q^P,m_q\}.
$$

These formulas and their expansion through quadratic total field order are given in Appendix A and Eq. (A.5) of `Alda:2026nrd` [@Alda:2026nrd].

The same Appendix A definition fixes the complementary building block
through quadratic total field order. With
$\bm\chi=\bm\chi_S+i\bm\chi_P$ and
$\mathcal A_\phi[A]=\{\phi,A\}$, Eq. (A.7) gives
$$
\chi_-=2i\cos(\mathcal A_\phi)[\chi_P]
-2i\sin(\mathcal A_\phi)[\chi_S]
=2i\left(\chi_P-\{\phi,\chi_S\}
-\frac12\{\phi,\{\phi,\chi_P\}\}+\cdots\right).
$$
For one ALP and at most one meson, use
$\chi_0=2B_0m_q$,
$\chi_S=\chi_0+i(a/f_a)[c_q^S,\chi_0]$ and
$\chi_P=(a/f_a)\{c_q^P,\chi_0\}$. Thus
$$
\chi_-=\frac{2ia}{f_a}\{c_q^P,\chi_0\}
-2i\{\phi,\chi_0\}
+\frac{2a}{f_a}\{\phi,[c_q^S,\chi_0]\}
+O(\text{fields}^3,a^2).
$$
The sign of the linear-meson term is therefore explicitly negative in
$\chi_-$; this is the source of the positive $i\chi_-$ meson mass insertion
used in the odd-parity operators. [@Alda:2026nrd, App. A, Eqs. (A.1), (A.5)--(A.7)]

For one-ALP amplitudes, the paper prescribes the field redefinition

$$
\kappa=-c_q^A-\mathbb M_a^{-1}(C),\qquad \delta=-c_q^S,
$$

which gives

$$
\widehat c_q^A=-\mathbb M_a^{-1}(C),\qquad
\widehat c_q^P=-c_q^{AP}-\mathbb M_a^{-1}(C),\qquad
\widehat c_q^V=c_q^{VS},\qquad
\widehat c_q^S=0,
$$

and

$$
C=-B_0\{m_q,c_q^{AP}\}-\frac{m_0^2}{N_f}c_{GA}\mathbf 1.
$$

This is the no-ALP-meson-mixing prescription of Eqs. (III.30), (III.33), and Sec. III E of `Alda:2026nrd` [@Alda:2026nrd].

The invariant vector/scalar combination appearing here is defined in the
archived source, lines 264--274 and 290, by

$$
\bm{c_q^{VS}}=\bm{c_q^V}-\bm{c_q^S},\qquad
\bm{c_q^V}\longrightarrow\bm{c_q^V}+\bm\delta,\qquad
\bm{c_q^S}\longrightarrow\bm{c_q^S}+\bm\delta.
$$

The choice $\bm\delta=-\bm{c_q^S}$ and the resulting
$\bm{\widehat c_q^V}=\bm{c_q^{VS}}$, $\bm{\widehat c_q^S}=0$ are also
explicit at source lines 812--814. These are literature definitions, not a
new inference from the two-pseudoscalar amplitudes. [@Alda:2026nrd]
