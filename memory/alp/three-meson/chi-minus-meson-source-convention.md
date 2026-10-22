---
title: Appendix A derivation of the outgoing-meson chi-minus source sign
status: derived
created: 2026-10-22
updated: 2026-10-22
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Assumptions and inputs

This record settles the sign of the one-meson mass insertion directly from
Alda:2026nrd Appendix A. The source definitions are recorded in the
[Appendix A extract](../../extracts/Alda:2026nrd/appendix-a-expansion.md)
[@Alda:2026nrd, App. A, Eqs. (A.1), (A.5)--(A.7)]. The calculation is
linear in the ALP, through at most one meson field, and uses
$\chi_0=2B_0m_q$. RChT fields obey
$u=\exp(i\Phi_R/(\sqrt2F))$ and $U=u^2$; $F=f_\pi$ at leading order.
No equations of motion, quark-mass approximation, or field rotation beyond
the no-mixing source definitions is used.

The executable derivation is
[`scripts/derive_chi_minus_appendix_a.py`](../../../scripts/derive_chi_minus_appendix_a.py):

```text
.venv/bin/python scripts/derive_chi_minus_appendix_a.py
.venv/bin/python scripts/derive_chi_minus_appendix_a.py --test
```

# Appendix A expansion

The paper defines $\phi=\Phi/(2f_\pi)$,
$\mathcal A_\phi[A]=\{\phi,A\}$,
$\chi=\chi_S+i\chi_P$, and
$$
\chi_\pm=u^\dagger\chi u^\dagger\pm u\chi^\dagger u,
\qquad u=e^{i\phi}.
$$
Its Eq. (A.7) gives
$$
\begin{aligned}
\chi_+&=2\cos(\mathcal A_\phi)[\chi_S]
+2\sin(\mathcal A_\phi)[\chi_P],\\
\chi_-&=2i\cos(\mathcal A_\phi)[\chi_P]
-2i\sin(\mathcal A_\phi)[\chi_S]\\
&=2i\left(\chi_P-\{\phi,\chi_S\}
-\frac12\{\phi,\{\phi,\chi_P\}\}+\cdots\right).
\end{aligned}
$$
At the retained order,
$$
\chi_S=\chi_0+i\frac{a}{f_a}[c_q^S,\chi_0],
\qquad
\chi_P=\frac{a}{f_a}\{c_q^P,\chi_0\}.
$$
Substitution and field-order expansion yield
$$
\boxed{\chi_-=
\frac{2ia}{f_a}\{c_q^P,\chi_0\}
-2i\{\phi,\chi_0\}
+\frac{2a}{f_a}\{\phi,[c_q^S,\chi_0]\}
+O(\text{fields}^3,a^2).}
$$
The first and second displayed terms are the one-ALP and one-meson pieces.
The mixed ALP-meson term is quadratic in fields and vanishes in the
no-mixing convention $c_q^S=0$.

# RChT outgoing-meson source and discrepancy resolution

Use Alda:2026nrd's pseudoscalar field and generator convention,
$$
\Phi=\sum_P\lambda_P P,\qquad
\operatorname{tr}(\lambda_i\lambda_j)=2\delta_{ij},\qquad
u=\exp\!\left(\frac{i\Phi}{2f_\pi}\right),\qquad
f_\pi=0.093\,{\rm GeV}.
$$
For an outgoing pseudoscalar with coefficient $\lambda^\dagger$, the
coefficient of that field in $\chi_-$ is
$$
\chi_-|_P=-\frac{i}{f_\pi}\{\lambda^\dagger,\chi_0\},
\qquad
\boxed{i\chi_-|_P=+\frac{1}{f_\pi}
\{\lambda^\dagger,\chi_0\}.}
$$
This settles the discrepancy: the positive sign is authoritative for the
mass-source variable $\mathcal Y_P\equiv i\chi_-|_P$ used in the odd-parity
operator contractions. The earlier note mixed Alda's $f_\pi$ with the
RChT generator $T_P=\lambda_P/\sqrt2$, which has
$\operatorname{tr}(T_iT_j)=\delta_{ij}$. In that alternate normalization the
same formula reads $i\chi_-|_P=+\sqrt2\{T_P^\dagger,\chi_0\}/F$ with
$F=f_\pi$; the $\sqrt2$ belongs to that generator conversion and is absent
when using Alda's $\lambda_P$. The outgoing-meson source sign error in the
earlier $a\to VPP$ calculation is a separate error and has also been
corrected. Its ALP source is different and remains
$$
i\chi_-|_a=-\frac{2}{f_a}\{c_q^P,\chi_0\},
$$
so this correction does not change the ALP insertion or the $s$-channel
ALP-source numerator. It reverses the crossed pseudoscalar-source $d_2$
pieces in the $a\to VPP$ amplitude. The corrected sign is used in both
[the strong $a\to VPP$ record](alp-vpp-amplitude.md) and [the radiative
$a\to\gamma PP$ record](alp-photon-pp-amplitude.md).

# Verification

The SymPy script multiplies the noncommuting exponential series in the
definition of $\chi_\pm$, extracts total field degrees zero through two,
substitutes the Alda:2026nrd spurions, and maps the one-meson coefficient to
RChT normalization. Three checks pass: agreement with the Appendix A sine
and cosine expansions; an independent commuting exact-exponential limit;
and the no-mixing-basis meson-sign check. The analytic sign also agrees with
the source-backed RChT expansion of $\chi_-$ used in the
[radiative SM vertex derivation](../../rcht/odd-parity/radiative-sm-amplitude.md).
