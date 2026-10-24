---
title: Two-vector pion TFF constraints and odd-coefficient slope map
status: derived
created: 2026-10-24
updated: 2026-10-24
bibtex_keys: [Alda:2026nrd, Knecht:2001qf, Ruiz-Femenia:2003jdx]
provenance: analytical-derivation
---

# Scope and inputs

This record derives a two-vector-pole meromorphic diagnostic for the chiral
$\pi^0$ transition form factor. It is a candidate improvement to the original
one-vector tree TFF, not a matched replacement for the full physical
$\pi^0,\eta,\eta'$ system or the ALP form factor. The existing canonical
one-vector amplitude, phase, and original-coefficient dependence are in
[the normalized TFF record](radiative-tff-amplitude.md). The sourced pion
short-distance targets and the limitation of LMD versus LMD+V are in the
[Knecht extract](../../extracts/Knecht:2001qf/pion-tff-short-distance.md)
[@Knecht:2001qf]. The odd-operator basis is sourced in the linked amplitude
record [@Ruiz-Femenia:2003jdx].

Use natural units, $\epsilon^{0123}=+1$, the same external-pion phase as the
canonical amplitude, and $F_{\pi}(0,0)<0$. Following Alda:2026nrd, identify
the chiral-limit constant $F$ with $f_\pi=93\,\mathrm{MeV}=0.093\,\mathrm{GeV}$
in this convention [@Alda:2026nrd; see the source-backed
[neutral-meson input record](alda-neutral-amplitude-inputs.md)]. $N_C$ is the
number of colors, and $M_1,M_2$ are two vector-pole masses. The polynomial
numerator is truncated to the displayed terms; this does not define its
behavior beyond the stated tree meromorphic model.

# Executable derivation

Run [`scripts/derive_improved_pion_tff.py`](../../../scripts/derive_improved_pion_tff.py):

```text
.venv/bin/python scripts/derive_improved_pion_tff.py
```

SymPy constructs the Bose-symmetric ansatz

$$
F_\pi(s_1,s_2)=\frac{f_\pi\,[s_1s_2(s_1+s_2)+h_1(s_1+s_2)^2+h_2s_1s_2+h_5(s_1+s_2)+h_7]}{3\prod_{i=1}^2(M_i^2-s_1)(M_i^2-s_2)}.
$$

The anomaly fixes $h_7$ by requiring
$F_\pi(0,0)=-N_C/(12\pi^2f_\pi)$, which is
$-0.27237\,\mathrm{GeV}^{-1}$ at $N_C=3$ and $f_\pi=0.093\,\mathrm{GeV}$.
Expanding the single-virtual path
$(s_1,s_2)=(-Q^2,0)$, the absence of a constant requires $h_1=0$. Matching
the leading $1/Q^2$ coefficient to a symbolic source-specific input
$C_\pi^{\rm single}$ then gives

$$
h_5=-\frac{3C_\pi^{\rm single}M_1^2M_2^2}{f_\pi},\qquad
h_7=-\frac{N_CM_1^4M_2^4}{4\pi^2f_\pi^2}.
$$

The repository's sourced symmetric pion limit is reproduced:

$$
\lim_{Q^2\to\infty}Q^2F_\pi(-Q^2,-Q^2)=-\frac{2f_\pi}{3}
=-0.062\,\mathrm{GeV}.
$$

For unequal fixed hard virtualities $(-rQ^2,-Q^2)$, $r>0$, the same ansatz
instead predicts

$$
\lim_{Q^2\to\infty}Q^2F_\pi(-rQ^2,-Q^2)
=-\frac{f_\pi(r+1)}{3r}.
$$

This ray dependence is a model prediction, not a sourced factorization
coefficient. The unequal-ratio leading-twist input must be derived before
claiming full hard matching. The symmetric limit alone does not establish
that matching.

# Low-energy map to the original odd coefficients

Expanding the canonical original-eleven-operator expression at the origin,
the executable derivation obtains

$$
\left.\frac{\partial F_\pi}{\partial s_1}\right|_{0,0}
=-\frac{8F_V}{3f_\pi M_V^4}
\left[F_Vd_3-\sqrt2M_V(c_5-c_6)\right],
$$

and

$$
\left.\frac{\partial^2F_\pi}{\partial s_1\partial s_2}\right|_{0,0}
=-\frac{8F_V}{3f_\pi M_V^6}
\left[2F_Vd_3-\sqrt2M_V(c_1-c_2+c_5)\right].
$$

These derivatives are the two original-operator combinations that a
low-energy slope and mixed-curvature measurement can constrain in this
chiral pion channel. One may equate them to the corresponding derivatives
of the two-vector expression to preserve these local observables. That
produces a conditional slope equation and a value of $h_2$; the latter is
an optional low-energy matching assumption. Neither equation follows from
the leading OPE. The hard constraints leave $h_2$ free, while the anomaly
and single-virtual coefficient set $h_7$ and $h_5$ respectively.

# Interpretation and limits

This calculation replaces the incorrect constant single-virtual limit of
the original one-vector expression and reproduces the sourced symmetric
pion coefficient. It also makes explicit which local low-energy derivatives
carry the original $c_i,d_i$ combinations. It does not show that a
two-vector ansatz is sufficient: its unequal-hard-ray coefficient needs an
independent QCD comparison, and higher virtuality dependence may require
more terms or spectral strength.

The value of $C_\pi^{\rm single}$ is kept symbolic rather than copied from
an unstated convention. No physical $\eta,\eta'$ axial-decay-constant
matrix, singlet matching/evolution, or source-specific ALP hard coefficient
is supplied by this pion calculation. Do not transfer the pion inputs to
those sources. This distinction is required before fitting shared odd
coefficients with the proposed improved TFFs; see the
[TFF null-lifting recommendations](vector-tff-null-lifting-recommendations.md).
