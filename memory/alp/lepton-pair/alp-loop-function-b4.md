---
title: Triangle loop function B_4 for a pseudoscalar coupled to two photons through a fermion loop
status: derived
created: 2026-11-15
updated: 2026-11-15
bibtex_keys: [Bauer:2021mvw, Alda:2024cxn, Alda:2026nrd]
provenance: analytical-derivation
---

# Definition and notation

For a fermion of mass $m_h$ in the loop, an ALP of mass $m_a$ and photon virtualities $q_1^2,q_2^2$ (momenta $q_1,q_2$ with $q_1+q_2=p$, $p^2=m_a^2$), set
$$
\tau_a=\frac{4m_h^2}{m_a^2},\qquad \tau_1=\frac{4m_h^2}{q_1^2},\qquad \tau_2=\frac{4m_h^2}{q_2^2},
$$
with $\tau_i\to\infty$ for a real photon ($q_i^2=0$). The function $B_4$, defined here, is
$$
\boxed{B_4(\tau_a,\tau_1,\tau_2)=1-\int_0^1\!dx\int_0^{1-x}\!dy\;
\frac{2m_h^2}{m_h^2-xy\,m_a^2-x(1-x-y)\,q_2^2-y(1-x-y)\,q_1^2-i0}.}
$$
The ALP comes first and the photon virtualities follow, as in the one-real-one-virtual function $B_3(\tau_1,\tau_2)$ of
[@Bauer:2021mvw, Eq. (2.82)] and [@Alda:2024cxn, Eq. (A1)] (arguments $4m_f^2/m_a^2$ and $4m_f^2/\text{virtuality}$). The letter $B_3$ is already used there, so the function with two virtual
photons is called $B_4$; $B_1$ is the real-photon function of [@Bauer:2021mvw, Eq. (2.71)] (also the $B_1$ in the $c_{gg}$ formula of [@Alda:2026nrd]).
The prescription $-i0$ is $m_h^2\to m_h^2-i0$, from the Feynman propagator.

It enters the heavy-quark contribution to the two-photon form factor
(see the [triangle record](alp-heavy-quark-triangle.md)):
$$
F_h=\frac{2N_cQ_h^2}{4\pi^2f_a}\Bigl[c_{hh}^P+\bigl(c_{hh}^A-c_{hh}^P\bigr)B_4(\tau_a,\tau_1,\tau_2)\Bigr],
$$
identical to the form with the integral written out (checked symbolically). The [script](../../../scripts/derive_alp_heavy_quark_triangle.py) (Part 8) executes every property below.

# Properties

1. **Full symmetry.** $B_4$ is symmetric under the exchange of any two of its arguments. With $z=1-x-y$ the three invariants pair with the three parameter products $xy$, $xz$, $yz$;
   relabeling the parameters ($x\leftrightarrow y$ for $\tau_1\leftrightarrow\tau_2$, $x\leftrightarrow z$ for $\tau_a\leftrightarrow\tau_1$, $y\leftrightarrow z$ for $\tau_a\leftrightarrow\tau_2$) leaves the integrand
   unchanged, with unit Jacobian and the same domain. Numerically all orderings agree to $10^{-11}$ (spacelike and timelike below-threshold values) and to $10^{-10}$ with one invariant above threshold (complex mass).
2. **The function $f$.**
$$
f(\tau)=\begin{cases}\arcsin\dfrac1{\sqrt\tau},&\tau\ge1,\\[2mm]\dfrac\pi2+\dfrac i2\ln\dfrac{1+\sqrt{1-\tau}}{1-\sqrt{1-\tau}},&0<\tau<1,\\[2mm]\dfrac i2\ln\dfrac{\sqrt{1-\tau}+1}{\sqrt{1-\tau}-1},&\tau<0,\end{cases}
$$
   the first two as in [@Bauer:2021mvw, Eq. (2.71)]. The form for $\tau<0$ (spacelike) is derived in the script from the closed form for one spacelike invariant; $f^2$ is real there. The printed expression for $0<\tau<1$ continued
   to $\tau<0$ gives $-i\pi$ instead of $+i\pi$ from $m_h^2\to m_h^2-i0$ and then $f^2$ is not real, so the explicit form for $\tau<0$ must be used. Only $f^2$ enters, so the overall sign of $f$ is a convention.
3. **Reductions to the literature.**
$$
B_4(\tau_a,\tau_1,\infty)=B_3(\tau_a,\tau_1)=1+\frac{\tau_a\tau_1}{\tau_a-\tau_1}\bigl[f^2(\tau_a)-f^2(\tau_1)\bigr],\qquad
B_4(\tau_a,\infty,\infty)=B_1(\tau_a)=1-\tau_af^2(\tau_a).
$$
   The reduction to $B_3$ is proven: with $x=(1-y)\xi$ and the substitution $y\to1-y$ the integrand becomes $m_h^2/\{m_h^2-y(1-y)[\xi m_a^2+(1-\xi)q_1^2]\}$, which is the integrand of the one-real-photon
   expression of [@Bauer:2021mvw, Eq. (2.73)], and the integral over $\xi$ gives $B_3$. It is also checked numerically against the two-dimensional integral and Eq. (2.73) (real mass: $10^{-10}$; complex mass and a
   one-dimensional form above threshold: $10^{-9}$ and $10^{-12}$). With symmetry, a single finite argument in any slot gives $B_1$ of that argument.
4. **Closed forms for spacelike photons ($\tau<0$).**
$$
B_4(\infty,\tau,\tau)=1+\frac{\tau}{2\sqrt{1-\tau}}\ln\frac{\sqrt{1-\tau}+1}{\sqrt{1-\tau}-1},\qquad
B_4(\infty,\tau,\infty)=B_1(\tau)=1+\frac\tau4\ln^2\frac{\sqrt{1-\tau}+1}{\sqrt{1-\tau}-1},
$$
   for $m_a\to0$ with equal virtualities, and for $m_a\to0$ with one real photon; these equal the closed forms in $|q^2|$ of the [triangle record](alp-heavy-quark-triangle.md).
5. **Limits.** Heavy fermion ($m_h\to\infty$):
$$
B_4=-\frac{m_a^2+q_1^2+q_2^2}{12m_h^2}-\frac{m_a^4+q_1^4+q_2^4+m_a^2q_1^2+m_a^2q_2^2+q_1^2q_2^2}{90\,m_h^4}+\dots
=-\frac13\sum_i\frac1{\tau_i}-\frac8{45}\Bigl(\sum_i\frac1{\tau_i^2}+\sum_{i<j}\frac1{\tau_i\tau_j}\Bigr)+\dots
$$
   (sums over $i\in\{a,1,2\}$). Light fermion or large virtuality ($m_h\to0$ or $|q^2|\to\infty$ at fixed $m_h$): $B_4\to1$, with
   $B_4-1\simeq-\frac{2m_h^2}{|q^2|}\ln\frac{|q^2|}{m_h^2}$ for $q_1^2=q_2^2=-|q^2|$ and $m_a^2=0$;
   $-\frac{m_h^2}{|q^2|}\ln^2\frac{|q^2|}{m_h^2}$ for $q_1^2=-|q^2|$, $q_2^2=0$, $m_a^2=0$; and a term $-\frac{m_a^2}{3|q^2|}$ from the first order in $m_a^2$ at equal virtualities (this expansion requires $m_a^2<4m_h^2$).
   In the single-heavy-flavor decoupling language: $B_4\to0$ for $m_h\to\infty$ (the derivative coupling decouples) and $B_4\to1$ for $m_h\to0$ (the anomaly constant).
6. **Imaginary part.** $\operatorname{Im}B_4=-\operatorname{Im}(\text{integral})\le0$ above threshold. For real photons with $0<\tau_a<1$ and $\beta=\sqrt{1-\tau_a}$,
$$
B_4=1-\tau_a\Bigl(\frac\pi2+\frac i2\ln\frac{1+\beta}{1-\beta}\Bigr)^2,\quad
\operatorname{Re}B_4=1-\frac{\tau_a}4\Bigl(\pi^2-\ln^2\frac{1+\beta}{1-\beta}\Bigr),\quad
\operatorname{Im}B_4=-\frac{\pi\tau_a}2\ln\frac{1+\beta}{1-\beta}.
$$
   It is continuous at $m_a=2m_h$ with a square-root behavior of the derivative. For spacelike photons with $m_a>2m_h$ the imaginary part is minus the one-dimensional closed form
   of the [triangle record](alp-heavy-quark-triangle.md) (Part 7 of the script).

# Limits of what is verified

- With all three invariants nonzero and one of them above threshold there is no closed form here; only the symmetry is checked.
- At real mass above threshold only the reductions to $B_3$ and $B_1$ are checked directly (one-dimensional integrals); the two-dimensional integral is checked there at complex mass.
- The large-$|q^2|$ behavior with $m_a\neq0$ is supported numerically ($|q^2|=10^2$ and $10^4$ GeV$^2$), not proven in general. The heavy-mass expansion is tested to $10^{-7}$ at $m_h=10$ GeV.
- The equal-argument form $B_3(\tau,\tau)=1-\tau f(\tau)/\sqrt{\tau-1}$ for $\tau>1$ is checked numerically only.

# Verification

Run
```text
.venv/bin/python scripts/derive_alp_heavy_quark_triangle.py
.venv/bin/python scripts/derive_alp_heavy_quark_triangle.py --test
```
Thirty-eight unit tests pass in about 140 s; nine of them cover this function (equivalence with $F_h$, symmetry, $f$ in the three ranges, the reductions, the closed forms, the limits, and the imaginary part).
