---
title: Matrix Taylor coefficients of radial and phase dependent quark masses
status: derived
created: 2026-11-18
updated: 2026-11-18
bibtex_keys: [DiLuzio:2026poh]
provenance: analytical-derivation
---

# Matrix Taylor coefficients of radial and phase dependent quark masses

## Inputs and conventions

The mass spurion and polar field parametrization are those of DiLuzio:2026poh, Eqs. (7) and (9) [@DiLuzio:2026poh]. Rescale the dimensionless complex field by $\Phi=f_\varphi\varphi$ and write the dimensionful radial and phase coordinates as
$$
\Phi=\frac{f_\varphi+\sigma}{\sqrt2}e^{ia/f_\varphi},
$$
so $a$ and $\sigma$ have mass dimension one. Define diagonal matrices
$$
m_q=\operatorname{diag}(m_u,m_d,m_s),\qquad N=\operatorname{diag}(n_u,n_d,n_s).
$$
In the paper's spurion convention,
$$
M_q(a,\sigma)=m_q\left(\mathbf1+\frac{\sigma}{f_\varphi}\right)^N e^{iNa/f_\varphi}=S(a,\sigma)+iP(a,\sigma),
$$
where matrix powers are defined entrywise. Therefore
$$
S=m_q\left(\mathbf1+\frac{\sigma}{f_\varphi}\right)^N\cos\!\left(\frac{Na}{f_\varphi}\right),\qquad
P=m_q\left(\mathbf1+\frac{\sigma}{f_\varphi}\right)^N\sin\!\left(\frac{Na}{f_\varphi}\right).
$$

## Taylor coefficients

Define ordinary Taylor coefficients by
$$
S(a,\sigma)=\sum_{r,s\ge0}S_{rs}\sigma^r a^s,\qquad P(a,\sigma)=\sum_{r,s\ge0}P_{rs}\sigma^r a^s.
$$
Then
$$
\boxed{\begin{aligned}
S_{rs}&=\frac{m_q}{f_\varphi^{r+s}}\binom Nr N^s\frac{\cos(\pi s/2)}{s!},\\
P_{rs}&=\frac{m_q}{f_\varphi^{r+s}}\binom Nr N^s\frac{\sin(\pi s/2)}{s!},
\end{aligned}}
$$
for nonnegative integers $r,s$, where $\binom Nr=N(N-\mathbf1)\cdots(N-(r-1)\mathbf1)/r!$ and $\binom N0=\mathbf1$. These are matrix expressions; no division by the mass matrix occurs. If coefficients are defined as derivatives at the origin rather than monomial coefficients, multiply $S_{rs},P_{rs}$ by $r!s!$.

For illustration,
$$
\begin{aligned}
S={}&m_q+\frac{m_qN}{f_\varphi}\sigma+\frac{m_q\binom N2}{f_\varphi^2}\sigma^2-\frac{m_qN^2}{2f_\varphi^2}a^2+\cdots,\\
P={}&\frac{m_qN}{f_\varphi}a+\frac{m_qN^2}{f_\varphi^2}\sigma a+\frac{m_q\binom N2N}{f_\varphi^3}\sigma^2a-\frac{m_qN^3}{6f_\varphi^3}a^3+\cdots.
\end{aligned}
$$

## Derivation and scope

For each diagonal entry, expand $(1+\sigma/f_\varphi)^{n_q}$ by the binomial theorem and $e^{in_qa/f_\varphi}$ by its exponential series. Real and imaginary parts select respectively $\cos(\pi s/2)/s!$ and $\sin(\pi s/2)/s!$ in the coefficient of $a^s$. Since $m_q$ and $N$ are diagonal, the matrix factors commute. Generalized binomial coefficients give the Taylor series for negative integer charges as well, locally around the nonzero vacuum.

These coefficients describe the paper's quark-mass spurion before combining it with the chiral field $U$. In the alternate partonic convention $X$ defined by $-\bar q_LX^\dagger q_R+\mathrm{h.c.}$, $X=M_q^\dagger=S-iP$, so the pseudoscalar part has the opposite sign.
