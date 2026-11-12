---
title: Closed form for a LO meson pole far from the other unmixed meson masses, for any $f_\pi/f_a$
status: mixed
created: 2026-11-12
updated: 2026-11-12
bibtex_keys: [Alda:2026nrd]
provenance: assumption
---

# Closed form for a LO meson pole far from the other unmixed meson masses, for any $f_\pi/f_a$

## Assumptions, inputs and notation

- Exact meson pole equation, Step 5 of the [exact pole-equation record](alp-lo-exact-pole-equations.md), in the mass basis of the [arrow-matrix record](alp-lo-arrow-matrix-spectrum.md): unmixed masses $m_{n,0}^2$, ALP entry $M^2_{\hat a\hat a,\rm can}$, $\epsilon=f_\pi/f_a$. Conventions of Alda:2026nrd [@Alda:2026nrd]. Meson $n$ non-degenerate with $\operatorname{tr}(\lambda_n\mathcal C)\neq0$.
- Notation (approved): the ALP mass seen by the meson $n$, which includes the effect of the other mesons and depends on $s$,
$$
\widetilde M_n(s)=M^2_{\hat a\hat a,\rm can}+\epsilon^2\sum_{n'\neq n}\frac{\operatorname{tr}(\lambda_{n'}\mathcal C)^2}{s-m_{n',0}^2},\qquad \widetilde M_n^{(0)}=\widetilde M_n(m_{n,0}^2).
$$
  Mesons with $\operatorname{tr}(\lambda_{n'}\mathcal C)=0$ do not enter.
- **Assumption (approximation, stated explicitly).** The meson $n$ is far from the other unmixed meson masses: $|s-m_{n,0}^2|\ll|m_{n,0}^2-m_{n',0}^2|$ for the pole $s$ of interest and $\epsilon^2\operatorname{tr}(\lambda_{n'}\mathcal C)^2\ll(m_{n,0}^2-m_{n',0}^2)^2$ for all $n'\neq n$. No assumption on the size of $\epsilon$ and none on the distance between $m_{n,0}^2$ and $M^2_{\hat a\hat a,\rm can}$.

## Step 1. Exact equation

Step 5 of the exact record reads, with no approximation,
$$
(s-m_{n,0}^2)\,\big[s-\widetilde M_n(s)\big]=\epsilon^2\operatorname{tr}(\lambda_n\mathcal C)^2 .
$$

## Step 2. Frozen closed form

Replace $\widetilde M_n(s)$ by its value $\widetilde M_n^{(0)}$ at $s=m_{n,0}^2$ (approximation under the assumption above). The equation becomes a quadratic with roots
$$
\boxed{\,s_\pm=\tfrac12\big(m_{n,0}^2+\widetilde M_n^{(0)}\big)\pm\sqrt{\tfrac14\big(\widetilde M_n^{(0)}-m_{n,0}^2\big)^2+\epsilon^2\operatorname{tr}(\lambda_n\mathcal C)^2}\,}
$$
the two levels formed by the meson $n$ and the ALP. The meson-like root is the one on the side of $m_{n,0}^2$ away from $\widetilde M_n^{(0)}$ ($s_-$ if $\widetilde M_n^{(0)}>m_{n,0}^2$, $s_+$ otherwise); at $\widetilde M_n^{(0)}=m_{n,0}^2$ both roots are equal mixtures, split by $2\epsilon|\operatorname{tr}(\lambda_n\mathcal C)|$. If the ALP entry is also far from the other mesons, the other root approximates the ALP-like pole: this is the two-level ALP-meson form. Dropping the sum in $\widetilde M_n$ gives the two-level formula of the [near-degeneracy record](alp-lo-near-degeneracy-two-level.md); the sum is the correction from the other mesons.

## Step 3. Error bound and validity

For fixed $b$ and $\epsilon\operatorname{tr}(\lambda_n\mathcal C)\neq0$, the roots of $(s-m_{n,0}^2)(s-b)=\epsilon^2\operatorname{tr}(\lambda_n\mathcal C)^2$ lie outside the interval between $m_{n,0}^2$ and $b$, so $s-m_{n,0}^2$ and $s-b$ have the same sign, and along each root
$$
\frac{ds}{db}=\frac{s-m_{n,0}^2}{(s-m_{n,0}^2)+(s-b)}\in(0,1).
$$
An exact pole $s_*$ is a root with $b=\widetilde M_n(s_*)$; the frozen roots have $b=\widetilde M_n^{(0)}$. Hence every exact pole lies within $|\widetilde M_n(s_*)-\widetilde M_n^{(0)}|$ of $s_+$ or $s_-$ (exact statement), with
$$
\widetilde M_n(s_*)-\widetilde M_n^{(0)}=-\epsilon^2\,(s_*-m_{n,0}^2)\sum_{n'\neq n}\frac{\operatorname{tr}(\lambda_{n'}\mathcal C)^2}{(s_*-m_{n',0}^2)(m_{n,0}^2-m_{n',0}^2)} .
$$
Relative to the distance of the pole from the unmixed mass, the error is therefore at most $\epsilon^2\sum_{n'\neq n}\operatorname{tr}(\lambda_{n'}\mathcal C)^2/\big(|s_*-m_{n',0}^2|\,|m_{n,0}^2-m_{n',0}^2|\big)$, which under the assumption is about $\epsilon^2\sum_{n'\neq n}\operatorname{tr}(\lambda_{n'}\mathcal C)^2/(m_{n,0}^2-m_{n',0}^2)^2$. At fixed couplings the error falls as the inverse square of the spacing to the other unmixed masses. At fixed spacing it grows with $\epsilon$, through $\epsilon^2\operatorname{tr}(\lambda_{n'}\mathcal C)^2$ and through the distance $|s_*-m_{n,0}^2|$. What controls the approximation is $\epsilon\operatorname{tr}(\lambda_{n'}\mathcal C)/(m_{n,0}^2-m_{n',0}^2)$ and the ratio of $|s_*-m_{n,0}^2|$ to the spacing, not $\epsilon$ alone.

## Verification

[`scripts/derive_alp_lo_arrow_matrix_poles.py`](../../../scripts/derive_alp_lo_arrow_matrix_poles.py), exact poles from the 50-digit eigenvalues of the arrow matrix.

Arrow model (test values, GeV$^2$): $m_{n,0}^2=1$, $M^2_{\hat a\hat a,\rm can}=1.1$, other mesons at $1+\Delta\,(1,\tfrac32,2,3)$, all $\operatorname{tr}(\lambda_{n'}\mathcal C)=1$. The two lowest poles $\mu_1$, $\mu_2$ are the meson-ALP pair. Some poles are negative for $\epsilon\ge1$; the model only tests the algebra.

| $\epsilon$ | $\Delta$ | $\lvert s_--\mu_1\rvert$ | $\lvert s_+-\mu_2\rvert$ |
|---|---|---|---|
| 0.1 | 2 | $8.6\times10^{-5}$ | $5.1\times10^{-4}$ |
| 0.1 | 64 | $7.6\times10^{-8}$ | $5.2\times10^{-7}$ |
| 0.5 | 2 | $3.2\times10^{-2}$ | $2.0\times10^{-2}$ |
| 0.5 | 64 | $2.3\times10^{-5}$ | $3.3\times10^{-5}$ |
| 1 | 2 | $0.30$ | $6.7\times10^{-2}$ |
| 1 | 64 | $2.0\times10^{-4}$ | $2.4\times10^{-4}$ |
| 2 | 2 | $2.2$ | $0.13$ |
| 2 | 16 | $3.2\times10^{-2}$ | $2.3\times10^{-2}$ |
| 2 | 64 | $1.8\times10^{-3}$ | $1.8\times10^{-3}$ |

Fitted exponent of the error in $\Delta$ ($\Delta=16,32,64$): $-2.00$, $-1.97$, $-1.95$, $-2.09$ for $\epsilon=0.1,0.5,1,2$. The bound of Step 3 holds at all 24 points.

LO blocks (neutral test sources of the exact record, $m_{a,0}^2=3/5$ GeV$^2$): the bound holds for every coupled meson at $\epsilon=1/2,1,2$. Error of the meson-like root at $\epsilon=1/2$: $10^{-13}$ GeV$^2$ for the lightest meson ($m_{n,0}^2=0.020$ GeV$^2$), $2.4\times10^{-4}$ GeV$^2$ for the coupled $K^0$ combination, $1.4\times10^{-2}$ GeV$^2$ for the heaviest ($1.03$ GeV$^2$, close to $M^2_{\hat a\hat a,\rm can}=1.08$ GeV$^2$). For the meson at $0.257$ GeV$^2$, only $0.0095$ GeV$^2$ from $m_{K^0,0}^2$, the assumption fails: the bound is larger than $10$ GeV$^2$ and the error is $3.6\times10^{-2}$ GeV$^2$.

```bash
.venv/bin/python scripts/derive_alp_lo_arrow_matrix_poles.py
.venv/bin/python -m unittest discover -s scripts -p test_derive_alp_lo_arrow_matrix_poles.py
```

All checks pass.
