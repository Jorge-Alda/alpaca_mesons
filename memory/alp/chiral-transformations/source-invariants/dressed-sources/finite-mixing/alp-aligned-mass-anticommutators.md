---
title: Aligned mass-source derivatives as nested anticommutators
status: derived
created: 2026-11-20
updated: 2026-11-20
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Aligned mass-source derivatives

## Inputs and assumptions

Use natural units and $x=a/f_a$. The mass source has mass dimension one; the axial generator is dimensionless. The [finite-mixing three-pseudoscalar derivation](../../../../three-pseudoscalar/finite-mixing/alp-finite-mixing-three-pseudoscalar-amplitudes.md) defines $\widehat M(x)$ in derivative-source-free coordinates and the constant alignment generator $\widehat c_q^A=-\mathbb M_{m_a^2}^{-1}(\mathcal C)$ at the exact ALP pole. Its aligned mass source is $M'(x)=e^{ix\widehat c_q^A}\widehat M(x)e^{ix\widehat c_q^A}$. Assume an analytic source near $x=0$ and a regular isolated pole. The anticommutator map $\mathcal A_X[Y]=\{X,Y\}$ follows the [Appendix A extract](../../../../../extracts/Alda:2026nrd/appendix-a-expansion.md) of `Alda:2026nrd` [@Alda:2026nrd].

## Derivation

1. Differentiation of both exponentials gives
$$
\frac{dM'}{dx}=i\{\widehat c_q^A,M'\}+e^{ix\widehat c_q^A}\frac{d\widehat M}{dx}e^{ix\widehat c_q^A}.
$$
2. Apply the product rule repeatedly and evaluate at $x=0$. With $\widehat M_n=d^n\widehat M/dx^n|_0$ and $M'_n=d^nM'/dx^n|_0$,
$$
M'_n=\sum_{k=0}^n\binom nk i^k\mathcal A_{\widehat c_q^A}^{\,k}[\widehat M_{n-k}].
$$
For example,
$$
M'_2=\widehat M_2+2i\{\widehat c_q^A,\widehat M_1\}-\{\widehat c_q^A,\{\widehat c_q^A,\widehat M_0\}\}.
$$
This is exact in $f_\pi/f_a$. The independent higher derivatives $\widehat M_n$ cannot generally be determined from first-order source invariants alone.

3. Multiplying the formal Taylor series of $\widehat M(x)$ by the operator exponential gives
$$
\sum_{n=0}^{\infty}\frac{x^n}{n!}M'_n
=\exp\!\left(ix\mathcal A_{\widehat c_q^A}\right)[\widehat M(x)].
$$
The exponential denotes its formal power series as a linear map on matrices. Extracting its $x^n$ coefficient recovers the finite binomial sum above. At fixed $n$, the exponential cannot act on $\widehat M_n$ alone: it also acts on the lower derivatives of $\widehat M$. This identity is formal and needs no convergence assumption; analytic sources give its local Taylor interpretation.

4. Since the alignment generator is independent of $x$, the derivative and anticommutator operations commute. Repeated differentiation therefore also gives the compact expression
$$
M'_n=\left.\left(\frac{d}{dx}+i\mathcal A_{\widehat c_q^A}\right)^n[\widehat M(x)]\right|_{x=0}.
$$
This expression retains all derivatives of the source; it is an operator power, not an operator acting on $\widehat M_n$ alone. If $\widehat M$ is constant, it reduces to $i^n\mathcal A_{\widehat c_q^A}^{\,n}[\widehat M_0]$. This reduction is not valid for a general source.

5. Write $\widehat M=\widehat S+i\widehat P$ with Hermitian scalar and pseudoscalar matrices, as in the [second-order invariant derivation](../alp-dressed-mass-second-order-invariants.md). The constant Hermitian generator makes $\mathcal A_{\widehat c_q^A}$ preserve Hermiticity. Separating even and odd powers of the operator exponential yields
$$
S'=\cos(x\mathcal A_{\widehat c_q^A})[\widehat S]-\sin(x\mathcal A_{\widehat c_q^A})[\widehat P],\quad
P'=\sin(x\mathcal A_{\widehat c_q^A})[\widehat S]+\cos(x\mathcal A_{\widehat c_q^A})[\widehat P].
$$
Applying the product rule, for each nonnegative integer $n$,
$$
S'_n=\sum_{k=0}^n\binom nk\mathcal A_{\widehat c_q^A}^{k}[\cos(k\pi/2)\widehat S_{n-k}-\sin(k\pi/2)\widehat P_{n-k}],
$$
$$
P'_n=\sum_{k=0}^n\binom nk\mathcal A_{\widehat c_q^A}^{k}[\sin(k\pi/2)\widehat S_{n-k}+\cos(k\pi/2)\widehat P_{n-k}].
$$
The trigonometric functions of the operator are formal series; no commutation of the flavour matrices is assumed. The script derives these series to fourth order and the unit test verifies both formulas at those orders. The all-order result follows from the formal series and product rule.

## Reproduction and checks

[`derive_aligned_mass_anticommutators.py`](../../../../../../scripts/derive_aligned_mass_anticommutators.py) differentiates both matrix exponentials and an arbitrary noncommuting mass series through fourth order, then compares with the nested anticommutators. The general formula follows from the product rule above; the finite-order test does not by itself prove all orders.

```sh
.venv/bin/python scripts/derive_aligned_mass_anticommutators.py
.venv/bin/python -m unittest discover -s scripts -p test_derive_aligned_mass_anticommutators.py
```
