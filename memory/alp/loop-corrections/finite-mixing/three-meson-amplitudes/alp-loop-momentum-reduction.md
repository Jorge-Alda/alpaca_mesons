---
title: Momentum reduction for the three-meson one-loop estimate
status: derived
created: 2026-11-22
updated: 2026-11-22
bibtex_keys: [Alda:2026nrd, Manohar:1983md]
provenance: analytical-derivation
---

# Momentum reduction for the three-meson one-loop estimate

The [three-meson estimate](alp-three-meson-one-loop-estimate.md) contains kinetic vertices with two derivatives. This record reduces their loop-momentum dependence to scalar moments without integrating them. The vertex normalization and the [loop-size prescription](../power-counting-and-validity/loop-counting/alp-loop-naturalness-assumption.md) are inputs. Momenta and masses have units GeV and squared momenta have units GeV$^2$. The derivation uses dimensional regularization in Euclidean signature. Continue the result to the physical kinematics with the same prescription for every denominator.

## Parameter shift

For a bubble use the first two denominators below. For a triangle use all three. Let $q_1,q_2$ be the sums of external momenta on the internal lines and $m_{n_1},m_{n_2},m_{n_3}$ their pole masses. Feynman parameters $z_1,z_2$ lie in the simplex $z_1\geq0$, $z_2\geq0$, $z_1+z_2\leq1$. The bubble is obtained by setting $z_2=0$. Completing the square gives

$$
\begin{aligned}
&(1-z_1-z_2)(k^2+m_{n_1}^2)
+z_1\big((k+q_1)^2+m_{n_2}^2\big)
+z_2\big((k+q_2)^2+m_{n_3}^2\big)\\
&=(k+z_1q_1+z_2q_2)^2
+(1-z_1-z_2)m_{n_1}^2+z_1m_{n_2}^2+z_2m_{n_3}^2\\
&\quad+z_1(1-z_1)q_1^2+z_2(1-z_2)q_2^2-2z_1 z_2\,q_1\!\cdot q_2 .
\end{aligned}
$$

The denominator powers are two for a bubble and three for a triangle, with the usual factors $1$ and $2$ multiplying their parameter integrals. The square completion is an identity before integration, checked in the [SymPy derivation](../../../../../scripts/derive_alp_loop_tensor_estimate.py).

## Angular average

Shift the loop momentum to the variable appearing in the square. For any scalar function of its square, rotational invariance in $d$ dimensions gives

$$
\int d^d k\,k_\mu k_\nu F(k^2)
=\frac{g_{\mu\nu}}{d}\int d^d k\,k^2 F(k^2),
$$
$$
\int d^d k\,k_\mu k_\nu k_\rho k_\sigma F(k^2)
=\frac{g_{\mu\nu}g_{\rho\sigma}+g_{\mu\rho}g_{\nu\sigma}+g_{\mu\sigma}g_{\nu\rho}}
{d(d+2)}
\int d^d k\,(k^2)^2 F(k^2).
$$

Odd moments vanish. Each kinetic quartic vertex has a numerator at most quadratic in the loop momentum, and the quadratic part is proportional to its square because the two derivatives are Lorentz contracted. Consequently the product of two quartic numerators requires scalar moments through $(k^2)^2$. A triangle with two potential cubic vertices and one quartic vertex requires moments through $k^2$. For their product, SymPy verifies the average of two arbitrary polynomials with scalar, linear and $k^2$ parts. The scalar product of the two linear coefficients appears with $k^2/d$. The scalar terms in the shifted numerator retain the external momenta and the Feynman parameters.

## Scalar moments and derivatives

Write $\Delta$ for the part of the combined denominator independent of the shifted loop momentum. For nonnegative integer $t$ and positive integer $r$, numerator cancellation gives

$$
\frac{(k^2)^t}{(k^2+\Delta)^r}
=\sum_{j=0}^{t}\binom{t}{j}(-\Delta)^{t-j}(k^2+\Delta)^{j-r}.
$$

Terms with $j\geq r$ are scaleless polynomial integrals in dimensional regularization. Combine the coefficients of each remaining scalar moment before estimating its magnitude. For the external momentum derivative needed in the residue, differentiate the shifted numerator and use

$$
\partial_s(k^2+\Delta)^{-r}
=-r(\partial_s\Delta)(k^2+\Delta)^{-r-1}.
$$

Apply this derivative before assigning any loop size. The same rule with a mass derivative handles double propagators from the mixing expansion. Both operations lower the mass dimension of the integral by two.

At a positive reference squared scale $m^2$ for the loop, the size assigned to a scalar moment with $r$ remaining denominators is

$$
\frac{(m^2)^{2-r}}{16\pi^2}.
$$

For this estimate, factor that size out of the parameter integral and retain the polynomial parameter dependence in the coefficients. For a triangle the normalized simplex measure gives

$$
2\int_0^1 dz_1\int_0^{1-z_1}dz_2\,z_1^i z_2^j (1-z_1-z_2)^k
=\frac{2i!j!k!}{(i+j+k+2)!}.
$$

These are polynomial weights, with the scalar loop functions unevaluated. They are a prescription for estimating the magnitude when the loop scales are comparable and there is no enhancement from a threshold, an infrared singularity or a large logarithm. Relative signs and phases of different moments cannot be obtained from this prescription.

## Verification

The executable script expands the three-denominator square completion and the product of two generic vertex numerators. It verifies the angular averages, moment reduction, differentiation, and simplex weights with symbolic inputs.

```sh
.venv/bin/python scripts/derive_alp_loop_tensor_estimate.py
PYTHONPATH=scripts .venv/bin/python -m unittest derive_alp_loop_tensor_estimate
```
