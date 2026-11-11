---
title: Axial-vector nonet field, flavor and phase conventions
status: derived
created: 2026-11-10
updated: 2026-11-10
bibtex_keys: [Ecker:1988te, Alda:2026nrd, ParticleDataGroup:2024cfk]
provenance: analytical-derivation
---

# Axial-vector nonet field, flavor and phase conventions

## Inputs and assumptions

- The field is an antisymmetric Lorentz tensor $A_{\mu\nu}$ transforming as $A_{\mu\nu}\to hA_{\mu\nu}h^\dagger$, with $h$ the CCWZ compensator. This is the same choice as the vector nonet in [RChT basics](../foundations/rcht-basics.md) [@Ecker:1988te]. The Proca alternative is not used here; its equivalence map is a later task.
- Generators: $\langle\lambda_a\lambda_b\rangle=2\delta_{ab}$ with $T_a=\lambda_a/\sqrt2$ and $T_0=\mathbf 1/\sqrt3$, so $\langle T_iT_j\rangle=\delta_{ij}$ [@Ecker:1988te].
- Metric and Levi-Civita: $(+,-,-,-)$ with $\epsilon^{0123}=+1$, as fixed in the [Levi-Civita convention](../odd-parity/levi-civita-convention.md). No odd-parity content is used in this record.
- Strange phase: the ideal $s\bar s$ state enters with a positive sign, as in the [vector nonideal mass record](../foundations/vector-nonideal-mass-action.md).
- Assumption: the physical state labelled $f_1(1285)$ is the one that is mostly nonstrange at small mixing angle. This is a labelling choice. The angle is not fitted here.

## Nonet matrix

Writing $A_{\mu\nu}=\sum_i T_i\,A^i_{\mu\nu}$ with ideal-basis fields $a_1^0$, $a_1^\pm$, $f_1^N=(u\bar u+d\bar d)/\sqrt2$, $f_1^S=s\bar s$ and the $K_1$ pair $K_1^{+},K_1^{0},K_1^{-},\bar K_1^{0}$, the matrix has the same form as the vector matrix in RChT basics with $\rho\to a_1$, $\omega\to f_1^N$, $\phi\to f_1^S$, $K^*\to K_1$:

$$
A_{\mu\nu}=
\begin{pmatrix}
\dfrac{a_1^0+f_1^N}{\sqrt2}&a_1^+&K_1^+\\
a_1^-&\dfrac{-a_1^0+f_1^N}{\sqrt2}&K_1^0\\
K_1^-&\bar K_1^0&f_1^S
\end{pmatrix}_{\mu\nu}.
$$

Charged and octet components relate to the fields as $A^1\mp iA^2=\sqrt2\,a_1^\pm$, $A^8=(f_1^N-\sqrt2 f_1^S)/\sqrt3$ and $A^0=(\sqrt2 f_1^N+f_1^S)/\sqrt3$. The same relations hold for the $K_1$ pairs.

## Nonideal mixing of the neutral isoscalars

Define the physical states by $\mathbf A_{\rm phys}=(f_1(1285),f_1(1420))^T$ and $\mathbf A_I=(f_1^N,f_1^S)^T$, related by
$\mathbf A_I=R\,\mathbf A_{\rm phys}$ with $R=\begin{pmatrix}\cos\theta&\sin\theta\\-\sin\theta&\cos\theta\end{pmatrix}$. The physical generators are

$$
T_{f_1(1285)}=\operatorname{diag}\!\left(\tfrac{\cos\theta}{\sqrt2},\tfrac{\cos\theta}{\sqrt2},-\sin\theta\right),\qquad
T_{f_1(1420)}=\operatorname{diag}\!\left(\tfrac{\sin\theta}{\sqrt2},\tfrac{\sin\theta}{\sqrt2},\cos\theta\right).
$$

The sign of $\theta$ follows the vector-record convention. The angle is not fitted here.

## Executable check

Run:

```text
.venv/bin/python scripts/derive_axial_nonet_conventions.py
.venv/bin/python -m unittest scripts.implementation.test_axial_nonet_conventions
```

The [script](../../../scripts/derive_axial_nonet_conventions.py) builds the generators, the nonet matrix and the rotated generators symbolically with SymPy, and prints them as LaTeX. The [tests](../../../scripts/implementation/test_axial_nonet_conventions.py) check:

- $\langle T_iT_j\rangle=\delta_{ij}$ for all nine generators;
- the $U(3)$ completeness relation $\sum_i (T_i)_{ab}(T_i)_{cd}=\delta_{ad}\delta_{bc}$;
- orthonormality of the ideal-basis matrices and the explicit entries above;
- orthogonality of $R$ and preservation of normalization under rotation;
- the dictionary to Alda:2026nrd (below).

Last run: 9 tests, all passing; the script ends with `Assertions passed`.

## Dictionary to Alda:2026nrd

Alda:2026nrd uses $\Phi_{\rm Alda}=\sum_P\lambda_P P$, $u=\exp(i\Phi_{\rm Alda}/(2f_\pi))$, $\langle\lambda_i\lambda_j\rangle=2\delta_{ij}$ [@Alda:2026nrd, see the [neutral generator extract](../../extracts/Alda:2026nrd/neutral-meson-generators.md)]. The RChT convention uses $\Phi_R=\sum_P T_PP$ and $u=\exp(i\Phi_R/(\sqrt2F))$. With $T_P=\lambda_P/\sqrt2$ and $F=f_\pi$ the two exponents agree, since $\Phi_{\rm Alda}=\sqrt2\,\Phi_R$. The same map is used for the axial nonet. Alda:2026nrd does not define an axial nonet, so this is a dictionary for the pseudoscalar fields only.

## Evidence status

- Established by the derivation and the tests: normalization, ideal-basis orthonormality, the explicit matrix, the rotation, and the pseudoscalar dictionary.
- Assumption: the $f_1(1285)$ labelling and the sign of $\theta$.
- Not checked: the axial nonet content against literature values. Kampf:2011ty and Cirigliano:2004ue are recorded in the [extract index](../../extracts/INDEX.md) but are not used for any number here.
