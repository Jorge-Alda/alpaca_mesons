---
title: Charge-conserving U(3) vector redefinition of the two-meson vertex
status: derived
created: 2026-11-11
updated: 2026-11-11
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Charge-conserving $U(3)$ vector redefinition of the two-meson vertex

## Inputs and assumptions

Use the conventions, outgoing-state matrices, source terms and leading-order masses in the [flavour-changing two-pseudoscalar derivation](../two-pseudoscalar/alp-flavour-violating-two-pseudoscalar-amplitudes.md), and the normalized mixed poles in the [finite-mixing record](finite-alp-meson-mixing.md), following `Alda:2026nrd` [@Alda:2026nrd]. Work in natural units, at tree level and leading chiral order, about a stationary parity-even vacuum. The ALP is electrically neutral, the outgoing charged mesons do not mix with it, and the pure-meson leading-order even-parity Lagrangian has no cubic terms. The meson masses and quark masses satisfy the same leading-order mass relation. The calculation is exact in the incoming mixed pole vector; it uses only the cubic vertex, which is the first derivative in the ALP coordinate at the vacuum. No small-mixing expansion is made.

## Derivation

1. In the quark basis $(u,d,s)$, choose the Hermitian vector generator $(|d\rangle\langle s|+|s\rangle\langle d|)/2$. It commutes with the electromagnetic charge matrix. Applying $i[T,\Phi]$ to the charged $\pi$ and $K$ fields mixes the two positive-charge fields and, separately, the two negative-charge fields. Thus the following formal change of variables preserves electric charge.

2. Substitute the inverse field map through first order in the ALP coordinate into the free charged-meson Lagrangian. The momentum-space contact coefficient for the outgoing $K^+\pi^-$ channel is, apart from the vector-rotation parameter,
$$
\frac{i}{2}\left[(p_K^2-m_K^2)-(p_\pi^2-m_\pi^2)\right].
$$
It is generally nonzero off shell and vanishes when both charged legs are on shell. The script derives it from explicit $3\times3$ matrices and the transformed kinetic and mass terms, rather than assigning the masses independently of the channel.

3. Contract the original cubic vertex with the incoming pole. Since the two outgoing mesons are charged, their pole vectors have no neutral ALP component. The pure-meson cubic vertex vanishes under the stated assumptions. The result is the ALP component of the normalized incoming pole times the known bare-ALP vertex:
$$
 i\mathcal M(a_{\rm phys}\to K^+\pi^-)
 =-\frac{B_0(m_s-m_d)(c_q^{VS})_{ds}}{f_a\sqrt{K_{\hat a\hat a}+2K_{i\hat a}w_i^{(a)}+K_{ij}w_i^{(a)}w_j^{(a)}}}.
$$
Here $w_{\hat a}^{(a)}=1$, and the expression is valid when that normalization exists. The factor under the square root is positive for a regular positive-norm pole. This is generally nonzero. The vector-rotation contact in step 2 vanishes on shell, so the same amplitude is obtained in the rotated coordinates. The exact physical pole vector and kinetic matrix are inputs, not solved here.

## Verification and limitations

The script checks the charge commutator, the charged field map, the off-shell contact and the on-shell comparison, and uses the existing channel extraction to check the invariant combination $c_q^{VS}=c_q^V-c_q^S$. The contact check and the physical amplitude comparison are both restricted to the three-point vertex. The three-meson amplitude and its exchange diagrams require a separate calculation.

```sh
.venv/bin/python scripts/derive_alp_charged_u3_vector_redefinition.py
.venv/bin/python -m unittest discover -s scripts -p test_derive_alp_charged_u3_vector_redefinition.py
```
