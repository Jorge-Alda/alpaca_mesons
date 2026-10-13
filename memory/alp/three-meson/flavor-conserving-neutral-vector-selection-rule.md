---
title: Neutral-vector selection rule in flavor-conserving ALP three-pseudoscalar exchange
status: derived
created: 2026-09-24
updated: 2026-09-24
bibtex_keys: [Ecker:1988te, Ecker:1989yg, Alda:2026nrd]
provenance: analytical-derivation
---

# Neutral-Vector Selection Rule for Flavor-Conserving $a\to3P$

## Scope and inputs

This record audits the absence of $\rho^0$, $\omega$, and $\phi$ exchange in the diagonal-$C_A$ part of [`vector_exchange_3P.md`](../../../vector_exchange_3P.md). The classified leading-order, even-parity RChT topology has the flavor factor

$$
\operatorname{Tr}\!\left(T_V[T_{\rm spectator},C_A]\right)
\operatorname{Tr}\!\left(T_{\bar V}[T_{{\rm pair},1},T_{{\rm pair},2}]\right),
$$

as derived in [the RChT flavor expansion](../../rcht/foundations/rcht-vpp-vpa-photon-pa-expansion.md) from the standard $G_V$ operator [@Ecker:1988te; @Ecker:1989yg] and the stated ALP-source expansion [@Alda:2026nrd]. The executable audit is [`scripts/audit_neutral_vector_classification.py`](../../../scripts/audit_neutral_vector_classification.py).

Assumptions are an ordinary finite-dimensional flavor trace, a flavor-conserving diagonal spurion

$$
C_A=\operatorname{diag}(c_u,c_d,c_s),
$$

and the commutator-form $VPa$ vertex used by the classifier. No ideal $\omega$--$\phi$ mixing, generator normalization, pseudoscalar mixing angle, mass degeneracy, or numerical coupling values are required.

## Derivation

For any spectator generator $T_P$, cyclicity of the trace gives

$$
\begin{aligned}
\operatorname{Tr}\!\left(T_V[T_P,C_A]\right)
&=\operatorname{Tr}(T_VT_PC_A-T_VC_AT_P)\\
&=\operatorname{Tr}(C_AT_VT_P-T_VC_AT_P)\\
&=\operatorname{Tr}\!\left([C_A,T_V]T_P\right).
\end{aligned}
$$

Consequently, every vector direction in the centralizer of $C_A$ obeys

$$
[C_A,T_V]=0
\quad\Longrightarrow\quad
\operatorname{Tr}\!\left(T_V[T_P,C_A]\right)=0
$$

for arbitrary $T_P$. The physical $\rho^0$, $\omega$, and $\phi$ generators are linear combinations of diagonal quark-flavor generators, even without ideal neutral-vector mixing. They therefore commute with diagonal $C_A$, so the first vertex vanishes before the $VPP$ vertex or propagator is considered.

Equivalently, in the diagonal quark basis $[T_P,C_A]$ has zero diagonal entries, while multiplication by any diagonal $T_V$ cannot generate a diagonal entry. Its trace is therefore zero. The audit constructs a fully generic $3\times3$ matrix $T_P$ in SymPy, obtains this zero directly, and then verifies a $3\times9$ zero table for $V=(\rho^0,\omega,\phi)$ and every pseudoscalar generator used by the classifier.

This is a centralizer or flavor-flow selection rule, not a prohibition on the second vertex. A neutral vector can have a nonzero $VPP$ coupling such as $\rho^0\to\pi^+\pi^-$, but the diagonal spurion cannot produce that vector through the classified $aPV$ commutator vertex. Flavor-off-diagonal vectors instead give differences such as

$$
\operatorname{Tr}\!\left(T_{\rho^-}[T_{\pi^+},C_A]\right)=c_d-c_u,
$$

and need not vanish. Electrically neutral but flavor-off-diagonal $K^{*0}$ and $\bar K^{*0}$ are likewise not covered by the rule.

## Scope limits and evidence status

- **Derived result:** the absence of $\rho^0$, $\omega$, and $\phi$ exchange in the flavor-conserving table is correct for the classified LO commutator topology.
- **Stronger check:** the zero is independent of ideal $\omega$--$\phi$ mixing; any mixing restricted to diagonal neutral-vector directions preserves it.
- **Counterexample:** an off-diagonal flavor-changing spurion need not commute with diagonal vectors, consistent with the nonzero neutral-vector entries in the $C_A=c_{sd}E_{sd}$ table.
- **Additional zero:** if $C_A\propto\mathbf1$, the complete classified vector-exchange topology vanishes because $[T_P,C_A]=0$ for every $T_P$.
- **Explicit limitation:** this is not a theorem for the full hadronic amplitude. Higher-order operators, loops, odd-parity interactions, extra symmetry-breaking spurions, or other topologies can evade the commutator argument. The direct Goldstone-current $WPa$ term recorded in the RChT expansion also vanishes for diagonal $C_A$ because it is proportional to $[Q,C_A]$.
- **Literature inputs:** the RChT operator basis and tensor conventions [@Ecker:1988te; @Ecker:1989yg] and the ALP-source expansion [@Alda:2026nrd].
- **Analytical evidence:** the trace identity and all concrete classifier checks are produced by the executable SymPy audit; no new literature source was consulted.

Run with

```text
.venv/bin/python scripts/audit_neutral_vector_classification.py
```

All assertions pass.
