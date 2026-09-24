---
title: U(3) vector-exchange flavor classifier for ALP three-pseudoscalar decays
status: derived
created: 2026-09-24
updated: 2026-09-24
bibtex_keys: [Ecker:1988te, Ecker:1989yg]
provenance: analytical-derivation
---

# Vector-Exchange Flavor Classifier for $a\to3P$

## Scope and inputs

The executable derivation is [`scripts/classify_a3p_vector_exchange.py`](../scripts/classify_a3p_vector_exchange.py). It uses the $U(3)$ component and trace conventions recorded in [the RChT basics](rcht-basics.md) and the two vertex structures derived in [the RChT flavor expansion](rcht-vpp-vpa-photon-pa-expansion.md):

$$
\operatorname{Tr}\!\left(T_V[T_{\rm spectator},C_A]\right),\qquad
\operatorname{Tr}\!\left(T_{\bar V}[T_{{\rm pair},1},T_{{\rm pair},2}]\right).
$$

The calculation uses natural units, $U=\exp(i\sqrt2\Phi/f_\pi)$, and explicitly sets

$$
f_\pi=\frac{93}{1000}\ \mathrm{GeV}.
$$

The pseudoscalar generators are those of $\pi^\pm,\pi^0,K^\pm,K^0,\bar K^0$ and

$$
T_\eta=\cos\theta_P T_{\eta_8}-\sin\theta_P T_{\eta_0},\qquad
T_{\eta'}=\sin\theta_P T_{\eta_8}+\cos\theta_P T_{\eta_0}.
$$

The vector generators are $\rho^\pm,\rho^0,K^{*\pm},K^{*0},\bar K^{*0}$ and the ideal-mixing approximation $T_\omega=\operatorname{diag}(1,1,0)/\sqrt2$, $T_\phi=\operatorname{diag}(0,0,1)$. Ideal vector mixing is an explicit approximation; no value is chosen for $\theta_P$.

## Method

1. Construct every generator as an explicit SymPy $3\times3$ matrix.
2. Enumerate every size-three pseudoscalar multiset with electric charge zero and total strangeness zero or $+1$.
3. For each distinct spectator choice and every vector, calculate both traces and their product. The second trace uses $T_{\bar V}=T_V^T$, including charged internal vectors.
4. Repeat for $C_A=\operatorname{diag}(c_u,c_d,c_s)$ and for $C_A=c_{sd}E_{sd}$. In the matrix convention, $E_{sd}$ has strangeness $-1$ and therefore selects final-state strangeness $+1$.
5. Factor the exact symbolic products so zeros and conditional factors involving $c_u,c_d,c_s$ and $\theta_P$ remain visible.

Pair ordering follows the fixed pseudoscalar list in the script. Reversing the pair changes the sign of the second trace, but not whether a vector contribution is present.

## Outputs and checks

The script prints all enumerated channels, using `--` when every product vanishes, and displays each nonzero spectator/vector/pair product in LaTeX. Repeated identical spectator choices are compressed with an explicit multiplicity. It asserts:

- zero total electric charge and allowed total strangeness for every enumerated channel;
- charge and strangeness flow independently at both trace vertices;
- equality of every classified product to the direct product of the two SymPy matrix traces;
- use of the transposed generator for every conjugate internal vector;
- the global $S=0$ selection for diagonal $C_A$ and $S=+1$ selection for $C_A=E_{sd}$.

Run with

```text
.venv/bin/python scripts/classify_a3p_vector_exchange.py
```

The $S=-1$ table follows by Hermitian conjugation together with $E_{sd}\to E_{ds}$.
