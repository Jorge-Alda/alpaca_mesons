---
title: Tadpole integral and subtraction used in the one-loop two-point functions
status: derived
created: 2026-11-06
updated: 2026-11-06
bibtex_keys: [Gasser:1983yg, Gasser:1984gg, Kaiser:1998ds, Alda:2026nrd]
provenance: analytical-derivation
---

# Tadpole integral and subtraction

## Assumptions and inputs

- Only tadpole graphs exist in $\Sigma_{ij}(s)$ and $\Sigma_{i\hat a}(s)$, see [vertices-and-topologies.md](vertices-and-topologies.md).
- The subtraction, $\lambda(\mu)$ and $L_i=L_i^r(\mu)+\Gamma_i\lambda(\mu)$ are those of [subtraction-scheme-and-scales.md](../counting/subtraction-scheme-and-scales.md) [@Gasser:1983yg; @Kaiser:1998ds]. No new symbol is introduced: the finite part is written out.
- Units: masses and $\mu$ in GeV, $s$ in GeV$^2$.

## Step 1: coincident-point propagator

For a field $\phi_k$ with LO mass $m_k$ (eigenstate of the LO mass matrix, including the $\pi^0$–$\eta$–$\eta'$ mixing), in dimensional regularization [@Gasser:1983yg]

$$
\langle\phi_k\phi_k\rangle=2m_k^2\,\lambda(\mu)+\frac{m_k^2}{16\pi^2}\log\frac{m_k^2}{\mu^2}.
$$

The first term is the pole, removed by the counterterms, and the second is the finite part after the subtraction with the $+1$ in $\lambda(\mu)$.

## Step 2: derivative contractions

In dimensional regularization $\int d^dk\,k^2/(k^2-m^2)=m^2\int d^dk\,1/(k^2-m^2)$, and an integrand odd in $k$ vanishes. Hence

$$
\langle\partial_\mu\phi_k\,\partial^\mu\phi_k\rangle=m_k^2\langle\phi_k\phi_k\rangle,\qquad
\langle\phi_k\,\partial_\mu\phi_k\rangle=0 .
$$

For the nonet, $\langle\phi_c\phi_d\rangle=\sum_kO_{ck}O_{dk}\langle\phi_k\phi_k\rangle$ with $O$ the orthogonal matrix that diagonalizes the LO mass matrix of [neutral-meson-matrices.md](../neutral-meson-matrices.md); the charged and strange states are already mass eigenstates.

## Step 3: consequence

Since the loop integrals do not depend on the external momentum, $\Sigma_{ij}(s)$ is linear in $s$ and real. The $\mu$ dependence of the finite part, $\mu\,\partial_\mu\big[\tfrac{m_k^2}{16\pi^2}\log(m_k^2/\mu^2)\big]=-\tfrac{m_k^2}{8\pi^2}$, cancels that of $2m_k^2\lambda(\mu)$ ($\mu\,d\lambda/d\mu=-1/(16\pi^2)$ at $d=4$ in the sense of the running of $L_i^r$); the full cancellation against the local terms is in the next record.

## Verification

`--meson` of [`scripts/derive_nlo_loop_self_energies.py`](../../../scripts/derive_nlo_loop_self_energies.py) asserts the contraction against an explicit Wick expansion of the vertex.
