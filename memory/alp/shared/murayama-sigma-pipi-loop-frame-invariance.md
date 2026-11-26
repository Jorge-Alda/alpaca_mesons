---
title: Constant axial redefinition and the traces of the one-loop sigma to pi+ pi- estimate
status: derived
created: 2026-11-26
updated: 2026-11-27
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: analytical-derivation
---

# Frame invariance of the traces of the $\sigma\to\pi^+\pi^-$ estimate

Task 4.2 of the OpenSpec change `estimate-sigma-pipi-one-loop-murayama`. No loop is evaluated.

## Setting and assumptions

- The source $M(\sigma,a)=\sum_{r,s}M_{r,s}\,y^rx^s$, with $y=\sigma/f_a$ and $x=a/f_a$, follows [the radial frame record](murayama-radial-frame-invariance.md). The coefficients have mass dimension one.
- Redefinition: $g_R=g_L^\dagger=e^{ix\kappa}$, with a constant Hermitian $\kappa$ and $[\kappa,Q]=0$. $\sigma$ is not rotated. Starting from the connection-free frame, the new frame has $M'=g_RMg_L^\dagger$, $C'_R=\kappa$, $C'_L=-\kappa$ and $\Theta'=\Theta-x\langle\kappa\rangle$.
- Dressing: transport along $x$ at fixed $y$, $\partial_xW=iCW$ with $W(0)=1$. Then $\widehat M=W_R^\dagger MW_L$ and $\widehat\Theta=\Theta+\int_0^x\langle(C_R-C_L)/2\rangle$.
- Redefinitions that depend on $\sigma$ are not covered.

## Inputs

The kept traces $T_\Gamma$ are taken from the [1PI trace record](murayama-sigma-pipi-loop-traces-1pi.md) and the [external-leg record](murayama-sigma-pipi-loop-traces-external-legs.md): all 42 kept graphs (26 1PI graphs, 6 $\sigma$ and 5 + 5 $\pi^\pm$ self-energies), with the D4 terms removed exactly, the per-term loop scales (D3) and the termwise routing maximum (D12). The masses, $N(m_a^2)$ and the conditions of the Piecewise coefficients (inequalities between squared masses) do not contain the sources, so only the products of traces are transformed. The sources that appear in them are $\widehat S_{0,0}$, $\widehat S_{0,1}$, $\widehat S_{0,2}$, $\widehat S_{1,0}$, $\widehat S_{1,1}$, $\widehat S_{1,2}$, $\widehat S_{2,0}$, $\widehat S_{2,1}$, $\widehat S_{3,0}$, $\widehat P_{0,1}$, $\widehat P_{1,1}$ and $\widehat P_{2,0}$, read from the enumeration. All have $r\le3$ and $s\le3$, the range covered by check 1; the script checks this.

## Checks and result

1. **Dressed sources.** For noncommuting symbols, all 32 dressed coefficients $\widehat S_{r,s}$, $\widehat P_{r,s}$ with $r\le3$ and $s\le3$ are identical in the two frames. This includes $r=3$, which the [radial frame record](murayama-radial-frame-invariance.md) did not cover.
2. **Traces.** $\mathbb M^{-1}_{m_a^2}(\mathcal C)$, $N(m_a^2)$ and $\lambda_{\pi^\pm}$ are built from $\widehat S_{0,0}$, $\widehat P_{0,1}$ and $\widehat\Theta_1$ only. Every kept $T_\Gamma$, written as a polynomial in the dressed symbols with these matrices opaque, is therefore unchanged. The script checks this for all 42 kept graphs by inserting the dressed coefficients of the new frame. The check is on the traces only; no loop enters.
3. **Negative control.** Inserting the raw coefficients of the new frame, to first order in $\kappa$, changes exactly those kept $T_\Gamma$ that contain a source with $s\ge1$ (29 of the 42). The thirteen that contain only sources with $s=0$ (graphs 1, 4, 7, 10, 16, 22, σ1, and π1, π2, π3 for each charged leg) are not changed, because the $x^0$ coefficients are not affected by the redefinition. The $\delta Z$ of π1 and π2 contain only $\mathbb M^{-1}_{m_a^2}(\mathcal C)$ and $\lambda_{\pi^\pm}$ (held opaque), with the per-term scales of D3. The test checks that the set of changed graphs is exactly the set with a source of $s\ge1$.
4. **$\widehat\Theta_1$.** $\widehat\Theta_{0,1}$ is unchanged, while the raw $\Theta_{0,1}$ shifts by $-\langle\kappa\rangle$. The assumption $\widehat\Theta_1=0$ is therefore frame independent.

The estimate is invariant under a constant axial redefinition for the same physical legs, at the level of its traces.

## Reproduction

```sh
PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_loop_frame_invariance.py
PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_loop_frame_invariance
```

[Script](../../../scripts/derive_sigma_pipi_loop_frame_invariance.py) (an independent reimplementation of the transport of [`derive_murayama_radial_frame_invariance.py`](../../../scripts/derive_murayama_radial_frame_invariance.py), extended to $r\le3$), [tests](../../../scripts/test_derive_sigma_pipi_loop_frame_invariance.py).
