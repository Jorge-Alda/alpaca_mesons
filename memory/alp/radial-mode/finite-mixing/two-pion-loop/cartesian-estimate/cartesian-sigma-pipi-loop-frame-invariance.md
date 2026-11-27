---
title: Cartesian frame - constant axial redefinition and the traces of the one-loop sigma to pi+ pi- estimate
status: derived
created: 2026-11-27
updated: 2026-11-27
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: analytical-derivation
---

# Cartesian frame: invariance of the traces under a constant axial redefinition

Task 4.2 of the OpenSpec change `estimate-sigma-pipi-one-loop-murayama-cartesian`. **No loop integral, scalar loop function, loop coefficient or loop amplitude is computed**; the check is made on the products of traces of the estimate.

**Caveat.** The estimate whose traces are checked here is an order of magnitude per graph class. Graphs with the same $\epsilon$ scaling can cancel in the exact loop, differently in the Cartesian and in the polar parametrization. Agreement with the polar estimate, including its up-only prefactor $-21/(64\pi^2)$, is not expected, and neither estimate is authoritative.

## Setting and assumptions

- Taylor variables $y=s/f_a$, $x=p/f_a$. The source is $\widehat M(s,p)=\sum_{r,s}(\widehat S^C_{r,s}+i\widehat P^C_{r,s})\,y^rx^s$ ([Cartesian frame and tree record](cartesian-sigma-pipi-loop-frame-and-tree.md)) and $\Theta(s,p)=\sum_{r,s}\Theta_{r,s}y^rx^s$ is the coefficient function of the anomalous singlet term. The starting frame has no derivative coupling of $s$ or $p$ to the mesons (assumption 4 of the frame record), i.e. vanishing connections $C_R=C_L=0$.
- **Redefinition (assumption on its form).** $g_R=g_L^\dagger=e^{ix\kappa}$, with a constant Hermitian $3\times3$ matrix $\kappa$ that commutes with $Q$; $s$ is not rotated. The new frame has the raw source $\widehat M'=g_R\widehat Mg_L^\dagger$, the connections $C'_R=\kappa$, $C'_L=-\kappa$ along $x$, and $\Theta'=\Theta-x\langle\kappa\rangle$.
- **Dressing.** Transport along $x$ at fixed $s$: $\partial_xW=iCW$ with $W(0)=\mathbf 1$, then $\widehat M_{\rm dressed}=W_R^\dagger\widehat MW_L$ and $\widehat\Theta=\Theta+\int_0^x\langle(C_R-C_L)/2\rangle$. In the starting frame the dressed coefficients equal the raw ones.
- **Not covered.** Redefinitions that depend on $s$ (for which the transport along $x$ at fixed $s$ is not enough) are not covered by this check.

## Inputs

- The termwise routing maxima $X_\Gamma^{\max}$ and $\delta Z^{\max}$ of the 37 kept graphs (26 1PI graphs, 6 $s$ self-energies, 5 $\pi^\pm$ self-energies; the $\pi^-$ self-energies equal the $\pi^+$ ones graph by graph), with D4 terms removed and the D3 scales: [Cartesian 1PI trace record](cartesian-sigma-pipi-loop-traces-1pi.md), [Cartesian external-leg record](cartesian-sigma-pipi-loop-traces-external-legs.md).
- Pole data: $A=B_0(\langle\widehat S^C_{1,0}\rangle-2\langle\widehat S^C_{0,2}\rangle)$, $\kappa_\sigma=B_0(\langle\widehat S^C_{1,0}\rangle-2\langle\widehat S^C_{2,0}\rangle)$, $\mathcal C=B_0\widehat P^C_{0,1}$ and $\mathbb M_s$ built from $\widehat S^C_{0,0}$: [Cartesian frame and tree record](cartesian-sigma-pipi-loop-frame-and-tree.md), Steps 2 and 3.

## Method

Every product is computed for noncommuting symbols (no matrices). The series in $x$ are cut at $x^3$ and $y^3$; this only limits the range of coefficients that are compared ($r,s\le3$) and does not enter the estimate.

1. $\widehat M'$, $W_R$, $W_L$ and the dressed coefficients are computed in both frames. $\widehat S^C$ is the Hermitian part and $\widehat P^C$ the anti-Hermitian part divided by $i$ of each Taylor coefficient.
2. Every source letter of each kept $T_\Gamma$ is replaced by the dressed coefficient of the new frame. The masses, $N(m_a^2)$, $\mathbb M^{-1}_{m_a^2}(\mathcal C)$, $\lambda_{\pi^\pm}$ and the Piecewise conditions are held fixed. The script checks that no coefficient of a kept $T_\Gamma$ contains a source letter, and lists the dressed coefficients on which the pole data depend, which are covered by step 1.
3. Negative control: the raw coefficients of the new frame are inserted to first order in $\kappa$.
4. $\Theta$: the raw and dressed $\Theta_{0,1}$ are compared.

## Result

1. **Dressed coefficients.** All 32 dressed coefficients $\widehat S^C_{r,s}$, $\widehat P^C_{r,s}$ with $r\le3$, $s\le3$ are identical in the two frames. The raw coefficients change at first order in $\kappa$ exactly for $s\ge1$: $\delta\widehat S^C_{r,s}=-\{\kappa,\widehat P^C_{r,s-1}\}$, $\delta\widehat P^C_{r,s}=\{\kappa,\widehat S^C_{r,s-1}\}$; the coefficients with $s=0$ do not change.
2. **Traces.** The source letters that enter the kept traces, $M_0$ and the pole data are $\widehat S^C_{0,0}$, $\widehat S^C_{0,1}$, $\widehat S^C_{0,2}$, $\widehat S^C_{1,0}$, $\widehat S^C_{1,1}$, $\widehat S^C_{1,2}$, $\widehat S^C_{2,0}$, $\widehat S^C_{2,1}$, $\widehat S^C_{3,0}$, $\widehat P^C_{0,1}$, $\widehat P^C_{1,1}$, $\widehat P^C_{2,0}$, all with $r\le3$ and $s\le3$. The pole data depend on $\widehat S^C_{0,0}$, $\widehat S^C_{0,2}$, $\widehat S^C_{1,0}$, $\widehat S^C_{2,0}$ and $\widehat P^C_{0,1}$. Every kept $T_\Gamma$ (all 37) is unchanged when the dressed coefficients of the new frame are inserted.
3. **Negative control.** With the raw coefficients of the new frame, to first order in $\kappa$, 20 of the 37 kept $T_\Gamma$ change: $C2$, $C3$, $C11$–$C15$, $C41$–$C48$, $C\sigma5$–$C\sigma9$. There are 27 kept graphs with a source of $s\ge1$. The seven that do not change ($C7$, $C9$, $C19$, $C21$, $C39$, $C\pi7$, $C\pi9$) contain only $\widehat S^C_{0,1}$ and $\widehat S^C_{1,1}$ among the sources with $s\ge1$. Their raw shifts $-\{\kappa,\widehat P^C_{0,0}\}$ and $-\{\kappa,\widehat P^C_{1,0}\}$ vanish because $\widehat P^C_{0,0}=\widehat P^C_{1,0}=0$ (assumptions 1 and 2) hold in both frames. The script checks that the set of changed graphs is exactly the set of graphs that contain a source of $s\ge1$ with a nonvanishing shift.
4. **$\widehat\Theta$.** $\widehat\Theta_{0,1}$ is unchanged, while the raw $\Theta_{0,1}$ shifts by $-\langle\kappa\rangle$; all other $\widehat\Theta_{r,s}$ are unchanged. The Cartesian form of $\widehat\Theta_1=0$ (assumption 6 of the frame record) is therefore frame independent.

The traces of the estimate are invariant under a constant axial redefinition for the same physical legs. The check covers the traces only; the size rule itself is not involved.

## Verification

Exact tests with noncommuting symbols: the transport of a constant connection is the exponential; $e^{ix\kappa}e^{-ix\kappa}=\mathbf 1$ order by order; the dagger of words; check 1 (32 coefficients, raw shifts equal to the anticommutators above, zero for $s=0$); the entering letters lie in the covered range; the coefficients of the traces contain no source letter; check 2 for all 37 graphs; the negative control (20 changed graphs, 27 with a source of $s\ge1$, the seven unchanged graphs and the letters $\widehat S^C_{0,1}$, $\widehat S^C_{1,1}$ with vanishing shift); check 4.

## Reproduction

```sh
PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_cartesian_loop_frame_invariance.py
PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_cartesian_loop_frame_invariance
```

[Script](../../../../../../scripts/derive_sigma_pipi_cartesian_loop_frame_invariance.py), [tests](../../../../../../scripts/test_derive_sigma_pipi_cartesian_loop_frame_invariance.py). The transport is written anew for the Cartesian variables; the procedure is the one of the [polar frame-invariance record](../one-loop-estimate/murayama-sigma-pipi-loop-frame-invariance.md), used for its rules only.
