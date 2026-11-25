---
title: Constant axial redefinition of the radial-mode loop-size inputs
status: derived
created: 2026-11-23
updated: 2026-11-23
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: analytical-derivation
---

# Frame invariance of the radial-mode loop-size inputs

This record is task 3.4 of the OpenSpec change `phi-dependent-loop-tree-crossing-murayama`. It extends the [ALP-only frame-invariance record](alp-loop-estimate-frame-invariance.md) to the dressed radial invariants and to the contractions that the radial counting record uses. No loop value is computed.

## Setting

- A constant Hermitian axial generator $\kappa$ acts on the spurion $M(\sigma,a)$ of [the radial coefficient record](alp-murayama-radial-phase-mass-coefficients.md).
- The redefinition depends on $a$ only. The transport is integrated along $x=a/f_a$ at fixed $y=\sigma/f_a$. The radial mode is not rotated, and no transport in the $\sigma$ direction is needed for this redefinition. Redefinitions that depend on $\sigma$ are not covered.

## Checks

The [executable check](../../../scripts/derive_murayama_radial_frame_invariance.py) and its [unit tests](../../../scripts/test_derive_murayama_radial_frame_invariance.py) verify:

1. **Dressed invariants.** The dressed $\widehat S_{rs}$, $\widehat P_{rs}$ and $\widehat\Theta_{rs}$ are unchanged under the transport, for general noncommuting sources (all 36 residuals for $r\le2$, $s\le3$, zero) and in the connection-free polar representative ($r\le2$, $s\le6$, all 63 residuals zero).
2. **Pole data.** $\mathcal C$, $\mathbb M_{m_i^2}^{-1}(\mathcal C)$, $N(m_i^2)$ and $m_\sigma^2$ are unchanged, in the diagonal sector.
3. **Estimate contractions.** The $\sigma$-PP coupling $g_{ij}$, the $\sigma aP$ vertex and contraction, and the class-R2 traces $\operatorname{tr}\widehat S_{20}\{\Omega_i,\Omega_j\}$, $\operatorname{tr}\widehat P_{21}(\Omega_i+\Omega_j)$, $\operatorname{tr}\widehat S_{22}$ and $\operatorname{tr}\widehat S_{20}$ are unchanged, in the diagonal sector.
4. **No raw inputs.** The estimate uses only the dressed invariants and the pole data, so it is not built from the raw sources or connections of any frame. The untransported coefficients are not used anywhere (they are frame dependent, which is why they are excluded).

## Limits

- Items 2 and 3 are checked in the flavour-diagonal sector with diagonal $\kappa$, where $\mathbb M_s^{-1}$ can be written out. With commuting entries the equality is close to automatic. The substantive statement for noncommuting sources is item 1.
- The class-R2 vertex is not stated as one assembled expression in the notes. The check uses its separate traces.
- The singlet mass appears as $m_0^2/N_f$ in the ALP-only pole record and as $m_1^2/3$ in the radial record. The script uses one symbol for both.
- The invariance of the contractions for the other classes follows from the common dressed inputs, by the same argument as the ALP-only record, step 4. It is not checked symbolically here.

## Status

Verified symbolically for the items above, with the limits stated. The loop-size estimate is therefore invariant under a constant axial redefinition at the level of its inputs.
