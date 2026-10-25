---
title: Isoscalar vector to three-pion partial wave for the ALP TFF
status: derived
created: 2026-10-25
updated: 2026-10-25
bibtex_keys: [Alda:2026nrd, Niecknig:2012sj, Hoferichter:2023bjm, Hoferichter:2025lcz]
provenance: analytical-derivation
---

# Scope and conventions

This record identifies the on-shell $V_i\to\pi^+\pi^-\pi^0$ decay function
needed by the isoscalar ALP TFF pole residue, for $V_i=\omega,\phi$. Use
$g_{\mu\nu}=(+---)$, $\epsilon^{0123}=+1$, and the no-mixing ALP source
convention of Alda:2026nrd [@Alda:2026nrd]. The related strict spin-one
factorization and its production numerator are derived in the
[pole-factorization record](alp-three-pion-vector-pole-factorization.md).

The dispersive decay amplitude is parameterized as
$$
\mathcal M(V_i\to3\pi)=i\epsilon_{\mu\nu\alpha\beta}
n^\mu p_+^\nu p_-^\alpha p_0^\beta
\mathcal F_i(s,t,u),\qquad
\mathcal F_i(s,t,u)=\mathcal F_i(s)+\mathcal F_i(t)+\mathcal F_i(u),
$$
where $n^\mu$ is the vector polarization and $s=(p_-+p_0)^2$,
$t=(p_++p_0)^2$, $u=(p_++p_-)^2$. Thus the scalar decay function is symmetric
under pion-pair permutations. The dominant two-pion partial wave in any one
pair channel is
$$
f_{1,i}(s)=\mathcal F_i(s)+\widehat{\mathcal F}_i(s),
$$
with $\widehat{\mathcal F}_i$ the crossed-channel angular average. The same
function is used in the transverse spin-one current
$$
J^\beta_{3\pi}=\epsilon^{\beta\lambda\kappa\delta}
p_{+\lambda}p_{-\kappa}p_{0\delta},\qquad
\mathcal J^\beta_{V_i\to3\pi}=F_i^{\rm on}(t,u)J^\beta_{3\pi},
$$
where $F_i^{\rm on}$ denotes the full on-shell Bose-symmetric scalar
amplitude, including its normalization in the selected S-matrix convention.
The notations $f_{1,i}$ and $F_i^{\rm on}$ refer to the same decay solution
after the partial-wave reconstruction, not to the electromagnetic-current
partial wave's fitted normalization.

# Dispersive equation and result

Elastic $\pi\pi$ P-wave rescattering gives the single-variable equation
[@Niecknig:2012sj]:
$$
\operatorname{disc}\mathcal F_i(s)=2i\bigl(\mathcal F_i(s)
+\widehat{\mathcal F}_i(s)\bigr)\sin\delta_1^1(s)e^{-i\delta_1^1(s)},
$$
and the once-subtracted Omnès representation
$$
\mathcal F_i(s)=\Omega(s)\left[a_i+\frac{s}{\pi}\int_{4m_\pi^2}^\infty
\frac{ds'}{s'}\frac{\sin\delta_1^1(s')\widehat{\mathcal F}_i(s')}
{|\Omega(s')|(s'-s-i0)}\right],
\quad
\Omega(s)=\exp\left\{\frac{s}{\pi}\int_{4m_\pi^2}^\infty
\frac{\delta_1^1(s')}{s'(s'-s-i0)}ds'\right\}.
$$
The crossed angular average couples the three pair channels and must be solved
self-consistently. A twice-subtracted form introduces an additional complex
subtraction parameter; the sum rule and fit choice are part of the amplitude
prescription, not fixed by pole factorization.

The factorized ALP production residue multiplies $F_i^{\rm on}(t,u)$ and
the transverse tensor contraction given in the linked pole-factorization
record. Its relative phase $\zeta_i$ remains explicit until the on-shell
$a\to V_i\gamma$ and $V_i\to3\pi$ amplitudes are matched in one convention.
This fixes the on-shell pole residue only. It supplies no regular off-shell
continuation or direct $V3\pi$ contact term.

# Inputs and provenance

- The P-wave phase $\delta_1^1$ and crossed-channel equation are literature
  inputs from [Niecknig:2012sj](../../extracts/Niecknig:2012sj/dispersive-three-pion-analysis.md).
- The current decay-function decomposition and the fact that its fitted
  normalization is distinct from the partial-wave shape follow
  [Hoferichter:2023bjm](../../extracts/Hoferichter:2023bjm/three-pion-amplitude-and-normalization.md).
- Updated current-fit choices and limits on available complex grids are in
  [Hoferichter:2025lcz](../../extracts/Hoferichter:2025lcz/gamma-star-to-three-pion-fit.md).
- The ALP-side pole numerator and tensor projection are in the
  [factorization record](alp-three-pion-vector-pole-factorization.md).

The SM $\gamma^*\to3\pi$ current solution is not an ALP production amplitude.
Likewise, its fitted normalization cannot replace the source-dependent
$a\gamma^*\to3\pi$ side of the unitarity product.

# Verification and limits

The approved invocation
`.venv/bin/python scripts/derive_alp_three_pion_vector_pole.py` completed.
SymPy returned zero for the inverse-Hessian decomposition residual, the local
tensor-projector term on the transverse $k\wedge J_{3\pi}$ source, the three-
pion current divergence, both production Ward contractions, and the spin-one
projector reduction residual. It derived the ordered VJP/VVP numerator and
isolated the strict pole coefficient $N_i(M_i^2,y)/(M_i^2-s)$. The script
does not solve the Omnès integral equation or evaluate a numerical decay
function.

Neither the 2023/2025 current-fit records nor the cited decay analysis as
represented in this repository provides a ready complex $F_i^{\rm on}(t,u)$
grid with correlated fit/phase replicas. Consequently the established result
here is the partial-wave equation and its insertion into the strict ALP pole
residue; a numerical spectral integral remains out of scope until that grid,
its normalization, and the matched relative phase are supplied or reproduced.
