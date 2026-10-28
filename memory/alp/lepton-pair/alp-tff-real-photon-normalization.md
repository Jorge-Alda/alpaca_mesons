---
title: Real-photon normalization of the ALP transition form factor
status: derived
created: 2026-10-28
updated: 2026-10-28
bibtex_keys: [Alda:2026nrd, Bauer:2020jbp, Ecker:1988te, Ecker:1989yg, Messerli:2025rnv, Ruiz-Femenia:2003jdx]
provenance: analytical-derivation
---

# Convention and calculated real-photon coupling

Use the scalar convention of the [ALP tree TFF](alp-doubly-virtual-tff.md)
and Messerli:2025rnv,
$$
i\mathcal M=ie^2R F_{a\gamma^*\gamma^*}(s_1,s_2),\qquad
R=\epsilon_{\mu\nu\rho\sigma}\varepsilon_1^\mu\varepsilon_2^\nu
q_1^\rho q_2^\sigma,
$$
with metric $(+---)$ and $\epsilon^{0123}=+1$. Define
$$
F_{a\gamma\gamma}=F_{a\gamma^*\gamma^*}(0,0).
$$
This is the unnormalized TFF used throughout the memory and in the
[dilepton-loop extract](../../extracts/Messerli:2025rnv/dilepton-loop-formalism.md).

The source/WZW calculation gives the local photon coefficient
$C_{a\gamma\gamma}(m_a)$. The separate executable vector-photon calculation
gives the real-photon vector correction
$$
\Delta C_{a\gamma\gamma}^{\rm vector}
=-4\pi^2 f_a\bigl(G_{\rm VJP}+G_{\rm VVP}\bigr),
\qquad
C_{a\gamma\gamma}^{\rm eff}
=C_{a\gamma\gamma}+\Delta C_{a\gamma\gamma}^{\rm vector},
$$
where $G_{\rm VJP}$ and $G_{\rm VVP}$ are the explicitly calculated
zero-virtuality scalars in the $i e^2E G$ convention of the
[vector-photon record](../photon/alp-photon-vector-mixing.md). Since $E=-R$,
the corresponding vector contribution to this note's $F$ convention is
$$
F_{\rm vector}(0,0)=-(G_{\rm VJP}+G_{\rm VVP})
=\frac{\Delta C_{a\gamma\gamma}^{\rm vector}}{4\pi^2f_a}.
$$
The local contribution is
$$
\frac{C_{a\gamma\gamma}}{4\pi^2f_a}.
$$
Therefore the computed real-photon normalization, including the vector
correction, is
$$
\boxed{F_{a\gamma\gamma}^{\rm calc}
=\frac{C_{a\gamma\gamma}^{\rm eff}}{4\pi^2f_a}
=\frac{C_{a\gamma\gamma}+\Delta C_{a\gamma\gamma}^{\rm vector}}
{4\pi^2f_a}.}
$$
The meson-pole WZW contribution is already included in
$C_{a\gamma\gamma}$ in the no-mixing basis; do not add it again
[@Alda:2026nrd; @Bauer:2020jbp]. The vector correction above is the
additional VJP/VVP photon-conversion result from the existing calculation.

# Effective-pole normalization

Write the completed TFF as
$$
F_{a\gamma^*\gamma^*}=F_{\rm local}+F_{\rm low}^{I=1}
+F_{\rm low}^{I=0}+F_{\rm eff}+F_{\rm asym},
$$
where $F_{\rm low}^{I=1}$ and $F_{\rm low}^{I=0}$ are the explicitly
included low-energy spectral contributions. The local term is already
fixed, $F_{\rm local}(0,0)=C_{a\gamma\gamma}/(4\pi^2f_a)$. If omitted
strength is represented by the symmetric pole
$$
F_{\rm eff}(s_1,s_2)=g_{\rm eff}
\frac{M_{\rm eff}^4}
{(M_{\rm eff}^2-s_1)(M_{\rm eff}^2-s_2)},
$$
then $F_{\rm eff}(0,0)=g_{\rm eff}$. The normalization condition
$F_{a\gamma\gamma}=F_{a\gamma\gamma}^{\rm calc}$ gives
$$
\boxed{g_{\rm eff}
=\frac{\Delta C_{a\gamma\gamma}^{\rm vector}}{4\pi^2f_a}
-F_{\rm low}^{I=1}(0,0)-F_{\rm low}^{I=0}(0,0)-F_{\rm asym}(0,0).}
$$
Here the known contact term has canceled between the full normalization
$F_{a\gamma\gamma}^{\rm calc}=(C_{a\gamma\gamma}
+\Delta C_{a\gamma\gamma}^{\rm vector})/(4\pi^2f_a)$ and
$F_{\rm local}(0,0)$. Thus the pole residue supplies the vector correction
and any other origin strength required by the chosen spectral completion,
after subtracting the low-energy and asymptotic terms already included. If
the low-energy spectra represent some or all of the vector strength, that
overlap is subtracted in the displayed terms and counted only once.

This fixes the pole residue under the stated one-pole model; it does not fix
$M_{\rm eff}$ or establish that one pole accurately represents the omitted
spectrum. The existing vector correction uses the documented provisional
odd-coupling fit, so its numerical value inherits that fit's assumptions and
does not include omitted operators, widths, or loops
[@Ruiz-Femenia:2003jdx; @Ecker:1988te; @Ecker:1989yg].

# Relation to the dilepton cut

The two-real-photon cut constrains the imaginary part of the reduced
dilepton-loop amplitude after normalizing the TFF by $F_{a\gamma\gamma}$,
when nonzero [@Messerli:2025rnv]. It does not determine
$F_{\rm low}^{I=1}(0,0)$, $F_{\rm low}^{I=0}(0,0)$, or the pole residue. Here the residue follows from the
real-photon normalization sum above, not from the dilepton cut. If
$F_{a\gamma\gamma}=0$, the normalized loop representation is undefined;
the unnormalized TFF and loop remain the appropriate quantities.

# Derivation and inputs

The relation between the vector correction and the TFF follows directly
from $E=-R$: $i e^2E(G_{\rm VJP}+G_{\rm VVP})$ equals
$i e^2R[-(G_{\rm VJP}+G_{\rm VVP})]$. Substituting the vector record's
definition of $\Delta C_{a\gamma\gamma}^{\rm vector}$ gives the displayed
factor $1/(4\pi^2f_a)$. At zero virtuality, the effective-pole denominator
equals one, so its value is $g_{\rm eff}$; requiring the sum to equal
$F_{a\gamma\gamma}^{\rm calc}$ yields the boxed residue.

Inputs are the [source/WZW photon coefficient](../photon/alp-photon-invariant-from-mixing.md),
the [calculated vector-photon correction](../photon/alp-photon-vector-mixing.md),
the [canonical tensor conversion](../photon/alp-photon-odd-action-phase-matching.md),
and the [tree TFF convention](alp-doubly-virtual-tff.md). The two low-energy
spectral contributions and the ALP asymptotic term are not yet calculated.
The next step is to calculate the two spectral contributions and the
source-specific asymptotic term, then use the residue equation with a
pole/continuum separation that counts each vector contribution once.
