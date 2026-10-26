---
title: Gauged Wess-Zumino-Witten gamma three-pion contribution
status: derived
created: 2026-10-26
updated: 2026-10-26
bibtex_keys: [Wess:1971yu, Ruiz-Femenia:2003jdx, Ecker:1988te, Alda:2026nrd]
provenance: analytical-derivation
---

# Inputs and conventions

This record expands the source-linear, cubic-Goldstone part of the full
gauged Wess-Zumino functional and its photon--vector-mixing contribution.
The Wess--Zumino functional is the anomaly-matching term in
[@Wess:1971yu]. The RChT setup identifies it as $Z_{\rm WZ}[v,a]$ and
separates it from chirally invariant resonance operators
[@Ruiz-Femenia:2003jdx].

Use $U=\exp(i\Phi/f_\pi)$, the charged-pion matrices and source convention
of [the Alda:2026nrd expansion](../../extracts/Alda:2026nrd/appendix-a-expansion.md),
natural units, and the repository's canonical Levi-Civita convention
$\epsilon^{0123}=+1$. Write the electromagnetic source as
$v_\mu=-eQA_\mu$, with
$Q=\operatorname{diag}(2/3,-1/3,-1/3)$. The sign of the Wess--Zumino
functional's overall orientation is left explicit; reversing it reverses
all anomalous vertices together.

# The $\gamma3\pi$ term in the full functional

The full gauged functional contains a term linear in the vector source whose
leading three-Goldstone part can be represented by
$$
Z_{\rm WZ}[v,a]\big|_{v\Phi^3}
=\frac{N_c}{48\pi^2}\,\mathrm{orientation}
\int d^4x\,\epsilon^{\mu\nu\rho\sigma}
\left\langle v_\mu(R_\nu R_\rho R_\sigma+L_\nu L_\rho L_\sigma)\right\rangle,
\qquad
R_\mu=(\partial_\mu U)U^\dagger,\quad
L_\mu=U^\dagger\partial_\mu U.
$$
At cubic order, $R_\mu=L_\mu=i\partial_\mu\Phi/f_\pi+O(\Phi^2)$.
The charged-pion exterior-product flavor trace is
$$
\sum_{\pi^+,\pi^-,\pi^0\ {\rm permutations}}
\operatorname{sgn}(\text{permutation})
\operatorname{Tr}(Q\,T_{\text{permuted }+}
T_{\text{permuted }-}T_{\text{permuted }0})=2.
$$
Consequently, in the script's field and source normalization,
$$
\mathcal L_{\gamma3\pi}
=i\,\mathrm{orientation}\,
\frac{N_c e}{12\pi^2f_\pi^3}
\epsilon^{\mu\nu\rho\sigma}A_\mu
\partial_\nu\pi^+\partial_\rho\pi^-\partial_\sigma\pi^0.
$$
For $N_c=3$, its magnitude is the familiar anomaly coefficient
$e/(4\pi^2f_\pi^3)$. This is a direct anomalous $\gamma3\pi$ vertex. The
previous version of the [WZW extract](../../extracts/Ruiz-Femenia:2003jdx/wzw-anomaly-contribution.md)
incorrectly stated that this vertex was absent.

# What photon--vector diagonalization induces

In antisymmetric-tensor RChT, $Z_{\rm WZ}[v,a]$ is a functional of external
sources; the vector resonance is a separate field. The even-parity mixing
term and first-order photon diagonalization are given in the
[photon--vector record](../foundations/rcht-photon-vector-diagonalization.md).
At first order in the mixing, the source photon contains a vector-field
piece proportional to $g_i/M_i=\sqrt2eR_i/M_i$. Substituting this into
$\mathcal L_{\gamma3\pi}$ induces a local $V_i3\pi$ term whose coefficient
is proportional to
$$
\frac{N_c e}{12\pi^2f_\pi^3}\frac{g_i}{M_i}
=\frac{N_c\sqrt2e^2R_i}{12\pi^2M_if_\pi^3},
$$
up to the common orientation and field-redefinition sign. It is of order
$e^2$ in this source convention. It is therefore an electromagnetic,
field-basis-induced contribution, not the leading strong $V_i\to3\pi$
amplitude. It must not be counted again if the photon-mediated contribution
is retained explicitly. The mixing residue $R_i$ can be nonzero for ideal
$\omega$ and $\phi$. This photon--vector-mixing contribution is distinct
from strong local odd-parity resonance interactions; see the separate
[vector three-pion contact note](vector-three-pion-contact.md).

# Reproducibility

The executable
[derive_wzw_gamma_three_pion_and_vector_contact.py](../../../scripts/derive_wzw_gamma_three_pion_and_vector_contact.py)
computes the charged flavor trace and anomaly coefficient. It also contains
the separate contact-operator checks documented in the
[vector three-pion contact note](vector-three-pion-contact.md). Run:

    .venv/bin/python scripts/derive_wzw_gamma_three_pion_and_vector_contact.py

Six focused checks pass. The WZW-specific checks reproduce the charged
flavor trace and anomaly coefficient; they do not enumerate the full odd
sector.
