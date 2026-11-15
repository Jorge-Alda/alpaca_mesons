---
title: Light-cone ALP transition form factor at leading twist
status: derived
created: 2026-11-15
updated: 2026-11-15
bibtex_keys: [Alda:2026nrd, Hoferichter:2020lap, Ludtke:2020moa]
provenance: analytical-derivation
---

# Inputs and assumptions

Use the tensor, metric and no-mixing conventions in [the conventions record](alp-high-energy-tff-conventions.md), natural units, and $Q_i^2=-q_i^2>0$. The light axial-current matrix elements at the low scale follow from [the LO chiral derivation](alp-axial-current-matrix-element.md); their singlet evolution and thresholds are as in [the symmetric result](alp-symmetric-high-energy-tff.md). The hard coefficient below is at tree level in $\alpha_s$ and leading twist, with $m_a^2/Q_i^2$ corrections omitted. The nonlocal axial matrix element includes a straight Wilson line. Its local limit defines the normalized distribution amplitude. This uses the leading-twist factorization in [Hoferichter:2020lap](../../extracts/Hoferichter:2020lap/pseudoscalar-light-cone-tff.md), not its meson-specific coefficient.

# Derivation

For each light flavor define its distribution amplitude by the nonlocal axial matrix element, normalized so its integral is one and its local limit is $\langle0|\bar q_f\gamma^\mu\gamma_5q_f|a(p)\rangle$. Contracting the two electromagnetic currents with a free hard quark gives denominators from the two photon orderings. For charge-conjugation-even distributions the orderings combine to give
$$
F_{a\gamma^*\gamma^*}\big|_{\rm light,twist\,2}=
\frac{4f_\pi^2}{f_a}\sum_{f=u,d,s}Q_f^2[\mathbb M_a^{-1}(\bm{\mathcal C})]_{ff}(\mu)
\int_0^1du\frac{\phi_f(u,\mu)}{uQ_1^2+(1-u)Q_2^2}.
$$
The diagonal entries stand for the renormalized current matrix elements; the chiral expression supplies their low-scale approximation. For general current matrix elements, replace $2f_\pi^2[\mathbb M_a^{-1}(\bm{\mathcal C})]_{ff}/f_a$ by $p_\mu\langle0|\bar q_f\gamma^\mu\gamma_5q_f|a\rangle/(i m_a^2)$. This expression is linear in the ALP source matrix elements and does not use a meson decay constant as ALP input. Photon contact and heavy-quark triangle contributions are the separate terms in [the symmetric assembly](alp-symmetric-high-energy-tff.md); do not add an overlapping heavy-current contribution to the explicit triangle.

As a **shape assumption**, use independently for each flavor
$$
\phi_f(u,\mu)=6u(1-u)\left[1+a_2^f(\mu)C_2^{3/2}(2u-1)\right].
$$
The named Gegenbauer moment $a_2^f$ is an input, not fixed by the local current or by LO ChPT. Setting it to zero gives the asymptotic-shape assumption. The truncated ansatz is illustrative; higher moments and singlet quark-gluon evolution are omitted. No positivity constraint is imposed.

The script integrates the assumed distribution and gives, for each flavor,
$$
\int_0^1du\,\phi_f=1,\qquad
\int_0^1du\frac{\phi_f}{uQ^2+(1-u)Q^2}=\frac1{Q^2},\qquad
\int_0^1du\frac{\phi_f}{uQ^2}=\frac{3(1+a_2^f)}{Q^2}.
$$
Thus equal spacelike virtualities reproduce the light axial-current term of task 3.1 at common tree order; its $1-\alpha_s/\pi$ correction is not the unequal-virtuality NLO coefficient. The singly virtual result depends on $a_2^f$, whereas the symmetric result does not. The other photon ordering follows from $u\leftrightarrow1-u$.

# Domain and verification

The fixed-ratio expression applies when both spacelike virtualities are hard with a finite nonzero ratio. Its one-real-photon endpoint is the leading-twist factorization limit, conditional on endpoint convergence and factorization; it is not a local OPE at fixed soft virtuality. Near thresholds the named current matrix elements and heavy contributions retain the limitations of the symmetric record. The renormalized currents use the $\overline{\rm MS}$ convention and active-flavor intervals described there; distribution amplitudes and their singlet mixing require corresponding scale-dependent inputs. The explicit contact, massive triangles and higher-twist corrections are not described by the convolution.

Run `.venv/bin/python scripts/derive_alp_light_cone_tff.py` and `.venv/bin/python scripts/derive_alp_light_cone_tff.py --test`. The script integrates the Gegenbauer ansatz, verifies normalization, exchange symmetry, equal and singly virtual limits, and a Ward contraction. The tensor is antisymmetric, so both Ward identities follow by the same argument.
