---
title: Pion-source check of the light-cone form factor
status: derived
created: 2026-11-15
updated: 2026-11-15
bibtex_keys: [Alda:2026nrd, Hoferichter:2020lap]
provenance: analytical-derivation
---

# Inputs and assumptions

Use natural units and the current and tensor conventions of the [ALP light-cone calculation](alp-light-cone-high-energy-tff.md). This check uses the pion state phase and $F_\pi$ normalization in [Hoferichter:2020lap](../../extracts/Hoferichter:2020lap/pseudoscalar-light-cone-tff.md), with $\langle0|\bar q\gamma_\mu\gamma_5\lambda^3 q/2|\pi^0(p)\rangle=ip_\mu F_\pi$. Assume isospin symmetry, vanishing isoscalar pion matrix elements, tree-level leading twist, and the symmetric normalized distribution in the ALP calculation. The singly virtual asymptotic limit additionally assumes $a_2=0$ and vanishing higher moments. The overall sign changes if the pion state phase is reversed.

# Calculation

The pion axial-current coefficients for $u,d,s$ are $F_\pi,-F_\pi,0$. Inserting these into the same charge-weighted current and convolution used in the ALP calculation gives
$$
\sum_f Q_f^2\frac{\langle0|\bar q_f\gamma_\mu\gamma_5q_f|\pi^0(p)\rangle}{ip_\mu}=\frac{F_\pi}{3}.
$$
The division by $p_\mu$ denotes the scalar coefficient of the four-vector matrix element. For $Q^2=-q^2>0$, the script produces
$$
F_{\pi^0\gamma^*\gamma^*}(-Q^2,-Q^2)=\frac{2F_\pi}{3Q^2},\qquad
F_{\pi^0\gamma^*\gamma}(-Q^2,0)=\frac{2F_\pi(1+a_2)}{Q^2}.
$$
At $a_2=0$ the second expression is $2F_\pi/Q^2$. Both reproduce the corresponding formulas of Hoferichter:2020lap in the same state phase. The singly virtual result assumes leading-twist endpoint factorization; finite-mass and higher-twist corrections are omitted. The pion inputs are used only in this check: the ALP coefficient remains derived from its own matrix elements and sources.

# Verification

Run `.venv/bin/python scripts/derive_pion_light_cone_crosscheck.py` and `.venv/bin/python scripts/derive_pion_light_cone_crosscheck.py --test`. The script imports the distribution and integrals from the ALP derivation, computes the flavor trace, and tests the two limits and shape dependence.
