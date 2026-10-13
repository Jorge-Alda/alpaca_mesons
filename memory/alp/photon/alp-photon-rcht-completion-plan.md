---
title: Complete the resonance-corrected physical ALP two-photon amplitude
status: mixed
created: 2026-10-13
updated: 2026-10-13
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx, Ecker:1988te, Ecker:1989yg, Bai:2025fvl]
provenance: assumption
---

# Physical two-photon ALP completion plan

## Starting evidence and independent inputs

The [no-vector WZW/contact invariant](alp-photon-invariant-from-mixing.md)
already derives

$$
C_{a\gamma\gamma}=c_\gamma^{\rm inv}
+2N_c\operatorname{tr}(Q^2\mathbb M_a^{-1}(\mathcal C)),\qquad
\mathcal L\supset\frac{\alpha_{\rm em}}{4\pi f_a}
C_{a\gamma\gamma}aF\widetilde F
$$

under its stated photon-anomaly normalization. The [two-source insertions](../../rcht/odd-parity/two-source-alp-rcht-odd-parity-extension.md)
and [operator kernel](../../rcht/odd-parity/odd-parity-rcht-numerical-kernel.md)
do not add up to a normalized resonance-corrected ALP predictor. The
[real-photon sensitivity audit](alp-photon-rcht-blind-direction-audit.md)
is zero along all three current structural blind directions but is stripped
of the overall physical normalization.

Use Alda arXiv-v1 $X=\mathbb M_a^{-1}(\mathcal C)$ and $Y=c_q^{AP}+X$
independently. The WZW/contact term, derivative resonance term and pseudoscalar
mass-source term do not all have the same effective flavor matrix. UV ALP
coefficients and any independently matched invariant local photon term remain
external inputs rather than SM-fit outputs.

## Ordered calculation

1. Declare the retained tree/single-trace action and leading $e^2/f_a$ order.
   Fix the dual-field convention and the invariant UV/contact coefficient.
   Verify Alda's source signs against the archived arXiv formulas; do not use
   the local manuscript or infer a singlet matching coefficient from the anomaly.
2. In `scripts/derive_physical_alp_gamma_gamma.py`, independently derive the
   VJP term with each photon placement and the VVP term with both vector
   assignments, using actual Fourier/tensor contractions. Insert $W(0)$ into
   each internal vector slot and retain the explicit photon $Q$ in VJP.
   Mass traces use $Y$ while derivative traces use $X$.
3. Assemble the WZW/contact invariant with these resonance pieces in one
   field basis. Compare direct-plus-meson-pole LSZ with the no-mixing basis;
   never add a meson-pole anomaly contribution already absorbed into the
   transformed photon coefficient. Retain any necessary induced local terms,
   but avoid adding resonance-saturated local ChPT terms a second time.
4. Derive the signed coefficient of
   $\epsilon_{\alpha\beta\rho\sigma}\varepsilon_1^\alpha
   \varepsilon_2^\beta q_1^\rho q_2^\sigma$ including the explicit $e^2$
   convention. Check both Ward identities, Bose symmetry, field-redefinition
   invariance and the zero-resonance limit of the existing invariant.
5. Derive the identical-photon polarization sum and width. Specify whether
   $e^2$ and $1/f_a$ are inside or outside the scalar form factor; reproduce
   the direct $aF\widetilde F$ normalization in an independent limiting test.
6. Project the completed resonance gradient onto the SM blind columns without
   OPE substitutions, for independent $X,Y$ and every retained vector/source
   hypothesis. The current zero projection must persist in the declared
   completion; a failure requires identifying the additional operator or
   missing term rather than setting a blind coefficient by hand.

## Deliverables and later extensions

Pair the executable with a derivation record and numerical predictor; include
Ward/Bose, zero-resonance, reparametrization, normalization and null-shift tests.
State the pointwise ALP–meson near-degeneracy domain from Alda, and replace the
perturbative inverse with a coupled treatment for retained violating points.
Virtual photons, loops, further singlet/OZI terms and high-energy matching
are separate extensions, not implied by the real-photon tree calculation.
