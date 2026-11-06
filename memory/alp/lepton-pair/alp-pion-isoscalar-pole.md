---
title: Isoscalar pole residue and phase in the two-pion production wave
status: derived
created: 2026-11-06
updated: 2026-11-06
bibtex_keys: [Alda:2026nrd, Messerli:2025rnv, Watson:1952ji, Omnes:1958hv]
provenance: analytical-derivation
---

# Scope and conventions

This record derives the isoscalar photon-pole content of the charged two-pion production $P$ wave, its elastic phase relation, and its subtraction representation. It uses the source-resolved production numerators and current convention in [the pion production record](alp-pion-vector-and-production-form-factors.md), the no-mixing sources of Alda:2026nrd, and the single-cut normalization of Messerli:2025rnv [@Alda:2026nrd; @Messerli:2025rnv]. Use $S=1+iT$, metric $(+---)$, $\epsilon^{0123}=+1$, the real-epsilon crossed-current convention, and $\mathcal F_{a\pi\pi\gamma^*}=f_1/e$.

# Isoscalar pole continuation and residue

The retained photon response is $W_\gamma(y)=\sum_i R_iT_i/(M_i^2-y-i0)$. Since the VJP/VVP and VPPP production terms are linear in this response, write the projected production wave as $f_1(s,y)=f_{1,\mathrm{reg}}(s,y)+\sum_{I=\omega,\phi} C_I(s,y)/(M_I^2-y-i0)$, where $C_I$ is obtained by replacing $W_\gamma$ by $R_IT_I$ in the indexed matrix traces, applying the $3/4\int_{-1}^1(1-z^2)dz$ projector, and retaining the VJP/VVP and VPPP terms. The WZW and explicit-current pieces have no isoscalar photon pole. This coefficient extraction preserves arbitrary noncommuting source matrices; no relation between axial and pseudoscalar sources is used.

For a dressed isoscalar response, continue $y$ through the physical cut to the sheet containing the resonance pole $s_I$, defined by $M_I^2-s_I-\Delta_I^{\mathrm{self}}(s_I)=0$. If its numerator is $C_I(s,y)$ and the inverse propagator has a simple zero, the Laurent residue of the charge-stripped wave in the $(y-s_I)^{-1}$ convention is
$$
\operatorname*{Res}_{y=s_I}\mathcal F_{a\pi\pi\gamma^*}(s,y)
=-\frac{C_I(s,s_I)}{e[1+\Delta_I^{\mathrm{self}\prime}(s_I)]}.
$$
The minus sign follows from differentiating $M_I^2-y-\Delta_I^{\mathrm{self}}(y)$. The tree result follows with $\Delta_I^{\mathrm{self}}=0$ and $s_I=M_I^2-i0$. The complex residue factorizes as the analytically continued signed production coupling times the compatible $I\to\pi\pi$ decay/current coupling, divided by the inverse-propagator derivative; the same LSZ phases must be used in both factors. A simultaneous external-state rephasing cancels in their product.

# Elastic phase and subtractions

For the two-pion $P$ wave, write $f_\pm=f_1(s\pm i0,y)$ and let $\widehat f_1$ be the crossed-channel angular inhomogeneity. Elastic unitarity gives
$$
f_+-f_-=2i(f_++\widehat f_1)\sin\delta_1^1 e^{-i\delta_1^1},\qquad
f_+=e^{2i\delta_1^1}f_-+2i e^{i\delta_1^1}\sin\delta_1^1\widehat f_1.
$$
Thus Watson's phase is the homogeneous result only when the crossed inhomogeneity vanishes (or is accounted for in the corresponding Omnès solution). With $\Omega_+/\Omega_-=e^{2i\delta_1^1}$ and $|\Omega|=|\Omega_+|$, the quotient jump is $\operatorname{disc}(f_1/\Omega)=2i\widehat f_1\sin\delta_1^1/|\Omega|$. For $n$ subtractions at $s=0$, the solution is
$$
f_1(s,y)=\Omega(s)\left[P_{n-1}(s,y)+\frac{s^n}{\pi}\int_{4M_\pi^2}^{\infty}\frac{dx}{x^n}\frac{\widehat f_1(x,y)\sin\delta_1^1(x)}{|\Omega(x)|(x-s-i0)}\right],
$$
where $P_{n-1}$ contains the $n$ subtraction functions in $y$ and any explicitly retained left-hand-cut contribution must be added through its own dispersive integral. The indexed production amplitude has crossed $t/u$ terms, so their inhomogeneity is not set to zero. This formula specifies the phase and subtraction structure without assigning numerical subtraction constants.

For a simple pole term with residue $\mathcal R_I(s)$, subtraction at $y=y_0$ gives $\mathcal R_I[(y-s_I)^{-1}-(y_0-s_I)^{-1}]$ for one subtraction and subtracts additionally $(y-y_0)\partial_y[\mathcal R_I/(y-s_I)]_{y_0}$ for two. These identities preserve the pole residue. The remaining regular subtraction functions are independent analytic input functions and are retained symbolically.

# Reproducibility and limits

Run `.venv/bin/python scripts/derive_alp_pion_isoscalar_pole.py` and `.venv/bin/python scripts/derive_alp_pion_isoscalar_pole.py --test`. The SymPy script derives the response coefficient, pole derivative, Laurent and subtraction identities, and crossed-inhomogeneity Omnès jump. It does not solve the physical three-body/elastic integral equation or assign fitted subtraction functions. The physical amplitude is therefore specified analytically as a continued, source-resolved Omnès representation with a matched complex production/decay factorization; numerical evaluation is not part of this derivation.
