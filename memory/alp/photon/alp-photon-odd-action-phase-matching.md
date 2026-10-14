---
title: Signed photon-dual dictionary and bounded odd-action assembly
status: mixed
created: 2026-10-14
updated: 2026-10-14
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx, Bai:2025fvl]
provenance: analytical-derivation
---

# Inputs and evidence status

Inputs are Ruiz's [negative WZW derivative vertex](../../extracts/Ruiz-Femenia:2003jdx/pseudoscalar-two-photon-form-factor.md),
the [normalized SM TFF](../../rcht/odd-parity/radiative-tff-amplitude.md),
the [bounded neutral cubic gradients](../../rcht/odd-parity/odd-power-counted-channel-maps.md),
and Alda's [source/meson LSZ replacement](../../rcht/odd-parity/odd-alp-source-lsz-assembly.md).
The positive field conversion is $\Phi_A=\sqrt2\Phi_R$, with $f_\pi=F$.
The [inherited photon record](alp-photon-invariant-from-mixing.md) assumes
a positive pseudoscalar coefficient of $F\widetilde F$ and a photon Jacobian
shift. Its displayed dual formula originally did **not** specify
the orientation of $\epsilon^{0123}$. Completing that orientation below is
a convention choice required by its assumed coefficient and the source
action; it is not a photon-anomaly result supplied by Alda, which excludes WZW.

The [Bai digest](../../extracts/Bai:2025fvl/background-fields-wzw.md) establishes
the need for the full gauged WZW functional and its gauge counterterm.
This calculation matches the existing photon-only assumption to the tensor
action. It does not import a background-vector model or match new resonance
anomaly coefficients. A UV/local photon coefficient remains an input in its
explicitly declared dual convention.

# Tensor calculation

Run the deterministic SymPy derivation:

```text
.venv/bin/python scripts/derive_odd_photon_phase_matching.py
```

Use $g_{\mu\nu}=\operatorname{diag}(1,-1,-1,-1)$ and
$\epsilon^{\rm src}_{0123}=+1$, hence
$\epsilon^{0123}_{\rm src}=-1$ after metric raising. All momenta are incoming,
$\partial\to iq$. Define the source-basis tensors

$$
E=\epsilon_{\mu\nu\rho\sigma}
q_1^\mu\varepsilon_1^\nu q_2^\rho\varepsilon_2^\sigma,
\qquad R=\epsilon_{\mu\nu\rho\sigma}
\varepsilon_1^\mu\varepsilon_2^\nu q_1^\rho q_2^\sigma,
\qquad E=-R.
$$

Let $\epsilon^{\mu\nu\rho\sigma}_\gamma
=\eta_\gamma\epsilon^{\mu\nu\rho\sigma}_{\rm src}$ define the photon dual.
Explicit metric lowering and both photon Wick assignments give

$$
\frac{\langle F_{\mu\nu}\widetilde F_\gamma^{\mu\nu}\rangle}{E}
=-4\eta_\gamma,
\qquad
\frac{\langle\epsilon_{\mu\nu\rho\sigma}
\partial^\mu A^\nu\partial^\rho A^\sigma\rangle}{E}=-2.
$$

Matching Ruiz's coefficient $-\sqrt2N_c/(8\pi^2F)$ of the second operator
to the inherited positive $\sqrt2N_c/(16\pi^2F)$ coefficient of the first
solves

$$
\eta_\gamma=-1,\qquad \epsilon^{0123}_\gamma=+1.
$$

The inherited photon dual is therefore oppositely oriented to the
metric-raised source tensor. Their two epsilon symbols cannot be related
by metric raising. Equivalently, in the common raised-source dual,

$$
\widetilde F_\gamma=-\widetilde F_{\rm src},\qquad
c_\gamma^{\rm src}=-c_\gamma,\qquad
C_{a\gamma\gamma}^{\rm src}=-C_{a\gamma\gamma}.
$$

The inherited shift is $c_\gamma'=c_\gamma-2N_c\operatorname{tr}(Q^2\kappa)$;
the common-dual coefficient accordingly shifts with the opposite sign.
This translates photon coefficients without changing Alda's $c_G$, axial
and mass-source definitions, or the physical $c_i,d_i$.

# Common amplitude and LSZ phase

The bounded gradients $G_j$ multiply $i e^2E$ in the source field basis.
The legacy positive SM scalar $F_P$ equals that scalar coefficient: its
$i e^2R F_P$ definition includes the opposite external-meson phase.
Using that rephased meson amplitude in an ALP LSZ sum requires the opposite
meson-overlap phase as well. Flipping only WZW changes physical interference.

For $S=\mathbb M_a^{-1}(\mathcal C)$ and $H=c_q^A+S$, the source-meson
replacement $\Phi_R/a=\sqrt2F H/f_a$ gives

$$
G_{\rm WZW/mix}=\frac{2N_c\operatorname{tr}(Q^2H)}{4\pi^2f_a},
\qquad C_{a\gamma\gamma}=c_\gamma+2N_c\operatorname{tr}(Q^2H).
$$

The signed tree amplitude under the retained photon assumption is

$$
G_a=\frac{C_{a\gamma\gamma}}{4\pi^2f_a}
+\sum_{j=1}^{241}g_jG_j(U_a,Y_a;W_1,W_2),
\qquad i\mathcal M=i e^2E\,G_a=-i e^2R\,G_a,
$$

$$
U_a=-2S/f_a,\qquad Y_a=2\{c_q^{AP}+S,\chi\}/f_a.
$$

The 241 raw columns are a redundant presentation of the bounded 228
operator quotient, rather than 241 independent fit parameters. Both VJP
photon placements and VVP assignments remain. Swap the local momentum
polynomial before inserting $W_2(q_2^2)$; swapping its already-inserted
pole denominator would put it on the wrong photon. Responses stay signed
and are not expanded across poles.

Direct-plus-meson insertion equals the no-mixing result for every raw photon
parent. The contact becomes $\widehat c_\gamma=C_{a\gamma\gamma}$.
Adding a meson anomaly or resonance contribution to the no-mixing result
again double counts its already-included LSZ overlaps.

# Verification and limits

The [exact artifact](../../../data/vector_sector/odd_photon_phase_matching.json)
and [SymPy transcript](../../../data/vector_sector/odd_photon_phase_matching.md)
record the tensor dictionary, zero-resonance meson LSZ limit, full original
eleven-coupling SM recovery without OPE, 241 photon-parent LSZ identities,
241 photon-parent Bose identities, two contact Ward identities, coefficient
rotation and a nonzero anomaly–resonance interference example. Flipping both
pieces preserves its squared amplitude; flipping only the anomaly changes
it by the derived nonzero interference difference.

This closes the shared signed-convention gate. The photon normalization and
leading nonet WZW input remain assumptions. The tree/one-ALP/bounded-action
limits, LO resolvent pole exclusions and perturbative ALP–meson overlap
requirement remain. No extra gauged-WZW matching, completed width, numerical
predictor, measurement rank or fit is inferred.
