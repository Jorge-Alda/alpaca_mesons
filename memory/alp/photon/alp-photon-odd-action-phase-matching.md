---
title: Signed photon-dual dictionary and bounded odd-action assembly
status: mixed
created: 2026-10-14
updated: 2026-10-15
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx, Bai:2025fvl, Bauer:2020jbp]
provenance: analytical-derivation
---

# Inputs and evidence status

Inputs are the canonically adapted [WZW derivative vertex](../../extracts/Ruiz-Femenia:2003jdx/pseudoscalar-two-photon-form-factor.md),
the [normalized SM TFF](../../rcht/odd-parity/radiative-tff-amplitude.md),
the [bounded neutral cubic gradients](../../rcht/odd-parity/odd-power-counted-channel-maps.md),
and Alda:2026nrd's [source/meson LSZ replacement](../../rcht/odd-parity/odd-alp-source-lsz-assembly.md).
The positive field conversion is $\Phi_A=\sqrt2\Phi_R$, with $f_\pi=F$.
The [photon record](alp-photon-invariant-from-mixing.md) uses the source-backed
positive photon-operator normalization of [@Bauer:2020jbp, Eqs. (39)/(54)]
and a separately declared photon Jacobian shift and leading-nonet WZW input.
The canonical convention is the single metric-related tensor
$\epsilon^{0123}=+1$, following [@Bauer:2020jbp, after Eq. (1)]. The
[conversion audit](../../rcht/odd-parity/levi-civita-convention.md) identifies
the older oppositely oriented source tensor in the saved reference-backend
transcripts. This is a convention conversion, not a photon-anomaly result
supplied by Alda:2026nrd, which excludes WZW.

The [Bai digest](../../extracts/Bai:2025fvl/background-fields-wzw.md) establishes
the need for the full gauged WZW functional and its gauge counterterm.
This calculation converts the photon-only operator convention to the tensor
action. It does not import a background-vector model or match new resonance
anomaly coefficients. A UV/local photon coefficient remains an input in its
explicitly declared dual convention.

# Tensor calculation

Run the deterministic SymPy derivation:

```text
.venv/bin/python scripts/derive_odd_photon_phase_matching.py
```

Use $g_{\mu\nu}=\operatorname{diag}(1,-1,-1,-1)$ and
$\epsilon^{0123}=+1$, hence $\epsilon_{0123}=-1$ after metric lowering.
The original executable above is retained as a reference-convention backend;
the executable canonical conversion is linked in the audit. All momenta are incoming,
$\partial\to iq$. Define the source-basis tensors

$$
E=\epsilon_{\mu\nu\rho\sigma}
q_1^\mu\varepsilon_1^\nu q_2^\rho\varepsilon_2^\sigma,
\qquad R=\epsilon_{\mu\nu\rho\sigma}
\varepsilon_1^\mu\varepsilon_2^\nu q_1^\rho q_2^\sigma,
\qquad E=-R.
$$

Use the same tensor in the photon dual. Explicit metric lowering and both
photon Wick assignments give, in the canonical convention,

$$
\frac{\langle F_{\mu\nu}\widetilde F^{\mu\nu}\rangle}{E}
=-4,
\qquad
\frac{\langle\epsilon_{\mu\nu\rho\sigma}
\partial^\mu A^\nu\partial^\rho A^\sigma\rangle}{E}=-2.
$$

The canonically adapted WZW derivative coefficient is
$+\sqrt2N_c/(8\pi^2F)$; its positive dual coefficient remains
$\sqrt2N_c/(16\pi^2F)$. The canonical photon dictionary is

$$
\epsilon^{0123}=+1,\qquad \epsilon_{0123}=-1.
$$

The prior reference-backend calculation used a lower-index source tensor
with positive orientation and a separately oppositely oriented photon dual.
Their former $\eta_\gamma=-1$ is conversion history only. The canonical
convention uses one tensor; $c_\gamma$, $C_{a\gamma\gamma}$ and their
shift $c_\gamma'=c_\gamma-2N_c\operatorname{tr}(Q^2\kappa)$ remain unchanged.
The odd resonance coefficients reverse relative to that backend,
$c_i=-c_i^R$, $d_i=-d_i^R$, as documented in the audit.

# Common amplitude and LSZ phase

The canonical bounded gradients $G_j$ multiply $i e^2E$ in the source field basis.
Their normalized column scalars satisfy $G_j=G_j^R$, since both the
operator contraction and its epsilon normalization reverse; their canonical
coefficients satisfy $g_j=-g_j^R$. Convert the coefficients once, rather
than reversing both the normalized gradient and its coefficient.
The canonical signed SM scalar $F_P$ equals that scalar coefficient: its
$i e^2R F_P$ definition includes the opposite external-meson phase.
Using that rephased meson amplitude in an ALP LSZ sum requires the opposite
meson-overlap phase as well. Flipping only WZW changes physical interference.

For $S=\mathbb M_a^{-1}(\mathcal C)$ and $H=c_q^A+S$, the source-meson
replacement $\Phi_R/a=\sqrt2F H/f_a$ gives

$$
G_{\rm WZW/mix}=-\frac{2N_c\operatorname{tr}(Q^2H)}{4\pi^2f_a},
\qquad C_{a\gamma\gamma}=c_\gamma+2N_c\operatorname{tr}(Q^2H).
$$

The signed tree amplitude with this photon-operator normalization and the
retained leading-nonet WZW input is

$$
G_a=-\frac{C_{a\gamma\gamma}}{4\pi^2f_a}
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
record the original reference-backend tensor dictionary, zero-resonance meson LSZ limit, full original
eleven-coupling SM recovery without OPE, 241 photon-parent LSZ identities,
241 photon-parent Bose identities, two contact Ward identities, coefficient
rotation and a nonzero anomaly–resonance interference example. Flipping both
pieces preserves its squared amplitude; flipping only the anomaly changes
it by the derived nonzero interference difference.

The canonical conversion artifact verifies action and signed-interference
preservation; old transcripts are not silently relabeled as regenerated outputs.
This closes the shared signed-convention gate. The photon-operator normalization
is literature-backed; the leading nonet WZW treatment remains a stated
approximation, with independent enlarged-action local/anomaly matching still
uncomputed. The tree/one-ALP/bounded-action
limits, LO resolvent pole exclusions and perturbative ALP–meson overlap
requirement remain. No extra gauged-WZW matching, completed width, numerical
predictor, measurement rank or fit is inferred.
