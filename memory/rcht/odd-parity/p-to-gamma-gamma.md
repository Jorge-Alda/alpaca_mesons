---
title: Real-photon pseudoscalar decays and coupling sensitivity
status: derived
created: 2026-10-13
updated: 2026-10-13
bibtex_keys: [Ruiz-Femenia:2003jdx]
provenance: analytical-derivation
---

# Real-photon pseudoscalar decays

## Inputs and assumptions

Use the action-derived [TFF amplitude](radiative-tff-amplitude.md), its
explicit $F_Q,F_H$ definitions, field normalization and common $M_V,F_V,F$.
The [source digest](../../extracts/Ruiz-Femenia:2003jdx/pseudoscalar-two-photon-form-factor.md)
provides the published pion amplitudes and width, Eqs. (34)–(35).
Assume tree level, leading electromagnetic order, neutral diagonal flavor
matrices, one vector multiplet and the recorded single-trace action.
The nonet extension remains an assumption; singlet/OZI matching is absent.
No ALP spurion is inserted and no neutral state is renormalized.
Natural units, masses in GeV, $\alpha_{\rm em}=e^2/(4\pi)$.

## Operator audit and method

The complete two-photon content within this action consists of the WZW
contact, VJP with both photon placements and photon–vector mixing, and VVP
with both vector assignments. Neutral $[T_P,Q]=0$ eliminates $c_4$;
$c_7,d_4$ vanish in the tensor contraction. Setting $s_1=s_2=0$ also removes
$c_6$. Covariant-derivative photon terms give no additional neutral vertex.
Thus the real-photon process is a limit of the existing TFF, rather than a
new operator calculation. This completeness statement does not exclude
additional singlet/OZI terms, higher-order local terms or loops in an enlarged
action.

The executable [derivation](../../../scripts/derive_p_to_gamma_gamma.py)
calls the existing tensor contractions, substitutes both virtualities to
zero, explicitly sums the two transverse polarization bases in the rest
frame, and integrates the massless two-body phase space with $1/2!$.
It then substitutes the five OPE relations and projects the coupling gradient
onto their homogeneous nullspace. Run:

```text
.venv/bin/python scripts/derive_p_to_gamma_gamma.py
```

## Results

For the signed amplitude convention in the TFF record,

$$
\begin{aligned}
F_P(0,0)={}&\frac{\sqrt2N_CF_Q}{8\pi^2F}
+\frac{8F_V}{FM_V^3}\left[F_Qm_P^2(c_1+c_2-c_5)+4c_3F_H\right]\\
&-\frac{4\sqrt2F_V^2}{FM_V^4}
\left[F_Qm_P^2(d_1-d_3)+4d_2F_H\right].
\end{aligned}
$$

The polarization sum of the squared Lorentz tensor is $m_P^4/2$, yielding

$$
\boxed{\Gamma(P\to\gamma\gamma)
=\frac{\pi\alpha_{\rm em}^2m_P^3}{4}|F_P(0,0)|^2.}
$$

After imposing the OPE relations, define $\Delta_H=F_H-2m_P^2F_Q$:

$$
F_P(0,0)=\frac{\sqrt2N_CF_Q}{8\pi^2F}
-\frac{\sqrt2FF_Qm_P^2}{2M_V^4}
+\Delta_H\left[\frac{32F_Vc_3}{FM_V^3}
-\frac{16\sqrt2F_V^2d_2}{FM_V^4}\right].
$$

For the pion, $F_Q=\sqrt2/3$, $F_H=2m_\pi^2F_Q$, $N_C=3$ give

$$
\Gamma(\pi^0\to\gamma\gamma)=
\frac{\alpha_{\rm em}^2m_\pi^3}{64\pi^3F^2}
\left(1-\frac{4\pi^2F^2m_\pi^2}{3M_V^4}\right)^2,
$$

reproducing the source result. Equal-spurion real-photon rows add no rank to
the five OPE rows. With independent spurions, every channel's projected row
is proportional to the same $c_3,d_2$ direction when $F_V/M_V$ is common;
the maximum combined rank is 6 (attained for nonzero $\Delta_H$). Widths
cannot determine the overall sign, and their local sensitivity can vanish
at zero amplitude. This is amplitude rank, not a statistical information
claim. Off-shell TFFs can provide a second direction, as recorded separately.

## Verification and limitations

The script passed the explicit width normalization, published pion correction,
equal-spurion nullspace cancellation and independent-spurion rank assertions.
Physical $\eta,\eta'$ require the field/decay-constant and spurion inputs
identified in [next steps](next-steps.md); their anomalous masses do not imply
$F_H=2m_P^2F_Q$. No PDG extraction, experimental comparison or fit was run.
