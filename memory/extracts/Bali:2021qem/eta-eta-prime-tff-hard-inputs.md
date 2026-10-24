---
title: Eta and eta-prime TFF asymptotic inputs from lattice decay constants
status: literature
created: 2026-10-24
updated: 2026-10-24
bibtex_keys: [Bali:2021qem]
provenance: literature
---

# Axial-current normalization and scheme

The decay constants are defined by
$\langle 0|A_\mu^a|P(p)\rangle=iF_P^a p_\mu$ in Eq. (3). In the octet and
singlet basis, there are four independent constants. The singlet current
and $F_P^0$ depend on the QCD renormalization scale. The paper reports
$F_P^8$ and $F_P^0$ in the $\overline{\mathrm{MS}}$ scheme for $N_f=3$, and
describes evolution and matching to $N_f=4$ in Sec. 6.1 before evaluating
the hard TFF input [@Bali:2021qem].

The reported angle representation at the physical point gives
$F^8=115.0(2.8)\,\mathrm{MeV}$, $\theta_8=-25.8(2.3)^\circ$,
$\theta_0=-8.1(1.8)^\circ$, and
$F^0(\mu=2\,\mathrm{GeV})=100.1(3.0)\,\mathrm{MeV}$ in $N_f=3$.
The definitions are

$$
F^8_\eta=F^8\cos\theta_8,\quad
F^8_{\eta'}=F^8\sin\theta_8,\quad
F^0_\eta=-F^0\sin\theta_0,\quad
F^0_{\eta'}=F^0\cos\theta_0.
$$

These are the source's lattice inputs, not inputs derived from the ALP
no-mixing source or from Alda:2026nrd meson generators.

# Single-virtual hard limits

Equation (5.80) gives the asymptotic single-virtual coefficient in the
source's TFF scalar convention:

$$
\lim_{Q^2\to\infty}Q^2F_{\gamma\gamma^*\to P}(Q^2)
=\frac{2}{\sqrt3}\left[F_P^8+2\sqrt2F_P^0(\mu=\infty,N_f=4)\right],
\qquad P=\eta,\eta'.
$$

After the stated singlet evolution and flavor-threshold matching, the paper
reports $160.5(10.0)\,\mathrm{MeV}$ for $\eta$ and
$230.5(10.1)\,\mathrm{MeV}$ for $\eta'$. It stresses that neglecting the
scale evolution overestimates the $\eta'$ coefficient because its octet
and singlet contributions have opposite signs. The repository's fixed-phase
scalar uses the opposite overall sign; apply the common tensor-convention
conversion when inserting these magnitudes in canonical formulas.

The paper discusses singly virtual hard limits and gives light-cone sum-rule
predictions at finite virtuality. It does not by itself establish a
source-matched doubly-virtual, all-ratio TFF representation for the RChT
odd coefficients. Finite-virtuality predictions also assume specific
light-cone distribution amplitudes and a two-gluon input; those assumptions
must be retained if using the plotted finite-range predictions.

Source locations: Eq. (3), Sec. 6.1, Eq. (5.80), and the paragraph following
Eq. (5.80), in `memory/sources/Bali:2021snd/etaetapr.tex`. The INSPIRE record
is stored unmodified in `bibliography.bib`.
