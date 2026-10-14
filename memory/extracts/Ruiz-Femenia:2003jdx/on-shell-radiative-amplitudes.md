---
title: On-shell vector-to-pseudoscalar-photon amplitudes and widths
status: literature
created: 2026-09-30
updated: 2026-10-15
bibtex_keys: [Ruiz-Femenia:2003jdx]
provenance: literature
---

Displayed tensors and odd couplings use $\epsilon^{0123}=+1$ and
$\epsilon_{0123}=-1$ for $(+---)$. The canonical coefficients are
$c_i=-c_i^R$, $d_i=-d_i^R$ relative to the opposite source convention.
The [conversion audit](../../rcht/odd-parity/levi-civita-convention.md)
derives the action-preserving map; archived source text remains literal.


# Exact published formulas

The local auditable source is `memory/sources/Ruiz-Femenia:2003jdx/vvp.tex`; source
line numbers below refer to that file. The antisymmetric tensor convention is
$\epsilon_{0123}^R=+1$ (lines 352--354); the displayed formulas have been
converted to the canonical tensor as explained above. The octet matrices are normalized as
`Phi=lambda^a phi^a/sqrt(2)` and `V=lambda^a V^a/sqrt(2)` (Eqs. (4),(8), lines
206--219 and 291--304). Ideal mixing is
$$
|\omega\rangle=\sqrt{2/3}|\omega_1\rangle+\sqrt{1/3}|\omega_8\rangle,
\qquad
|\phi\rangle=-\sqrt{1/3}|\omega_1\rangle+\sqrt{2/3}|\omega_8\rangle
$$
(lines 784--797). The singlet is added as `omega_1/sqrt(3)` to the octet
matrix. These mixing statements are needed for any singlet-containing channel.

The even-parity tensor normalization is
$$
\mathcal L_2(V)=\frac{F_V}{2\sqrt2}\langle V_{\mu\nu}f_+^{\mu\nu}\rangle
+\frac{iG_V}{\sqrt2}\langle V_{\mu\nu}u^\mu u^\nu\rangle,
$$
with `f_+` defined from the left/right external field strengths (Eq. (7),
lines 278--290). The odd terms use
$$
\mathcal L_{VJP}=\sum_{a=1}^7\frac{c_a}{M_V}\mathcal O^a_{VJP},
\qquad
\mathcal L_{VVP}=\sum_{a=1}^4d_a\mathcal O^a_{VVP}
$$
(Eq. (12), lines 415--430). Thus $c_a$ are dimensionless while $d_a$ carry
the normalization of the displayed operators. The source does not print a
separate $V$--photon mixing equation or an explicit `v_mu=... A_mu`
substitution; the factors of $e$ in the following published amplitudes are
the operative photon normalization.

For $\omega(q)\to\pi(k)\gamma$, the uneliminated direct and $\rho$-exchange
amplitudes are (Eq. (27), source lines 810--831):
$$
i\mathcal M^{\rm direct}=i\epsilon_{\alpha\beta\rho\sigma}
\epsilon_\omega^\alpha\epsilon_\gamma^\beta q^\rho k^\sigma
\frac{2\sqrt2e}{M_\omega M_VF}
\left[(c_2-c_1+c_5-2c_6)M_\omega^2
+(c_1+c_2+8c_3-c_5)m_\pi^2\right],
$$
$$
i\mathcal M^\rho=-i\epsilon_{\alpha\beta\rho\sigma}
\epsilon_\omega^\alpha\epsilon_\gamma^\beta q^\rho k^\sigma
\frac{4e}{M_V^2M_\omega}\frac{F_V}{F}
\left[d_3M_\omega^2+(d_1+8d_2-d_3)m_\pi^2\right].
$$
After the OPE constraints, the published on-shell result is (Eq. (28),
source lines 832--843):
$$
i\mathcal M_{\omega\to\pi\gamma}=-i\epsilon_{\alpha\beta\rho\sigma}
\epsilon_\omega^\alpha\epsilon_\gamma^\beta q^\rho k^\sigma
\frac{e}{F_V}\left[\frac{N_C}{8\pi^2}\frac{M_\omega}{F}
-\frac F2\frac{M_\omega}{M_V^2}\left(1+\frac{m_\pi^2}{M_\omega^2}\right)\right].
$$
Using $F_V=\sqrt2F$, its width is (Eq. (29), source lines 856--868):
$$
\Gamma(\omega\to\pi\gamma)=\frac{\alpha}{192}M_\omega
\left(1-\frac{m_\pi^2}{M_\omega^2}\right)^3
\left[\frac{N_C}{4\pi^2}\frac{M_\omega^2}{F^2}
-\frac{M_\omega^2}{M_V^2}\left(1+\frac{m_\pi^2}{M_\omega^2}\right)\right]^2.
$$

## Scope limit: this is the only on-shell vector radiative decay in the source

The paper contains exactly three phenomenology subsections: $\omega\to\pi\gamma$
(line 778), $\omega\to\pi^+\pi^-\pi^0$ (line 891), and $\pi\to\gamma\gamma$
(line 994). **No general $V\to P\gamma$ formula is given**, and no other
radiative vector decay is treated.

On extending to other channels the authors write (source lines 880--889):

> The extension of our analysis to other decay channels (e.g.
> $K^*\to K\gamma$, $\phi\to\eta\gamma$) requires that exact
> SU(3)$_{\mathrm{V}}$ symmetry is left aside in order not to lose the
> predictive power shown in $\omega\to\pi\gamma$. This study would require to
> consider the OPE expansion in the asymptotic regime keeping distinct masses
> for each quark flavour, a rather non trivial task.

Two readings must be distinguished, because an earlier version of this extract
conflated them:

1. **For the numerical OPE-eliminated prediction (Eq. (28) and its width
   Eq. (29)) the extension is genuinely non-trivial.** It requires abandoning
   exact SU(3)$_{\mathrm V}$ and redoing the short-distance matching with
   distinct quark masses. The two channels named above are given as *examples*
   ("e.g."), **not** as an exhaustive list of what is excluded. Any other channel,
   including $\rho\to\pi^0\gamma$ and $\omega\to\eta\gamma$, falls under the
   same statement.
2. **For the uneliminated couplings $c_a,d_a$ the requirement is weaker.** The
   fit uses the uneliminated amplitude, so no new OPE matching is required to
   *evaluate* a channel; what is required is a derivation of that channel's
   operator-level amplitude, including its flavor trace and singlet/octet
   content. That is genuine EFT work rather than transcription, but it is
   tractable and is **not performed in this source**.

A claim that the displayed coupling combinations simply "generalize after
replacing the flavor trace, external masses, vector pole, and mixing factors" is
an inference, not something this source states, and it is unsafe: for an ALP
source $\mathbb M_a^{-1}(C)$ is a general Hermitian matrix rather than an
SU(3)-diagonal spurion, so the singlet/octet decomposition that would justify
such a recipe does not hold automatically. Treat any extension to a new
radiative channel as an open derivation.
