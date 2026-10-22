---
title: Widths for an ALP to a photon and two pseudoscalars
status: derived
created: 2026-10-22
updated: 2026-10-22
bibtex_keys: [Alda:2026nrd, Ecker:1988te, Ecker:1989yg, Ruiz-Femenia:2003jdx]
provenance: analytical-derivation
---

# Assumptions and inputs

This record reduces the odd tree resonance amplitude in [the $a\to\gamma PP$
amplitude record](alp-photon-pp-amplitude.md) to an invariant width. It uses
natural units, metric $(+---)$, $\epsilon^{0123}=+1$, a massless real photon,
and the Alda:2026nrd pseudoscalar normalization and source conventions
[@Alda:2026nrd]. The channels below are those selected by the flavor traces
in the parent amplitude; the photon is flavor diagonal and does not change
the strangeness selection rule.

The width derived here is $\Gamma_{\rm odd}$ for the retained odd resonance
amplitude. It is not the complete physical width: local anomalous terms and
the even scalar/vector-source radiative amplitude are not included. Internal
resonance denominators and all exchange contributions are kept coherently in
the complex coefficient $\mathcal A_{\gamma PP}$; any finite-width prescription
must be specified there. The executable symbolic derivation is
[`scripts/derive_alp_photon_two_pseudoscalar_width.py`](../../../scripts/derive_alp_photon_two_pseudoscalar_width.py):

```text
.venv/bin/python scripts/derive_alp_photon_two_pseudoscalar_width.py
```

# Polarization sum

For
$$
\mathcal M=\epsilon_{\mu\nu\rho\sigma}
\varepsilon_\gamma^\mu p_1^\nu p_2^\rho k^\sigma\,
\mathcal A_{\gamma PP}(s,t,u),
$$
the tensor multiplying $\varepsilon_\gamma^\mu$ is transverse to $k$.
Consequently the gauge-dependent terms in the real-photon polarization sum
drop out, and the Levi-Civita contraction gives
$$
\sum_{\rm pol}|\mathcal M|^2
=\Delta_G(s,t,u)\,|\mathcal A_{\gamma PP}(s,t,u)|^2,
\qquad
\Delta_G=\det\begin{pmatrix}
0 & k\!\cdot p_1 & k\!\cdot p_2\\
k\!\cdot p_1 & m_1^2 & p_1\!\cdot p_2\\
k\!\cdot p_2 & p_1\!\cdot p_2 & m_2^2
\end{pmatrix}.
$$
With $s=(p_1+p_2)^2$, $t=(k+p_2)^2$, and
$u=(k+p_1)^2=m_a^2+m_1^2+m_2^2-s-t$, this is equivalently
$$
\Delta_G(s,t)=\frac14\left[-m_2^2(m_a^2+m_2^2-s-t)^2
+(m_a^2+m_2^2-s-t)(t-m_2^2)(s-m_1^2-m_2^2)
-m_1^2(t-m_2^2)^2\right].
$$
The Källén function used in the limits below is
$$
\lambda(x,y,z)=x^2+y^2+z^2-2xy-2xz-2yz.
$$

# Dalitz integral and channel widths

The differential and integrated widths are
$$
\frac{d^2\Gamma_{\rm odd}}{ds\,dt}
=\frac{\Delta_G(s,t)|\mathcal A_{\gamma PP}(s,t,u)|^2}
{256\pi^3m_a^3 S_{PP}},
\qquad
\Gamma_{\rm odd}(a\to\gamma P_1P_2)
=\frac{1}{256\pi^3m_a^3 S_{PP}}
\int_{(m_1+m_2)^2}^{m_a^2}ds
\int_{t_-(s)}^{t_+(s)}dt\,
\Delta_G(s,t)|\mathcal A_{\gamma PP}(s,t,u)|^2,
$$
where $S_{PP}=2!$ for identical pseudoscalars and $1$ otherwise, and
$$
t_\pm(s)=m_2^2+
\frac{(m_a^2-s)(s+m_2^2-m_1^2)
\pm(m_a^2-s)\sqrt{\lambda(s,m_1^2,m_2^2)}}{2s}.
$$
The integral is zero if $m_a<m_1+m_2$. Every channel selected by the flavor
rules in this sector has distinguishable pseudoscalars, hence $S_{PP}=1$:
$$
\begin{aligned}
\Gamma_{\rm odd}(a\to\gamma\pi^+\pi^-),\quad
\Gamma_{\rm odd}(a\to\gamma K^+K^-),\quad
\Gamma_{\rm odd}(a\to\gamma K^0\bar K^0)
&\quad\text{for diagonal }c_q^A\text{ and/or }c_q^P,\\[2pt]
\Gamma_{\rm odd}(a\to\gamma K^0\pi^0),\quad
\Gamma_{\rm odd}(a\to\gamma K^0\eta),\quad
\Gamma_{\rm odd}(a\to\gamma K^0\eta'),\quad
\Gamma_{\rm odd}(a\to\gamma K^+\pi^-)
&\quad\text{for }c_{ds},\\[2pt]
\Gamma_{\rm odd}(a\to\gamma\bar K^0\pi^0),\quad
\Gamma_{\rm odd}(a\to\gamma\bar K^0\eta),\quad
\Gamma_{\rm odd}(a\to\gamma\bar K^0\eta'),\quad
\Gamma_{\rm odd}(a\to\gamma K^-\pi^+)
&\quad\text{for }c_{sd}.
\end{aligned}
$$
For each listed mode, insert its masses and the corresponding coherent
$\mathcal A_{\gamma PP}$ from the amplitude record into the same integral.
If both Hermitian-conjugate flavor entries are present, use their coherent
sum in each final state's amplitude before taking its absolute square.

# Verification and limitations

The executable derivation constructs the Gram determinant from the external
momenta, eliminates $u$ in favor of $s$ and $t$, derives the physical $t$
limits, and obtains the three-body phase-space normalization. It checks the
massless-photon limit of the general invariant phase-space expression. The
result is a tree-level width formula for the truncation stated above; it does not include loops,
rescattering, additional odd operators, local anomalous contacts, or the even
radiative amplitude.
