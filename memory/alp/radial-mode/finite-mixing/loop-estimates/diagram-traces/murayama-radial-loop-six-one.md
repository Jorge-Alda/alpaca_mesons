---
title: Radial bubble diagram 6.1 vertex product
status: derived
created: 2026-11-26
updated: 2026-11-26
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: analytical-derivation
---

# Radial bubble diagram 6.1

Use the conventions and inputs of [diagram 6.7](murayama-radial-loop-six-seven.md), including the quartic radial potential and $V_{\sigma^3}=-6\lambda f_a+6B_0f_\pi^2\langle\widehat S_{30}\rangle/f_a^3$, and the compact vertex rules of the [generator-sum record](murayama-loop-generator-trace-products.md). The [script](../../../../../../scripts/derive_radial_loop_six_seven.py) prints $V(\sigma,\sigma;x,y)$ for general pole legs $x,y$ and the product of diagram 6.1: the $\sigma^3$ vertex with the charged-pion tree exchange, times $V(\sigma,\sigma;a,\pi^0)$. With incoming momenta,

$$
\begin{aligned}
\mathcal V_{6.1}={}&\frac{3B_0(\lambda-B_0f_\pi^2\langle\widehat S_{30}\rangle/f_a^4)}{f_a^4(m_\sigma^2-s_{+-})}\langle\widehat S_{10}\{\lambda_{\pi^+},\lambda_{\pi^-}\}\rangle\Big[2B_0 f_af_\pi\big(\xi_a\langle\widehat P_{21}\lambda_{\pi^0}\rangle+\xi_{\pi^0}\langle\widehat P_{21}\lambda_a\rangle\big)\\
&-B_0f_a^2\langle\widehat S_{20}\{\lambda_a,\lambda_{\pi^0}\}\rangle+4B_0f_\pi^2\xi_a\xi_{\pi^0}\langle\widehat S_{22}\rangle-2f_a^2\xi_a\xi_{\pi^0}\,k_a\cdot k_{\pi^0}\Big].
\end{aligned}
$$

The graph symmetry factor $1/2$ and loop factor $1/(16\pi^2)$ are not included. The expression holds for general dressed sources under the stated potential and pole assumptions. In the up-only model $\widehat S_{2s}=\widehat P_{2s}=0$, and only the last term of the bracket remains. No loop integral is evaluated.

Run `PYTHONPATH=scripts .venv/bin/python scripts/derive_radial_loop_six_seven.py` and `PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_radial_loop_six_seven`. The test checks exactly, with general noncommuting sources and $\langle\widehat S_{30}\rangle$ kept, that $V(\sigma,\sigma;x,y)$ and the product of 6.1 equal the direct expansion of the Lagrangian of the [generator-sum record](murayama-loop-generator-trace-products.md); that the product equals the display above; and that setting $\widehat S_{2s}=\widehat P_{2s}=0$ leaves only the kinetic term of the bracket.
