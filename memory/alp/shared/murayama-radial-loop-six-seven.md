---
title: Radial self-energy diagram 6.7 vertex product
status: derived
created: 2026-11-26
updated: 2026-11-26
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: analytical-derivation
---

# Radial self-energy diagram 6.7

## Inputs and assumptions

Use the polar-field conventions, incoming momenta and restricted graph set of the [estimate record](murayama-radial-up-only-estimate-result.md), and the compact vertex rules of the [generator-sum record](murayama-loop-generator-trace-products.md). The radial potential is the quartic potential of DiLuzio:2026poh, with $\mathcal L\supset-\lambda f_a\sigma^3$; conventions follow Alda:2026nrd. Units are GeV. The diagram contains two radial lines in the loop and two radial tree exchanges, both with denominator $m_\sigma^2-s_{+-}$.

## Derivation

The [script](../../../scripts/derive_radial_loop_six_seven.py) prints the cubic vertex $V(\sigma;x,y)$ for general pole legs $x,y$ in compact traces, and the three-radial vertex
$$
V_{\sigma^3}=-6\lambda f_a+\frac{6B_0f_\pi^2}{f_a^3}\langle\widehat S_{30}\rangle ,
$$
whose second term is the potential part with $r=3$. Diagram 6.7 is the product of two $\sigma^3$ vertices, the two tree exchanges and the two cubic vertices:
$$
\mathcal V_{6.7}=\frac{36f_a^2(\lambda-B_0f_\pi^2\langle\widehat S_{30}\rangle/f_a^4)^2 V(\sigma;a,\pi^0)V(\sigma;\pi^+,\pi^-)}{(m_\sigma^2-s_{+-})^2}.
$$
After inserting both cubic vertices and the zero angular components of the charged pions,
$$
\begin{aligned}
\mathcal V_{6.7}={}&-\frac{9B_0(\lambda-B_0f_\pi^2\langle\widehat S_{30}\rangle/f_a^4)^2}{f_a^2(m_\sigma^2-s_{+-})^2}\langle\widehat S_{10}\{\lambda_{\pi^+},\lambda_{\pi^-}\}\rangle\Big[2B_0 f_af_\pi\big(\xi_a\langle\widehat P_{11}\lambda_{\pi^0}\rangle+\xi_{\pi^0}\langle\widehat P_{11}\lambda_a\rangle\big)\\
&-B_0f_a^2\langle\widehat S_{10}\{\lambda_a,\lambda_{\pi^0}\}\rangle+4B_0f_\pi^2\xi_a\xi_{\pi^0}\langle\widehat S_{12}\rangle-4f_a^2\xi_a\xi_{\pi^0}\,k_a\cdot k_{\pi^0}\Big].
\end{aligned}
$$
The symmetry factor $1/2$ and the loop factor $1/(16\pi^2)$ are not included. The script prints the traces in full, with no local symbols standing for traces. This is a dimensional estimate; no loop integral is evaluated. Diagram 6.1 is in its [own record](murayama-radial-loop-six-one.md).

## Reproduction and verification

Run `PYTHONPATH=scripts .venv/bin/python scripts/derive_radial_loop_six_seven.py` and `PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_radial_loop_six_seven`. The exact checks, with general noncommuting sources and $\langle\widehat S_{30}\rangle$ kept, are: $V_{\sigma^3}$ and $V(\sigma;x,y)$ equal the direct expansion of the Lagrangian of the [generator-sum record](murayama-loop-generator-trace-products.md); $V_{\sigma^3}$ equals the closed form above; the product of 6.7 equals the direct expansion with both radial lines inserted; and it equals the second display above.
