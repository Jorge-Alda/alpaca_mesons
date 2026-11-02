---
title: NLO quadratic meson and meson–ALP Lagrangian
status: derived
created: 2026-11-02
updated: 2026-11-02
bibtex_keys: [Kaiser:1998ds, Alda:2026nrd]
provenance: analytical-derivation
---

# NLO quadratic terms

## Assumptions and inputs

The NLO basis is the displayed $\mathcal{L}^{(1)}$ and $\mathcal{L}^{(2)}$ in [the NLO Lagrangian record](lagrangian.md). The field and spurion expansions use [Appendix A of Alda:2026nrd](../extracts/Alda:2026nrd/appendix-a-expansion.md), with $\phi=\Phi/(2f_\pi)$, $\chi_0=2B_0m_q$, and $X_{TC}$ from [the topological source record](../extracts/Alda:2026nrd/lo-topological-action.md). The literature inputs are [@Kaiser:1998ds; @Alda:2026nrd].

Only terms of total field degree two are retained: $\Phi^2$ and $a\Phi$. Terms proportional to $a^2$, loop contributions, and the indicated $\mathcal{O}(N_c p^6)$ terms are omitted. Flavor matrices are arbitrary noncommuting matrices; traces are cyclic. Lorentz indices are contracted.

## Expansion used

The linear vielbein needed in products of two vielbeins is

$$
u_\mu^{(1)}=-\frac{\partial_\mu\Phi}{f_\pi}+2\frac{\partial_\mu a}{f_a}c_q^A.
$$

The relevant spurion pieces are

$$
\chi_+^{(0)}=2\chi_0,\qquad
\chi_+^{(\Phi^2)}=-\{\phi,\{\phi,\chi_0\}\},\qquad
\chi_+^{(a\Phi)}=2\frac{a}{f_a}\{\phi,\{c_q^P,\chi_0\}\},
$$

$$
\chi_-^{(1)}=2i\frac{a}{f_a}\{c_q^P,\chi_0\}-2i\{\phi,\chi_0\},\qquad
\chi_-^{(a\Phi)}=2\frac{a}{f_a}\{\phi,[c_q^S,\chi_0]\},
$$

and

$$
X_{TC}^{(1)}=-\frac{\operatorname{tr}\Phi}{f_\pi}-2c_G\frac{a}{f_a}.
$$

The $c_q^V$ term in $u_\mu$ starts at degree two and therefore cannot contribute to these quadratic products. The linear-$a$ pieces involving $c_q^S$ have degree one, while the $a\Phi$ piece of $\chi_-$ is multiplied by another field in the NLO operators considered here. Any single-trace commutators vanish by cyclicity. Thus neither $c_q^S$ nor $c_q^V$ occurs in the retained quadratic terms.

## Result: two pseudoscalar mesons

The $\Phi^2$ part is

$$
\begin{aligned}
\mathcal{L}_{\mathrm{NLO}}^{(\Phi^2)}={}&
L_5\frac{2\operatorname{tr}(\chi_0\partial_\mu\Phi\partial^\mu\Phi)}{f_\pi^2}
-L_8\frac{2[\operatorname{tr}(\Phi^2\chi_0^2)+\operatorname{tr}(\Phi\chi_0\Phi\chi_0)]}{f_\pi^2}
+\Lambda_1\frac{\operatorname{tr}(\partial_\mu\Phi)\operatorname{tr}(\partial^\mu\Phi)}{12}
-\Lambda_2\frac{\operatorname{tr}(\Phi\chi_0)\operatorname{tr}(\Phi)}{6}\\
&+L_4\frac{2\operatorname{tr}(\partial_\mu\Phi\partial^\mu\Phi)\operatorname{tr}(\chi_0)}{f_\pi^2}
-L_6\frac{4\operatorname{tr}(\Phi^2\chi_0)\operatorname{tr}(\chi_0)}{f_\pi^2}
-L_7\frac{4\operatorname{tr}(\Phi\chi_0)^2}{f_\pi^2}\\
&+L_{18}\frac{2\operatorname{tr}(\chi_0\partial_\mu\Phi)\operatorname{tr}(\partial^\mu\Phi)}{f_\pi^2}
+L_{25}\frac{4\operatorname{tr}(\Phi\chi_0^2)\operatorname{tr}(\Phi)}{f_\pi^2}.
\end{aligned}
$$

## Result: one meson and one ALP

The $a\Phi$ part is

$$
\begin{aligned}
\mathcal{L}_{\mathrm{NLO}}^{(a\Phi)}={}&
-L_5\frac{4\partial_\mu a}{f_af_\pi}\left[\operatorname{tr}(c_q^A\chi_0\partial^\mu\Phi)+\operatorname{tr}(c_q^A\partial^\mu\Phi\chi_0)\right]
+L_8\frac{4a}{f_af_\pi}\left[\operatorname{tr}(\Phi\chi_0^2c_q^P)+\operatorname{tr}(\Phi c_q^P\chi_0^2)+2\operatorname{tr}(\Phi\chi_0c_q^P\chi_0)\right]\\
&-\Lambda_1\frac{f_\pi\partial_\mu a}{3f_a}\operatorname{tr}(c_q^A)\operatorname{tr}(\partial^\mu\Phi)
+\Lambda_2\frac{af_\pi}{3f_a}\left[-c_G\operatorname{tr}(\Phi\chi_0)+\operatorname{tr}(c_q^P\chi_0)\operatorname{tr}(\Phi)\right]\\
&-L_4\frac{8\partial_\mu a}{f_af_\pi}\operatorname{tr}(c_q^A\partial^\mu\Phi)\operatorname{tr}(\chi_0)
+L_6\frac{8a}{f_af_\pi}\left[\operatorname{tr}(\Phi c_q^P\chi_0)+\operatorname{tr}(\Phi\chi_0c_q^P)\right]\operatorname{tr}(\chi_0)\\
&+L_7\frac{16a}{f_af_\pi}\operatorname{tr}(\Phi\chi_0)\operatorname{tr}(c_q^P\chi_0)
-L_{18}\frac{4\partial_\mu a}{f_af_\pi}\left[\operatorname{tr}(c_q^A\chi_0)\operatorname{tr}(\partial^\mu\Phi)+\operatorname{tr}(\chi_0\partial^\mu\Phi)\operatorname{tr}(c_q^A)\right]\\
&+L_{25}\frac{8a}{f_af_\pi}\left[c_G\operatorname{tr}(\Phi\chi_0^2)-\operatorname{tr}(c_q^P\chi_0^2)\operatorname{tr}(\Phi)\right].
\end{aligned}
$$

## Method and reproducibility

The SymPy derivation expands the Appendix A expressions with noncommuting flavor matrices, extracts the coefficients of the scalar markers for $\Phi^2$ and $a\Phi$ operator by operator, and prints the results in LaTeX. A minimal cyclic-rotation helper is used because the installed SymPy trace object does not collect cyclically equivalent words. Run the script with .venv/bin/python scripts/derive_nlo_quadratic_chipt.py.

The expressions above are the script output. Inputs are the NLO operator basis and the linked Alda:2026nrd conventions; no numerical parameters or phenomenological approximations enter.
