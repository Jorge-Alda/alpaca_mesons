---
title: Full NLO Chiral perturbation Lagrangian
status: derived
created: 2026-11-02
updated: 2026-11-02
bibtex_keys: [Kaiser:1998ds,Alda:2026nrd]
provenance: analytical-derivation
---

## Notation

The relations needed to rewrite the NLO Lagrangian in the CCWZ formalism in the notation of Alda:2026nrd are the following

$$
\bm{U} = \bm{u}^2
$$

$$
D_\mu\bm{U} = -i \bm{u}\bm{u_\mu}\bm{u}
$$

$$
(D_\mu\bm{U})^\dagger = i \bm{u}^\dagger\bm{u_\mu}\bm{u}^\dagger
$$

$$
\bm{\chi}=\frac{1}{2}\bm{u}\bm{\chi_+}\bm{u} + \frac{1}{2}\bm{u} \bm{\chi_-}\bm{u}
$$

$$
\bm{\chi}^\dagger=\frac{1}{2}\bm{u}^\dagger\bm{\chi_+}\bm{u}^\dagger - \frac{1}{2}\bm{u}^\dagger \bm{\chi_-}\bm{u}^\dagger
$$

$$
\psi+\theta = -X_{TC}
$$

$$
D_\mu\psi = -\langle \bm{u_\mu}\rangle
$$

## NLO Lagrangian

Kaiser:1998ds organizes the NLO terms according to their power counting in the combined $p$, $m_q$ and $1/N_c$ expansion, with expansion parameter $\delta$,

$$\partial_\mu \sim \mathcal{O}(\sqrt{\delta}),\qquad m_q\sim\mathcal{O}(\delta),\qquad 1/N_c\sim\mathcal{O}(\delta)$$

The $\mathcal{O}(\delta^0)$ Lagrangian is already covered in Alda:2026nrd. The next term includes contributions of $\mathcal{O}(N_c p^4)$, as well as OZI-violating terms of $\mathcal{O}(p^2)$,

$$
\mathcal{L}^{(1)} = L_5 \langle\bm{u_\mu}\bm{u^\mu}\bm{\chi_+}\rangle + \frac{L_8}{2}\langle\bm{\chi_+}^2 + \bm{\chi_-}^2\rangle +\frac{f_\pi^2}{12}\Lambda_1 \langle\bm{u_\mu}\rangle\langle\bm{u^\mu}\rangle +i\frac{f_\pi^2}{12}\Lambda_2 X_{TC} \langle\bm{\chi_-}\rangle.
$$

while the rest of the $\mathcal{O}(p^4)$ terms belong to $\mathcal{O}(\delta^2)$,

$$
\mathcal{L}^{(2)} =L_4 \langle\bm{u_\mu}\bm{u^\mu}\rangle\langle\bm{\chi_+}\rangle+L_6\langle\bm{\chi_+}\rangle^2+L_7\langle\bm{\chi_-}\rangle^2+L_{18}\langle\bm{u_\mu}\rangle\langle\bm{u^\mu}\bm{\chi_+}\rangle-iL_{25} X_{TC}\langle\bm{\chi_+}\bm{\chi_-}\rangle + \mathcal{O}(N_c p^6).
$$