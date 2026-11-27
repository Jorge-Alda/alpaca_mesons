---
title: Sign convention for external neutral-meson legs at finite ALP-meson mixing
status: mixed
created: 2026-11-14
updated: 2026-11-14
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Sign convention for external neutral-meson legs at finite ALP-meson mixing

## Scope

An external neutral meson pole $j$ in a decay amplitude enters through a leg $L_j$, a flavour matrix in the unmixed generator basis $U_k$. This note fixes its sign and the value used where the general formula is singular. It applies to every amplitude that uses such legs (the WZW $a\to\gamma P_iP_j$ record, and any future amplitude with neutral external mesons). It supersedes the sentence "a neutral daughter leg tends to its unmixed generator" in [the WZW record](../../../../photon-two-pseudoscalar/finite-mixing/alp-lo-finite-mixing-wzw-pp-amplitude.md), which holds only when $\operatorname{sign}(\operatorname{tr}(U_j\mathcal C)D_j)=+1$.

## Notation

- $\epsilon=f_\pi/f_a$, $\mathcal C=B_0\widehat P_1-\frac{m_0^2}{N_f}\widehat\Theta_1\mathbf 1$ (see [the exact pole-equation record](../../../../meson-mixing/finite-mixing/spectrum-and-poles/pole-equations/alp-lo-exact-pole-equations.md)).
- Unmixed masses $\mu_k$ and generators $U_k$ of the neutral meson states, with $\operatorname{tr}(U_kU_l)=2\delta_{kl}$.
- $t_k=\operatorname{tr}(U_k\mathcal C)$.
- $s_j$ the physical squared pole mass, $m_a^2$ the physical ALP pole.
- $D_j=\mathbb M^{\rm can}_{aa}-s_j-\epsilon^2\sum_{k\neq j}\dfrac{t_k^2}{\mu_k-s_j}$, with $\mathbb M^{\rm can}_{aa}=m_a^2+\epsilon^2\sum_k\dfrac{t_k^2}{\mu_k-m_a^2}$.
- The exact meson pole equation gives $(\mu_j-s_j)D_j=\epsilon^2t_j^2$ (Step 5 of [the exact pole-equation record](../../../../meson-mixing/finite-mixing/spectrum-and-poles/pole-equations/alp-lo-exact-pole-equations.md)).

## 1. Analytic structure of the sign (two-level case)

For one neutral meson coupled to the ALP, the canonical mass matrix is $\begin{pmatrix}\mu&g\\g&m_a^2\end{pmatrix}$ with $g=-\epsilon t$. Put $\delta=m_a^2-\mu$. The eigenvalues are $s_\pm=\mu+\delta/2\pm R$, $R=\sqrt{\delta^2/4+g^2}$. Then $s_+>\mu$ and $s_-<\mu$ for every $\delta$. The eigenvector ratio $r_\pm=x_\pm/y_\pm=\epsilon t/(\mu-s_\pm)$ therefore has a fixed sign on each branch:
$$
\operatorname{sign}(r_+)=-\operatorname{sign}(t),\qquad \operatorname{sign}(r_-)=+\operatorname{sign}(t),\qquad\text{for all }\delta.
$$
The meson weight $r^2/(1+r^2)$ selects the branch with the smaller $|\mu-s|$. That branch is $+$ for $\delta<0$ and $-$ for $\delta>0$. So a pole labelled by largest overlap changes branch at $\delta=0$, i.e. when the ALP mass equals the unmixed meson mass, and its sign flips there. No eigenvector passes through zero along a branch.

**Rule.** Signs of neutral external legs are fixed per branch, not per overlap label. A branch keeps its sign for every ALP mass and every $\epsilon$. A label switch is not a physical sign change of the state.

## 2. Finite-mixing leg

The leg of pole $j$ in the boxed WZW record, written so that no difference $\mu_j-s_j$ appears, is
$$
L_j=\operatorname{sign}(t_j)\;
\frac{D_jU_j+\epsilon^2t_j\sum_{k\neq j}\dfrac{t_k}{\mu_k-s_j}\,U_k}
{\sqrt{D_j^2+\epsilon^2t_j^2+\sum_{k\neq j}\Big(\dfrac{\epsilon^2t_jt_k}{\mu_k-s_j}\Big)^2}} ,
$$
with $\operatorname{sign}(t_j)$ the sign of the real number $t_j$. This is exact for every $\epsilon$. Its weak-mixing limit is
$$
L_j\to\operatorname{sign}\big(t_j(\mu_j-s_j)\big)\,U_j=\operatorname{sign}(t_jD_j)\,U_j .
$$
The implementation is `neutral_pole_leg` in [`scripts/implementation/alp_photon_pp_wzw_finite_mixing.py`](../../../../../../scripts/implementation/alp_photon_pp_wzw_finite_mixing.py).

## 3. The singular case $t_j=0$

When $t_j=0$:

- The pole equation gives $(\mu_j-s_j)D_j=0$. The pole is exactly unmixed, $s_j=\mu_j$, with no coupling to the ALP.
- The boxed form $\frac{2\epsilon X_j}{\sqrt{1+2\epsilon^2\operatorname{tr}X_j^2}}$ is $0/0$: the numerator $t_jU_j/(2(\mu_j-s_j))$ and the normalization both involve $t_j/(\mu_j-s_j)$, whose limit depends on the path.
- The finite form above is discontinuous: $L_j\to\operatorname{sign}(t_j)\operatorname{sign}(D_j)U_j$ as $t_j\to0$, with magnitude $1$ on both sides.

**Convention.** At $t_j=0$ the leg is
$$
L_j=\sigma_j\,U_j,\qquad \sigma_j=\operatorname{sign}(D_j)\ \text{with}\ \operatorname{sign}(0)\equiv+1,
$$
where $D_j$ is evaluated at $s_j=\mu_j$. This is a stated convention: the discontinuous value is kept, not the limit. The magnitude $1$ is exact. It is not a threshold, and no exception is raised.

The value of $\sigma_j$ is a convention of the external state. It must be the same in every diagram that uses this state. Because the sign of $t_j$ is not defined at $t_j=0$, $\sigma_j$ may differ from the value obtained by approaching $t_j=0$ from one side. That is allowed. Only the relative sign between diagrams sharing the leg matters, and that is unaffected.

## 4. What is convention-free

- **Internal exchanges.** Propagators $D^{-1}(s)$ and residues $v^{(n)}v^{(n)T}$ are invariant under $v^{(n)}\to-v^{(n)}$. They carry no sign ambiguity.
- **Flavour states.** Legs for $K^0$, $\bar K^0$, $K_S$, $K_L$ built from the real kaon poles are invariant under the sign of each pole's leg, because $\operatorname{tr}(\lambda^\dagger d_k)\,L_k/\sqrt{\operatorname{tr}(d_k^2)/2}$ is even in $L_k$. Their relative phases $p,q$ are fixed by $\epsilon_K$ and are physical.
- **Rates.** A single external leg enters $|\mathcal A|^2$ only up to a global sign. Relative signs matter only between diagrams that share all external legs but differ in internal structure, or between different flavour combinations of the same particle.

## 5. Rules for implementation

1. Compute the external legs once per external state, and reuse them for every diagram and every resonance term.
2. Use the finite form of section 2 for $t_j\neq0$.
3. Use the value of section 3 for $t_j=0$, with $\operatorname{sign}(0)=+1$.
4. Do not use overlap labels to pick the sign of an external leg. The sign is a property of the branch, as in section 1. Where the implementation labels poles by overlap, it must be checked that a label switch does not flip the sign of a leg used in several diagrams.
5. Keep internal exchanges as full propagators or residues, never as sign-fixed legs.

## Validity

Derived for the leading-order neutral system with stable simple poles. The two-level case is exact. For several coupled mesons, section 1 relies on the interlacing of the poles with the unmixed masses ([arrow-matrix record](../../../../meson-mixing/finite-mixing/pole-limits/alp-lo-arrow-matrix-spectrum.md)); the switch condition for the overlap label has not been derived in closed form. Degenerate poles and the $K^0$ doublet ([doublet record](../../../../meson-mixing/finite-mixing/spectrum-and-poles/pole-equations/alp-lo-k0-doublet-exact-pole.md)) need the corresponding subspace treatment and are not covered here.
