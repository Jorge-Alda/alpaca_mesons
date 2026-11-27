---
title: Tree-level partial-wave unitarity of ALP scattering in the only-up model
status: derived
created: 2026-11-21
updated: 2026-11-21
bibtex_keys: [Murayama:2026ioh, DiLuzio:2026poh, Lee:1977eg, Alda:2026nrd]
provenance: analytical-derivation
---

# Tree-level partial-wave unitarity of ALP scattering in the only-up model

## Inputs

- Model, charges $(n_u,n_d,n_s)=(1,0,0)$ and the choice of the quartic coupling: [Murayama:2026ioh model extract](../../../extracts/Murayama:2026ioh/model-definition.md) and DiLuzio:2026poh [@Murayama:2026ioh; @DiLuzio:2026poh].
- Polar field $\Phi=(f_a+\sigma)e^{ia/f_a}/\sqrt2$, $f_a=f_\varphi$, and the quark-mass source: [radial and phase coefficient record](alp-murayama-radial-phase-mass-coefficients.md).
- Potential at $U=\mathbf 1$, vacuum condition $m_\varphi^2=B_0m_uf_\pi^2/f_a^2-\lambda f_a^2$, $m_\sigma^2=B_0m_uf_\pi^2/f_a^2+2\lambda f_a^2$, Cartesian vertex $\partial^3V/\partial s\,\partial p^2=2\lambda f_a$, $\sigma P_iP_j$ couplings and the connection-free source $\widehat M=\operatorname{diag}[m_u(1+\sigma/f_a)e^{-ia/f_a},m_d,m_s]$: [radial-mode decay record](finite-mixing/murayama-sigma-two-pseudoscalar-decays.md) (Cartesian section).
- Equality of on-shell polar and Cartesian amplitudes: [polar–Cartesian record](../../decay-amplitudes/vertices-and-kinematics/tree-vertices/finite-mixing/alp-polar-cartesian-amplitude-equivalence.md).
- Partial-wave normalization $(a_J)_{fi}=\frac{N_iN_f}{32\pi}\int_{-1}^1d\cos\theta\,P_J\,\mathcal M_{fi}$, $N=1/\sqrt2$ for identical particles; primary condition $|\operatorname{Re}a_0|\le\frac12$, secondary $|a_0|\le1$; meson vertex formula $\mathcal V_{a^nP_{i_1}\cdots P_{i_k}}$: [LO partial-wave unitarity record](../../decay-amplitudes/vertices-and-kinematics/kinematics-and-validity/finite-mixing/alp-lo-partial-wave-unitarity.md), generalizing Lee:1977eg [@Lee:1977eg].
- Conventions of Alda:2026nrd for $U=\exp(i\Phi/f_\pi)$, $\langle\lambda_i\lambda_j\rangle=2\delta_{ij}$ [@Alda:2026nrd].

Natural units; masses, $B_0$, $f_\pi$, $f_a$ in GeV; $\lambda$ dimensionless. Feynman rule $i\times$vertex with all momenta incoming, $i\mathcal M=$ sum of diagrams, so a scalar of mass $m$ exchanged in the channel $x$ contributes $\mathcal V_1\mathcal V_2/(m^2-x)$. The combination $B_0m_uf_\pi^2/f_a^2$ is written out: the symbol $m_0$ is the topological mass in this repository and is not used for it.

## Assumptions and approximations

1. **Tree level, LO chiral Lagrangian,** $\theta=0$, CP-conserving vacuum, $m_\varphi^2>0$, $\lambda>0$.
2. **ALP–meson mass mixing neglected (approximation, Parts A–E).** The external ALP is the unmixed field $a$ (polar) or $p$ (Cartesian), with mass squared $\partial^2V/\partial a^2|_0=B_0m_uf_\pi^2/f_a^2$. Validity: the mixing entry $M^2_{3a}=B_0m_uf_\pi/f_a$ (from $\frac{B_0f_\pi}{f_a}a\langle\widehat P_{01}\Phi\rangle$, analogously for $\eta_8,\eta_0$) must satisfy $|M^2_{ia}|\ll|m_a^2-m_{P_i}^2|$.
3. **High-energy kinematics:** $4m^2\ll s$ wherever the phase-space factor is set to one; results quoted for $s\to\infty$ are at fixed angle or after the exact S-wave projection.
4. **Coupled channels with $\sigma$ (Part C):** $\sigma$ is treated as an asymptotic state although $\sigma\to aa$ is open when $m_\sigma>2m_a$. As in Lee:1977eg, only the contact terms at $s\gg m_\sigma^2$ are used; at finite $s$ the $u$-channel $p$ pole of $\sigma p\to\sigma p$ lies in the physical region and the wave acquires an imaginary part.
5. **Part E** keeps only the constant $aP_iP_jP_k$ contact; the growing four-meson amplitudes of the coupled meson system are not included.

## Method

1. **Vacuum.** From $V=m_\varphi^2|\Phi|^2+\lambda|\Phi|^4-B_0m_uf_\pi^2(1+\sigma/f_a)\cos(a/f_a)$ with $|\Phi|^2=(f_a+\sigma)^2/2$, stationarity in $\sigma$ and the second derivatives reproduce $m_\varphi^2$ and $m_\sigma^2$ of the inputs and $m_a^2=B_0m_uf_\pi^2/f_a^2$.
2. **Polar vertices.** The kinetic term $\frac12(1+\sigma/f_a)^2(\partial a)^2$ contains $\frac{\sigma}{f_a}(\partial a)^2$, vertex $-2p_1\!\cdot p_2/f_a$. The cosine term gives $-\partial^3V/\partial\sigma\partial a^2=-B_0m_uf_\pi^2/f_a^3$ and the contact $-\partial^4V/\partial a^4=+B_0m_uf_\pi^2/f_a^4$. With on-shell ALPs ($2p_1\!\cdot p_2=q^2-2m_a^2$) the off-shell $\sigma aa$ vertex is
$$
\mathcal V_{\sigma aa}(q^2)=-\frac{1}{f_a}\left(q^2-\frac{B_0m_uf_\pi^2}{f_a^2}\right).
$$
3. **Cartesian vertices.** From $V(s,p)$ with $\lambda|\Phi|^4=\frac\lambda4[(f_a+s)^2+p^2]^2$ and the linear source $-B_0m_uf_\pi^2(1+s/f_a)$: $\mathcal V_{spp}=-2\lambda f_a$, $\mathcal V_{sss}=-6\lambda f_a$, $\mathcal V_{pppp}=\mathcal V_{ssss}=-6\lambda$, $\mathcal V_{sspp}=-2\lambda$, no $p^3$ or $s^2p$ vertex.
4. **Amplitudes** $aa\to aa$ in both coordinates, $s+t+u=4m_a^2$.
5. **S wave** of $aa\to aa$ with $t,u=-(s-4m_a^2)(1\mp\cos\theta)/2$, integrated symbolically.
6. **Coupled channels** from the Cartesian contacts; exact S waves of the exchange parts with general two-body kinematics, evaluated at large $s$.
7. **Phase-only potential** $V_{\rm ext}(a)$ (even) added to the polar potential; fixed-angle $s\to\infty$ limit.
8. **Meson vertices** from $f_\pi^2B_0\operatorname{Re}\langle\widehat MU^\dagger\rangle$ by direct expansion; $aP_3\to aP_3$ in both coordinates.

## Results

**A. Amplitude.** In polar coordinates
$$
\mathcal M(aa\to aa)=\frac{B_0m_uf_\pi^2}{f_a^4}+\sum_{x=s,t,u}\frac{\big(x-B_0m_uf_\pi^2/f_a^2\big)^2}{f_a^2\,(m_\sigma^2-x)}
=-6\lambda+4\lambda^2f_a^2\sum_{x=s,t,u}\frac{1}{m_\sigma^2-x},
$$
and the right-hand side is the Cartesian amplitude ($-6\lambda$ contact, exchange with $\mathcal V_{spp}=-2\lambda f_a$). The identity uses $s+t+u=4B_0m_uf_\pi^2/f_a^2$: the ALP mass in the kinematics must be the unmixed pole of the same Lagrangian. If the mass squared in the kinematics and in $2p_1\cdot p_2$ is shifted to $B_0m_uf_\pi^2/f_a^2+\delta$ while the vertices are kept (inconsistent at order $\delta$, which comes from mixing), then
$$
\mathcal M^{\rm polar}-\mathcal M^{\rm Cart}=\frac{8\delta}{f_a^2}+\sum_{x}\frac{(2\lambda f_a^2-2\delta)^2-4\lambda^2f_a^4}{f_a^2(m_\sigma^2-x)}\ \xrightarrow{s\to\infty}\ \frac{8\delta}{f_a^2},
$$
so a physical mass inserted by hand spoils the cancellation of the growing terms by a constant; a consistent treatment includes the mixing in the vertices as well.

**B. Energy dependence.** Exact S wave ($N_iN_f=\frac12$):
$$
a_0(aa\to aa)=\frac{1}{64\pi}\left[-12\lambda+\frac{8\lambda^2f_a^2}{m_\sigma^2-s}+\frac{16\lambda^2f_a^2}{s-4m_a^2}\ln\!\left(1+\frac{s-4m_a^2}{m_\sigma^2}\right)\right],\qquad m_a^2=\frac{B_0m_uf_\pi^2}{f_a^2}.
$$
For $4m_a^2\ll s\ll m_\sigma^2$ the polar form gives $\mathcal M\simeq\sum_xx^2/(f_a^2m_\sigma^2)$ and $a_0\simeq 5s^2/(96\pi f_a^2m_\sigma^2)$, which reaches $\frac12$ at $s^2=\frac{48\pi}{5}f_a^2m_\sigma^2$. This scale lies below $m_\sigma^2$ only if $\lambda>\frac{24\pi}{5}-\frac{B_0m_uf_\pi^2}{2f_a^4}$, i.e. only for a strongly coupled radial sector (it is relevant when $\sigma$ is integrated out). For $s\gg m_\sigma^2$ the $s$-channel term falls like $1/s$, the $t,u$ terms like $\ln s/s$, and
$$
a_0(aa\to aa)\to-\frac{3\lambda}{16\pi}\quad\Rightarrow\quad\lambda\le\frac{8\pi}{3}\ \text{(primary)},\qquad\lambda\le\frac{16\pi}{3}\ \text{(secondary)}.
$$

**C. Coupled channels** ($s\gg m_\sigma^2$, Cartesian fields). On $\{aa/\sqrt2,\sigma\sigma/\sqrt2\}$ and the separate channel $a\sigma$ (no mixing with the others since $p\to-p$ is a symmetry):
$$
a_0=\frac{\lambda}{16\pi}\begin{pmatrix}-3&-1\\-1&-3\end{pmatrix},\qquad a_0(a\sigma\to a\sigma)=-\frac{\lambda}{8\pi}.
$$
The eigenvalues are $-\lambda/(4\pi)$ and $-\lambda/(8\pi)$. All exchange parts of the four entries vanish asymptotically (checked numerically at $s=10^4,10^6,10^8$ GeV$^2$, decreasing and below $50\ln s/s$). Bounds, also written with $2\lambda f_a^2=m_\sigma^2-B_0m_uf_\pi^2/f_a^2$:

| channel | primary | secondary | primary as mass condition |
|---|---|---|---|
| $aa\to aa$ | $\lambda\le8\pi/3$ | $16\pi/3$ | $m_\sigma^2-B_0m_uf_\pi^2/f_a^2\le\frac{16\pi}{3}f_a^2$ |
| $\{aa,\sigma\sigma\}$ | $\lambda\le2\pi$ | $4\pi$ | $m_\sigma^2-B_0m_uf_\pi^2/f_a^2\le4\pi f_a^2$ |
| $a\sigma$ | $\lambda\le4\pi$ | $8\pi$ | $m_\sigma^2-B_0m_uf_\pi^2/f_a^2\le8\pi f_a^2$ |

The strongest condition is $m_\sigma^2\le4\pi f_a^2+B_0m_uf_\pi^2/f_a^2$.

**D. Phase-dependent potential.** In polar coordinates $|\Phi|^2=(f_a+\sigma)^2/2$, so $\lambda|\Phi|^4$ (any function of $|\Phi|$) gives no vertex with ALPs and no $a^4$ contact. The $a^4$ contact is $-\partial^4V_{\rm tot}/\partial a^4|_0$, from the QCD-induced cosine term ($+B_0m_uf_\pi^2/f_a^4$) and from any phase-only potential ($-V_{\rm ext}''''(0)$). The contact alone does not fix the high-energy limit: with $V_{\rm ext}$ included (and $m_a^2=B_0m_uf_\pi^2/f_a^2+V_{\rm ext}''(0)$),
$$
\mathcal M(aa\to aa)\xrightarrow{s\to\infty}-6\lambda-V_{\rm ext}''''(0)+\frac{8V_{\rm ext}''(0)}{f_a^2},
$$
because the growing $\sigma aa$ exchange no longer cancels exactly. Without $V_{\rm ext}$ only $-6\lambda$ remains.

**E. Mesons.** The $\sigma P_iP_j$ vertices from $\widehat S_{10}$ are momentum independent, so exchanges involving only them fall with $s$. In the polar frame the $t$-channel $\sigma$ exchange in $aP_3\to aP_3$ contains the growing $\sigma aa$ vertex; it cancels the polar contact $B_0m_u/f_a^2$ (from $\widehat S_{02}$) exactly:
$$
\mathcal M(aP_3\to aP_3)=\frac{2\lambda B_0m_u}{m_\sigma^2-t}\quad\text{(both coordinates)},
$$
which tends to $B_0m_u/f_a^2$ for $|t|\ll m_\sigma^2$ and $B_0m_uf_\pi^2/f_a^2\ll m_\sigma^2$, the contact of the frame without $\sigma$, and falls like $\ln s/s$ in the S wave for $s\gg m_\sigma^2$. The constant contact vertices are
$$
\mathcal V_{aP_3P_3P_3}=\frac{B_0m_u}{f_af_\pi},\qquad \mathcal V_{aP_3P_1P_1}=\mathcal V_{aP_3P_2P_2}=\frac{B_0m_u}{3f_af_\pi}=\mathcal V_{aP_3\pi^+\pi^-},
$$
giving $a_0(aP_3\to P_3P_3)=\sqrt2B_0m_u/(32\pi f_af_\pi)$ and $a_0(aP_3\to\pi^+\pi^-)=B_0m_u/(48\pi f_af_\pi)$, hence $f_\pi/f_a\le8\sqrt2\pi f_\pi^2/(B_0m_u)\simeq52.5$ and $f_\pi/f_a\le24\pi f_\pi^2/(B_0m_u)\simeq111$ (inputs as in the LO partial-wave record: $m_u=2.16$ MeV, $B_0=m_{\pi^0}^2/(m_u+m_d)=2.656$ GeV, $f_\pi=92.07$ MeV). These are weak; $P_3$ is the $\lambda_3$ field, not the physical $\pi^0$. Any $\sigma$ exchange in $aP_3\to P_3P_3$ uses only momentum-independent vertices.

**Choice in Murayama:2026ioh.** The quartic term is allowed but "mostly ignored", i.e. $\lambda=0$ with $a$ and $\sigma$ degenerate at $m_\varphi$ [@Murayama:2026ioh]. Then $\mathcal M(aa\to aa)=0$ and $\mathcal M(aP_3\to aP_3)=0$ at tree level with mixing neglected, and all bounds above are satisfied. With $\lambda>0$ the condition $m_\sigma^2-m_a^2\le4\pi f_a^2$ (unmixed $m_a^2$) applies.

## Verification

[`derive_murayama_alp_unitarity.py`](../../../../scripts/derive_murayama_alp_unitarity.py), tests in [`test_derive_murayama_alp_unitarity.py`](../../../../scripts/test_derive_murayama_alp_unitarity.py) (16 tests, all pass):

- vacuum, $m_\sigma^2$ and $m_a^2$ agree with the radial-mode decay record; polar vertices derived from the Lagrangian; $\lambda|\Phi|^4$ independent of $a$; Cartesian vertices;
- polar amplitude $=$ hypothesis form $=$ closed form $=$ Cartesian amplitude (exact symbolic); shifted-mass difference and its limit $8\delta/f_a^2$;
- integrated S wave equals the closed form at five values of $s$ (real, to $10^{-12}$); pole-integral helper against direct integration; $s\to\infty$ limit; $s^2$ term of the exact amplitude at $B_0\to0$ equals $\sum_xx^2/(f_a^2m_\sigma^2)$; intermediate wave and scale;
- contact matrix, eigenvalue bounds, mass condition; exchange parts below $50\ln s/s$ and decreasing, full matrices within $10^{-7}$ of the contact values at $s=10^8$ GeV$^2$ (sample point $f_a=0.5$ GeV, $\lambda=1.3$);
- $V_{\rm ext}$ limit; $aP_3\to aP_3$ polar $=$ Cartesian and vanishing at large $|t|$; three-meson vertices and bounds.

## Invocation

```sh
.venv/bin/python scripts/derive_murayama_alp_unitarity.py
.venv/bin/python -m unittest scripts.test_derive_murayama_alp_unitarity
```
