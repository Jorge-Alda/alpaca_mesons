---
title: Tree-level partial-wave unitarity of the LO ALP-meson theory at finite $f_\pi/f_a$
status: mixed
created: 2026-11-21
updated: 2026-11-21
bibtex_keys: [Alda:2026nrd, Chivukula:1992gi, Lee:1977eg]
provenance: analytical-derivation
---

# Tree-level partial-wave unitarity of the LO ALP-meson theory at finite $f_\pi/f_a$

## Scope, conventions and assumptions

LO chiral Lagrangian with an ALP, tree level, in the frame without ALP current sources of the [derivative-coupling removal record](../../../../chiral-transformations/source-invariants/dressed-sources/finite-mixing/alp-dressed-frame-derivative-coupling-removal.md) (on-shell equivalent to the original frame, exactly in $f_\pi/f_a$). $U=\exp(i\Phi/f_\pi)$, $\operatorname{tr}(\lambda_i\lambda_j)=2\delta_{ij}$, $\chi=2B_0\widehat M(x)$, $x=a/f_a$, $\widehat M_n=d^n\widehat M/dx^n|_0$ [@Alda:2026nrd]. Pure QCD sources only: no photons, no weak spurions. Sections 1-2 include only the bare ALP mass term outside the chiral Lagrangian; Section 3 adds a general external ALP potential $V_{\rm ext}(a)$. Natural units, GeV.

Partial waves for $2\to2$ scattering of spinless particles:
$$
(a_J)_{fi}=\frac{N_iN_f}{32\pi}\int_{-1}^{1}d\cos\theta\,P_J(\cos\theta)\,\mathcal M_{fi},\qquad N=\tfrac1{\sqrt2}\ \text{for two identical particles, else }1 .
$$
For $\pi^0\pi^0$ or the isospin channels this is the $1/(64\pi)$ normalization of Chivukula:1992gi, Eq. (2.4) ([extract](../../../../../extracts/Chivukula:1992gi/su-n-goldstone-unitarity-bound.md)). Elastic unitarity $\operatorname{Im}a_J=\sigma|a_J|^2$, $\sigma=\sqrt{1-4m^2/s}\to1$, gives the primary condition $|\operatorname{Re}a_0|\le\frac12$ for $s\gg4m^2$; the secondary condition $|a_0|\le1$ is the criterion used in Chivukula:1992gi. The argument generalizes Lee:1977eg [@Lee:1977eg; @Chivukula:1992gi]. $f$ of Chivukula:1992gi equals $f_\pi$ here: $\Sigma=\exp(2i\pi^aT^a/f)$ with $T^a=\lambda_a/2$ is $\exp(i\Phi/f)$, and both use $\frac{f^2}4\langle\partial U\partial U^\dagger\rangle$.

Validity (applies to all results): LO only; away from $s$-channel poles; external states are stable poles; the derived scales are where tree-level LO amplitudes stop being consistent with unitarity, not predictions of where new states appear.

## 1. Growth with energy (H3)

**Claim.** In this frame the only vertices that grow with momenta are those of $\frac{f_\pi^2}4\langle u_\mu u^\mu\rangle$ with mesons only; cubic vertices are momentum independent, so exchange diagrams are $O(s^0)$. The four-meson vertex is
$$
\mathcal V_{i_1i_2i_3i_4}=-\frac1{24f_\pi^2}\sum_{\sigma\in S_4}\Big[\langle\lambda_{\sigma1}\lambda_{\sigma2}\lambda_{\sigma3}\lambda_{\sigma4}\rangle-\langle\lambda_{\sigma1}\lambda_{\sigma3}\lambda_{\sigma2}\lambda_{\sigma4}\rangle\Big]\,p_{\sigma2}\cdot p_{\sigma4},
$$
from $\mathcal L_4=\langle[\Phi,\partial_\mu\Phi][\Phi,\partial^\mu\Phi]\rangle/(48f_\pi^2)$ (checked against the CCWZ generator). For $s,|t|\gg$ all masses, the amplitude among poles $n_1..n_4$ is this vertex contracted with the meson components $v^{(n)}_i$ of the pole vectors. With $K=\mathbf 1$ the pole vectors form an orthogonal matrix $O$ on $(P_i,a)$.

**Eigenvalue invariance.** Let $P_{IJKL}$ be the coefficient of $s$ in the S-wave projection ($\langle t\rangle=\langle u\rangle=-s/2$), padded with zeros whenever an index is the ALP. On normalized symmetric two-particle states the $a_0$ matrix is $s/(32\pi)$ times the restriction of $P$ to the symmetric subspace, which transforms with $O\otimes O$. Its eigenvalues are therefore those of the meson-only matrix padded with zeros, for any $O$: the leading growth, and hence $\Lambda_{\rm unit}$, does not depend on $f_\pi/f_a$ or on the ALP couplings.

**Largest eigenvalue.** In units of $s/(32\pi f_\pi^2)$: U(3) (45 states) and SU(3) (36 states) give largest eigenvalue $3$, i.e. $a_0^{\max}=N s/(32\pi f_\pi^2)$ with $N=3$ (the $\eta_0$ is free in the chiral limit); SU(2) gives 2. With a common meson mass $m$ the singlet state $\sum_a|aa\rangle/\sqrt{N^2-1}$ gives exactly
$$
a_0^{\rm singlet}=\frac{Ns}{32\pi f_\pi^2}-\frac{m^2}{16\pi Nf_\pi^2}\qquad(N=2,3),
$$
reproducing Chivukula:1992gi, Eq. (2.5), $a_{\Delta0}$; for $N=2$ it is the $\pi\pi$ $I=0$ wave $(2s-m^2)/(32\pi f_\pi^2)$. For U(3) with a common mass the mass correction of the largest eigenvalue was not computed (the trace state $\sum_{a=0}^8|aa\rangle/3$ is not an eigenchannel; it gives $s/(12\pi f_\pi^2)+m^2/(48\pi f_\pi^2)$).

**Unitarity scale.** From $a_0^{\rm singlet}=\frac12$ (primary) and $=1$ (secondary):
$$
\Lambda_{\rm unit}^2=\frac{16\pi f_\pi^2}{N}+\frac{2m^2}{N^2},\qquad
\Lambda_{\rm unit}'^2=\frac{32\pi f_\pi^2}{N}+\frac{2m^2}{N^2}.
$$
The mass term assumes degenerate mesons (approximation, as in Chivukula:1992gi). Input $f_\pi=92.07$ MeV; $m=m_{\pi^0}=134.977$ MeV for the illustrative mass term:

| channel | $\vert\operatorname{Re}a_0\vert\le\frac12$ | $\vert a_0\vert\le1$ |
|---|---|---|
| U(3)/SU(3) singlet, $N=3$, $m=0$ | 376.9 MeV | 533.0 MeV |
| same, $m=m_{\pi^0}$ | 382.2 MeV | 536.8 MeV |
| $\pi\pi$ $I=0$, $N=2$, $m=0$ | 461.6 MeV | 652.8 MeV |
| same, $m=m_{\pi^0}$ | 471.3 MeV | 659.7 MeV |

The $f_\pi/f_a$-independence holds for the $O(s)$ term only; the constant terms below depend on $f_\pi/f_a$.

## 2. Constant contact terms (H4)

**General vertices.** Since $\frac{f_\pi^2}4\langle\chi_+\rangle=f_\pi^2B_0\operatorname{Re}\langle\widehat M(x)U^\dagger\rangle$, the vertex with $n$ ALPs and mesons $P_{i_1}\dots P_{i_k}$ is
$$
\mathcal V_{a^nP_{i_1}\cdots P_{i_k}}=\frac{f_\pi^2B_0}{f_\pi^kf_a^n}\operatorname{Re}\Big[\frac{(-i)^k}{k!}\sum_{\sigma\in S_k}\big\langle\widehat M_n\lambda_{i_{\sigma1}}\cdots\lambda_{i_{\sigma k}}\big\rangle\Big],
$$
plus, from $-\frac{f_\pi^2m_0^2}{4N_f}\big(\langle\Phi\rangle/f_\pi+2\widehat\Theta(x)\big)^2$ (convention of the [joint-vertex record](../../tree-vertices/finite-mixing/alp-lo-joint-vertex-expansion.md)), terms with $k\le2$ only. Examples:
$$
\mathcal V_{aaaa}=\frac{f_\pi^2B_0}{f_a^4}\operatorname{Re}\operatorname{tr}\widehat M_4,\qquad
\mathcal V_{aaP_iP_j}=-\frac{B_0}{2f_a^2}\operatorname{Re}\operatorname{tr}\big(\widehat M_2\{\lambda_i,\lambda_j\}\big).
$$
For $n=0,k=2$ this reproduces $-M^2_{ij}$. The topological term adds no vertex with three or more legs when $\widehat\Theta$ is linear in $x$, which is the case for constant $c_q^A$ and $\Theta=c_Gx$ ($\widehat\Theta=c_{GA}x$). With non-zero $\widehat\Theta_2,\widehat\Theta_3$ it adds e.g. $\mathcal V_{aaaa}\supset-\frac{f_\pi^2m_0^2}{N_ff_a^4}\sum_r\binom4r\widehat\Theta_r\widehat\Theta_{4-r}$.

**Cubic vertices.** For CP-conserving sources (real symmetric $c_q^{A,P}$, real $m_q$), $\overline{\widehat M(x)}=\widehat M(-x)$, so $\widehat M_n$ is real for even $n$ and imaginary for odd $n$. For real symmetric generators ($\lambda_{0,3,8,6}$) the vertex vanishes unless $n+k$ is even. Hence all cubic vertices vanish (checked including a $d$-$s$ entry), and in the diagonal neutral sector the tree amplitudes in this frame have no exchange diagrams.

**S waves with mass mixing neglected (approximation).** Treating the external ALP and meson poles as the fields $a$ and $P_i$:
$$
a_0(aa\to aa)=\frac{f_\pi^2B_0\operatorname{Re}\operatorname{tr}\widehat M_4}{32\pi f_a^4},\quad
a_0(aP_i\to aP_i)=-\frac{B_0\operatorname{Re}\operatorname{tr}(\widehat M_2\lambda_i^2)}{16\pi f_a^2},\quad
a_0(P_iP_i\to aa)=-\frac{B_0\operatorname{Re}\operatorname{tr}(\widehat M_2\lambda_i^2)}{32\pi f_a^2}.
$$
These are energy independent, so $|\operatorname{Re}a_0|\le\frac12$ gives (for $V_{\rm ext}''''(0)=0$; with an external quartic the $aa\to aa$ bound is modified as in Section 3, the other two are unchanged)
$$
\Big(\frac{f_\pi}{f_a}\Big)^4\le\frac{16\pi f_\pi^2}{B_0|\operatorname{Re}\operatorname{tr}\widehat M_4|},\qquad
\Big(\frac{f_\pi}{f_a}\Big)^2\le\frac{8\pi f_\pi^2}{B_0|\operatorname{Re}\operatorname{tr}(\widehat M_2\lambda_i^2)|}\ (aP_i),\qquad
\Big(\frac{f_\pi}{f_a}\Big)^2\le\frac{16\pi f_\pi^2}{B_0|\operatorname{Re}\operatorname{tr}(\widehat M_2\lambda_i^2)|}\ (P_iP_i\to aa).
$$
When the approximation holds: the dressed-frame mixing $M^2_{ia}=\frac{f_\pi}{2f_a}\operatorname{tr}(\lambda_i\{c_q^{AP},\chi_0\})+\frac{f_\pi m_0^2c_{GA}}{N_ff_a}\operatorname{tr}\lambda_i$ (for $[c_q^{VS},m_q]=0$) must satisfy $|M^2_{ia}|\ll|M^2_{aa}-M^2_{ii}|$, i.e. small $f_\pi/f_a$ or an ALP mass far from the meson masses; in addition $4m^2\ll s\lesssim\Lambda_{\rm unit}^2$ so that $\sigma\simeq1$ and the LO theory applies. Single channels only; the coupled-channel matrix including the growing $PP\to PP$ entries was not analysed.

**Linear-spurion specialization.** For $M(x)=m_q+ix\{c_q^P,m_q\}$, constant $c_q^A$, all diagonal ($c_q^A=\operatorname{diag}(a_q)$, $c_q^P=\operatorname{diag}(p_q)$):
$$
\widehat M_n=m_q(-2i)^n\,a_q^{\,n-1}(a_q-n\,p_q),
$$
so $\operatorname{tr}\widehat M_4=16\sum_qm_qa_q^3(a_q-4p_q)$ and $\operatorname{tr}(\widehat M_2\lambda_3^2)=-4\sum_{q=u,d}m_qa_q(a_q-2p_q)$. **Finding:** for $n\ge2$, $\widehat M_n$ depends on $c_q^A$ and $c_q^P$ separately, not only on $c_q^{AP}=c_q^A-c_q^P$ (the two linear-spurion Lagrangians with equal $c_q^{AP}$ differ at $O(x^2)$; the redefinition relating them makes the spurion non-linear). Expressing the result through $c_q^{AP}=\operatorname{diag}(c_u,c_d,c_s)$ alone requires choosing a representative. For $c_q^P=0$ (only derivative couplings, $c_q^{AP}=c_q^A$), $\widehat M_n=m_q(-2ic_q)^n$ and
$$
a_0(aa\to aa)=\frac{B_0}{2\pi f_\pi^2}\Big(\frac{f_\pi}{f_a}\Big)^4\sum_qm_qc_q^4,\quad
a_0(a\pi_3\to a\pi_3)=\frac{B_0(m_uc_u^2+m_dc_d^2)}{4\pi f_\pi^2}\Big(\frac{f_\pi}{f_a}\Big)^2,\quad
a_0(\pi_3\pi_3\to aa)=\tfrac12a_0(a\pi_3\to a\pi_3),
$$
$$
\frac{f_\pi}{f_a}\le\Big(\frac{\pi f_\pi^2}{B_0\sum_qm_qc_q^4}\Big)^{1/4},\qquad
\frac{f_\pi}{f_a}\le\Big(\frac{2\pi f_\pi^2}{B_0(m_uc_u^2+m_dc_d^2)}\Big)^{1/2},\qquad
\frac{f_\pi}{f_a}\le\Big(\frac{4\pi f_\pi^2}{B_0(m_uc_u^2+m_dc_d^2)}\Big)^{1/2}.
$$
For $c_q^A=0$ (only the linear pseudoscalar spurion) $\widehat M_n=0$ for $n\ge2$ and these contact bounds are absent at this order.

**Illustration** (inputs: $c_u=1$, $c_d=c_s=0$, $c_q^P=0$; $m_u=2.16$ MeV, $m_d=4.70$ MeV, $B_0=m_{\pi^0}^2/(m_u+m_d)=2.656$ GeV, $f_\pi=92.07$ MeV; no external ALP potential beyond the mass term, $V_{\rm ext}''''(0)=0$, see Section 3): $f_\pi/f_a\le1.47$ ($aa\to aa$), $3.05$ ($a\pi_3\to a\pi_3$), $4.31$ ($\pi_3\pi_3\to aa$). Here $P_3$ is the $\lambda_3$ field, not the physical $\pi^0$. At $f_\pi/f_a\simeq1.5$ the mass mixing $|M^2_{3a}|\sim(f_\pi/f_a)\,2B_0m_u\approx0.017$ GeV$^2\approx m_\pi^2$, so the no-mixing approximation requires $|m_a^2-m_\pi^2|\gg0.02$ GeV$^2$; choosing $c_G=-1$ ($c_{GA}=0$) removes the topological mixing with $\eta_0$.

## 3. External ALP potential (H5)

**Setup (assumption, stated).** Besides the LO chiral Lagrangian, the ALP may have a potential $V_{\rm ext}(a)$ not induced by QCD, a general function of $a$ only (no mesons). It contains the bare mass term $\frac12m_{a,0}^2a^2$ used elsewhere in this project. The total ALP potential is
$$
V_{\rm tot}(a)=V_{\rm QCD}(a)+V_{\rm ext}(a),\qquad
V_{\rm QCD}(a)=-f_\pi^2B_0\operatorname{Re}\operatorname{tr}\widehat M(a/f_a)+\frac{f_\pi^2m_0^2}{N_f}\,\widehat\Theta(a/f_a)^2,
$$
i.e. $V_{\rm QCD}=-\mathcal L_{\rm LO}|_{P=0}$ in the dressed frame ($\frac{f_\pi^2}4\langle\chi_+\rangle=f_\pi^2B_0\operatorname{Re}\langle\widehat MU^\dagger\rangle$ and $X_{TC}=-2\widehat\Theta$ at $P=0$). The vacuum is assumed at $a=0$: $V_{\rm tot}'(0)=0$ (for CP-conserving sources $V_{\rm QCD}'(0)=0$, so this requires $V_{\rm ext}'(0)=0$). Mass mixing with mesons is neglected (same approximation as Section 2).

**Feynman rules and normalization.** Pure-ALP vertices are $\mathcal V_{a^n}=\partial^n\mathcal L/\partial a^n|_0=-V_{\rm tot}^{(n)}(0)$; the Feynman rule is $i\mathcal V$, the propagator $i/(q^2-m_a^2)$ and $i\mathcal M$ is the sum of diagrams. The S wave uses $(a_0)_{fi}=\frac{N_iN_f}{32\pi}\int_{-1}^1d\cos\theta\,\mathcal M_{fi}$ with $N_i=N_f=1/\sqrt2$ for the $aa$ states.

**Contact term.** With $V_{\rm tot}'''(0)=0$: $\mathcal M(aa\to aa)=-V_{\rm tot}''''(0)$ and
$$
a_0(aa\to aa)=\frac12\cdot\frac{2\,\mathcal M}{32\pi}=-\frac{1}{32\pi}\frac{d^4V_{\rm tot}}{da^4}\Big|_0,\qquad
|\operatorname{Re}a_0|\le\tfrac12\iff\Big|\frac{d^4V_{\rm tot}}{da^4}\Big|_0\Big|\le16\pi .
$$
**Confirmed** (sign and identical-particle factor). For $V_{\rm ext}''''=0$ and linear $\widehat\Theta$, $-V_{\rm QCD}''''(0)=f_\pi^2B_0\operatorname{Re}\operatorname{tr}\widehat M_4/f_a^4=\mathcal V_{aaaa}$, reproducing $a_0(aa\to aa)$ of Section 2. For non-linear $\widehat\Theta$ the script gives $-V_{\rm QCD}''''(0)=\frac{f_\pi^2B_0}{f_a^4}\operatorname{Re}\operatorname{tr}\widehat M_4-\frac{f_\pi^2m_0^2}{N_ff_a^4}\sum_{r=0}^4\binom4r\widehat\Theta_r\widehat\Theta_{4-r}$, in agreement with Section 2.

**ALP mass.** $M^2_{aa}=V_{\rm QCD}''(0)+V_{\rm ext}''(0)$; with $V_{\rm ext}=\frac12m_{a,0}^2a^2$ this is the $M^2_{aa}$ of `quadratic_blocks` (checked). Only $V_{\rm ext}''''(0)$ enters the contact bound; $V_{\rm ext}''(0)$ changes the ALP mass and, through it, the mixing, but not the contact vertex. Higher derivatives $V^{(n\ge5)}$ do not enter $2\to2$ at tree level.

**Cubic term and ALP exchange.** A cubic $V_{\rm tot}'''(0)\neq0$ requires CP violation (the ALP is CP odd) or a $V_{\rm ext}$ without definite parity in $a$; for CP-conserving QCD sources $V_{\rm QCD}'''(0)=0$, so then $V_{\rm tot}'''=V_{\rm ext}'''$. It adds $s$-, $t$-, $u$-channel ALP exchange:
$$
\mathcal M(aa\to aa)=-V_{\rm tot}''''-\big(V_{\rm tot}'''\big)^2\Big[\frac1{s-m_a^2}+\frac1{t-m_a^2}+\frac1{u-m_a^2}\Big],\quad u=4m_a^2-s-t,
$$
$$
a_0(aa\to aa)=\frac1{32\pi}\Big[-V_{\rm tot}''''-\frac{(V_{\rm tot}''')^2}{s-m_a^2}+\frac{2(V_{\rm tot}''')^2}{s-4m_a^2}\ln\frac{s-3m_a^2}{m_a^2}\Big]\qquad(s>4m_a^2),
$$
(derivatives at $a=0$). The exchange part decreases like $\ln s/s$, so for $s\gg m_a^2$ and $s\gg(V_{\rm tot}''')^2/|V_{\rm tot}''''|$ the bound $|V_{\rm tot}''''|\le16\pi$ is unchanged; at $s=4m_a^2$ the wave is $\frac1{32\pi}[-V_{\rm tot}''''+\frac53(V_{\rm tot}''')^2/m_a^2]$, but there $\sigma\neq1$ and the bound $|\operatorname{Re}a_0|\le\frac12$ is not the right condition. The bounds below assume $V_{\rm tot}'''(0)=0$ or the large-$s$ regime.

**Mixed channels.** $V_{\rm ext}$ contains no mesons, so the contact vertices $aP\to aP$ and $PP\to aa$ and their S waves are unchanged. Exchange diagrams with $V_{\rm ext}'''$ in these channels need an $aaP$ or $PPa$ vertex, which vanishes for CP-conserving QCD sources (Section 2). With mass mixing included, $V_{\rm ext}''$ changes the pole vectors and thus all channels; not analysed.

**General bound.** For linear $\widehat\Theta$ and $V_{\rm tot}'''=0$:
$$
\Big|\frac{f_\pi^2B_0\operatorname{Re}\operatorname{tr}\widehat M_4}{f_a^4}-V_{\rm ext}''''(0)\Big|\le16\pi,
$$
with $V_{\rm ext}''''$ dimensionless.

**Specialization $c_q^P=0$, diagonal $c_q^A=\operatorname{diag}(c_u,c_d,c_s)$.** $\operatorname{tr}\widehat M_4=16\sum_qm_qc_q^4\ge0$, so the condition reads $\big|16B_0\sum_qm_qc_q^4\,(f_\pi/f_a)^4/f_\pi^2-V_{\rm ext}''''(0)\big|\le16\pi$. Solving for $f_\pi/f_a>0$:
- $V_{\rm ext}''''(0)<-16\pi$: no solution (the external quartic alone violates the bound);
- $|V_{\rm ext}''''(0)|\le16\pi$: $\displaystyle\frac{f_\pi}{f_a}\le\Big(\frac{f_\pi^2\,[16\pi+V_{\rm ext}''''(0)]}{16B_0\sum_qm_qc_q^4}\Big)^{1/4}$, which for $V_{\rm ext}''''=0$ is the bound $(\pi f_\pi^2/(B_0\sum_qm_qc_q^4))^{1/4}$ of Section 2;
- $V_{\rm ext}''''(0)>16\pi$: $\displaystyle\Big(\frac{f_\pi^2\,[V_{\rm ext}''''(0)-16\pi]}{16B_0\sum_qm_qc_q^4}\Big)^{1/4}\le\frac{f_\pi}{f_a}\le\Big(\frac{f_\pi^2\,[V_{\rm ext}''''(0)+16\pi]}{16B_0\sum_qm_qc_q^4}\Big)^{1/4}$.

The last case is a cancellation between the QCD quartic ($-V_{\rm QCD}''''>0$ here) and a positive $V_{\rm ext}''''$; the allowed window in $(f_\pi/f_a)^4$ has relative width $32\pi/V_{\rm ext}''''(0)$, so large $f_\pi/f_a$ at large $V_{\rm ext}''''$ requires tuning $V_{\rm ext}''''$ against the QCD term. A cosine potential $\Lambda^4[1-\cos(a/f)]$ has $V''''(0)=-\Lambda^4/f^4<0$ and adds to the QCD term instead of cancelling it (the QCD potential for $c_q^P=0$, $-f_\pi^2B_0\sum_qm_q\cos(2c_qa/f_a)$, is of the same type). Only a positive quartic (for example a bare $\lambda_aa^4/4!$ with $\lambda_a>0$) can compensate.

**a-independent terms.** Terms of $V_{\rm ext}$ that do not depend on $a$ do not contribute. Example: if the ALP is the phase of a complex field, $\phi=(v+\sigma)e^{ia/v}/\sqrt2$, then $-\lambda|\phi|^4=-\lambda(v+\sigma)^4/4$ and $\partial_a^4$ of it vanishes at fixed $\sigma$ (checked symbolically). The exchange of $\sigma$ through the couplings of the kinetic term $\frac12(v+\sigma)^2(\partial a)^2/v^2$ is a separate, momentum-dependent effect, treated in the Murayama-model record `memory/alp/radial-mode/couplings-and-decays/murayama-alp-unitarity.md` (being written in parallel).

**Illustration** (inputs of Section 2, $c_u=1$, $c_q^P=0$): allowed $f_\pi/f_a$ is $[0,1.468]$ for $V_{\rm ext}''''=0$; $[0,1.624]$ for $8\pi$; $[0,1.234]$ for $-8\pi$; $\{0\}$ for $-16\pi$; $[1.464,1.930]$ for $V_{\rm ext}''''=100$ (cancellation window). These are illustrative numbers, not predictions.

## Verification

Script [`derive_alp_dressed_frame_unitarity.py`](../../../../../../scripts/derive_alp_dressed_frame_unitarity.py); tests `H3HighEnergy` and `H4ContactTerms` in `scripts/test_derive_alp_dressed_frame_unitarity.py`: closed four-meson vertex vs CCWZ generator; momentum independence of ALP vertices; largest eigenvalues (U(3), SU(3): 3; SU(2): 2) to $10^{-10}$; exact singlet waves with common mass for SU(2), SU(3); eigenvalue invariance for an exact rational Cayley orthogonal $O$ and for the dressed pole vectors at $f_\pi/f_a=1,3$; full pole amplitude at $s=2000,8000$ GeV$^2$ minus the four-meson contraction stays constant (the $O(s^0)$ part); closed-form constant and topological vertices vs CCWZ generator (non-diagonal CP-conserving sources); vanishing cubic vertices; pole amplitudes without mixing equal the contact formulas; bound formulas and numbers; dependence of $\widehat M_2$ on $c_q^A,c_q^P$ separately. Test class `H5ExternalPotential`: $-\partial_a^nV_{\rm QCD}(0)$ ($n=3,4$, non-linear $\widehat\Theta$, non-diagonal sources) equals the $n$-ALP vertex of the joint-vertex code, and $M^2_{aa}=V_{\rm tot}''(0)$ with $V_{\rm ext}=\frac12m_{a,0}^2a^2$; the closed-form amplitude with exchange equals `projected_amplitude` for a one-field theory with vertices $-V_{\rm tot}''',-V_{\rm tot}''''$ (sign of the Feynman rule); the closed-form S wave equals the numerical projection with $N=1/\sqrt2$ (to $10^{-25}$) and tends to $-V_{\rm tot}''''/(32\pi)$; the bound interval for both signs of the QCD coefficient, the cancellation window, and $1.468$ at $V_{\rm ext}''''=0$; $\partial_a^4(-\lambda|\phi|^4)=0$ in the phase model. All pass.

```sh
.venv/bin/python scripts/derive_alp_dressed_frame_unitarity.py
.venv/bin/python -m unittest discover -s scripts -p test_derive_alp_dressed_frame_unitarity.py
```
