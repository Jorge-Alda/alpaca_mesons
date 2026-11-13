---
title: DiLuzio:2026poh model Lagrangian, field content and field normalization
status: literature
created: 2026-11-13
updated: 2026-11-13
bibtex_keys: [DiLuzio:2026poh]
provenance: literature
---

# Model Lagrangian, field content and field normalization

Equation numbers follow the order of numbered displayed equations in the arXiv source (counted from the TeX file, not from a compiled PDF).

## Field content and symmetry

The model adds one complex scalar $\varphi$, a Standard Model singlet with unit charge under a global $\mathrm{U}(1)_{\rm PQ}$, to QCD with light quarks $q=(u,d,s)$ of PQ charges $n_q=(n_u,n_d,n_s)$. Below the electroweak scale the Lagrangian is [@DiLuzio:2026poh, Eq. (3)]

$$
\mathscr{L}=-\tfrac14 G^a_{\mu\nu}G^{a,\mu\nu}+\theta\frac{g_s^2}{32\pi^2}G^a_{\mu\nu}\tilde G^{a,\mu\nu}+i\sum_q\bar q\slashed D q
+f_\varphi^2\Big[|\partial_\mu\varphi|^2-m_\varphi^2|\varphi|^2-\lambda f_\varphi^2|\varphi|^4\Big]-\sum_q m_q(\sqrt2\varphi)^{n_q}\bar q_Lq_R+\text{h.c.}
$$

- $\varphi$ is introduced as a dimensionless complex field; $f_\varphi$ carries the mass dimension. The factor $f_\varphi^2$ in front of the bracket makes the kinetic term $f_\varphi^2|\partial\varphi|^2$ ("for later convenience", Sec. 2.1).
- $m_q$ is made real and positive by field re-phasing. If $n_q=0$, $m_q$ is the ordinary quark mass; if $n_q\neq0$ the quark couples to $\varphi$.
- $|n_q|>1$ gives a non-renormalizable interaction (with $\varphi\to\varphi^*$ if $n_q<0$). The model of Murayama (arXiv 2601.02465, cited in the paper, no BibTeX key in this repository) is $(n_u,n_d,n_s)=(1,0,0)$.
- $\mathrm{U}(1)_{\rm PQ}$ is anomalous when $\sum_qn_q\ne0$ [@DiLuzio:2026poh, Eq. (5)]:
$$
\frac EN=\frac23\frac{4n_u+n_d+n_s}{\sum_qn_q},\qquad N_{\rm DW}=2N=\Big|\sum_qn_q\Big|.
$$

## Symmetry breaking

The quark condensate $\langle\bar q_{Li}q_{Rj}\rangle=-\delta_{ij}B_0f_\pi^2/2$ induces a vacuum expectation value of $\varphi$. The parameter $m_\varphi^2$ is fixed so that $\langle\varphi\rangle=1/\sqrt2$ (Sec. 2.1), and the PQ-charged quarks acquire the mass [@DiLuzio:2026poh, Eq. (4)]
$$
m_q=\kappa_q\frac{f_\varphi}{\sqrt2}.
$$
The paper assumes $m_\varphi^2>0$ and $\lambda>0$, so PQ is broken only after confinement. A negative $m_\varphi^2$ would give a DFSZ-like model (footnote in Sec. 2.1).

## Chiral Lagrangian at leading order

With $U=\exp(i\sqrt2\,\Pi/f_\pi)\to LUR^\dagger$ under $\mathrm{U}(3)_L\otimes\mathrm{U}(3)_R$ and $f_\pi\approx92$ MeV [@DiLuzio:2026poh, Eq. (6)]:
$$
\mathscr{L}_{\pi\varphi}=f_\varphi^2\Big[|\partial_\mu\varphi|^2-m_\varphi^2|\varphi|^2-\lambda f_\varphi^2|\varphi|^4\Big]
+\frac{f_\pi^2}4\mathrm{tr}[D_\mu UD^\mu U^\dagger]+\frac{f_\pi^2B_0}2\mathrm{tr}[M_qU^\dagger+M_q^\dagger U]
-\frac{m_1^2f_\pi^2}{12}\Big[\Big(-\tfrac i2\ln\det U+\text{h.c.}\Big)-\theta\Big]^2 .
$$
The last term is the QCD anomaly at leading order in $1/N$. The quark mass matrix is promoted to a field-dependent spurion, $M_q\to LM_qR^\dagger$ [@DiLuzio:2026poh, Eq. (7)]:
$$
M_q=\mathrm{diag}\big[m_u(\sqrt2\varphi)^{n_u},\,m_d(\sqrt2\varphi)^{n_d},\,m_s(\sqrt2\varphi)^{n_s}\big].
$$
Conventions are those of Pich:1991fq (cited in the paper, no BibTeX key in this repository), with $B_0=m_{\pi^\pm}^2/(m_u+m_d)$. The $\eta'$ is kept in $U$ and $f_{\eta'}\approx1.3f_\pi$ is replaced by $f_\pi$ in the leading-order chiral Lagrangian.

## Meson field normalization (before canonical rescaling)

The meson matrix is written for the dimensionless combination $\Pi/f_\pi$ [@DiLuzio:2026poh, Eq. (8)]:
$$
\frac\Pi{f_\pi}=\begin{pmatrix}
\frac{\pi^0}{\sqrt2}+\frac{\eta}{\sqrt6}+\frac{\eta'}{\sqrt3}&\pi^+&K^+\\
\pi^-&-\frac{\pi^0}{\sqrt2}+\frac{\eta}{\sqrt6}+\frac{\eta'}{\sqrt3}&K^0\\
K^-&\bar K^0&-\frac{2\eta}{\sqrt6}+\frac{\eta'}{\sqrt3}
\end{pmatrix}.
$$
The fields $\pi^0,\eta,\eta',\dots$ in this matrix are dimensionless and "related to the canonically normalized fields by $\pi^0_{\rm can}=f_\pi\pi^0$, etc." (text below Eq. (8)). The same convention applies to $\varphi$ and $a$; see [angular-mode-parametrization.md](angular-mode-parametrization.md).

## Ambiguities and assumptions

- The $\theta$ term in Eq. (6) expands to $\frac{m_1^2f_\pi^2}2(\eta'-\theta/\sqrt6)^2$ for $\det U=e^{i\sqrt6\eta'}$, whereas Eq. (14) prints $(\eta'-\theta)^2$. This is an apparent mismatch of a factor $\sqrt6$ in the $\theta$ term (checked by the extractor from Eq. (6)); it does not affect results at $\theta=0$.
- Equations (6)-(7) assume the leading-order chiral Lagrangian, $f_{\eta'}\to f_\pi$ and the large-$N$ anomaly term; the paper states the chiral description is valid only for $E\lesssim$ GeV.
- Electromagnetic corrections are neglected.
