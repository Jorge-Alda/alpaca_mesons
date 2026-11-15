---
title: Leading-twist light-cone result for pseudoscalar transition form factors and its convention map
status: literature
created: 2026-11-15
updated: 2026-11-15
bibtex_keys: [Hoferichter:2020lap]
provenance: literature
---

# Source and scope

Read from `memory/sources/Hoferichter:2020lap/src/sections/{Introduction,BTT,BL}.tex`
(arXiv 2004.06127) [@Hoferichter:2020lap]. Labels are LaTeX labels. Order in $\alpha_s$:
leading order for the coefficient function; the source explicitly restricts to
leading order and cites NLO work (Phys. Rev. D 28 (1983) 524,
delAguila:1981nk) for $\alpha_s$ corrections (Introduction). The
light-cone factorization is credited to Lepage:1979zb and Phys. Rev. D 22 (1980) 2157
Phys. Rev. D 22 (1980) 2157.

# Conventions (BTT.tex)

- $\epsilon^{0123}=+1$ (line after label `TFF_pion`).
- Matrix element $\mathcal M^{\mu\nu}(p\to q_1,q_2)=i\int d^4x\,e^{iq_1\cdot x}\langle0|T\{j^\mu_{\rm em}(x)j^\nu_{\rm em}(0)\}|P(p)\rangle$
  and (label `TFF_pion`) $\mathcal M_{\mu\nu}=\epsilon_{\mu\nu\alpha\beta}q_1^\alpha q_2^\beta F_{P\gamma^*\gamma^*}(q_1^2,q_2^2)$
  [@Hoferichter:2020lap]. This is the project tensor convention.
- $j_{\rm em}^\mu=\bar q\hat Q\gamma^\mu q$, $\hat Q={\rm diag}(2,-1,-1)/3$, $q=(u,d,s)$.
  Normalization: $|F(0,0)|^2=4\Gamma_{\gamma\gamma}/(\pi\alpha^2m_P^3)$.

# Decay constants and wave functions (BL.tex)

$\langle0|\bar q\gamma_\mu\gamma_5\frac{\lambda^a}{2}q|P(p)\rangle=ip_\mu F_P^a$ with
$\lambda_0=\sqrt{2/3}\,\mathbb 1$; the wave function $\phi_P^a(u)$ is defined by
the nonlocal matrix element (label `wave_function`), normalized to 1.
Asymptotic form (label `phi_asym`): $\phi(u)=6u(1-u)$. Scale dependence
(labels `wave_function_expansion` and the following equations):
$$
\phi(u,\mu)=6u(1-u)\sum_{n}a_n(\mu)C_n^{3/2}(2u-1),\quad
a_n(\mu)=a_n(\mu_0)\Bigl(\frac{\alpha_s(\mu)}{\alpha_s(\mu_0)}\Bigr)^{\gamma_n^{(0)}/\beta_0},
$$
$\gamma_n^{(0)}=C_F\bigl(1-\frac{2}{(n+1)(n+2)}+4\sum_{m=2}^{n+1}\frac1m\bigr)$,
$\beta_0=\frac{11}{3}N_c-\frac23N_f$, $a_0=1$. For $P=P^c$ and $a\in\{0,1,3,4,6,8\}$ the
odd coefficients vanish [@Hoferichter:2020lap]. The source notes that the path-ordered
gauge factor is omitted. It gives no singlet anomalous dimension; the
scale dependence above is the flavor-nonsinglet one.

# Result (labels `weights`, `pion_TFF_result`, `eq:pQCD`, `eq:OPE_limit`, `eq:BL_limit`)

Flavor weights $C_a=\tfrac12{\rm Tr}(\hat Q^2\lambda^a)$:
$C_3=1/6$, $C_8=1/(6\sqrt3)$, $C_0=2/(3\sqrt6)$. Leading twist:
$$
F_{P\gamma^*\gamma^*}(q_1^2,q_2^2)=-4\sum_aC_aF_P^a\int_0^1du\,
\frac{\phi(u)}{uq_1^2+(1-u)q_2^2-u(1-u)m_P^2}.
$$
The source warns that the $m_P^2$ term is not a consistent treatment of mass
corrections (they need subleading light-cone terms). For the pion:
$F(q_1^2,q_2^2)=-\frac{2F_\pi}{3}\int du\,\phi_\pi/(uq_1^2+(1-u)q_2^2)+\mathcal O(q_i^{-4})$,
symmetric limit $F(q^2,q^2)=-2F_\pi/(3q^2)$, singly virtual $F(q^2,0)=-2F_\pi/q^2$
for the asymptotic $\phi$ [@Hoferichter:2020lap]. The source adds that the singly-virtual
Brodsky-Lepage limit goes beyond a strict OPE since the wave function resums
higher-order OPE terms (citing Gorsky:1987, Manohar:1990hu), and that the
formula is strictly justified from an OPE only when both virtualities are large.
Here $F_\pi=92.28(19)$ MeV.

# Map to the project hard-current formula (interpretation, by hand)

Using $\sum_fQ_f^2\bar q_f\gamma_\sigma\gamma_5q_f=2\sum_aC_a\,\bar q\gamma_\sigma\gamma_5\frac{\lambda^a}{2}q$
(unscripted check of the trace identity $\hat Q^2=\sum_aC_a\lambda^a$), the project formula
$F\supset\frac{2i}{(q^2+i0)m_a^2}\sum_fQ_f^2p^\sigma\langle0|\bar q_f\gamma_\sigma\gamma_5q_f|a\rangle$
gives $F=-4\sum_aC_aF_P^a/q^2$ at $q_1^2=q_2^2=q^2$, $m_P\to0$, which equals the
source's `pion_TFF_result` in the symmetric limit including sign. The
sign and the factor are therefore consistent, provided the same phase of the
meson state is used. The integral gives $1/q^2$ at equal virtualities and $3/q^2$
at $q_2^2=0$ for $\phi=6u(1-u)$, which reproduces the $2/\sqrt3\,F^8$ and
$4\sqrt{2/3}\,F^0$ coefficients of $Q^2F$ in
[Bali:2021qem](../Bali:2021qem/eta-eta-prime-tff-hard-inputs.md).

# Flags

- Sign: the source pion limit is $-2F_\pi/(3q^2)$, positive for spacelike
  $q^2$, whereas the project's fixed-phase pion scalar is negative (see
  [Knecht:2001qf](../Knecht:2001qf/pion-tff-short-distance.md)); the difference is
  the external-meson phase. For an ALP there is no such phase freedom: use the ALP
  matrix element defined in the project sources.
- No singlet scale dependence, no gluon operator, no heavy-flavor treatment and
  no $\alpha_s$ correction in this source.
